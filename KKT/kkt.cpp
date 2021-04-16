//
//  kkt.cpp
//  KKT
//
//  Created by Cheng Lu on 4/2/20.
//  Copyright © 2020 Cheng Lu. All rights reserved.
//

//  Implementation of kkt.hpp

#include <chrono>
#include <cmath>
#include <iostream>
#include "kkt.hpp"
#include "utils.hpp"

// Auxiliary function.
// Overload function of pow, to avoid the pow(0, 0) domain error.
// In our implementation, we always assume pow(*, 0) = 1.
static inline data_type Pow(data_type base, data_type exponent) {
    if (exponent == 0) return 1.0;
    return pow(base, exponent);
}

// For piecewise functions.
// Use to determine which piece a variable belongs to.
// Find the interval [\lambda_l, lambda_r) that contains x.
// Implying right-subgradients.
static inline int getPQIndex(int pwDeg, const data_type* pw, int bkpNum, int stIndex,
                      const data_type& x) {
    assert(pwDeg == 1 || pwDeg == 2);  // Only support piecewise l1 and l2 for now.
    assert(bkpNum >= 0 && stIndex >= 0);
    if (bkpNum == 0) return 0;
    int head = 0, tail = bkpNum;
    while (head < tail) {
        int mid = (head + tail) / 2;
        data_type lambda = pw[stIndex + pwDeg + (pwDeg + 1) * mid];
        if (x < lambda) {
            if (mid == 0 || x >= pw[stIndex + pwDeg + (pwDeg + 1) * (mid - 1)]) {
                return mid;
            } else {
                tail = mid;
            }
        } else {
            if (mid == bkpNum - 1 || x < pw[stIndex + pwDeg + (pwDeg + 1) * (mid + 1)]) {
                return mid + 1;
            } else {
                head = mid + 1;
            }
        }
    }
    return head;
}

static inline bool pwValid(int pwDeg, const data_type* pw, int bkpNum, int stIndex,
                    const data_type& x, int pwIndex) {
    bool b1 = (pwIndex >= 0 && pwIndex <= bkpNum);
    bool b2 = true;
    if (bkpNum) {
        if (pwIndex == 0) {
            b2 = x < pw[stIndex + pwDeg];
        } else if (pwIndex == bkpNum) {
            b2 = x >= pw[stIndex + pwDeg + (pwDeg + 1) * (pwIndex - 1)];
        } else {
            b2 = (x >= pw[stIndex + pwDeg + (pwDeg + 1) * (pwIndex - 1)] &&
                  x < pw[stIndex + pwDeg + (pwDeg + 1) * pwIndex]);
        }
    }
    return b1 && b2;
}

void KKTSolver::compObj(const InputData& inputData, OutputData* outputData) {
    assert(outputData != NULL);
    outputData->_objVal = 0;
    if (inputData._deviationType == InputData::LP) {
        for (int i = 0; i < inputData._n; ++i) {
            outputData->_objVal += (1.0 / inputData._p) * inputData._cDev[i] *
                Pow(fabs(outputData->_x[i] - inputData._aDev[i]), inputData._p);
        }
    } else if (inputData._deviationType == InputData::PIECEWISE_LP) {
        assert(inputData._p == inputData._pwDeg);
        if (inputData._p == 2) {
            outputData->_objVal = pqTV(inputData._n, inputData._pw, inputData._bkpNums, inputData._cSep, outputData->_x);
            return;
        } else if (inputData._p == 1) {
            outputData->_objVal = plTV(inputData._n, inputData._pw, inputData._bkpNums, inputData._cSep, outputData->_x);
            return;
        }
    } else if (inputData._deviationType == InputData::HUBER_D) {
        for (int i = 0; i < inputData._n; ++i) {
            outputData->_objVal += inputData._cDev[i] *
                huberObj(outputData->_x[i] - inputData._aDev[i], inputData._huberD[i]);
        }
    }

    if (inputData._separationType == InputData::LQ) {
        for (int i = 0; i < inputData._n - 1; ++i) {
            outputData->_objVal += (1.0 / inputData._q) * inputData._cSep[i] *
                Pow(fabs(outputData->_x[i] - outputData->_x[i + 1]), inputData._q);
        }
    } else if (inputData._separationType == InputData::HUBER_S) {
        for (int i = 0; i < inputData._n - 1; ++i) {
            outputData->_objVal += inputData._cSep[i] *
                huberObj(outputData->_x[i] - outputData->_x[i + 1], inputData._huberS[i]);
        }
    }
}

void KKTSolver::compDrvt(const InputData& inputData, const OutputData& outputData,
                 int index, bool inclPrev, data_type* out_fDrvtValue) {
    assert(out_fDrvtValue != NULL);
    assert(index >= 0 && index < inputData._n);
    if (inclPrev) {
        if (index == 0) {
            inclPrev = false;
        }
    }
    assert(inputData._p >= 1 && inputData._q >= 1);
    int p = inputData._p;
    int q = inputData._q;
    data_type fDrvtValue = 0;
    if (inputData._deviationType == InputData::LP) {
        fDrvtValue = inputData._cDev[index] *
            Pow(outputData._x[index] - inputData._aDev[index], p - 1);
        if (p % 2 == 1 && outputData._x[index] - inputData._aDev[index] < 0) {
            fDrvtValue = -fDrvtValue;
        }
    } else if (inputData._deviationType == InputData::PIECEWISE_LP) {
        // Piecewise cases.
        assert(inputData._p == inputData._pwDeg);
        int stIndex = outputData._stIndex;
        int pwIndex = getPQIndex(inputData._pwDeg, inputData._pw, inputData._bkpNums[index], stIndex, outputData._x[index]);
        assert(pwValid(inputData._pwDeg, inputData._pw, inputData._bkpNums[index], stIndex, outputData._x[index], pwIndex));
        data_type a = inputData._pw[stIndex + (inputData._pwDeg + 1) * pwIndex];
        fDrvtValue = a * Pow(outputData._x[index], inputData._pwDeg - 1);
        for (int i = 1; i < inputData._pwDeg; ++i) {
            data_type b = inputData._pw[stIndex + (inputData._pwDeg + 1) * pwIndex + i];
            fDrvtValue -= b * Pow(outputData._x[index], inputData._pwDeg - i - 1);
        }
    } else if (inputData._deviationType == InputData::HUBER_D) {
        fDrvtValue = inputData._cDev[index] *
            huberDrvt(outputData._x[index] - inputData._aDev[index], inputData._huberD[index]);
    }

    if (inclPrev) {
        if (inputData._separationType == InputData::LQ) {
            if (q % 2 == 0) {
                fDrvtValue -= inputData._cSep[index - 1] *
                    Pow(outputData._x[index - 1] - outputData._x[index], q-1);
            } else {
                if (q == 1) {
                    // Take right sub-gradients
                    if (outputData._x[index - 1] - outputData._x[index] >= 0) {
                        fDrvtValue -= inputData._cSep[index - 1];
                    } else {
                        fDrvtValue += inputData._cSep[index - 1];
                    }
                } else {
                    if (outputData._x[index - 1] - outputData._x[index] >= 0) {
                        fDrvtValue -= inputData._cSep[index - 1] *
                            Pow(outputData._x[index - 1] - outputData._x[index],
                                q-1);
                    } else {
                        fDrvtValue += inputData._cSep[index - 1] *
                            Pow(outputData._x[index] - outputData._x[index - 1],
                                q-1);
                    }
                }
            }
        } else if (inputData._separationType == InputData::HUBER_S) {
            fDrvtValue += inputData._cSep[index - 1] *
                huberDrvt(outputData._x[index] - outputData._x[index - 1],
                          inputData._huberS[index - 1]);
        }
    }
    *out_fDrvtValue = fDrvtValue;
}

data_type KKTSolver::compSepInv(const InputData& inputData, const data_type& fDrvtValue,
                     int index) {
    data_type z = 0;
    if (inputData._separationType == InputData::LQ) {
        int q = inputData._q;
        if (fDrvtValue >= 0) {
            if (q > 1) {
                z = Pow(fDrvtValue / inputData._cSep[index], 1.0/(q - 1));
            } else {
                // TV-l1: Right end of the inverse.
                if (fDrvtValue < inputData._cSep[index]) {
                    z = 0;
                } else {
                    z = inputData._infinity;
                }
            }
        } else {
            // fDevValue < 0
            if (q > 1) {
                z = -Pow(-fDrvtValue / inputData._cSep[index], 1.0/(q - 1));
            } else {
                // TV-l1: Right end of the inverse.
                if (-fDrvtValue <= inputData._cSep[index]) {
                    z = 0;
                } else {
                    z = -inputData._infinity;
                }
            }
        }
    } else if (inputData._separationType == InputData::HUBER_S) {
        data_type delta = inputData._huberS[index];
        delta = delta * inputData._cSep[index];
        if (fDrvtValue > -delta && fDrvtValue < delta) {
            z = fDrvtValue / inputData._cSep[index];
        } else if (fDrvtValue >= delta) {
            z = inputData._infinity;
        } else if (fDrvtValue < -delta) {
            z = -inputData._infinity;
        } else {
            z = -inputData._huberS[index];
        }
    }
    return z;
}

int KKTSolver::propagate(const InputData& inputData, OutputData* outputData,
                    int index, data_type* out_fDrvtValue) {
    assert(outputData != NULL && index >= 0);
    assert(out_fDrvtValue != NULL);
    int n = inputData._n;
    compDrvt(inputData, *outputData, index, true, out_fDrvtValue);
    data_type& fDrvtValue = *out_fDrvtValue;

    // Start the propagation
    for (int i = index; i < n - 1; ++i) {
        data_type z = compSepInv(inputData, fDrvtValue, i);
        
        outputData->_x[i + 1] = outputData->_x[i] + z;
        // Check whether we need to continue
        if (outputData->_x[i + 1] < outputData->_bounds[i+1][0]) {
            for (int j = index; j < i + 1; ++j) {
                assert(outputData->_x[j] >= outputData->_bounds[j][0]);
                outputData->_bounds[j][0] = outputData->_x[j];
            }
            return -1;
        }

        if (outputData->_x[i + 1] > outputData->_bounds[i+1][1]) {
            for (int j = index; j < i + 1; ++j) {
                assert(outputData->_x[j] <= outputData->_bounds[j][1]);
                outputData->_bounds[j][1] = outputData->_x[j];
            }
            return 1;
        }

        data_type drvtDelta = 0;
        if (inputData._deviationType == InputData::PIECEWISE_LP) {
            // Move to the next piecewise deviation function.
            outputData->_stIndex += (inputData._pwDeg + 1) *
                                        inputData._bkpNums[i]
                                    + inputData._pwDeg;
        }
        compDrvt(inputData, *outputData, i + 1, false, &drvtDelta);
        fDrvtValue += drvtDelta;
    }

    if (fDrvtValue > 0) {
        // New upper divergence bound
        for (int i = index; i < n; ++i) {
            assert(outputData->_x[i] <= outputData->_bounds[i][1]);
            outputData->_bounds[i][1] = outputData->_x[i];
        }
    } else if (fDrvtValue < 0) {
        // New lower divergence bound
        for (int i = index; i < n; ++i) {
            assert(outputData->_x[i] >= outputData->_bounds[i][0]);
            outputData->_bounds[i][0] = outputData->_x[i];
        }
    }

    return 0;
}

void KKTSolver::solve(const InputData& inputData, OutputData* result) {
    assert(inputData._n >= 1 && inputData._p >= 1 && inputData._q >= 1);
    assert(result != NULL);

    for (int i = 0; i < inputData._n; ++i) {
        data_type l, u;
        l = result->_bounds[i][0];
        u = result->_bounds[i][1];

        if (u - l < inputData._solEsp) {
            result->_x[i] = (u + l) / 2;
            if (inputData._deviationType == InputData::PIECEWISE_LP) {
                // Move to the next piecewise deviation function.
                result->_stIndex += (inputData._pwDeg + 1) *
                                        inputData._bkpNums[i]
                                    + inputData._pwDeg;
            }
            continue;
        }
        result->_x[i] = (l + u) / 2;
        int stIndex = result->_stIndex;
        data_type fDrvtValue;
        int state = propagate(inputData, result, i, &fDrvtValue);
        while (u - l >= inputData._solEsp) {
            if (state < 0) {
                // Go up.
                l = result->_x[i];
            } else if (state > 0) {
                // Go down.
                u = result->_x[i];
            } else {
                if (fabs(fDrvtValue) < inputData._drvtEsp) {
                    return;
                } else if (fDrvtValue < 0) {
                    // Go up.
                    l = result->_x[i];
                } else {
                    // Go down.
                    u = result->_x[i];
                }
            }
            result->_x[i] = (l + u) / 2;
            result->_stIndex = stIndex;
            state = propagate(inputData, result, i, &fDrvtValue);
        }
        result->_stIndex = stIndex;
        if (inputData._deviationType == InputData::PIECEWISE_LP) {
            // Move to the next piecewise deviation function.
            result->_stIndex += (inputData._pwDeg + 1) * inputData._bkpNums[i]
                                + inputData._pwDeg;
        }
    }
}

static inline data_type l1Slope(data_type value, data_type anchor, data_type slope) {
    assert(slope >= 0);
    // Right sub-derivative.
    if (value > anchor) {
        return slope;
    } else {
        return -slope;
    }
}

static inline int getStIndex(const std::vector<int>& boundIndex) {
    assert(boundIndex.size() == 2);
    if (boundIndex[0] <= boundIndex[1]) {
        return 0;
    } else {
        return 1;
    }
}

void KKTSolver::fast_l2_l1(const InputData &inputData, OutputData *result) {
    assert(inputData._n >= 2 && inputData._p == 2 && inputData._q == 1);
    assert(result != NULL);

    int n = inputData._n;
    int i = 0;
    while (i < n) {
        std::vector<int> boundIndex(2, i);
        std::vector<data_type> accuDrvtCoeff(2, inputData._cDev[i]);
        std::vector<data_type> accuDrvtConst(2, -inputData._cDev[i] * inputData._aDev[i]);
        data_type l = result->_bounds[i][0];
        data_type u = result->_bounds[i][1];
        result->_x[i] = (l + u) / 2;
        while (u - l >= inputData._solEsp) {
            int binIndex = getStIndex(boundIndex);
            int stIndex = boundIndex[binIndex];
            data_type drvtCoeff = accuDrvtCoeff[binIndex];
            data_type drvtConst = accuDrvtConst[binIndex];
            data_type l1Const = 0;
            if (i > 0) {
                // Include previous slope
                l1Const = l1Slope(result->_x[i], result->_x[i - 1], inputData._cSep[i - 1]);
            }
            data_type drvtValue = drvtCoeff * result->_x[i] + drvtConst + l1Const;
            // +1: Go down; -1: Go up.
            int direction = drvtValue >= 0 ? 1 : -1;
            // Propagation
            while (stIndex < n - 1) {
                // Propagate
                if (drvtValue >= 0) {
                    if (drvtValue < inputData._cSep[stIndex]) {
                        // Propagate success
                        drvtCoeff += inputData._cDev[stIndex + 1];
                        drvtConst += -inputData._cDev[stIndex + 1] * inputData._aDev[stIndex + 1];
                        stIndex++;
                        drvtValue = drvtCoeff * result->_x[i] + drvtConst + l1Const;
                    } else {
                        // Propagate failed, upper bound exceeded
                        // Update the data structures
                        boundIndex[1] = stIndex;
                        accuDrvtCoeff[1] = drvtCoeff;
                        accuDrvtConst[1] = drvtConst;
                        direction = 1;
                        break;
                    }
                } else {
                    if (-drvtValue <= inputData._cSep[stIndex]) {
                        // Propagate success
                        drvtCoeff += inputData._cDev[stIndex + 1];
                        drvtConst += -inputData._cDev[stIndex + 1] * inputData._aDev[stIndex + 1];
                        stIndex++;
                        drvtValue = drvtCoeff * result->_x[i] + drvtConst + l1Const;
                    } else {
                        // Propagate failed, lower bound exceeded
                        // Update the data structures
                        boundIndex[0] = stIndex;
                        accuDrvtCoeff[0] = drvtCoeff;
                        accuDrvtConst[0] = drvtConst;
                        direction = -1;
                        break;
                    }
                }
            }

            if (stIndex == n - 1) {
                // Reach the end, no early stopping.
                direction = drvtValue >= 0 ? 1 : -1;
                // Update the data structures
                if (direction == 1) {
                    boundIndex[1] = stIndex;
                    accuDrvtCoeff[1] = drvtCoeff;
                    accuDrvtConst[1] = drvtConst;
                } else {
                    boundIndex[0] = stIndex;
                    accuDrvtCoeff[0] = drvtCoeff;
                    accuDrvtConst[0] = drvtConst;
                }
            }

            // Find the next search value.
            if (direction == -1) {
                l = result->_x[i];
            } else {
                u = result->_x[i];
            }
            result->_x[i] = (l + u) / 2;
        }
        int binIndex = getStIndex(boundIndex);
        int stIndex = boundIndex[binIndex];
        for (int j = i + 1; j <= stIndex; ++j) {
            result->_x[j] = result->_x[i];
        }
        i = stIndex + 1;
    }
}

void KKTSolver::fast_linear_l2(const InputData &inputData, OutputData *result) {
    assert(inputData._n >= 2 && inputData._p == 1 && inputData._q == 2);
    assert(result != NULL);
    data_type linCoeffSum = 0;
    int n = inputData._n;
    for (int i = 0; i < n; ++i) {
        linCoeffSum += inputData._cDev[i];
    }
    assert(fabs(linCoeffSum) < 1e-6);

    result->_x[0] = 0;
    linCoeffSum = 0;
    for (int i = 0; i < n - 1; ++i) {
        linCoeffSum += inputData._cDev[i];
        result->_x[i + 1] = result->_x[i] + linCoeffSum;
    }
}
