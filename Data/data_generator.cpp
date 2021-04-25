//
//  data_generator.cpp
//  KKT
//
//  Created by Cheng Lu on 4/11/20.
//  Copyright © 2020 Cheng Lu. All rights reserved.
//

// Implementation of functions in data_generator.hpp

#include <cmath>
#include "data_generator.hpp"
#include <iostream>
#include <random>

std::string toString(gen_data_type genDataType) {
    switch (genDataType) {
        case KKT_PWL1: return "KKT-PWL1";
        case KKT_PWL2: return "KKT-PWL2";
        case KKT_LP_LQ: return "KKT-LP-LQ";
        case MPO_NW_INPUT_SIZE: return "MPO-NW-INPUT-SIZE";
        case MPO_NW_LAMBDA: return "MPO-NW-LAMBDA";
        case MPO_W_INPUT_SIZE: return "MPO-W-INPUT-SIZE";
        case MPO_W_LAMBDA: return "MPO-W-LAMBDA";
        case CONDAT_WORST_CASE: return "CONDAT-WORST-CASE";
        case KKT_LINEAR_L2: return "KKT-Linear-L2";
        case KKT_HUBER: return "KKT-HUBER";
        default:
            return "";
    }
}

// Global random seed.
std::random_device rd;
std::mt19937 gen(rd());

std::vector<int> genPWBkpNums(int n) {
    assert(n >= 1);
    std::uniform_int_distribution<int>
        bkpNum_distribution(PW_BKPNUM_UNIF_LEFT, PW_BKPNUM_UNIF_RIGHT);
    std::vector<int> bkpNums(n, 0);
    for (int i = 0; i < n; ++i) {
        bkpNums[i] = bkpNum_distribution(gen);
    }
    return bkpNums;
}

std::vector<int> genPWBkpNums(int n, int bkpNum) {
    assert(n >= 1);
    assert(bkpNum >= 0);
    std::vector<int> bkpNums(n, 0);
    for (int i = 0; i < n; ++i) {
        bkpNums[i] = bkpNum;
    }
    return bkpNums;
}

std::vector<int> genPWBkpNums(int n, int bkpNumLB, int bkpNumUB) {
    assert(n >= 1);
    assert(bkpNumLB >= 0 && bkpNumUB >= 0 && bkpNumLB <= bkpNumUB);
    std::uniform_int_distribution<int>
        bkpNum_distribution(bkpNumLB, bkpNumUB);
    std::vector<int> bkpNums(n, 0);
    for (int i = 0; i < n; ++i) {
        bkpNums[i] = bkpNum_distribution(gen);
    }
    return bkpNums;
}

std::vector<data_type> genPWFuncs(int n, int pwDeg,
                                  const std::vector<int>& bkpNums) {
    assert(n >= 1 && n == bkpNums.size());
    assert(pwDeg == 1 || pwDeg == 2);  // Currently only support piecewise l1 and l2.
    int totalBkps = 0;
    for (int i = 0; i < n; ++i) {
        assert(bkpNums[i] >= 0);
        totalBkps += bkpNums[i];
    }
    std::vector<data_type> pw((pwDeg + 1) * totalBkps + pwDeg * n, 0);
    // Breakpoints
    std::uniform_real_distribution<data_type>
        bkp_distribution(PW_BKP_UNIF_LEFT, PW_BKP_UNIF_RIGHT);
    // Quadratic coefficients.
    // If pwDeg == 2, it must be non-negative.
    std::uniform_real_distribution<data_type>
        deg2_distribution(PW_DEG2_UNIF_LEFT, PW_DEG2_UNIF_RIGHT);
    // Linear coefficients
    // It is unrestricted in sign, even for pwDeg == 1.
    std::uniform_real_distribution<data_type>
        deg1_distribution(PW_DEG1_UNIF_LEFT, PW_DEG1_UNIF_RIGHT);
    // Incremental step, for all increments.
    // TODO: May want to have different inc random variable for different parameters.
    std::uniform_real_distribution<data_type>
        inc_distribution(PW_INC_UNIF_LEFT, PW_INC_UNIF_RIGHT);
    int pwIndex = 0;
    if (pwDeg == 2) {
        for (int i = 0; i < n; ++i) {
            // Initiate quadratic coefficient.
            data_type a = deg2_distribution(gen);
            pw[pwIndex] = a;
            // Initiate linear coefficient.
            data_type b = deg1_distribution(gen);
            pw[pwIndex + 1] = b;
            data_type lambda = bkp_distribution(gen);
            // Generate breakpoint list
            for (int j = 0; j < bkpNums[i]; ++j) {
                pw[pwIndex + 2 + 3 * j] = lambda;
                // Left quadratic gradient.
                data_type lGradient = a * lambda - b;
                // Inc for gradient.
                data_type inc = inc_distribution(gen);
                data_type rGradient = lGradient + inc;
                // Inc for quadratic coefficient.
                inc = inc_distribution(gen);
                a += inc;
                b = a * lambda - rGradient;
                pw[pwIndex + 2 + 3 * j + 1] = a;
                pw[pwIndex + 2 + 3 * j + 2] = b;
                // Inc for breakpoint.
                inc = inc_distribution(gen);
                lambda += inc;
            }
            pwIndex += 3 * bkpNums[i] + 2;
        }
    } else {
        // pwDeg == 1
        for (int i = 0; i < n; ++i) {
            // Initiate linear coefficient.
            data_type b = deg1_distribution(gen);
            pw[pwIndex] = b;
            data_type lambda = bkp_distribution(gen);
            // Generate breakpoint list
            for (int j = 0; j < bkpNums[i]; ++j) {
                pw[pwIndex + 1 + 2 * j] = lambda;
                // Inc for gradient
                data_type inc = inc_distribution(gen);
                b += inc;
                // Force the last half pieces to have positive gradients
                // in order to make the problem have bounded minimum.
                if (j >= bkpNums[i] / 2 && b <= 0) {
                    b = 1;
                }
                pw[pwIndex + 1 + 2 * j + 1] = b;
                // Inc for breakpoint.
                inc = inc_distribution(gen);
                lambda += inc;
            }
            pwIndex += 2 * bkpNums[i] + 1;
        }
    }
    return pw;
}

void fillSep(int n, InputData* inputData, data_type in_lambda,
             bool withSample) {
    assert(n >= 1);
    assert(inputData != NULL && inputData->_n == n);
    if (in_lambda >= 0) {
        if (!withSample) {
            // Uniform
            for (int i = 0; i < n - 1; ++i) {
                inputData->_cSep[i] = in_lambda;
            }
        } else {
            std::uniform_real_distribution<data_type>
                csep_distribution(0.5 * in_lambda, 1.5 * in_lambda);
            for (int i = 0; i < n - 1; ++i) {
                inputData->_cSep[i] = csep_distribution(gen);
            }
        }
    } else {
        // Weighted
        std::uniform_real_distribution<data_type>
            csep_distribution(PW_CSEP_UNIF_LEFT, PW_CSEP_UNIF_RIGHT);
        for (int i = 0; i < n - 1; ++i) {
            inputData->_cSep[i] = csep_distribution(gen);
        }
    }
}

// Shared implementation for cDev's and aDev's.
void genLpLqFuncInternal(int n, InputData* inputData, bool cDev1) {
    assert(n >= 1);
    assert(inputData != NULL && inputData->_n == n);
    std::uniform_real_distribution<data_type>
        adev_distribution(LPLQ_ADEV_UNIF_LEFT, LPLQ_ADEV_UNIF_RIGHT);
    for (int i = 0; i < n; ++i) {
        inputData->_aDev[i] = adev_distribution(gen);
    }

    if (cDev1) {
        for (int i = 0; i < n; ++i) {
            inputData->_cDev[i] = 1;
        }
    } else {
        std::uniform_real_distribution<data_type>
            cdev_distribution(LPLQ_CDEV_UNIF_LEFT, LPLQ_CDEV_UNIF_RIGHT);
        for (int i = 0; i < n; ++i) {
            inputData->_cDev[i] = cdev_distribution(gen);
        }
    }
}

void genLpLqFuncs(int n, InputData* inputData, bool cDev1) {
    assert(n >= 1);
    assert(inputData != NULL && inputData->_n == n);
    genLpLqFuncInternal(n, inputData, cDev1);
    // Generate cSep[i]s
    std::uniform_real_distribution<data_type>
        csep_distribution(LPLQ_CSEP_UNIF_LEFT, LPLQ_CSEP_UNIF_RIGHT);
    for (int i = 0; i < n - 1; ++i) {
        inputData->_cSep[i] = csep_distribution(gen);
    }
}

void genLpLqFuncs(int n, data_type lambda, InputData* inputData, bool resample, bool cDev1) {
    assert(n >= 1);
    assert(inputData != NULL && inputData->_n == n);
    assert(lambda >= 0);
    genLpLqFuncInternal(n, inputData, cDev1);
    if (!resample) {
        for (int i = 0; i < n - 1; ++i) {
            inputData->_cSep[i] = lambda;
        }
    } else {
        std::uniform_real_distribution<data_type>
            csep_distribution(0.5 * lambda, 1.5 * lambda);
        for (int i = 0; i < n - 1; ++i) {
            inputData->_cSep[i] = csep_distribution(gen);
        }
    }
}

void genLpLqFuncs(int n, data_type lambdaLB, data_type lambdaUB,
                  InputData* inputData, bool cDev1) {
    assert(n >= 1);
    assert(lambdaLB >= 0 && lambdaUB >= 0 && lambdaLB < lambdaUB);
    assert(inputData != NULL && inputData->_n == n);
    genLpLqFuncInternal(n, inputData, cDev1);
    // Generate cSep[i]s.
    std::uniform_real_distribution<data_type>
        csep_distribution(lambdaLB, lambdaUB);
    data_type lambda = csep_distribution(gen);
    for (int i = 0; i < n - 1; ++i) {
        inputData->_cSep[i] = lambda;
    }
}

void genLpLqFuncs(int n, InputData* inputData, gen_data_type genDataType,
                 data_type in_lambda, data_type* out_lambda) {
    assert(n >= 1 && inputData != NULL && inputData->_n == n);
    assert(in_lambda >= 0);
    switch (genDataType) {
        case MPO_NW_INPUT_SIZE: {
            std::uniform_real_distribution<data_type>
                lambda_distribution(0, 50);
            data_type lambda = lambda_distribution(gen);
            std::uniform_real_distribution<data_type>
                adev_distribution(-2*lambda, 2*lambda);
            for (int i = 0; i < n; ++i) {
                inputData->_cDev[i] = 1;
                inputData->_aDev[i] = adev_distribution(gen);
            }
            for (int i = 0; i < n - 1; ++i) {
                inputData->_cSep[i] = lambda;
            }
            break;
        }
        case MPO_NW_LAMBDA: {
            assert(in_lambda > 0);  // Otherwise the problem is boring.
            std::uniform_real_distribution<data_type>
                adev_distribution(-2, 2);
            for (int i = 0; i < n; ++i) {
                inputData->_cDev[i] = 1;
                inputData->_aDev[i] = adev_distribution(gen);
            }
            for (int i = 0; i < n - 1; ++i) {
                inputData->_cSep[i] = in_lambda;
            }
            break;
        }
        case MPO_W_INPUT_SIZE: {
            std::uniform_real_distribution<data_type> csep_distribution(0, 100);
            data_type cSepSum = 0;
            for (int i = 0; i < n - 1; ++i) {
                inputData->_cSep[i] = csep_distribution(gen);
                cSepSum += inputData->_cSep[i];
            }
            data_type lambda = cSepSum / (n - 1);
            if (out_lambda != NULL) {
                *out_lambda = lambda;
            }
            std::uniform_real_distribution<data_type>
                adev_distribution(-2*lambda, 2*lambda);
            for (int i = 0; i < n; ++i) {
                inputData->_cDev[i] = 1;
                inputData->_aDev[i] = adev_distribution(gen);
            }
            break;
        }
        case MPO_W_LAMBDA: {
            assert(in_lambda > 0);  // Otherwise the problem is boring.
            std::uniform_real_distribution<data_type>
                adev_distribution(-2, 2);
            for (int i = 0; i < n; ++i) {
                inputData->_cDev[i] = 1;
                inputData->_aDev[i] = adev_distribution(gen);
            }
            std::uniform_real_distribution<data_type>
                csep_distribution(0.5 * in_lambda, 1.5 * in_lambda);
            for (int i = 0; i < n - 1; ++i) {
                inputData->_cSep[i] = csep_distribution(gen);
            }
            break;
        }
        case CONDAT_WORST_CASE: {
            assert(n > 3);
            data_type alpha = 4.0 / ((n - 2) * (n - 3));
            for (int i = 0; i < n; ++i) {
                inputData->_cDev[i] = 1;
            }
            for (int i = 0; i < n - 1; ++i) {
                inputData->_cSep[i] = 1;
            }
            inputData->_aDev[0] = -2;
            for (int i = 2; i <= n - 1; ++i) {
                inputData->_aDev[i - 1] = alpha * (i - 2);
            }
            inputData->_aDev[n - 1] = alpha * (n - 3) + 2;
            break;
        }
        default:
            std::cout << "KKT_PWL1, KKT_PWL2: Call genPWBkpNums, genPWFuncs and fillSep\n"
                << "KKT_LP_LQ: Call other versions of genLpLqFuncs.\n";
    }
}

void genLinearL2Funcs(int n, InputData* inputData) {
    assert(n >= 1);
    assert(inputData != NULL && inputData->_n == n);

    std::uniform_real_distribution<data_type> cdev_distribution(-100, 100);
    data_type sum = 0;
    for (int i = 0; i < n - 1; ++i) {
        inputData->_cDev[i] = cdev_distribution(gen);
        sum += inputData->_cDev[i];
    }
    inputData->_cDev[n - 1] = -sum;
    for (int i = 0; i < n - 1; ++i) {
        inputData->_cSep[i] = 0.5;
    }
}

void genHuberFuncs(int n, const std::vector<data_type>& baselines,
                   InputData* inputData, bool isDev, double lRatio,
                   double rRatio) {
    assert(n >= 1 && baselines.size() == n);
    assert(lRatio > 0 && rRatio > 0 && lRatio < rRatio);
    assert((isDev && inputData->_deviationType == InputData::HUBER_D && n == inputData->_n)
           || (!isDev && inputData->_separationType == InputData::HUBER_S && n == inputData->_n - 1));

    std::uniform_real_distribution<data_type> ratio_distribution(lRatio, rRatio);
    if (isDev) {
        for (int i = 0; i < n; ++i) {
            inputData->_huberD[i] = ratio_distribution(gen) * fabs(baselines[i]);
        }
    } else {
        for (int i = 0; i < n; ++i) {
            inputData->_huberS[i] = ratio_distribution(gen) * fabs(baselines[i]);
        }
    }
}
