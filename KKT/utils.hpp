//
//  utils.hpp
//  KKT
//
//  Created by ChengLu on 4/15/21.
//  Copyright © 2021 Cheng Lu. All rights reserved.
//

#ifndef utils_h
#define utils_h

#include "kkt.hpp"

// Auxiliary functions
inline void outputSolution(const OutputData& result) {
    std::cout << "x = \n";
    for (int i = 0; i < result._n - 1; ++i) {
        std::cout << result._x[i] << ",";
    }
    std::cout << result._x[result._n - 1] << "\n\n";
}

inline void compareSolutions(const OutputData& result1, const OutputData& result2,
                      data_type* out_maxDiff = NULL,
                      data_type* out_minDiff = NULL,
                      data_type* out_meanDiff = NULL) {
    assert(result1._n == result2._n && result1._n > 0);
    data_type maxDiff, minDiff, meanDiff = 0;
    maxDiff = fabs(result1._x[0] - result2._x[0]);
    minDiff = maxDiff;
    meanDiff = maxDiff;
    int n = result1._n;
    for (int i = 1; i < n; ++i) {
        data_type diff = fabs(result1._x[i] - result2._x[i]);
        maxDiff = diff > maxDiff ? diff : maxDiff;
        minDiff = diff < minDiff ? diff : minDiff;
        meanDiff += diff;
    }
    meanDiff = meanDiff * 1.0 / n;
    if (out_maxDiff != NULL) {
        *out_maxDiff = maxDiff;
    }
    if (out_minDiff != NULL) {
        *out_minDiff = minDiff;
    }
    if (out_meanDiff != NULL) {
        *out_meanDiff = meanDiff;
    }
}

// Count the pairs of neighboring solution element that are different.
// Especially for fuse lasso verification.
inline int numSolChg(const OutputData& outputData) {
    assert(outputData._n > 1);
    int num = 0;
    for (int i = 0; i < outputData._n - 1; ++i) {
        if (fabs(outputData._x[i] - outputData._x[i + 1]) >= 1e-2) {
            ++num;
        }
    }
    return num;
}

// Huber functions
inline data_type huberObj(data_type x, data_type delta) {
    assert(delta > 0);
    if (fabs(x) <= delta) {
        return 0.5 * x * x;
    } else {
        return delta * (fabs(x) - 0.5 * delta);
    }
}

inline data_type huberDrvt(data_type x, data_type delta) {
    assert(delta > 0);
    if (fabs(x) <= delta) {
        return x;
    } else {
        if (x > 0) {
            return delta;
        } else {
            return -delta;
        }
    }
}

// Evaluation functions for piecewise linear and quadratic objective functions.

inline data_type plFunc(data_type* slopeAndBkps, int breakpointNum, data_type x) {
    if (breakpointNum == 0) return slopeAndBkps[0] * x;

    if (x <= slopeAndBkps[1]) return (x - slopeAndBkps[1]) * slopeAndBkps[0];

    data_type y = 0;

    for (slopeAndBkps += 2; breakpointNum > 1; slopeAndBkps += 2, breakpointNum--) {
        if (x <= slopeAndBkps[1]) break;
        y += (slopeAndBkps[1] - slopeAndBkps[-1]) * slopeAndBkps[0];
    }

    return y + (x - slopeAndBkps[-1]) * slopeAndBkps[0];
}

inline data_type plTV(int n, data_type* slopeAndBkps, int* breakpointNums, data_type* cSeps, double* x) {
    data_type cost = 0;

    // Deviation part.
    for (int i = 0; i < n; i++) {
        cost += plFunc(slopeAndBkps, breakpointNums[i], x[i]);
        slopeAndBkps += 2 * breakpointNums[i] + 1;
    }

    // Separation part.
    for (int i = 0; i < n-1; i++) {
        data_type diff = cSeps[i] * fabs(x[i+1] - x[i]);
        cost += diff;
    }

    return cost;
}

inline data_type quadraticFunc(data_type a, data_type b, data_type x0, data_type x) {
    return 0.5 * a * (x-x0) * (x+x0) - b * (x-x0);
}

inline data_type pqFunc(data_type* qpAndBkps, int breakpointNum, data_type x) {
    if (breakpointNum == 0) return 0.5 * qpAndBkps[0] * x * x - qpAndBkps[1] * x;

    if (x <= qpAndBkps[2]) return quadraticFunc(qpAndBkps[0], qpAndBkps[1], qpAndBkps[2], x);

    data_type y = 0;

    for (qpAndBkps += 3; breakpointNum > 1; qpAndBkps += 3, breakpointNum--) {
        if (x <= qpAndBkps[2]) break;
        y += quadraticFunc(qpAndBkps[0], qpAndBkps[1], qpAndBkps[-1], qpAndBkps[2]);
    }

    return y + quadraticFunc(qpAndBkps[0], qpAndBkps[1], qpAndBkps[-1], x);
}

inline data_type pqTV(int n, data_type* qpAndBkps, int* breakpointNums, data_type* cSeps, data_type* x) {
    data_type cost = 0;

    // Deviation part.
    for (int i = 0; i < n; i++) {
        cost += pqFunc(qpAndBkps, breakpointNums[i], x[i]);
        qpAndBkps += 3 * breakpointNums[i] + 2;
    }

    // Separation part.
    for (int i = 0; i < n-1; i++) {
        data_type diff = cSeps[i] * fabs(x[i+1] - x[i]);
        cost += diff;
    }

    return cost;
}

#endif /* utils_h */
