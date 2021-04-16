//
//  kkt.hpp
//  KKT
//
//  Created by Cheng Lu on 4/2/20.
//  Copyright © 2020 Cheng Lu. All rights reserved.
//

// Citation: Cheng Lu and Dorit S. Hochbaum. A unified approach for a 1D generalized total variation problem. Mathematical Programming, February 2021.
// https://doi.org/10.1007/s10107-021-01633-2

//  Header file for KKT method.


#ifndef kkt_h
#define kkt_h

#include <cassert>
#include <cstdlib>
#include <map>
#include <vector>

typedef signed long long time_ms_type;  // at least 64 bits.

typedef double data_type;
typedef std::map<int, data_type> l1_bound_type;

const data_type KKT_SOL_ESP = 1e-6;
const data_type KKT_DRVT_ESP = 1e-6;
const data_type KKT_INFINITY = 1e10;
const data_type KKT_LB = -1e4;  // Uniform solution bounds for all problems.
const data_type KKT_UB = 1e4;

// Input data for the generalized total variation model.
struct InputData {
    // Input parameters
    int _n;

    typedef enum DEVIATION_TYPE {
        LP = 0,
        PIECEWISE_LP = 1,  // Piecewise deviation functions.
        HUBER_D = 2,
        // Extensible to high order piecewise functions.
    } deviation_type;
    deviation_type _deviationType;

    typedef enum SEPARATIOIN_TYPE {
        LQ = 0,
        HUBER_S = 1,
    } separation_type;
    separation_type _separationType;

    // lp + lq.
    int _p, _q;
    data_type* _cDev = NULL;  // Appending a 1/_p in the front.
    data_type* _aDev = NULL;
    data_type* _cSep = NULL;  // Appending a 1/_q in the front.

    // Piecewise deviation functions + l1.
    // Extensible to high order piecewise functions.
    int _pwDeg = 0;  // Degree of composing pieces.
    data_type* _pw = NULL;
    int* _bkpNums = NULL;

    // Huber loss separation functions
    data_type* _huberD = NULL;
    data_type* _huberS = NULL;

    // Algorithm parameters
    data_type _lb, _ub;  // Solution lower and upper bounds.
    data_type _solEsp;  // Solution accuracy.
    data_type _drvtEsp;  // Derivative accuracy.
    data_type _infinity;  // For l1-TV.

    // By default, p = 2, q = 1.
    InputData(int n): _n(n) {
        assert(n >= 1);
        _deviationType = LP;
        _separationType = LQ;
        _p = 2;
        _q = 1;
        initParams();
    }

    InputData(int n, int p, int q): _n(n), _p(p), _q(q) {
        assert(n >= 1 && p >= 1 && q >= 1);
        _deviationType = LP;
        _separationType = LQ;
        initParams();
    }

    // Constructor for piecewise quadratic deviations.
    InputData(int n, int pwDeg, const std::vector<int>& bkpNums,
              const std::vector<data_type>& pw) {
        assert(n >= 1);
        assert(bkpNums.size() == n);
        assert(pwDeg == 1 || pwDeg == 2);  // Currently only support piecewise l1 and l2.
        _n = n;
        _pwDeg = pwDeg;
        _deviationType = PIECEWISE_LP;
        _separationType = LQ;
        _p = pwDeg;
        _q = 1;  // Currently only support L1-TV.
        _bkpNums = (int*)malloc(_n * sizeof(int));
        data_type bkp_lb = KKT_INFINITY;
        data_type bkp_ub = -KKT_INFINITY;
        int totalBkps = 0;
        for (int i = 0; i < _n; ++i) {
            assert(bkpNums[i] >= 0);
            _bkpNums[i] = bkpNums[i];
            totalBkps += bkpNums[i];
        }
        assert((_pwDeg + 1) * totalBkps + _pwDeg * _n == pw.size());
        _pw = (data_type*)malloc(pw.size() * sizeof(data_type));
        int pwIndex = 0;
        for (int i = 0; i < _n; ++i) {
            for (int j = 0; j < _pwDeg; ++j) {
                _pw[pwIndex] = pw[pwIndex];
                pwIndex++;
            }
            for (int j = 0; j < bkpNums[i]; ++j) {
                for (int k = 0; k < _pwDeg + 1; ++k) {
                    _pw[pwIndex] = pw[pwIndex];
                    if (k == 0) {
                        // Update lower and upper bound.
                        data_type bkp = _pw[pwIndex];
                        if (bkp < bkp_lb) bkp_lb = bkp;
                        if (bkp > bkp_ub) bkp_ub = bkp;
                    }
                    pwIndex++;
                }
            }
        }
        // Left separation parameter
        _cSep = (data_type*)calloc((_n - 1), sizeof(data_type));
        _lb = KKT_LB; // Use uniform lower and upper bounds for all problems.
        _ub = KKT_UB;
        _solEsp = KKT_SOL_ESP;
        _drvtEsp = KKT_DRVT_ESP;
        _infinity = KKT_INFINITY;
    }

    InputData(int n, int p, int q, deviation_type deviationType,
              separation_type separationType) {
        assert(n >= 1 && p >= 1 && q >= 1);
        _n = n;
        _p = p;
        _q = q;
        _deviationType = deviationType;
        _separationType = separationType;
        initParams();
        if (_deviationType == InputData::HUBER_D) {
            _p = 2;
            _huberD = (data_type*)calloc(_n, sizeof(data_type));
        }
        if (_separationType == InputData::HUBER_S) {
            _q = 2;
            _huberS = (data_type*)calloc((_n - 1), sizeof(data_type));
        }
    }

    ~InputData() {
        if (_cDev != NULL)
            free(_cDev);
        if (_aDev != NULL)
            free(_aDev);
        if (_cSep != NULL)
            free(_cSep);
        if (_pw != NULL)
            free(_pw);
        if (_bkpNums != NULL)
            free(_bkpNums);
        if (_huberD != NULL) {
            free(_huberD);
        }
        if (_huberS != NULL) {
            free(_huberS);
        }
    }

    void initParams() {
        _cDev = (data_type*)calloc(_n, sizeof(data_type));
        _aDev = (data_type*)calloc(_n, sizeof(data_type));
        _cSep = (data_type*)calloc((_n - 1), sizeof(data_type));
        _lb = KKT_LB;
        _ub = KKT_UB;
        _solEsp = KKT_SOL_ESP;
        _drvtEsp = KKT_DRVT_ESP;
        _infinity = KKT_INFINITY;
    }
};

// Output data for the generalized total variation model.
struct OutputData {
    int _n;
    data_type* _x = NULL;
    data_type _objVal;
    // Save divergent values.
    std::vector<std::vector<data_type>> _bounds;
    // For piecewise deviation functions.
    int _stIndex;

    OutputData() {}

    OutputData(const InputData& inputData) {
        _n = inputData._n;
        _x = (data_type*)calloc(_n, sizeof(data_type));
        _objVal = 0;
        // Rolling buffers.
        std::vector<data_type> item(2, 0);
        item[0] = inputData._lb;
        item[1] = inputData._ub;
        for (int i = 0; i < _n; ++i) {
            _bounds.push_back(item);
        }
        _stIndex = 0;
    }

    ~OutputData() {
        if (_x != NULL)
            free(_x);
    }

    void operator = (const OutputData& other) {
        if (_x != NULL) {
            free(_x);
        }

        _n = other._n;
        _x = (data_type*)calloc(_n, sizeof(data_type));
        for (int i = 0; i < _n; ++i) {
            _x[i] = other._x[i];
        }
        _objVal = 0;
        _bounds = other._bounds;
        _stIndex = other._stIndex;
    }
};

// KKT Solver
class KKTSolver {
public:
    // Main generic compute function.
    void solve(const InputData& inputData, OutputData* result);

    // Fast l2_l1 solver, working for both unweighted and weighted,
    // by adapting the general versions of the KKT algorithms.
    void fast_l2_l1(const InputData& inputData, OutputData* result);

    // Fast linear_l2 solver (1D graph Laplacian solver)
    // Problem:
    // min_{x_i} \sum_{i=1}^n c_ix_i + 0.5 * \sum_{i=1}^{n-1}(x_i - x_{i+1})^2.
    void fast_linear_l2(const InputData& inputData, OutputData* result);

    // Compute the objective value.
    virtual void compObj(const InputData& inputData, OutputData* outputData);

private:
    // Overridable for your specific fidelity/regularization functions.

    // Propagation function
    // Return type: +1: Go down; 0, depending on fDevValue; -1: Go up.
    int propagate(const InputData& inputData, OutputData* outputData,
                          int index, data_type* out_fDrvtValue);

    // Unified function to compute the derivative.
    // inclPrev: true: include the previous term for re-starting on sub-problems.
    // (Non-null) out_fDrvtValue: If to compute two derivatives, it represents the right one.
    void compDrvt(const InputData& inputData,
                          const OutputData& outputData, int index,
                          bool inclPrev, data_type* out_fDrvtValue);

    // Compute the inverse of the separation term.
    data_type compSepInv(const InputData& inputData,
                                 const data_type& fDrvtValue, int index);
};

#endif /* kkt_h */
