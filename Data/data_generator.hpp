//
//  data_generator.hpp
//  KKT
//
//  Created by Cheng Lu on 4/11/20.
//  Copyright © 2020 Cheng Lu. All rights reserved.
//

// Generate various types of synthetic data for algorithm comparison.


#ifndef data_generator_hpp
#define data_generator_hpp

#include "kkt.hpp"
#include <unordered_map>
#include <string>
#include <vector>

// Data types
typedef enum GEN_DATA_TYPE {
    KKT_PWL1 = 0,
    KKT_PWL2,
    KKT_LP_LQ,
    MPO_NW_INPUT_SIZE,
    MPO_NW_LAMBDA,
    MPO_W_INPUT_SIZE,
    MPO_W_LAMBDA,
    CONDAT_WORST_CASE,
    KKT_HUBER,
} gen_data_type;

// Map from data type to string for output.
std::string toString(gen_data_type genDataType);


//////////////////////////////////////////////////

// Generate parameters for ${n} piecewise functions.
// Each piece is a ${pqDeg} degree polynomial.
// Each piecewise function is convex, so one needs to guarantee that
//  -- Each piece is convex
//  -- (Sub-)Gradients are non-decreasing at breakpoints.
//
// Each p-degree piece is of the form (no constant terms):
//
// 1/p * a_p * x^p - \sum_{d=1}^{p-1} 1/(p-d) * a_{p-d} * x^{p-d}.
//
//
// Currently only support p = 1 and p = 2.
//
// TODO: Extend to arbitrary p >= 2.

// Tuning hyperparameters
const data_type PW_BKP_UNIF_LEFT = -15.0;
const data_type PW_BKP_UNIF_RIGHT = -1.0;
const data_type PW_DEG2_UNIF_LEFT = 0.5;
const data_type PW_DEG2_UNIF_RIGHT = 5.0;
const data_type PW_DEG1_UNIF_LEFT = -50.0;
const data_type PW_DEG1_UNIF_RIGHT = -25.0;
const data_type PW_INC_UNIF_LEFT = 0.1;
const data_type PW_INC_UNIF_RIGHT = 2.0;

const int PW_BKPNUM_UNIF_LEFT = 100;
const int PW_BKPNUM_UNIF_RIGHT = 200;

const data_type PW_CSEP_UNIF_LEFT = 0.0;
const data_type PW_CSEP_UNIF_RIGHT = 1.0;

// Generate the list of number of breakpoints for ${n} functions.
std::vector<int> genPWBkpNums(int n);
// Overload version for identical number of breakpoints.
std::vector<int> genPWBkpNums(int n, int bkpNum);
// Overload version for varied number of breakpoints.
std::vector<int> genPWBkpNums(int n, int bkpNumLB, int bkpNumUB);
// Based on the number of breakpoints for each function,
// generate the list of breakpoints and each piece's coefficients.
std::vector<data_type> genPWFuncs(int n, int pwDeg,
                                  const std::vector<int>& bkpNums);
// Fills in separation term coefficients
// If in_lambda >= 0, uniform version; o/w, weighted version.
// If withSample set, still weighted version, but the weights are sampled around in_lambda.
void fillSep(int n, InputData* inputData, data_type in_lambda = -1,
             bool withSample = false);

//////////////////////////////////////////////////////////

// Generate parameter for lp-lq problem of p, q >= 1.
// The problem is of the form:
//
// min_{i=1}^n f_i(x_i) + \sum_{i=1}^{n-1} h_i(x_i - x_{i+1})
//
// where
//
// f_i(x_i) = 1/p * cDev[i] * |x - aDev[i]|^p.
// h_i(x_i - x_{i+1}) = 1/q * cSep[i] * |x_i - x_{i + 1}|^q.
//
// For convexity, it must be cDev[i] >= 0 and cSep[i] >= 0 for all i.

// Tuning hyperparameters
const data_type LPLQ_ADEV_UNIF_LEFT = -1.0;
const data_type LPLQ_ADEV_UNIF_RIGHT = 1.0;
const data_type LPLQ_CDEV_UNIF_LEFT = 0.0;
const data_type LPLQ_CDEV_UNIF_RIGHT = 1.0;
const data_type LPLQ_CSEP_UNIF_LEFT = 0.0;
const data_type LPLQ_CSEP_UNIF_RIGHT = 1.0;
// Weighted version, cSep[i]s can be different.
void genLpLqFuncs(int n, InputData* inputData, bool cDev1 = false);
// Non-weighted and weighted version, cSep[i] = lambda for all i.
void genLpLqFuncs(int n, data_type lambda, InputData* inputData,
                  bool reSample = false, bool cDev1 = false);
// Non-weighted overload for sampling a uniform lambda coefficient from [lambdaLB, lambdaUB].
void genLpLqFuncs(int n, data_type lambdaLB, data_type lambdaUB,
                  InputData* inputData, bool cDev1 = false);

// Some ${gen_data_type} may generate lambda internally.
// In those cases, you do not need to pass the ${in_lambda} argument.
// The ${out_lambda} parameter is to take the average of
// the weighted coefficients as the penalty parameter for the uniform case.
void genLpLqFuncs(int n, InputData* inputData, gen_data_type genDataType,
                  data_type in_lambda = 0, data_type* out_lambda = NULL);

// Generate data for linear_l2 (graph Laplacian) problem.
void genLinearL2Funcs(int n, InputData* inputData);

// Generate Huber parameters.
void genHuberFuncs(int n, const std::vector<data_type>& baselines,
                   InputData* inputData, bool isDev = true,
                   double lRatio = 0.5, double rRatio = 1.0);

#endif /* data_generator_hpp */
