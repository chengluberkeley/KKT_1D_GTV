//
//  comparison_profiles_util.cpp
//  KKT
//
//  Created by Cheng Lu on 4/11/20.
//  Copyright © 2020 Cheng Lu. All rights reserved.
//

#include <chrono>
#include <cmath>
#include "comparison_profiles.hpp"
#include "data_generator.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include "utils.hpp"
#include <vector>

// Global KKT solver.
KKTSolver kktSolver;

// List of methods to compare for each problem type.
std::vector<std::vector<std::string>> cpAlgs = {
    {"KKT"}, //"Kolmogorov", "Kolmogorov-nloglogn"},
    {"KKT"}, //"Projected Newton", "Linearized Taut String", "Classic Taut String",
        //"Hybrid Taut String", "Condat", "Condat's Taut String", "Johnson", "Kolmogorov"},
    {"KKT"}, //"Projected Newton", "Taut String", "Kolmogorov"},
    {"KKT"}, //"Thomas Algorithm"},
    {"KKT"}, //"Kolmogorov", "Kolmogorov-nloglogn"},
    {"KKT"}, //"Kolmogorov"},
    {"KKT"}, //"ceres", "nlopt", "dlib"},
    {"KKT"},
    {"KKT"}, //"ceres", "nlopt", "dlib"},
};

// Tuning parameters fed from command line.
// The total number of scales for n or lambda to run the experiment.
int NUM_SCALES = 7;
// The fixed n size for varying lambda test.
int FIX_N = 100000;
// Initial lambda values for varying lambda experiment.
data_type INIT_LAMBDA = 1e-4;

// For solution validity check.
// For l1-l1, since it is not strictly convex,
// one may have vastly different optimal solutions.
// So one may need to relax the solution validity
// check criterion.
bool IS_COMP_OBJ = false;
data_type OBJ_ESP = 1e-2;
data_type SOL_ESP = 1e-4;

// For Kolmogorov's nloglogn algorithm
int NLOGLOGN_STOP_SCALE = 5;

std::string toString(problem_type problemType) {
    switch (problemType) {
        case L1_L1: return "L1-L1";
        case L2_L1_NW: return "L2-L1-NW";
        case L2_L1_W: return "L2-L1-W";
        case L2_L2: return "L2-L2";
        case PWL1_L1: return "PWL1-L1";
        case PWL2_L1: return "PWL2-L1";
        case LP_LQ: return "LP-LQ";
        case LINEAR_L2: return "Linear-L2";
        case HUBER: return "Huber";
        default:
            return "";
    }
}

bool solValid(const InputData& inputData, OutputData* kkt_outputData,
              OutputData* outputData) {
    bool b1 = false;
    if (IS_COMP_OBJ) {
        kktSolver.compObj(inputData, kkt_outputData);
        kktSolver.compObj(inputData, outputData);
        b1 = (fabs(kkt_outputData->_objVal - outputData->_objVal) < OBJ_ESP);
    }
    data_type maxDiff;
    compareSolutions(*kkt_outputData, *outputData, &maxDiff);
    bool b2 = (fabs(maxDiff) < SOL_ESP);
    if ((IS_COMP_OBJ && (!b1 && !b2)) || (!IS_COMP_OBJ && !b2)) {
        // Output for debugging.
        std::cout << "MaxDiff = " << maxDiff << std::endl;
        kktSolver.compObj(inputData, kkt_outputData);
        kktSolver.compObj(inputData, outputData);
        std::cout << "kkt obj = " << kkt_outputData->_objVal << "; other obj = " << outputData->_objVal << std::endl;
    }
    return b1 || b2;
}

void CSV::init(const std::vector<std::string>& cpAlgsList, int numScales) {
    assert(!cpAlgsList.empty() && numScales > 0);
    _colTitles.resize(numScales);
    size_t algNum = cpAlgsList.size();
    size_t plusItemCount = _plusItemSuffixes.size();
    size_t totalItemCount = 2 + plusItemCount;
    _rowTitles.resize(algNum * totalItemCount);
    _figures.resize(algNum * totalItemCount);
    for (int i = 0; i < algNum; ++i) {
        _rowTitles[i * totalItemCount] = cpAlgsList[i];
        _rowTitles[i * totalItemCount + 1] = cpAlgsList[i] + "-std";
        for (int j = 0; j < plusItemCount; ++j) {
            _rowTitles[i * totalItemCount + 2 + j] = cpAlgsList[i] + _plusItemSuffixes[j];
        }
        for (int j = 0; j < totalItemCount; ++j) {
            _figures[i * totalItemCount + j].resize(numScales);
        }
    }
}

void CSV::write(const std::string &filename) {
    // Check data format validity.
    assert(!_figures.empty() && !_figures[0].empty());
    size_t figureRows = _figures.size();
    size_t figureCols = _figures[0].size();
    assert(_colTitles.empty() || _colTitles.size() == figureCols);
    assert(_rowTitles.empty() || _rowTitles.size() == figureRows);
    assert(!filename.empty());
    std::ofstream file;
    file.open(filename, std::ofstream::out | std::ofstream::app);
    if (file.is_open()) {
        // Output problem type and data type first for record.
        file << toString(_problemType) << ","
            << toString(_genDataType) << "\n";
        // Output p and q values for record.
        file << "n," << _n << ",p," << _p << ",q," << _q << "\n";
        if (!_colTitles.empty()) {
            if (!_rowTitles.empty())
                file << ",";  // One space offset
            for (int i = 0; i < figureCols; ++i) {
                file << _colTitles[i] << ",";
            }
            file << "\n";
        }
        for (int i = 0; i < figureRows; ++i) {
            if (!_rowTitles.empty()) {
                file << _rowTitles[i] << ",";
            }
            for (int j = 0; j < figureCols; ++j) {
                if (j < _figures[i].size()) {
                    file << _figures[i][j] << ",";
                } else {
                    file << ",";  // Fill in empty slot.
                }
            }
            file << "\n";
        }
        file.close();
    } else {
        std::cerr << "File " << filename << " open error!\n";
    }
}
