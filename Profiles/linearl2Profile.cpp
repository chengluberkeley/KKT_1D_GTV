//
//  linearl2Profile.cpp
//  KKT
//
//  Created by Cheng Lu on 5/7/20.
//  Copyright © 2020 Cheng Lu. All rights reserved.
//

#include "comparison_profiles.hpp"
#include <iostream>

void compLinearL2Obj(const InputData& inputData, OutputData* outputData) {
    assert(outputData != NULL);
    outputData->_objVal = 0;
    for (int i = 0; i < inputData._n; ++i) {
        outputData->_objVal += inputData._cDev[i] * outputData->_x[i];
    }
    for (int i = 0; i <inputData._n - 1; ++i) {
        outputData->_objVal += 0.5 * (outputData->_x[i] - outputData->_x[i + 1]) *
                                (outputData->_x[i] - outputData->_x[i + 1]);
    }
}

void linearl2Profile(int rounds, const std::string& path) {
    assert(rounds > 0);
    int numScales = NUM_SCALES;

    int n = 1;

    for (int i = 0; i < numScales; ++i) {
        n *= 10;
        std::cout << "n = " << n << std::endl;
        for (int iter = 0; iter < rounds; ++iter) {
            InputData inputData(n, 1, 2);
            genLinearL2Funcs(n, &inputData);
            OutputData kkt_outputData(inputData);
            auto start = std::chrono::steady_clock::now();
            kktSolver.fast_linear_l2(inputData, &kkt_outputData);
            auto end = std::chrono::steady_clock::now();
            time_ms_type time = std::chrono::duration_cast
            <std::chrono::milliseconds>(end - start).count();
            std::cout << "Complete KKT in round " << iter
            << " in time " << time << " ms\n";
            compLinearL2Obj(inputData, &kkt_outputData);
            std::cout << "Objective = " << kkt_outputData._objVal << std::endl;
        }
    }
}
