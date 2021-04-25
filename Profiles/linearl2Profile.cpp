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
    std::vector<std::vector<time_ms_type>> runTimes;
    CSV csvData;
    csvData._problemType = LINEAR_L2;
    csvData._p = 1;
    csvData._q = 2;
    int numScales = NUM_SCALES;
    size_t algNum = cpAlgs[csvData._problemType].size();
    const std::vector<std::string>& cpAlgsList = cpAlgs[csvData._problemType];
    csvData.init(cpAlgsList, numScales);
    for (int i = 0; i < algNum; ++i) {
        runTimes.push_back(std::vector<time_ms_type>(rounds, 0));
    }
    csvData._genDataType = KKT_LINEAR_L2;

    int n = 1;

    std::cout << "Run " << toString(csvData._problemType) << " with data "
        << toString(csvData._genDataType)
        << " for varying n" << std::endl;

    for (int i = 0; i < numScales; ++i) {
        n *= 10;
        csvData._colTitles[i] = n;
        csvData._n = n;
        std::cout << "n = " << n << std::endl;
        for (int iter = 0; iter < rounds; ++iter) {
            InputData inputData(n, 1, 2);
            genLinearL2Funcs(n, &inputData);
            OutputData kkt_outputData(inputData);
            auto start = std::chrono::steady_clock::now();
            kktSolver.fast_linear_l2(inputData, &kkt_outputData);
            auto end = std::chrono::steady_clock::now();
            runTimes[0][iter] = std::chrono::duration_cast
            <std::chrono::milliseconds>(end - start).count();
            std::cout << "Complete KKT in round " << iter
            << " in time " << runTimes[0][iter] << " ms\n";
            compLinearL2Obj(inputData, &kkt_outputData);
            std::cout << "Objective = " << kkt_outputData._objVal << std::endl;

            std::cout << "****\n";
        }
        for (int j = 0; j < algNum; ++j) {
            double aveTime, stdTime;
            stat(runTimes[j], &aveTime, &stdTime);
            csvData._figures[j * 2][i] = aveTime;
            csvData._figures[j * 2 + 1][i] = stdTime;
        }
        std::cout << "===========\n";
    }
    std::string filename = path + "/out_" + toString(csvData._problemType)
        + "_" + toString(csvData._genDataType) + ".txt";
    csvData.write(filename);
    std::cout << "Written in file " << filename << std::endl;
    std::cout << "////////////////////\n";
}
