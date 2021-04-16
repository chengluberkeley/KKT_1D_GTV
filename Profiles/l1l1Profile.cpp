//
//  l1l1Profile.cpp
//  KKT
//
//  Created by Cheng Lu on 4/16/20.
//  Copyright © 2020 Cheng Lu. All rights reserved.
//

#include "comparison_profiles.hpp"
#include <iostream>

void l1l1Profile(int rounds, const std::string& path) {
    assert(rounds > 0);
    std::vector<std::vector<time_ms_type>> runTimes;
    CSV csvData;
    csvData._problemType = L1_L1;
    csvData._p = 1;
    csvData._q = 1;
    int numScales = NUM_SCALES;
    size_t algNum = cpAlgs[csvData._problemType].size();
    const std::vector<std::string>& cpAlgsList =
    cpAlgs[csvData._problemType];
    csvData.init(cpAlgsList, numScales);
    for (int i = 0; i < algNum; ++i) {
        runTimes.push_back(std::vector<time_ms_type>(rounds, 0));
    }
    int n;

    // Case 1: Varying input sizes
    std::vector<gen_data_type> inputSizeDataTypes = {KKT_LP_LQ};
    for (int dtIndex = 0; dtIndex < inputSizeDataTypes.size(); ++dtIndex)
    {
        csvData._genDataType = inputSizeDataTypes[dtIndex];
        std::cout << "Run " << toString(csvData._problemType) << " with data "
            << toString(csvData._genDataType)
            << " for varying n" << std::endl;
        n = 1;
        for (int i = 0; i < numScales; ++i) {
            n *= 10;
            csvData._colTitles[i] = n;
            csvData._n = n;
            std::cout << "n = " << n << std::endl;
            for (int iter = 0; iter < rounds; ++iter) {
                InputData inputData(n, 1, 1);
                switch (csvData._genDataType) {
                    case KKT_LP_LQ: {
                        genLpLqFuncs(n, &inputData);
                        inputData._lb = -1;
                        inputData._ub = 1;
                        break;
                    }
                    default: {
                        genLpLqFuncs(n, &inputData, csvData._genDataType);
                        inputData._lb = -200;
                        inputData._ub = 200;
                    }
                }

                OutputData kkt_outputData(inputData);
                // Function call to KKT
                auto start = std::chrono::steady_clock::now();
                kktSolver.solve(inputData, &kkt_outputData);
                auto end = std::chrono::steady_clock::now();
                runTimes[0][iter] = std::chrono::duration_cast
                    <std::chrono::milliseconds>(end - start).count();
                std::cout << "Complete KKT in round " << iter
                    << " in time " << runTimes[0][iter] << " ms\n";

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

    // Case 2: Fix input size, varying lambda.
    std::vector<gen_data_type> lambdaDataTypes = {KKT_LP_LQ};
    for (int dtIndex = 0; dtIndex < lambdaDataTypes.size(); ++dtIndex)
    {
        csvData._genDataType = lambdaDataTypes[dtIndex];
        std::cout << "Run " << toString(csvData._problemType) << " with data "
            << toString(csvData._genDataType)
            << " for varying lambdas" << std::endl;
        n = FIX_N;
        csvData._n = n;
        data_type in_lambda = INIT_LAMBDA;
        for (int i = 0; i < numScales; ++i) {
            in_lambda *= 10;
            csvData._colTitles[i] = in_lambda;
            std::cout << "lambda = " << in_lambda << std::endl;
            for (int iter = 0; iter < rounds; ++iter) {
                InputData inputData(n, 1, 1);
                switch (csvData._genDataType) {
                    case KKT_LP_LQ: {
                        genLpLqFuncs(n, in_lambda, &inputData, true);
                        inputData._lb = -1;
                        inputData._ub = 1;
                        break;
                    }
                    default: {
                        genLpLqFuncs(n, &inputData, csvData._genDataType,
                                     in_lambda);
                        inputData._lb = -2;
                        inputData._ub = 2;
                    }
                }

                OutputData kkt_outputData(inputData);
                // Function call to KKT
                auto start = std::chrono::steady_clock::now();
                kktSolver.solve(inputData, &kkt_outputData);
                auto end = std::chrono::steady_clock::now();
                runTimes[0][iter] = std::chrono::duration_cast
                    <std::chrono::milliseconds>(end - start).count();
                std::cout << "Complete KKT in round " << iter
                    << " in time " << runTimes[0][iter] << " ms\n";

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
}
