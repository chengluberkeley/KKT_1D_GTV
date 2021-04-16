//
//  huberProfile.cpp
//  KKT
//
//  Created by Cheng Lu on 7/19/20.
//  Copyright © 2020 Cheng Lu. All rights reserved.
//

#include <cmath>
#include "comparison_profiles.hpp"
#include <iostream>

void huberProfile(int rounds, const std::string& path) {
    assert(rounds > 0);
    std::vector<std::vector<time_ms_type>> runTimes;
    std::vector<std::vector<time_ms_type>> runTimesMax;
    std::vector<std::vector<data_type>> relErrors;
    std::vector<std::vector<data_type>> relErrorsMax;
    CSV csvData;
    csvData._problemType = HUBER;
    csvData._plusItemSuffixes = {"-max", "-max-std", "-rel-err", "-rel-err-max"};
    size_t totalItemCount = 2 + csvData._plusItemSuffixes.size();
    int numScales = NUM_SCALES;
    size_t algNum = cpAlgs[csvData._problemType].size();
    const std::vector<std::string>& cpAlgsList =
        cpAlgs[csvData._problemType];
    csvData.init(cpAlgsList, numScales);
    for (int i = 0; i < algNum; ++i) {
        runTimes.push_back(std::vector<time_ms_type>(rounds, 0));
        runTimesMax.push_back(std::vector<time_ms_type>(rounds, 0));
        relErrors.push_back(std::vector<data_type>(rounds, 0));
        relErrorsMax.push_back(std::vector<data_type>(rounds, 0));
    }
    int n;

    // Case 1: Huber + l1
    std::vector<gen_data_type> inputSizeDataTypes = {KKT_HUBER};
    for (int dtIndex = 0; dtIndex < inputSizeDataTypes.size(); ++dtIndex)
    {
        csvData._genDataType = inputSizeDataTypes[dtIndex];
        std::cout << "Run " << toString(csvData._problemType) << " with data "
            << toString(csvData._genDataType)
            << " for Huber-l1 and varying n" << std::endl;
        csvData._p = 2;
        csvData._q = 1;
        std::cout << "p = " << csvData._p << ", q = " << csvData._q << std::endl;
        n = 1;
        for (int i = 0; i < numScales; ++i) {
            n *= 10;
            csvData._colTitles[i] = n;
            csvData._n = n;
            std::cout << "n = " << n << std::endl;
            for (int iter = 0; iter < rounds; ++iter) {
                InputData inputData(n, csvData._p, csvData._q, InputData::HUBER_D, InputData::LQ);
                genLpLqFuncs(n, &inputData);
                std::vector<data_type> baselines(n, 0);
                for (int i = 0; i < n; ++i) {
                    baselines[i] = fabs(inputData._aDev[i]);
                }
                genHuberFuncs(n, baselines, &inputData, true);
                inputData._lb = -1;
                inputData._ub = 1;

                // KKT
                OutputData kkt_outputData(inputData);
                // Function call to KKT
                auto start = std::chrono::steady_clock::now();
                kktSolver.solve(inputData, &kkt_outputData);
                auto end = std::chrono::steady_clock::now();
                runTimes[0][iter] = std::chrono::duration_cast
                    <std::chrono::milliseconds>(end - start).count();
                kktSolver.compObj(inputData, &kkt_outputData);
                std::cout << "Complete KKT in round " << iter
                    << " in time " << runTimes[0][iter] << " ms,"
                    << " with objective value = " << kkt_outputData._objVal
                    << "\n";

                std::cout << "****\n";
            }
            for (int j = 0; j < algNum; ++j) {
                double aveTime, stdTime;
                stat(runTimes[j], &aveTime, &stdTime);
                csvData._figures[j * totalItemCount][i] = aveTime;
                csvData._figures[j * totalItemCount + 1][i] = stdTime;
                stat(runTimesMax[j], &aveTime, &stdTime);
                csvData._figures[j * totalItemCount + 2][i] = aveTime;
                csvData._figures[j * totalItemCount + 3][i] = stdTime;
                double aveRelError;
                stat(relErrors[j], &aveRelError);
                csvData._figures[j * totalItemCount + 4][i] = aveRelError;
                stat(relErrorsMax[j], &aveRelError);
                csvData._figures[j * totalItemCount + 5][i] = aveRelError;
            }
            std::cout << "===========\n";
        }
        std::string filename = path + "/out_" + toString(csvData._problemType)
        + "-l1_" + toString(csvData._genDataType) + ".txt";
        csvData.write(filename);
        std::cout << "Written in file " << filename << std::endl;
        std::cout << "////////////////////\n";
    }

    // Case 2: l2 + Huber
    for (int dtIndex = 0; dtIndex < inputSizeDataTypes.size(); ++dtIndex)
    {
        csvData._genDataType = inputSizeDataTypes[dtIndex];
        std::cout << "Run " << toString(csvData._problemType) << " with data "
            << toString(csvData._genDataType)
            << " for l2-Huber and varying n" << std::endl;
        csvData._p = 2;
        csvData._q = 2;
        std::cout << "p = " << csvData._p << ", q = " << csvData._q << std::endl;
        n = 1;
        for (int i = 0; i < numScales; ++i) {
            n *= 10;
            csvData._colTitles[i] = n;
            csvData._n = n;
            std::cout << "n = " << n << std::endl;
            for (int iter = 0; iter < rounds; ++iter) {
                InputData inputData(n, csvData._p, csvData._q, InputData::LP, InputData::HUBER_S);
                genLpLqFuncs(n, &inputData);
                std::vector<data_type> baselines(n - 1, 0);
                for (int i = 0; i < n - 1; ++i) {
                    baselines[i] = fabs(inputData._aDev[i] - inputData._aDev[i + 1]);
                }
                genHuberFuncs(n - 1, baselines, &inputData, false);
                inputData._lb = -1;
                inputData._ub = 1;

                // KKT
                OutputData kkt_outputData(inputData);
                // Function call to KKT
                auto start = std::chrono::steady_clock::now();
                kktSolver.solve(inputData, &kkt_outputData);
                auto end = std::chrono::steady_clock::now();
                runTimes[0][iter] = std::chrono::duration_cast
                    <std::chrono::milliseconds>(end - start).count();
                kktSolver.compObj(inputData, &kkt_outputData);
                std::cout << "Complete KKT in round " << iter
                    << " in time " << runTimes[0][iter] << " ms,"
                    << " with objective value = " << kkt_outputData._objVal
                    << "\n";

                std::cout << "****\n";
            }
            for (int j = 0; j < algNum; ++j) {
                double aveTime, stdTime;
                stat(runTimes[j], &aveTime, &stdTime);
                csvData._figures[j * totalItemCount][i] = aveTime;
                csvData._figures[j * totalItemCount + 1][i] = stdTime;
                stat(runTimesMax[j], &aveTime, &stdTime);
                csvData._figures[j * totalItemCount + 2][i] = aveTime;
                csvData._figures[j * totalItemCount + 3][i] = stdTime;
                double aveRelError;
                stat(relErrors[j], &aveRelError);
                csvData._figures[j * totalItemCount + 4][i] = aveRelError;
                stat(relErrorsMax[j], &aveRelError);
                csvData._figures[j * totalItemCount + 5][i] = aveRelError;
            }
            std::cout << "===========\n";
        }
        std::string filename = path + "/out_l2-" + toString(csvData._problemType)
        + "_" + toString(csvData._genDataType) + ".txt";
        csvData.write(filename);
        std::cout << "Written in file " << filename << std::endl;
        std::cout << "////////////////////\n";
    }
}
