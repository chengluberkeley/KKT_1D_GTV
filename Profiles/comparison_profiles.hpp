//
//  comparison_profiles.hpp
//  KKT
//
//  Created by Cheng Lu on 4/11/20.
//  Copyright © 2020 Cheng Lu. All rights reserved.
//

// Each problem type is handled by a separate function.

#ifndef comparison_profiles_hpp
#define comparison_profiles_hpp

#include "data_generator.hpp"
#include <string>
#include <vector>

// typedef signed long long time_ms_type;  // at least 64 bits.

// extern KKT solver
extern KKTSolver kktSolver;

// extern variables
extern std::vector<std::vector<std::string>> cpAlgs;
extern int NUM_SCALES;
extern int FIX_N;
extern data_type INIT_LAMBDA;

// For solution validity checks.
extern bool IS_COMP_OBJ;
extern data_type OBJ_ESP;
extern data_type SOL_ESP;

extern int NLOGLOGN_STOP_SCALE;

// The values are used as index for vector.
// DON'T CHANGE!!!
typedef enum PROBLEM_TYPE {
    L1_L1 = 0,
    L2_L1_NW,
    L2_L1_W,
    L2_L2,
    PWL1_L1,
    PWL2_L1,
    LP_LQ,
    LINEAR_L2,
    HUBER,
} problem_type;

// Map from problem type to string for output.
std::string toString(problem_type problemType);

// Functions for each comparison profile
// Parameter ${rounds} refers to the number of rounds to compute
// average run times.
void l1l1Profile(int rounds, const std::string& path);
void l2l1nwProfile(int rounds, const std::string& path);
void l2l1wProfile(int rounds, const std::string& path);
void l2l2Profile(int rounds, const std::string& path);
void pwl1Profile(int rounds, const std::string& path);
void pwl2Profile(int rounds, const std::string& path);
void lplqProfile(int rounds, const std::string& path);
void linearl2Profile(int rounds, const std::string& path);
void huberProfile(int rounds, const std::string& path);

// Utility functions
template <typename T>
void stat(const std::vector<T>& runTimes,
          double* aveTime, double* stdTime = NULL) {
    assert(aveTime != NULL);
    double ave = 0;
    for (int i = 0; i < runTimes.size(); ++i) {
        ave += runTimes[i];
    }
    ave = (ave * 1.0) / runTimes.size();
    *aveTime = ave;
    if (stdTime != NULL) {
        double stddev = 0;
        for (int i = 0; i < runTimes.size(); ++i) {
            stddev += (runTimes[i] - ave) * (runTimes[i] - ave);
        }
        stddev = (stddev * 1.0) / runTimes.size();
        stddev = sqrt(stddev);
        *stdTime = stddev;
    }
}

bool solValid(const InputData& inputData, OutputData* kkt_outputData,
              OutputData* outputData);

//////////////////////////////////////////////////////
// Data structure to write to csv.
struct CSV {
    // Column titles have one space offset if row titles are present.
    std::vector<std::string> _rowTitles;
    std::vector<data_type> _colTitles;
    std::vector<std::vector<data_type>> _figures;
    problem_type _problemType;
    gen_data_type _genDataType;
    int _n;
    int _p, _q;
    std::vector<std::string> _plusItemSuffixes;

    void init(const std::vector<std::string>& cpAlgsList, int numScales);
    void write(const std::string& filename);
    void clear() {
        _rowTitles.clear();
        _colTitles.clear();
        _figures.clear();
    }
};

#endif /* comparison_profiles_hpp */
