//
//  main.cpp
//  KKT
//
//  Created by Cheng Lu on 4/2/20.
//  Copyright © 2020 Cheng Lu. All rights reserved.
//

//  Compare performance of algorithms on various TV-based problems.
//  Problems to compare:
//  l1-l1
//     - Kolmogorov's algorithm
//  l2-l1 (weighted and non-weighted)
//  non-weighted:
//     - Condat's algorithm
//     - Variants of taut string algorithms
//     - Johnson's algorithm
//     - ProxTV's algorithms
//     - Kolmogorov's algorithm
//  weighted:
//     - Taut string algorithm
//     - ProxTV's algorithm
//     - Kolmogorov's algorithm
//  l2-l2:
//     - Thomas algorithm
//  Piecewise linear - l1
//     - Kolmogorov's algorithm
//  Piecewise quadratic - l1
//     - Kolmogorov's algorithm
//  lp-lq (p > 2 && q > 2):
//     - No other implementations to compare to.
//       Self-comparison with increasing p and q.
//  Huber-l1:
//     - ceres
//     - nlopt
//  l2-Huber:
//     - ceres
//     - nlopt

#include "comparison_profiles.hpp"
#include <iostream>
#include <string>

// Settable from command line arguments.
int ROUNDS = 5;
std::string PATH = "/tmp/";

// Output command line spec.
void help() {
    std::cout << "./main -profile <profile_name>\n"
        << "\t [-rounds <num_rounds> (default = 5)]\n"
        << "\t [-path <output_path> (default = /Users/user1/kkt_outs/)]\n"
        << "\t [-num_scales <num_scales> (default = 7)]\n"
        << "\t [-fix_n <fix_n> (default = 100000)]\n"
        << "\t [-init_lambda <init_lambda> (default = 1e-4)]\n"
        << "\t [--isCompObj (default = false)]\n"
        << "\t [-obj_esp <objective_esp> (default = 1e-2)]\n"
        << "\t [-sol_esp <solution_esp> (default = 1e-4)]\n"
        << "\t [-nloglogn_stop_scale <nloglogn_stop_scale> (default = 5)]\n";
    std::cout << "Available profiles: \n"
        << "1. l1-l1\n"
        << "2. l2-l1-nw\n"
        << "3. l2-l1-w\n"
        << "4. l2-l2\n"
        << "5. pwl1-l1\n"
        << "6. pwl2-l1\n"
        << "7. lp-lq\n"
        << "8. linear-l2\n"
        << "9. Huber\n";
}

void printParams() {
    std::cout << "Rounds: " << ROUNDS << " | "
        << "path: " << PATH << " | "
        << "num_scales: " << NUM_SCALES << " | "
        << "fix_n: " << FIX_N << " | "
        << "init_lambda: " << INIT_LAMBDA << " | "
        << "isCompObj: " << (IS_COMP_OBJ ? "true" : "false") << " | "
        << "obj_esp: " << OBJ_ESP << " | "
        << "sol_esp: " << SOL_ESP << " | "
        << "nloglogn_stop_scale: " << NLOGLOGN_STOP_SCALE << std::endl;
}

problem_type setProblemType(std::string problemTypeStr) {
    if (problemTypeStr.compare("l1-l1") == 0) {
        return L1_L1;
    }
    if (problemTypeStr.compare("l2-l1-nw") == 0) {
        return L2_L1_NW;
    }
    if (problemTypeStr.compare("l2-l1-w") == 0) {
        return L2_L1_W;
    }
    if (problemTypeStr.compare("l2-l2") == 0) {
        return L2_L2;
    }
    if (problemTypeStr.compare("pwl1-l1") == 0) {
        return PWL1_L1;
    }
    if (problemTypeStr.compare("pwl2-l1") == 0) {
        return PWL2_L1;
    }
    if (problemTypeStr.compare("lp-lq") == 0) {
        return LP_LQ;
    }
    if (problemTypeStr.compare("linear-l2") == 0) {
        return LINEAR_L2;
    }
    if (problemTypeStr.compare("Huber") == 0) {
        return HUBER;
    }
    return LP_LQ;  // Default profile.
}

// Compare performance of different algorithms.
int main(int argc, const char * argv[]) {
    // Parse input arguments
    if (argc < 3) {
        std::cout << "Not enough parameters!\n";
        help();
        return -1;
    }

    int arg_inc = 2;
    problem_type problemType = LP_LQ;
    for (int i = 1; i < argc; i += arg_inc) {
        if (strcmp(argv[i], "-profile") == 0) {
            problemType = setProblemType(argv[i + 1]);
            arg_inc = 2;
        } else if (strcmp(argv[i], "-rounds") == 0) {
            ROUNDS = std::stoi(argv[i + 1]);
            arg_inc = 2;
        } else if (strcmp(argv[i], "-path") == 0) {
            PATH = std::string(argv[i + 1]);
            arg_inc = 2;
        } else if (strcmp(argv[i], "-num_scales") == 0) {
            NUM_SCALES = std::stoi(argv[i + 1]);
            arg_inc = 2;
        } else if (strcmp(argv[i], "-fix_n") == 0) {
            FIX_N = std::stoi(argv[i + 1]);
            arg_inc = 2;
        } else if (strcmp(argv[i], "-init_lambda") == 0) {
            INIT_LAMBDA = std::stod(argv[i + 1]);
            arg_inc = 2;
        } else if (strcmp(argv[i], "--isCompObj") == 0) {
            IS_COMP_OBJ = true;
            arg_inc = 1;
        } else if (strcmp(argv[i], "-obj_esp") == 0) {
            OBJ_ESP = std::stod(argv[i + 1]);
            arg_inc = 2;
        } else if (strcmp(argv[i], "-sol_esp") == 0) {
            SOL_ESP = std::stod(argv[i + 1]);
            arg_inc = 2;
        } else if (strcmp(argv[i], "-nloglogn_stop_scale") == 0) {
            NLOGLOGN_STOP_SCALE = std::stoi(argv[i + 1]);
            arg_inc = 2;
        } else {
            std::cout << "Invalid flag.\n";
            help();
            return -1;
        }
    }

    std::cout << "Profile parameters:\n";
    printParams();

    switch (problemType) {
        case L1_L1: {
            std::cout << "Start l1-l1 profile:\n";
            l1l1Profile(ROUNDS, PATH);
            std::cout << "Complete l1-l1 profile.\n";
            break;
        }
        case L2_L1_NW: {
            std::cout << "Start l2-l1-nw profile:\n";
            l2l1nwProfile(ROUNDS, PATH);
            std::cout << "Complete l2-l1-nw profile.\n";
            break;
        }
        case L2_L1_W: {
            std::cout << "Start l2-l1-w profile:\n";
            l2l1wProfile(ROUNDS, PATH);
            std::cout << "Complete l2-l1-w profile.\n";
            break;
        }
        case L2_L2: {
            std::cout << "Start l2-l2 profile:\n";
            l2l2Profile(ROUNDS, PATH);
            std::cout << "Complete l2-l2 profile.\n";
            break;
        }
        case PWL1_L1: {
            std::cout << "Start pwl1-l1 profile:\n";
            pwl1Profile(ROUNDS, PATH);
            std::cout << "Complete pwl1-l1 profile.\n";
            break;
        }
        case PWL2_L1: {
            std::cout << "Start pwl2-l1 profile:\n";
            pwl2Profile(ROUNDS, PATH);
            std::cout << "Complete pwl2-l1 profile.\n";
            break;
        }
        case LINEAR_L2: {
            std::cout << "Start linear-l2 profile:\n";
            linearl2Profile(ROUNDS, PATH);
            std::cout << "Complete linear-l2 profile.\n";
            break;
        }
        case HUBER: {
            std::cout << "Start Huber profile:\n";
            huberProfile(ROUNDS, PATH);
            std::cout << "Complete Huber profile.\n";
            break;
        }
        case LP_LQ:
        default: {
            // Default to lp-lq.
            std::cout << "Start lp-lq profile:\n";
            lplqProfile(ROUNDS, PATH);
            std::cout << "Complete lp-lq profile.\n";
        }
    }

    return 0;
}
