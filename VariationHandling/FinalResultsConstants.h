#include "BASE.h"

#include <vector>

namespace AnalysisConstants_UL
{
    namespace FinalResultsConstans
    {
        std::vector<std::vector<float>> partonYAxisValuesNormalized = {
            {6E-6, 2},
            {1E-6, 2},
            {1E-3, 0.5},
            {2E-6, 1.},
            {2E-8, 1.},
            {0., 0.4},
            {0., 0.4},
            {0., 0.5},
            {0., 0.3},
            {0., 0.3}};
        
        std::vector<std::vector<float>> particleYAxisValuesNormalized = {
            {6E-6, 2},
            {1E-6, 2},
            {1E-3, 0.5},
            {2E-6, 1.},
            {2E-8, 1.},
            {0., 0.4},
            {0., 0.4},
            {0., 0.5},
            {0., 0.3},
            {0., 0.3}};

        std::vector<std::vector<float>> partonYAxisValues = {
            {6E-6, 1E-1},
            {1E-6, 1E-2},
            {1E-3, 2.5},
            {3E-6, 5E-2},
            {2E-8, 5E-2},
            {0., 2.},
            {0., 2.},
            {0., 2.},//chi 
            {0., 6.}, //cosTheta0
            {0., 6.}}; //cosTheta1;

        std::vector<std::vector<float>> particleYAxisValues = {
            {6E-6, 1E-2},
            {1E-6, 1E-2},
            {1E-3, 0.8},
            {3E-6, 5E-2},
            {2E-8, 5E-2},
            {0., 0.5},
            {0., 0.5},
            {0., 0.5},
            {0., 2.},
            {0., 2.}};

        std::vector<std::vector<float>> partonXAxisValues = {
            {450, 1500},
            {400, 1500},
            {-2.4, 2.4},
            {0, 1300},
            {1000, 5000},
            {0., 2.4},
            {0., 2.4},
            {1, 13},
            {0, 1},
            {0, 1}};

        std::vector<std::vector<float>> particleXAxisValues = {
            {450, 1500},
            {400, 1500},
            {-2.4, 2.4},
            {0, 1300},
            {1000, 5000},
            {0., 2.4},
            {0., 2.4},
            {1, 13},
            {0, 1},
            {0, 1}};

        std::vector<std::vector<float>> legendPositions = {
            {0.59, 0.59, 0.89, 0.89},
            {0.59, 0.59, 0.89, 0.89},
            {0.59, 0.59, 0.89, 0.89},
            {0.59, 0.59, 0.89, 0.89},
            {0.59, 0.59, 0.89, 0.89},
            {0.59, 0.59, 0.89, 0.89},
            {0.59, 0.59, 0.89, 0.89},
            {0.59, 0.59, 0.89, 0.89},
            {0.59, 0.59, 0.89, 0.89},
            {0.59, 0.59, 0.89, 0.89}};
    }
}