//
// Created by aabedraba on 6/11/19.
//

#ifndef INC_2_GENETIC_ALGORITHMS_UTILS_H
#define INC_2_GENETIC_ALGORITHMS_UTILS_H

#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

typedef vector<int> Chromosome;
struct Individual {
    Chromosome chromosome;
    unsigned int cost;
    bool evaluated;
};
typedef vector<Individual> Population;

namespace Utils {
    int solutionCost(const vector<int> &solution, const vector<vector<int>> &fluxMatrix,
                     const vector<vector<int>> &distanceMatrix, const bool simetric);

    void logPopulation(Population &population, int generation, string &log);

    void logCross(string position, const Chromosome *first, const Chromosome *second, string &log);

    void logVector(string position, const vector<int> *vector, string &log);

    void logElitesAndRabbles(Population *initialPopulation, Population *currentPopulation, vector<int> &elitePos,
                             vector<int> &rabblePos, string &log);

    void logBestSolution(const Chromosome *solution, int cost, string &log);

    int posMaxElement(const vector<pair<int, int>> &vector);

    int posMinElement(const vector<pair<int, int>> &vector);

    void writeInFile(int seed, const string type, const string airport, const string &log);

    void getResults(string type, string airportName, int bestCost, double executionTime, string &log);

//    int calculateSwapCost(vector<int> &solution, const int solutionCost, const pair<int, int> &swapPos, const vector<vector<int>> &fluxMatrix,
//                          const vector<vector<int>> &distanceMatrix, const bool simetric);
//    int partialSwapCost(const vector<int> &solution, const pair<int, int> swapPositions,
//                        const vector<vector<int>> &fluxMatrix, const vector<vector<int>> &distanceMatrix,
//                        const bool simetric);
//    const string getLog(string algorithmType, const string log, float time, string airportName, int sizeVectors, const vector<int> &solutionVector,
//                        int solutionCost);
//    string logSolution(string type, pair<int, int> movement, vector<int> &solutionVector, int cost, int iteration, int environmentIteration);

}


#endif //INC_2_GENETIC_ALGORITHMS_UTILS_H
