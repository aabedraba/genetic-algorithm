//
// Created by aabedraba on 6/11/19.
//

#include "Utils.h"

int Utils::posMaxElement(const vector<pair<int,int>> &vector) {
    int max = 0;
    for (int i = 0; i < vector.size(); ++i)
        if (vector[i].first > vector[max].first)
            max = i;
    return max;
}

int Utils::posMinElement(const vector<pair<int,int>> &vector) {
    int min = 0;
    for (int i = 0; i < vector.size(); ++i)
        if (vector[i].first < vector[min].first)
            min = i;
    return min;
}

int Utils::solutionCost(const vector<int> &solution, const vector<vector<int>> &fluxMatrix,
                        const vector<vector<int>> &distanceMatrix, const bool simetric) {
    int cost = 0;
    if (simetric)
        for (int i = 0; i < solution.size(); ++i)
            for (int j = i + 1; j < solution.size(); ++j)
                cost += 2 * (fluxMatrix[i][j] * distanceMatrix[solution[i]][solution[j]]);
    else
        for (int i = 0; i < solution.size(); ++i)
            for (int j = 0; j < solution.size(); ++j)
                if (i != j)
                    cost += fluxMatrix[i][j] * distanceMatrix[solution[i]][solution[j]];
    return cost;
}

void Utils::logPopulation(Population &population, int generation, string &log) {
    log += "====================Generation " + to_string(generation) + "====================\n";
    for (int i = 0; i < population.size(); ++i) {
        log += "Individual " + to_string(i) + ": ";
        for (int j = 0; j < population[i].chromosome.size(); ++j) {
            log += to_string(population[i].chromosome[j]) + " ";
        }
        log += "\n";
        log += "Cost: " + to_string(population[i].cost);
        log += "\n\n";
    }
}

void Utils::logCross(string position, const Chromosome *first, const Chromosome *second, string &log) {
    string log1 = position + " 1: ";
    string log2 = position + " 2: ";
    for (int i = 0; i < first->size(); ++i) {
        log1 += to_string(first->at(i)) + " ";
        log2 += to_string(second->at(i)) + " ";
    }
    log += log1 + "\n" + log2 + "\n";
}

void Utils::logVector(string position, const vector<int> *vector, string &log) {
    log += position + ": ";
    for (int i = 0; i < vector->size(); ++i) {
        log += to_string(vector->at(i)) + " ";
    }
    log += "\n";
}

void Utils::writeInFile(int seed, const string type, string airportName, const string &log) {
    ofstream fs;
    fs.open("../_logs/" + type + '-' + airportName + '-' + to_string(seed) + ".txt");
    if (fs.is_open()) {
        fs << log << endl;
    } else
        throw "File not properly opened!";
}

void Utils::logElitesAndRabbles(Population *initialPopulation, Population *currentPopulation, vector<int> &elitePos,
                                vector<int> &rabblePos, string &log) {
    string elites = "Elites: \n";
    string rabbles = "Rabbles: \n";
    for (int i = 0; i < elitePos.size(); ++i) {
        elites += "  " + to_string(i) + ":    ";
        rabbles += "  " + to_string(i) + ":    ";
        for (int j = 0; j < initialPopulation->at(0).chromosome.size(); ++j) {
            elites += to_string(initialPopulation->at(elitePos[i]).chromosome[j]) + " ";
            rabbles += to_string(currentPopulation->at(rabblePos[i]).chromosome[j]) + " ";
        }
        elites += "\n";
        elites += "  Cost: " + to_string(initialPopulation->at(elitePos[i]).cost) +"\n";
        rabbles += "\n";
        rabbles += "  Cost: " + to_string(currentPopulation->at(rabblePos[i]).cost) +"\n";
    }
    log += elites + "\n" + rabbles;
}

void Utils::logBestSolution(const Chromosome *solution, int cost, string &log) {
    log += "Best: ";
    for (int i = 0; i < solution->size(); ++i) {
        log += to_string(solution->at(i)) + " ";
    }
    log += "\n";
    log += "Cost: " + to_string(cost);
}

void Utils::getResults(string type, string airportName, int bestCost, double executionTime, string &log) {
    log += "-------------" + airportName + " " + type +  "-------------\n";
    log += "\n";
    log += "Best cost: " + to_string(bestCost) + "\n";
    log += "Execution time: " + to_string(executionTime) + "\n";
}

//int Utils::calculateSwapCost(vector<int> &solution, const int solutionCost, const pair<int, int> &swapPos,
//                             const vector<vector<int>> &fluxMatrix,
//                             const vector<vector<int>> &distanceMatrix, const bool simetric) {
//    int preSwapCost = 0, afterSwapCost = 0;
//    preSwapCost = partialSwapCost(solution, swapPos, fluxMatrix,
//                                  distanceMatrix,
//                                  simetric);
//    swap(solution[swapPos.first], solution[swapPos.second]);
//    afterSwapCost = partialSwapCost(solution, swapPos, fluxMatrix,
//                                    distanceMatrix,
//                                    simetric);
//    swap(solution[swapPos.first], solution[swapPos.second]);
//    return solutionCost + afterSwapCost - preSwapCost;
//}
//int Utils::partialSwapCost(const vector<int> &solution, const pair<int, int> swapPositions,
//                           const vector<vector<int>> &fluxMatrix,
//                           const vector<vector<int>> &distanceMatrix, const bool simetric) {
//    int cost = 0;
//    int swap1 = swapPositions.first;
//    int swap2 = swapPositions.second;
//    if (simetric) {
//        for (int i = 0; i < solution.size(); i++) {
//            if (i != swap1) cost += 2 * fluxMatrix[swap1][i] * distanceMatrix[solution[swap1]][solution[i]];
//            if (i != swap2) cost += 2 * fluxMatrix[swap2][i] * distanceMatrix[solution[swap2]][solution[i]];
//        }
//    } else {
//        for (int i = 0; i < solution.size(); i++) {
//            if (i != swap1) {
//                cost += fluxMatrix[swap1][i] * distanceMatrix[solution[swap1]][solution[i]];
//                cost += fluxMatrix[i][swap1] * distanceMatrix[solution[i]][solution[swap1]];
//            }
//            if (i != swap2) {
//                cost += fluxMatrix[swap2][i] * distanceMatrix[solution[swap2]][solution[i]];
//                cost += fluxMatrix[i][swap2] * distanceMatrix[solution[i]][solution[swap2]];
//            }
//        }
//    }
//    return cost;
//}
//const string Utils::getLog(string algorithmType, const string log,  float executionTime, string airportName, int sizeVectors, const vector<int> &solutionVector,
//                           int solutionCost) {
//    string solutionLog;
//    solutionLog += "---------" + algorithmType + " " + airportName + "----------\n";
//    solutionLog += "Solution:          ";
//    for (int i = 0; i < sizeVectors; ++i) {
//        solutionLog += to_string(solutionVector[i]) + " ";
//    }
//    solutionLog += "\n";
//    solutionLog += "Cost:              " + to_string(solutionCost) + '\n';
//    solutionLog += "Time:              " + to_string(executionTime) + " seconds" + "\n";
//    solutionLog += "------------------------------------------------------------------------\n\n";
//    return solutionLog + log;
//}
//
//string Utils::logSolution(string type, pair<int, int> movement, vector<int> &solutionVector, int cost, int iteration,
//                          int environmentIteration) {
//    string log;
//    string space((20 - type.length()), ' ');
//    if (type == "Initial Solution") {
//        log += type + space;
//        for (int i = 0; i < solutionVector.size(); ++i)
//            log += to_string(solutionVector[i]) + " ";
//        log += "\n";
//    } else
//        log += "Environment:        " + to_string(environmentIteration) + '\n';
//    if (movement.first != -1)
//        log += "Movement:           " + to_string(movement.first) + ' ' + to_string(movement.second) + '\n';
//    log += "Cost:               " + to_string(cost) + "\n";
//    log += "Iteration           " + to_string(iteration) + "\n";
//    log += "\n";
//    return log;
//}

