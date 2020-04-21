//
// Created by aabedraba on 6/11/19.
//

#ifndef INC_2_GENETIC_ALGORITHMS_GENETIC_H
#define INC_2_GENETIC_ALGORITHMS_GENETIC_H

#include <vector>
#include <tuple>
#include "Airport.h"
#include "Utils.h"

using namespace std;

class Genetic {
public:
    Genetic(const Airport *airport, const string crossoverType, const unsigned int eliteSize,
            const unsigned int numEvaluations, const unsigned int numTounarments);

    virtual ~Genetic();

    const string &getLog() const;

    const Airport *getAirport() const;

    const Chromosome &getBestSolution() const;

    double getExecutionTime() const;

    int getBestCost() const;

private:
    void generateSolution();

    void initializeIndividualRandomly(Individual &individual);

    void selection();

    void cross();

    void mutation();

    void selectElite(vector<int> &elitePositions);

    void selectRabble(vector<int> &rabblePositions, const int numRabbles);

    void reEval();

    void replacement(vector<int> &elitePos);

    bool checkIfSame(const Chromosome &first, const Chromosome &second);

    void ox2Cross(const Chromosome &first, Chromosome &second);

    void mocCross(const int cut, const Chromosome &first, Chromosome &second);

private:
    const Airport *_airport;
    const string _crossOverType;
    const int _eliteSize;
    const int _popSize;
    const int _chromosomeSize;

    Population _initialPopulation;
    Population _currentPopulation;
    Chromosome _bestSolution;

    int _numEvaluations;
    int _totalEvaluations;
    int _tournamentSize;
    string _log;
    double _executionTime;
    int _bestCost;
};


#endif //INC_2_GENETIC_ALGORITHMS_GENETIC_H
