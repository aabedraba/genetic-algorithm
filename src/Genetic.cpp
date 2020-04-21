//
// Created by aabedraba on 6/11/19.
//

#include <algorithm>
#include <chrono>
#include "Genetic.h"

Genetic::Genetic(const Airport *airport, const string crossover, const unsigned int eliteSize,
                 const unsigned int numEvaluations, const unsigned int numTournaments) :
        _airport(airport),
        _crossOverType(crossover),
        _eliteSize(eliteSize),
        _numEvaluations(0),
        _popSize(50),
        _initialPopulation(_popSize),
        _currentPopulation(_popSize),
        _totalEvaluations(numEvaluations),
        _tournamentSize(numTournaments),
        _chromosomeSize(airport->getNumDoors()) {
    for (int ind = 0; ind < _popSize; ++ind) {
        _initialPopulation[ind].chromosome.resize(_chromosomeSize);
    }
    generateSolution();

}

Genetic::~Genetic() {

}

void Genetic::generateSolution() {
    auto start = std::chrono::steady_clock::now();
    for (int ind = 0; ind < _popSize; ++ind)
        initializeIndividualRandomly(_initialPopulation[ind]);
    int generation = 1;
    while (_numEvaluations < _totalEvaluations) {
        vector<int> elitePos;
        Utils::logPopulation(_initialPopulation, generation, _log);
        selectElite(elitePos);
        selection();
        cross();
        mutation();
        reEval();
        replacement(elitePos);
        _initialPopulation = _currentPopulation;
        generation++;
    }
    int posBestCost = 0;
    for (int i = 1; i < _popSize; ++i) {
        if (_initialPopulation[i].cost < _initialPopulation[posBestCost].cost)
            posBestCost = i;
    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<float> diff = end - start;
    string log = "________Best individual________\n\n";
    Utils::logBestSolution(&_initialPopulation[posBestCost].chromosome, _initialPopulation[posBestCost].cost, log);
    log += "\n";
    log += "Time: " + to_string(diff.count());
    log += "\n----------------------------------\n\n";
    _log = log + _log;

    _executionTime = diff.count();
    _bestSolution = _initialPopulation[posBestCost].chromosome;
    _bestCost = _initialPopulation[posBestCost].cost;
}

void Genetic::cross() {
    _log += "________Crosses________\n\n";
    int chance;
    for (int i = 0; i < _popSize; i += 2) {
        chance = rand() % 100;
        if (chance < 70) {
            _log += "Cross individual " + to_string(i) + "\n";
            Utils::logCross("Parent", &_currentPopulation[i].chromosome, &_currentPopulation[i + 1].chromosome, _log);
            _currentPopulation[i].evaluated = _currentPopulation[i + 1].evaluated = false;
            if (_crossOverType == "ox2") {
                Chromosome copy = _currentPopulation[i + 1].chromosome;
                ox2Cross(_currentPopulation[i].chromosome, _currentPopulation[i + 1].chromosome);
                ox2Cross(copy, _currentPopulation[i].chromosome);
            } else if (_crossOverType == "moc") {
                int cut = _chromosomeSize - 2;
                Chromosome copy = _currentPopulation[i + 1].chromosome;
                mocCross(cut, _currentPopulation[i].chromosome, _currentPopulation[i + 1].chromosome);
                mocCross(cut, copy, _currentPopulation[i].chromosome);
            }
            Utils::logCross("Children", &_currentPopulation[i].chromosome, &_currentPopulation[i + 1].chromosome, _log);
            _log += "\n\n";
        }
    }
}

void Genetic::selection() {
    for (int i = 0; i < _popSize; ++i) {
        vector<int> rands;
        while (rands.size() < _tournamentSize) {
            int toInsert = rand() % _popSize;
            auto iter = find(rands.begin(), rands.end(), toInsert);
            if (iter == rands.end()) {
                rands.push_back(toInsert);
            }
        }
        auto iter = rands.begin();
        int best = *iter;
        while (iter != rands.end()) {
            if (_initialPopulation[*iter].cost > _initialPopulation[best].cost)
                best = *iter;
            iter++;
        }
        _currentPopulation[i] = _initialPopulation[best];
    }

}

void Genetic::initializeIndividualRandomly(Individual &individual) {
    for (int i = 0; i < individual.chromosome.size(); ++i)
        individual.chromosome[i] = i;

    for (int i = 0; i < individual.chromosome.size(); ++i) {
        int randInRange = rand() % individual.chromosome.size();
        swap(individual.chromosome[i], individual.chromosome[randInRange]);
    }
    individual.cost = Utils::solutionCost(individual.chromosome, _airport->getFluxMatrix(),
                                          _airport->getDistanceMatrix(),
                                          _airport->isSimetric());
    individual.evaluated = true;
}

void Genetic::ox2Cross(const Chromosome &first, Chromosome &second) {
    vector<int> firstSelectedGenes;
    vector<int> posSecondGenes;
    int chance;
    for (int i = 0; i < _chromosomeSize; ++i) {
        chance = rand() % 100;
        if (chance > 50)
            firstSelectedGenes.push_back(first[i]);
    }
    for (int i = 0; i < firstSelectedGenes.size(); ++i) {
        for (int j = 0; j < second.size(); ++j) {
            if (firstSelectedGenes[i] == second[j])
                posSecondGenes.push_back(j);
        }
    }
    sort(posSecondGenes.begin(), posSecondGenes.end());
    for (int i = 0; i < posSecondGenes.size(); ++i) {
        second[posSecondGenes[i]] = firstSelectedGenes[i];
    }
}

void Genetic::mocCross(const int cut, const Chromosome &first, Chromosome &second) {
    vector<int> firstSelectedGenes;
    vector<int> posSecondGenes;
    for (int i = cut; i < first.size(); ++i) {
        firstSelectedGenes.push_back(first[i]);
    }
    for (int i = 0; i < firstSelectedGenes.size(); ++i) {
        for (int j = 0; j < second.size(); ++j) {
            if (firstSelectedGenes[i] == second[j])
                posSecondGenes.push_back(j);
        }
    }
    sort(posSecondGenes.begin(), posSecondGenes.end());
    for (int i = 0; i < posSecondGenes.size(); ++i) {
        second[posSecondGenes[i]] = firstSelectedGenes[i];
    }
}

void Genetic::mutation() {
    _log += "________Mutation________\n\n";
    int chance, rand1, rand2;
    vector<int> posToChange;
    for (int i = 0; i < _popSize; ++i) {
        for (int j = 0; j < _chromosomeSize; ++j) {
            chance = rand() % 100;
            if (chance >= 95) {
                _log += "Mutation individual " + to_string(i) + "\n";
                _currentPopulation[i].evaluated = false;
                rand1 = rand2 = j;
                while (rand1 == rand2 && rand2 == j && rand1 == j) {
                    rand1 = rand() % _chromosomeSize;
                    rand2 = rand() % _chromosomeSize;
                }
                posToChange.insert(posToChange.end(), {j, rand1, rand2});
                Utils::logVector("Before", &_currentPopulation[i].chromosome, _log);
                sort(posToChange.begin(), posToChange.end());
                swap(_currentPopulation[i].chromosome[posToChange[0]],
                     _currentPopulation[i].chromosome[posToChange[1]]);
                swap(_currentPopulation[i].chromosome[posToChange[1]],
                     _currentPopulation[i].chromosome[posToChange[2]]);
                Utils::logVector("After", &_currentPopulation[i].chromosome, _log);
                _log += "\n\n";
                posToChange.clear();
            }
        }
    }
}


void Genetic::reEval() {
    for (int i = 0; i < _popSize; ++i) {
        if (_currentPopulation[i].evaluated == false) {
            _currentPopulation[i].cost = Utils::solutionCost(_currentPopulation[i].chromosome,
                                                             _airport->getFluxMatrix(),
                                                             _airport->getDistanceMatrix(),
                                                             _airport->isSimetric());
            _currentPopulation[i].evaluated = true;
            _numEvaluations++;
        }
    }
}

//TODO duplicate with elite.... probably implement something for both?
void Genetic::selectRabble(vector<int> &rabblePositions, const int numRabbles) {
    vector<pair<int, int>> rabbleCostAndPos;
    for (int i = 0; i < numRabbles; ++i) {
        rabbleCostAndPos.push_back(make_pair(_currentPopulation[i].cost, i));
    }
    int posBest = Utils::posMinElement(rabbleCostAndPos);
    for (int i = numRabbles; i < _popSize; ++i) {
        if (_currentPopulation[i].cost > rabbleCostAndPos[posBest].first) {
            rabbleCostAndPos[posBest] = make_pair(_currentPopulation[i].cost, i);
            posBest = Utils::posMinElement(rabbleCostAndPos);
        }
    }
    for (int i = 0; i < numRabbles; ++i) {
        rabblePositions.push_back(rabbleCostAndPos[i].second);
    }
}

void Genetic::selectElite(vector<int> &elitePositions) {
    vector<pair<int, int>> eliteCostAndPos;
    for (int i = 0; i < _eliteSize; ++i) {
        eliteCostAndPos.push_back(make_pair(_initialPopulation[i].cost, i));
    }
    int posWorse = Utils::posMaxElement(eliteCostAndPos);
    for (int i = _eliteSize; i < _popSize; ++i) {
        if (_initialPopulation[i].cost < eliteCostAndPos[posWorse].first) {
            bool alreadyIn = false;
            for (int j = 0; j < _eliteSize; ++j) {
                alreadyIn = checkIfSame(_initialPopulation[eliteCostAndPos[j].second].chromosome,
                                        _initialPopulation[i].chromosome);
                break;
            }
            if (alreadyIn) continue;
            eliteCostAndPos[posWorse] = make_pair(_initialPopulation[i].cost, i);
            posWorse = Utils::posMaxElement(eliteCostAndPos);
        }
    }
    for (int i = 0; i < _eliteSize; ++i) {
        elitePositions.push_back(eliteCostAndPos[i].second);
    }
}

void Genetic::replacement(vector<int> &elitePos) {
    _log += "________Elites and rabbles________\n\n";
    vector<int>::iterator iter;
    for (int i = 0; i < _popSize; ++i) {
        iter = elitePos.begin();
        while (iter != elitePos.end()) {
            if (_currentPopulation[i].cost == _initialPopulation[*iter].cost) {
                bool isSame = checkIfSame(_currentPopulation[i].chromosome,
                                          _initialPopulation[*iter].chromosome);
                if (isSame) {
                    iter = elitePos.erase(iter);
                } else
                    iter++;
            } else {
                iter++;
            }
        }
    }
    if (elitePos.empty() == true)
        return;
    vector<int> rabblePos;
    selectRabble(rabblePos, elitePos.size());
    Utils::logElitesAndRabbles(&_initialPopulation, &_currentPopulation, elitePos, rabblePos, _log);
    for (int i = 0; i < elitePos.size(); ++i) {
        _currentPopulation[rabblePos[i]] = _initialPopulation[elitePos[i]];
    }
    _log += "\n\n";
}

bool Genetic::checkIfSame(const Chromosome &first, const Chromosome &second) {
    for (int i = 0; i < first.size(); ++i) {
        if (first[i] != second[i])
            return false;
    }
    return true;
}

const string &Genetic::getLog() const {
    return _log;
}

const Airport *Genetic::getAirport() const {
    return _airport;
}

const Chromosome &Genetic::getBestSolution() const {
    return _bestSolution;
}

double Genetic::getExecutionTime() const {
    return _executionTime;
}

int Genetic::getBestCost() const {
    return _bestCost;
}



