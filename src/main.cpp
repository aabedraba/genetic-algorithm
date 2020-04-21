#include <iostream>
#include <string>
#include <vector>
#include <chrono>

#include "Utils.h"
#include "Airport.h"
#include "Genetic.h"

enum value {
    AIRPORT = 0,
    SEED = 1,
    CROSSOVER = 2,
    ELITE_SIZE = 3,
    NUM_EVALUATIONS = 4,
    NUM_TORNEO = 5
};

void readParametersFile(vector<string> &parameters) {
    ifstream file;
    file.open("../_parametros.txt");
    string linea;
    if (file.good()) {
        while (!file.eof()) {
            file.ignore(256, ' ');
            file.ignore(256, ' ');
            file >> linea;
            parameters.push_back(linea.c_str());
        }
    }
    file.close();
}


int main() {
    vector<string> parameters;
    readParametersFile(parameters);
    vector<string> airports = {"madrid01.dat","madrid02.dat","madrid03.dat","madrid04.dat",
                               "malaga01.dat","malaga02.dat","malaga03.dat","malaga04.dat"};
    vector<int> seeds = {23436383   };
    string log;
    for (int j = 0; j < seeds.size(); ++j) {
        for (int i = 0; i < airports.size(); ++i) {
            string file = "../_data/" + airports[i];
            srand(seeds[j]);
            Airport *airport = new Airport(file);
            int eliteSize = stoi(parameters[ELITE_SIZE], nullptr, 10);
            int numEvaluation = stoi(parameters[NUM_EVALUATIONS], nullptr, 10);
            int numTournaments = stoi(parameters[NUM_TORNEO], nullptr, 10);
            Genetic *genetic = new Genetic(airport, parameters[CROSSOVER], eliteSize, numEvaluation, numTournaments );
            string type = "Genetic-" + parameters[CROSSOVER] + "-" + parameters[ELITE_SIZE] + "-elites-" + to_string(seeds[j]);
            Utils::getResults(type, airport->getAirportName(), genetic->getBestCost(), genetic->getExecutionTime(), log);
            //Utils::writeInFile(seeds[j], type, airport->getAirportName(), genetic->getLog());
            delete airport;
            delete genetic;
        }
    }

    ofstream fs;
    fs.open("../_logs/results.txt");
    if (fs.is_open()) {
        fs << log << endl;
    } else
        throw "File not properly opened!";




    return 0;
}