#pragma once

#include <stdio.h>
#include <fstream>
#include <math.h>
#include <random>
#include <vector>
#include <set>
#include <iostream>

#include "Chromosome.h"

using namespace std;

class Individuals {
public:
	Individuals(Parameters);
	~Individuals();
	bool alive;
	bool reproduce;

	double w; //individual fitness

	Chromosome chromo;
	Chromosome neutral_genome;

	void initialise(double kk, Parameters para, std::gamma_distribution<> finds, std::uniform_real_distribution<>, std::normal_distribution<>);
	void benef_mutation(int, double, std::uniform_real_distribution<>);
	void back_mutation(int);
	void delet_mutation(int, double, std::uniform_real_distribution<>, std::gamma_distribution<>, std::uniform_real_distribution<>, Parameters para);
	void delet_mutation(int, std::uniform_real_distribution<>, double, double);
	void deleteInd(void);
	void SFS_sample(int g, int gg, std::ofstream* out);
	void neutral_genome_mut(int, std::uniform_real_distribution<>, std::uniform_real_distribution<>, Parameters);
	void neutralSFS_sample(int g, int gg, std::ofstream* out);

private:
};