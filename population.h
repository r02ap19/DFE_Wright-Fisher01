#pragma once

#include <stdio.h>
#include <vector>
#include <iostream>

#include "individuals.h"
#include "parameters.h"

using namespace std;

struct pop_muts {

	int count; //mutation occurrence to calculate frequency when output
	double s; //selection coefficient
	double h; //dominance coefficient
};

//map containing all mutations in the population: map<position, pop_muts> 
typedef std::map<double, pop_muts, std::less<double>> MapPopMuts;

class Population {
public:
	Population(int, int);
	~Population();
	int N, Noffs;
	int x, y;
	int DelMut; //nr. of deleterious mutations
	double meanDelMut;
	double DelMutSd; //standard deviation
	double W, Wsd; //mean population viability (calculated after survival) and standard deviation
	double Homoz, Hsd; //mean neutral homozygosity and standard deviation

	double Vneutral, VneutralSD; //neutral variance, sum of squares

	double Wmax; // maximum individual fitness in any given generation. Used as upper bound for sampling rep prob when loadEffect = 1 to drastically increase speed.

	vector<Individuals> inds;
	vector<Individuals> tmp_inds;

	MapPopMuts popMuts;

	void initialise_pop(double, double, Parameters, std::gamma_distribution<>, std::uniform_real_distribution<>, std::normal_distribution<>);
	void deleteAdults(void);
	void computeSums(Parameters, Individuals);
	void computeStats(Parameters);
	void set2zero(void);
	void outPop(int, int, std::ofstream*);
	void addMutation(
		double, //position
		double, //s
		double //h
	);
	void outMutations(int, int, int, std::ofstream*, Parameters);
};