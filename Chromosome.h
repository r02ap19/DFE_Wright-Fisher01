#pragma once

#include <map>
#include <set>
#include <stdio.h>
#include <iostream>
#include <fstream>

#include "parameters.h"

struct mutation {
	int homol; //0 = only on homologue 1; 1 = only on homologue 2; 2 = on both homologues (hence homozygote)
	double s; //selection coefficient
	double h; //dominance coefficient
	double pos; //position on continuous chromosome
};

//map containing all mutations: map<position, mutation> 
typedef std::map<double, mutation, std::less<double>> MapMuts;

class Chromosome {
public:
	Chromosome();
	~Chromosome();
	int nMut; //no. of deleterious mutations
	int Nho; //number of homozygote mutations

	double linkNeut[2]; //diploid neutral locus linked to deleterious mutations

	MapMuts mutations; //list of all deleterious mutations, their coefficients and whether they are homo or heterozygote 

	double addDelMutation(
		int, //homologue
		double, //position
		double, //s
		double //h
	);

	void add_neutral_mut(int, double, double, double); //analogous to addDelMutation 

	//Initialized neutral locus 
	void InitNeutral(
		double,
		double
	);
	void deleteChromo();

private:
};
