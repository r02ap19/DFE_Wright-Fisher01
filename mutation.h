#pragma once

#include <map>

struct mutation {
	int chrom; //chromosome
	double pos; //position on continuous chromosome
	double s; //selection coefficient
	double h; //dominance coefficient
};

//continuous chromosome containing deleterious mutations: map<position, mutation>
typedef std::map<double, mutation, std::less<double>> MapChrom;