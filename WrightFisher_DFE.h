#define CLUSTER 0

#include <stdio.h>
#include <stdlib.h>
#if CLUSTER 
#include <unistd.h>
#else
#include <tchar.h> 
#include <direct.h>
#include <io.h>
#endif
#include <vector>
#include <set>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <math.h>
#include <numeric>
#include <time.h>
#include <random>
#include <iterator>

#include "parameters.h"
#include "population.h"

using namespace std;

clock_t extime;

Parameters para;
string dir, dirOut;
int r, g; // counters for replicates and generations
ofstream par, pops, popmut, SFSsample, neutralSFSsample;
int count;

Population*** pop;

std::random_device rd;
std::mt19937 rdgen(rd());

std::uniform_int_distribution<> uniint(0, 1);
std::uniform_real_distribution<> unireal(0.0, 1.0);
std::bernoulli_distribution Bern(0.5);


//Add neutral locus linked to deleterious mutations (Keightley & Otto 2006)
//Normal distribution with Vm = 1
std::normal_distribution<> linked_neutral(0.0, 1.0);

//Distributions for sampling deleterious mutations
//sample no. of crossovers
std::poisson_distribution<> crossn(para.R);
//Deleterious mutation position (and crossovers positions)
std::uniform_real_distribution<> position(0.0, para.R);
//Selection coefficient for mildly deleterious mutations
std::gamma_distribution<> s_mild(1.0, double(para.mean_sd*2*para.assumed_Ne)); //(shape, mean / shape)
//Dominance coefficient for mildly deleterious mutations
double k_h = -log(2.0 * para.mean_hd) / para.mean_sd;

//Distributions for sh_dist
std::uniform_real_distribution<> uniform(0.0, 1.0);

// Functions declaration
const string Int2Str(const int x);
void RunModel(void);
void initialisation();
void reproduction_0(void); //offspring selection
void reproduction_1(void); //fertility selection
void inheritance_0(Individuals*, Individuals, Individuals); //centrometic neutral linked locus
void inheritance_1(Individuals*, Individuals, Individuals); //telomeric neutral linked locus
void neutral_inheritance(Individuals*, Individuals, Individuals); //inheritance of the neutral, unlinked genome
void outPop_header(void);
void outPopMut_header(void);
void outSFS_sample_header(void);
void housekeeping(void);
void outNeutral_SFS_sample_header(void);