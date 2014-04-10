#ifndef POP_H_INCLUDED
#define POP_H_INCLUDED

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <unistd.h>
#include <map>

#include "xorshift64.h"
#include "rexp.h"
#include "disperse.h"
#include "Individual.h"


#define foreach BOOST_FOREACH




class Population
{
private:
    int m_nMaxX;
    int m_nMaxY;
    int m_nPollen;
    int m_nOvule;
    int m_nMarkers;
    int m_nIndividuals;
    int m_nSalleles;
    int m_nAlleles;
    double m_dSigmaP;
    double m_dSigmaS;
	double m_pMut;
    double m_pDMut;
    double m_pSMut;
    int m_nMutCount;
	int m_nSample;
	xorshift64 m_myrand;
	Dispersal dist;
	std::ofstream & mout;
	std::vector<Individual> m_vPop1;
	std::vector<Individual> m_vPop2;
	std::vector<unsigned int> m_vWeights1;
	std::vector<unsigned int> m_vWeights2;
	std::vector<int> m_vtransIndex;
	std::vector<int> m_vtransDist;
	std::vector<int> m_DistCount;
	std::vector<double> m_vAvgIBD;
	double m_fAvgSig;
	void setMutCount();
	int disperse(int x, int y, double sigma);
	void pollenDispersal(int dad);
	void seedDispersal(int mom);
	void samplePop(int gen);
	void mutate(gamete &g);
	void mutCountDec();

public:
    Population(std::ofstream &o): mout(o) {};
    void initialize(int nMaxX, int nMaxY, int nPollen, int nOvule, int nMarkers, std::string si, std::string dist_name);
    void param(float dSigmaP, float dSigmaS, double dSMut, double dMMut, double dDMut, unsigned int seed);
	void evolve(int nGenerations, int nBurnIn, int nSample);
};

typedef std::vector<int> haplotype;
typedef std::vector<haplotype> genotype;

#endif // POP_H_INCLUDED
