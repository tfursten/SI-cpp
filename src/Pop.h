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
#include <boost/foreach.hpp>


#include "xorshift64.h"
#include "rexp.h"
#include "disperse.h"
#include "Individual.h"
#include "disk.h"


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
    int m_nLethal; //number of abortions
    double m_dSigmaP;
    double m_dSigmaS;
	double m_pMut;
    double m_pDMut;
    double m_pSMut;
    double m_pDel;
    int m_nMutCount;
	xorshift64 m_myrand;
	Dispersal dist;
	Disk pdisk;
	Disk sdisk;
	std::ofstream & pout;
	std::ofstream & dout;
	std::vector<Individual> m_vPop1;
	std::vector<Individual> m_vPop2;
	std::vector<unsigned int> m_vWeights1;
	std::vector<unsigned int> m_vWeights2;

	void setMutCount();
	int disperse(int x, int y, double sigma);
	int sDisperseDist(int x, int y);
	int pDisperseDist(int x, int y);
	int sDisperseDisk(int x, int y);
	int pDisperseDisk(int x, int y);
	inline int periodic(int x, int y);
	inline int absorbing(int x, int y);
	void pollenDispersal(int dad);
	void seedDispersal(int mom);
	void samplePop(int gen);
	void mutate(gamete &g);
	void mutCountDec();

protected:
    int(Population::*pDisperse)(int,int);
    int(Population::*sDisperse)(int,int);
    int(Population::*newXY)(int,int);

public:
    Population(std::ofstream &p, std::ofstream &d): pout(p), dout(d) {};
    void initialize(int nMaxX, int nMaxY, std::string bound, int nPollen, int nOvule, int nMarkers, float dSigmaP, float dSigmaS,std::string si, std::string dist_name);
    void param(double dSMut, double dMMut, double dDMut, double dPdel, unsigned int seed);
	void evolve(int nGenerations, int nBurnIn, int nSample);
};

typedef std::vector<unsigned int> haplotype;
typedef std::vector<haplotype> genotype;
typedef std::pair<int,int> xycoord;

#endif // POP_H_INCLUDED
