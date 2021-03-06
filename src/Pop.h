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
    std::string m_sBound; 
    int m_nPollen;
    int m_nOvule;
    int m_nMarkers;
    int m_nDel;
    int m_nLoci;
    int m_nIndividuals;
    int m_nSalleles;
    int m_nAlleles;
    int m_nFecundity;
    int m_nLethal; //number of abortions
    int m_nMutAlleles;
    double m_dSigmaP;
    double m_dSigmaS;
    double m_pMut;
    double m_pDMut;
    double m_pMMut;
    double m_pDel;
    int m_nMutCount;
	xorshift64 m_myrand;
	Dispersal pDisp;
    Dispersal sDisp;
	std::ofstream & pout;
	std::ofstream & dout;
	// TODO: See if using Individual **m_vPop1 will give better results
	std::vector<Individual> m_vPop1;
	std::vector<Individual> m_vPop2;
    std::vector<int> m_nPopSample;

	void setMutCount();
	void pollenDispersal(int dad);
	void seedDispersal(int mom);
	void samplePop(int gen);
	void mutate(Individual &I);
    void mutCountDec();
    int mutate_iam();
    int mutate_kam();

protected:
    int(Population::*mutated_allele)();


public:
    Population(std::ofstream &p, std::ofstream &d): pout(p), dout(d) {};
    void initialize(int nMaxX, int nMaxY, std::string bound, int nPollen, int nOvule, int nMarkers, int nDel, float dSigmaP, float dSigmaS,std::string si, std::string dist_name, float pp, float sp, bool fast, std::string mut_type, int nAlleles);
    void param(double dSMut, double dMMut, double dDMut, double dPdel, unsigned int seed);
    void evolve(int nGenerations, int nBurnIn, int nSample);
    
    ~Population()
      {

      }
    
};

//typedef std::vector<unsigned int> haplotype;
//typedef std::vector<haplotype> genotype;
typedef std::pair<int,int> xycoord;

#endif // POP_H_INCLUDED
