#ifndef INDIVIDUAL_H_INCLUDED
#define INDIVIDUAL_H_INCLUDED

#include <vector>
#include <boost/algorithm/string/predicate.hpp>
#include "disperse.h"



typedef std::vector<int> haplotype;
typedef std::vector<haplotype> genotype;
typedef std::vector<unsigned int> dominance;

struct genes
{
    genotype g;
    genotype par;
    genotype gpar;
    std::pair<int,int> del;
};

struct gamete
{
    haplotype h;
    haplotype par;
    haplotype gpar;
    int del;
};


class Individual
{

private:
    static std::string name;
    static dominance dRank;
    int m_nPosition;
    unsigned int m_nWeight;
    std::vector<gamete> m_vOvules;
    std::vector<unsigned int> m_vOWeights;
    genes m_genes;
    gamete m_gamete;
    genotype id;


protected:
    typedef bool(Individual::*fptr)(Individual);
public:
    static fptr op;
    Individual(int pos, int nOvules, int nMarkers)
    {
        m_nPosition = pos;
        m_nWeight = 0;
        m_vOWeights.assign(nOvules,0);
        haplotype p0;
        p0.assign(nMarkers,2*pos);
        id.push_back(p0);
        haplotype p1;
        p1.assign(nMarkers,2*pos+1);
        id.push_back(p1);
    }
    Individual(int pos, int nOvules, haplotype& h1, haplotype& h2)
    {
        m_nPosition = pos;
        m_nWeight = 1;
        m_vOWeights.assign(nOvules, 0);
        m_genes.g.push_back(h1);
        m_genes.g.push_back(h2);
        haplotype p0;
        p0.assign(sizeof(h1),2*pos);
        m_genes.par.push_back(p0);
        m_genes.gpar.push_back(p0);
        id.push_back(p0);
        haplotype p1;
        p1.assign(sizeof(h2),2*pos+1);
        m_genes.par.push_back(p1);
        m_genes.gpar.push_back(p1);
        id.push_back(p1);
        m_genes.del = std::make_pair(0,0);

    }
    void newIndividual(gamete& p1, gamete& p2)
    {
        m_genes.g = {p1.h,p2.h};
        m_genes.par = {p1.par,p2.par};
        m_genes.gpar = {p1.gpar,p2.gpar};
        m_genes.del = std::make_pair(p1.del,p2.del);
    }

    static void initDomRank(unsigned int rand, int nAlleles)
    {
        dRank.assign(rand,nAlleles);
    }

    static void addDomRank(unsigned int rand)
    {
        dRank.push_back(rand);
    }

    template<class A>
    static bool initialize(A &si) {
        static const char name_keys[][8]={
            "nsi","psi","gsi","ssi","bsi"
        };
        static const fptr si_ops[] = {
            &Individual::nsi,
            &Individual::psi,
            &Individual::gsi,
            &Individual::ssi,
            &Individual::bsi
        };
        int pos = key_switch(si,name_keys);
        if(pos == -1) {
            std::cerr << "ERROR: Invalid Self-Incompatibility System" << std::endl;
            return false;
        }
        name = std::string(name_keys[pos]);
        op = si_ops[pos];
        return true;
    }

    static inline std::string getName() {return name;}
    inline double operator()(Individual dad) {
        return (this->*op)(dad);
    }

    inline int position() {return m_nPosition;}
    inline genotype genes() {return m_genes.g;}
    inline unsigned int weight() {return m_nWeight;}
    inline gamete getGamete() {return m_gamete;}
    inline int ovuleWeight(int n) {return m_vOWeights[n];}
    inline gamete ovule(int n) {return m_vOvules[0];}
    void setWeight(unsigned int weight) {m_nWeight = weight;}
    void setOvuleWeight(unsigned int weight, int n) {m_vOWeights[n]=weight;}
    void setOvuleGenes(gamete &pollen, int n) {m_vOvules[n]= pollen;}
    void makeGamete(xorshift64& rand, int nMarkers);
    bool nsi(Individual dad);
    bool psi(Individual dad);
    bool gsi(Individual dad);
    bool ssi(Individual dad);
    bool bsi(Individual dad);

};


#endif // INDIVIDUAL_H_INCLUDED
