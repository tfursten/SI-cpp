#ifndef INDIVIDUAL_H_INCLUDED
#define INDIVIDUAL_H_INCLUDED

#include <vector>
#include <boost/algorithm/string/predicate.hpp>
#include "disperse.h"



typedef std::vector<unsigned int> haplotype;
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
    std::pair<unsigned int,unsigned int> s;
    gamete(){}
    gamete(int m)
    {
        h.assign((m+1),-1);
        par.assign((m+1),-1);
        gpar.assign((m+1),-1);
        del = 0;
        s = std::make_pair(0,0);
    }
};
void printGenotype(haplotype &h);

void initGenes(int m, genes &g, haplotype &dad, haplotype &mom);


class Individual
{

private:
    static std::string name;
    static dominance dRank;
    unsigned int m_nPosition;
    unsigned int m_nWeight;
    std::vector<gamete> m_vOvules;
    std::vector<unsigned int> m_vOWeights;
    genes m_genes;

protected:
    typedef bool(Individual::*fptr)(gamete&);
public:
    static fptr op;

    Individual(unsigned int pos, int nOvules, int nMarkers)
    {
        m_nPosition = pos;
        m_nWeight = 0;
        m_vOWeights.assign(nOvules,0);
        for(int i=0; i<nOvules; i++)
            m_vOvules.emplace_back(nMarkers);
        haplotype dad, mom;
        dad.assign((nMarkers+1),-1);
        mom.assign((nMarkers+1),-1);
        initGenes(nMarkers,dad,mom);
    }
    Individual(int pos, int nOvules, int nMarkers, haplotype& dad, haplotype& mom)
    {
        m_nPosition = pos;
        m_nWeight = 1;
        m_vOWeights.assign(nOvules, 0);
        for(int i=0; i<nOvules; i++)
            m_vOvules.emplace_back(nMarkers);
        initGenes(nMarkers,dad,mom);
    }
    void newIndividual(gamete& dad, gamete& mom, unsigned int weight)
    {
        m_nWeight = weight;
        m_genes.g[0].assign(dad.h.begin(),dad.h.end());
        m_genes.g[1].assign(mom.h.begin(),mom.h.end());
        m_genes.par[0].assign(dad.par.begin(),dad.par.end());
        m_genes.par[1].assign(mom.par.begin(),mom.par.end());
        m_genes.gpar.at(0) = dad.gpar;
        m_genes.gpar.at(1) = mom.gpar;
        m_genes.del = std::make_pair(dad.del,mom.del);
    }

    static void initDomRank(xorshift64& rand, int nAlleles)
    {
        for(int iii = 0; iii<nAlleles; iii++)
        {
            dRank.push_back(rand.get_uint64());
        }
    }

    static unsigned int getDomRank(int n)
    {
        return dRank.at(n);
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
    inline double operator()(gamete &dad) {
        return (this->*op)(dad);
    }

    inline int position() {return m_nPosition;}
    inline int dadGene(int n) {return m_genes.g[0][n];}
    inline int momGene(int n) {return m_genes.g[1][n];}
    inline int dadDel(){return m_genes.del.first;}
    inline int momDel(){return m_genes.del.second;}
    inline int dadGpar(int n) {return m_genes.gpar[0][n];}
    inline int momGpar(int n) {return m_genes.gpar[1][n];}
    inline int dadPos() {return m_genes.gpar[0][0];}
    inline int momPos() {return m_genes.gpar[1][0];}
    inline unsigned int weight() {return m_nWeight;}
    inline unsigned int ovuleWeight(int n) {return m_vOWeights[n];}
    inline gamete ovule(int n) {return m_vOvules[n];}
    void setWeight(unsigned int weight) {m_nWeight = weight;}
    void setOvuleWeight(unsigned int weight, int n) {m_vOWeights[n]=weight;}
    void setOvuleGenes(gamete &pollen, int n) {m_vOvules[n]= pollen;}
    gamete makeGamete(xorshift64& rand, int nMarkers);
    void initGenes(int m, haplotype &dad, haplotype &mom);
    void clearOvuleWeight(int nOvules);
    bool nsi(gamete &dad);
    bool psi(gamete &dad);
    bool gsi(gamete &dad);
    bool ssi(gamete &dad);
    bool bsi(gamete &dad);

};


#endif // INDIVIDUAL_H_INCLUDED
