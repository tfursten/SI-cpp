#ifndef INDIVIDUAL_H_INCLUDED
#define INDIVIDUAL_H_INCLUDED

#include <cstring>
#include <vector>
#include <boost/algorithm/string/predicate.hpp>
#include "disperse.h"


#define HAP_SIZE(m) (m+1)*sizeof(haplotype) 
#define NULL_GENE -1


typedef unsigned int haplotype;
typedef struct{
  haplotype *h1;
  haplotype *h2;
} genotype;

//typedef std::vector<haplotype*> genotype;
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
  haplotype *h;
  haplotype *par;
  haplotype *gpar;
  int del;
  int nLoci;
  std::pair<unsigned int,unsigned int> s;

  gamete(int m)
  {
    nLoci = m;
    h = new haplotype[m+1];
    par = new haplotype[m+1];
    gpar = new haplotype[m+1];
    std::fill_n(h, m+1, NULL_GENE);
    std::fill_n(par, m+1, NULL_GENE);
    std::fill_n(gpar, m+1, NULL_GENE);
    //memset(h, NULL_GENE, HAP_SIZE(m)); //TODO: Is this neccessary???
    //memset(par, NULL_GENE, HAP_SIZE(m));
    //memset(gpar, NULL_GENE, HAP_SIZE(m));
    del = 0;
    s = std::make_pair(0,0);
  }

  ~gamete()
  {
    delete[] h;
    delete[] par;
    delete[] gpar;
  }
};

void printGenotype(haplotype &h);

//void initGenes(int m, genes &g, haplotype &dad, haplotype &mom);



class Individual
{
  
 private:
  static std::string name;
  static dominance dRank;
  unsigned int m_nPosition;
  unsigned int m_nWeight;
  unsigned int m_nOvules;
  
  gamete **m_vOvules;
  //std::vector<gamete*> m_vOvules; // changed to pointer to prevent object copying.
  std::vector<unsigned int> m_vOWeights; // TODO: make this a int[] to run faster.
  genes m_genes;
  int m_nMarkers;
  int maxX;
  int maxY;
  int x_coord; // current position (x, y) on grid.
  int y_coord;
  
  void initVars(unsigned int pos, int nOvules, int nMarkers, unsigned int nWeight, unsigned int maxX, unsigned int maxY);
  gamete *makeGameteHelper(int, int, haplotype *, haplotype *, haplotype *, haplotype *);

 protected:
  typedef bool(Individual::*fptr)(gamete&);

 public:
  static fptr op;

  
  Individual(const Individual &i)
    {
      haplotype *dad = new haplotype[i.m_nMarkers+1];
      haplotype *mom = new haplotype[i.m_nMarkers+1];
      memcpy(dad, i.m_genes.g.h1, HAP_SIZE(i.m_nMarkers));
      memcpy(mom, i.m_genes.g.h2, HAP_SIZE(i.m_nMarkers));
   
      initVars(i.m_nPosition, i.m_nOvules, i.m_nMarkers,i.m_nWeight, i.maxX, i.maxY);
      initGenes(m_nMarkers, dad, mom);

    }

  Individual &operator=(const Individual &i) = default;
  

  Individual(unsigned int pos, int nOvules, int nMarkers, unsigned int maxX, unsigned int maxY)
    { 
      // create new set of known haplotypes
      haplotype *dad = new haplotype[nMarkers+1];
      haplotype *mom = new haplotype[nMarkers+1];
      memset(dad, NULL_GENE, HAP_SIZE(nMarkers));
      memset(mom, NULL_GENE, HAP_SIZE(nMarkers));
      std::fill_n(dad, nMarkers+1, NULL_GENE);
      std::fill_n(mom, nMarkers+1, NULL_GENE);

      initVars(pos, nOvules, nMarkers, 0, maxX, maxY);
      initGenes(nMarkers, dad, mom);
    }

  Individual(unsigned int pos, int nOvules, int nMarkers, unsigned int maxX, unsigned int maxY, haplotype *dad, haplotype *mom)
    {
      initVars(pos, nOvules, nMarkers, 1, maxX, maxY);
      initGenes(nMarkers, dad, mom);
    }

    /**
     * Updates the current individual's genotype, as well as his parent's and
     * grandparent's, based on mom & dad's gametes.
     */
    void newIndividual(gamete *dad, gamete *mom, unsigned int weight)
    {
      m_nWeight = weight;

      /*
      // Create a new genotype
      memcpy(m_genes.g.h1, &dad.h, HAP_SIZE(m_nMarkers));
      memcpy(m_genes.g.h2, &mom.h, HAP_SIZE(m_nMarkers));

      // Set parent's and grandparents genotype
      memcpy(m_genes.par.h1, &dad.par, HAP_SIZE(m_nMarkers));
      memcpy(m_genes.par.h2, &mom.par, HAP_SIZE(m_nMarkers));
      memcpy(m_genes.gpar.h1, &dad.gpar, HAP_SIZE(m_nMarkers));
      memcpy(m_genes.gpar.h2, &mom.gpar, HAP_SIZE(m_nMarkers));

      m_genes.del = std::make_pair(dad.del, mom.del);
      */

      
      // Create a new genotype
      memcpy(m_genes.g.h1, dad->h, HAP_SIZE(m_nMarkers));
      memcpy(m_genes.g.h2, mom->h, HAP_SIZE(m_nMarkers));

      // Set parent's and grandparents genotype
      memcpy(m_genes.par.h1, dad->par, HAP_SIZE(m_nMarkers));
      memcpy(m_genes.par.h2, mom->par, HAP_SIZE(m_nMarkers));
      memcpy(m_genes.gpar.h1, dad->gpar, HAP_SIZE(m_nMarkers));
      memcpy(m_genes.gpar.h2, mom->gpar, HAP_SIZE(m_nMarkers));

      m_genes.del = std::make_pair(dad->del, mom->del);
      
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

    ~Individual();

    static inline std::string getName() {return name;}
    inline double operator()(gamete &dad) {
        return (this->*op)(dad);
    }

    inline int coordX(){return x_coord;};
    inline int coordY(){return y_coord;};
    inline int position() {return m_nPosition;}
    inline int dadGene(int n) {return m_genes.g.h1[n];}
    inline int momGene(int n) {return m_genes.g.h2[n];}
    inline int dadDel(){return m_genes.del.first;}
    inline int momDel(){return m_genes.del.second;}
    inline int dadGpar(int n) {return m_genes.gpar.h1[n];}
    inline int momGpar(int n) {return m_genes.gpar.h2[n];}
    inline int dadPos() {return m_genes.par.h1[0];}
    inline int momPos() {return m_genes.par.h2[0];}
    inline unsigned int weight() {return m_nWeight;}
    inline unsigned int ovuleWeight(int n) 
    {return m_vOWeights[n];}
    unsigned int getMinOvuleWeight(){return *min_element(m_vOWeights.begin(),m_vOWeights.end());}
    int getMinOvuleID(){
      auto it = std::min_element(m_vOWeights.begin(),m_vOWeights.end());
      return it-m_vOWeights.begin();}  
    inline gamete *ovule(int n) {return m_vOvules[n];}
    void setDadGene(int m, int n) {m_genes.g.h1[m] = n;}
    void setMomGene(int m, int n) {m_genes.g.h2[m] = n;}
    void setDadDel(int n) {m_genes.del.first = n;}
    void setMomDel(int n) {m_genes.del.second = n;}
    void setCoordinates(int position, int nMaxX, int nMaxY);
    void setWeight(unsigned int weight) {m_nWeight = weight;}
    void setOvuleWeight(unsigned int weight, int n) {m_vOWeights[n]=weight;}
    void setOvuleGenes(gamete *pollen, int n) 
    {
      gamete *tmp = m_vOvules[n];
      m_vOvules[n] = pollen;
      delete tmp;
    }
    gamete *makeGamete(xorshift64& rand, int nMarkers);
    void initGenes(int m, haplotype *dad, haplotype *mom);
    void clearOvuleWeight(int nOvules);
    bool nsi(gamete &dad);
    bool psi(gamete &dad);
    bool gsi(gamete &dad);
    bool ssi(gamete &dad);
    bool bsi(gamete &dad);

};


#endif // INDIVIDUAL_H_INCLUDED
