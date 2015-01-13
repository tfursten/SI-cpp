#include <iostream>
#include "Individual.h"
using namespace std;



void printGenotype(haplotype *h, int hapSize)
{
  for(int i = 0; i < hapSize; i++)
    cout << h[i] << ",";
  cout << endl;
}

bool Individual::nsi(gamete &dad)
{
  return true;
}

bool Individual::psi(gamete &dad)
{
  if (m_nPosition == dad.par[0])
    return false;
  return true;
}

bool Individual::gsi(gamete &dad)
{
  if (m_genes.g.h1[0] == dad.h[0])
    return false;
  if (m_genes.g.h2[0] == dad.h[0])
    return false;
  return true;
}

bool Individual::ssi(gamete &dad)
{
  if (m_genes.g.h1[0] == dad.s.first)
    return false;
  if (m_genes.g.h2[0] == dad.s.first)
    return false;
  if (m_genes.g.h1[0] == dad.s.second)
    return false;
  if (m_genes.g.h2[0] == dad.s.second)
    return false;
  return true;

}

bool Individual::bsi(gamete &dad)
{
  unsigned int d1 = getDomRank(dad.s.first);
  unsigned int d2 = getDomRank(dad.s.second);
  unsigned int dominant;
  if (d1 >= d2)
    dominant = dad.s.first;
  else 
    dominant = dad.s.second;
  if (m_genes.g.h1[0] == dominant)
    return false;
  if (m_genes.g.h2[0] == dominant)
    return false;
  return true;
}

gamete *Individual::makeGameteHelper(int p, int m, haplotype *g0, haplotype *g1, haplotype *par0, haplotype *par1)
{
  // Helper class for makeGamete(). Replaced vector::assign() with cstring
  // memcpy (on array) for efficency.

  gamete *g = new gamete(m);

  // Take the first p alleles from g0 and par0, the remaining alleles
  // come from g1/par1 
  memcpy(g->h, g0, p*sizeof(haplotype));
  memcpy(g->gpar, par0, p*sizeof(haplotype));

  if(p<m+1)
    {
      // copy starting from the p'th index
      memcpy(&(g->h[p]), &(g1[p]), (m+1-p)*sizeof(haplotype));
      memcpy(&(g->gpar[p]), &(par1[p]), (m+1-p)*sizeof(haplotype));
    }    
  
  return g;
} 
		   

		   

gamete *Individual::makeGamete(xorshift64& rand, int m)
{

    unsigned int r = rand.get_uint32();
    unsigned int h0 = r & 1;
    //unsigned int h1 = (h0+1)&1;
    r >>= 1;
    int p = m+1;
    while ((r & 1) && (p >= 1))
    {
        r >>= 1;
        --p;
    }

    
    gamete *g;
    if(h0 == 0)
      {
	g = makeGameteHelper(p, m, m_genes.g.h1, m_genes.g.h2, m_genes.par.h1, m_genes.par.h2);
	g->s = make_pair(m_genes.g.h1[0], m_genes.g.h2[0]);
      }
    else
      {
	g = makeGameteHelper(p, m, m_genes.g.h2, m_genes.g.h1, m_genes.par.h2, m_genes.par.h1);
	g->s = make_pair(m_genes.g.h2[0], m_genes.g.h1[0]);
      }
    
    // Fill in the array 
    fill_n(g->par, m+1, m_nPosition);
    g->del = (r&1) ? m_genes.del.first : m_genes.del.second;
    
    return g;
}



void Individual::initVars(unsigned int pos, int nOvules, int nMarkers, unsigned int nWeight, unsigned int mX, unsigned int mY)
{
  // Added due to issue where one of the constructors was initializing a variable
  // but not the other. Initialize all Individual's variables/settings here
  m_nMarkers = nMarkers;
  m_nOvules = nOvules;
  m_nPosition = pos;
  m_nWeight = nWeight;
  m_vOWeights.assign(nOvules, 0);

  // store x and y position on grid at initialization instead of during the loop
  maxX = mX;
  maxY = mY;
  x_coord = pos/mX;
  y_coord = pos%mY;

  // allocate an array of gametes that will be filled in by 
  m_vOvules = new gamete*[nOvules];
  for(int i = 0; i < nOvules; i++)
    {
      gamete *gam = new gamete(nMarkers);
      m_vOvules[i] = gam;
    } 
}

void Individual::initGenes(int m, haplotype *dad, haplotype *mom)
{
  // genes are a reference to the genotype already created on the heap 
  m_genes.g.h1 = dad;
  m_genes.g.h2 = mom;

  
  // created holders for parent/grandparent genotype
  m_genes.par.h1 = new haplotype[m+1];
  m_genes.par.h2 = new haplotype[m+1];
  m_genes.gpar.h1 = new haplotype[m+1];
  m_genes.gpar.h2 = new haplotype[m+1];
  

  // The individual's parent and grandparent genome is unknown.
  // TODO: If NULL_GENE == 0 then we can use memset which is faster.
  fill_n(m_genes.par.h1, m+1, NULL_GENE);
  fill_n(m_genes.par.h2, m+1, NULL_GENE);
  fill_n(m_genes.gpar.h1, m+1, NULL_GENE);
  fill_n(m_genes.gpar.h2, m+1, NULL_GENE);
  
  // no deleterous allele yet
  m_genes.del = std::make_pair(0,0);
}

void Individual::setCoordinates(int position, int nMaxX, int nMaxY)
{
  // Since each individual doesn't change location on the grid, calculate
  // the x and y coordinates once instead of calling i2xy() every generation.
  // Could not find any condition where it would be neccessary to call function
  // but added it anyways.
  assert(0 <= position && position < nMaxX*nMaxY);
  x_coord = position/nMaxX;
  y_coord = position%nMaxY;
}


void Individual::clearOvuleWeight(int nOvules)
{
    m_vOWeights.assign(nOvules, 0);
}

dominance Individual::dRank;
string Individual::name;
typedef bool(Individual::*fptr)(gamete&);
fptr Individual::op;


Individual::~Individual()
{
  // This should only be called on initialization (when copy is put into pop
  // vector) and at the end of the program. 
  delete[] m_genes.g.h1;
  delete[] m_genes.g.h2;
  delete[] m_genes.par.h1;
  delete[] m_genes.par.h2;
  delete[] m_genes.gpar.h1;
  delete[] m_genes.gpar.h2;
  for(unsigned int a = 0; a < m_nOvules; a++)
    delete m_vOvules[a]; // delete every referenced object
  delete [] m_vOvules; // delete the array of reference
}
