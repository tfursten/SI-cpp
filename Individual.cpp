#include <iostream>
#include "Individual.h"
using namespace std;



bool Individual::nsi(Individual dad)
{
    return true;
}

bool Individual::psi(Individual dad)
{
    if (m_nPosition == dad.position())
        return false;
    return true;
}

bool Individual::gsi(Individual dad)
{
    if (m_genes.g[0][0] == dad.getGamete().h[0])
        return false;
    if (m_genes.g[1][0] == dad.getGamete().h[0])
        return false;
    return true;
}

bool Individual::ssi(Individual dad)
{
    if (m_genes.g[0][0] == dad.genes()[0][0])
        return false;
    if (m_genes.g[1][0] == dad.genes()[0][0])
        return false;
    if (m_genes.g[0][0] == dad.genes()[1][0])
        return false;
    if (m_genes.g[1][0] == dad.genes()[1][0])
        return false;
    return true;

}

bool Individual::bsi(Individual dad)
{
    int dominant = (dRank.at(dad.genes()[0][0])>dRank.at(dad.genes()[1][0])) ? dad.genes()[0][0] : dad.genes()[1][0];
    if (m_genes.g[0][0] == dominant)
        return false;
    if (m_genes.g[1][0] == dominant)
        return false;
    return true;
}

void Individual::makeGamete(xorshift64& rand, int nMarkers)
{
    gamete gam;
    int m = nMarkers + 1;
    unsigned int r = rand.get_uint32();
    unsigned int h0 = r & 1;
    unsigned int h1 = (h0+1)&1;
    r >>= 1;
    int p = 0;
    while (r & 1 && m > 0)
    {
        r >>= 1;
        ++p;
        --m;
    }

    gam.h.assign(m_genes.g[h0].begin(),m_genes.g[h0].end()-p);
    gam.par.assign(id[h0].begin(),id[h0].end()-p);
    gam.gpar.assign(m_genes.par[h0].begin(),m_genes.par[h0].end()-p);
    if(p>0)
        {
            gam.h.insert(gam.h.end(),m_genes.g[h1].begin()+p, m_genes.g[h1].end());
            gam.par.insert(gam.par.end(),id[h1].begin()+p, id[h1].end());
            gam.gpar.insert(gam.par.end(),m_genes.par[h1].begin()+p, m_genes.par[h1].end());
        }
    r>>= 1;
    gam.del = (r&1) ? m_genes.del.first : m_genes.del.second;
    m_gamete = gam;
}


dominance Individual::dRank;
string Individual::name;
typedef bool(Individual::*fptr)(Individual);
fptr Individual::op;
