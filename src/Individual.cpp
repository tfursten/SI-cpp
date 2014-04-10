#include <iostream>
#include "Individual.h"
using namespace std;


void printGenotype(haplotype &h)
{
    for(int i=0; i<h.size(); i++)
        cout << h.at(i) << ",";
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
    if (m_genes.g[0][0] == dad.h[0])
        return false;
    if (m_genes.g[1][0] == dad.h[0])
        return false;
    return true;
}

bool Individual::ssi(gamete &dad)
{
    if (m_genes.g[0][0] == dad.s.first)
        return false;
    if (m_genes.g[1][0] == dad.s.first)
        return false;
    if (m_genes.g[0][0] == dad.s.second)
        return false;
    if (m_genes.g[1][0] == dad.s.second)
        return false;
    return true;

}

bool Individual::bsi(gamete &dad)
{
    int dominant = max(getDomRank(dad.s.first), getDomRank(dad.s.second));
    if (m_genes.g[0][0] == dominant)
        return false;
    if (m_genes.g[1][0] == dominant)
        return false;
    return true;
}

gamete Individual::makeGamete(xorshift64& rand, int m)
{
    unsigned int r = rand.get_uint32();
    unsigned int h0 = r & 1;
    unsigned int h1 = (h0+1)&1;
    r >>= 1;
    int p = m+1;
    while ((r & 1) && (p >= 1))
    {
        r >>= 1;
        --p;
    }

    gamete g;
    g.h.assign(m_genes.g.at(h0).begin(),(m_genes.g.at(h0).begin()+p));
    g.gpar.assign(m_genes.par.at(h0).begin(),m_genes.par.at(h0).begin()+p);
    g.par.assign((m+1),m_nPosition);

    if(p<m+1)
    {
        g.h.insert(g.h.end(),m_genes.g.at(h1).begin()+p, m_genes.g.at(h1).end());
        g.gpar.insert(g.gpar.end(),m_genes.par.at(h1).begin()+p, m_genes.par.at(h1).end());
    }
    g.del = (r&1) ? m_genes.del.first : m_genes.del.second;
    g.s = make_pair(m_genes.g.at(h0).front(), m_genes.g.at(h1).front());
    return g;
}

void Individual::initGenes(int m, haplotype &dad, haplotype &mom)
{
    m_genes.g.push_back(dad);
    m_genes.g.push_back(mom);
    haplotype h;
    h.assign((m+1),-1);
    for(int i=0; i<2; ++i)
    {
        m_genes.par.push_back(h);
        m_genes.gpar.push_back(h);
    }
    m_genes.del = std::make_pair(0,0);
}

void Individual::clearOvuleWeight(int nOvules)
{
    m_vOWeights.assign(nOvules, 0);
}

dominance Individual::dRank;
string Individual::name;
typedef bool(Individual::*fptr)(gamete&);
fptr Individual::op;
