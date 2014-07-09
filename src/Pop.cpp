#include <fstream>
#include "Pop.h"

//urandom dev/urandom if it exists use it else use create random seed
using namespace std;

// random seed generator
inline unsigned int create_random_seed() {
	unsigned int v;
	ifstream urandom("/dev/urandom", ios::in|ios::binary);
	if(urandom.good()) {
		urandom >> v;
	} else {
		v = static_cast<unsigned int>(getpid());
		v += ((v << 15) + (v >> 3)) + 0x6ba658b3; // Spread 5-decimal PID over 32-bit number
		v^=(v<<17);
		v^=(v>>13);
		v^=(v<<5);
		v += static_cast<unsigned int>(time(NULL));
		v^=(v<<17);
		v^=(v>>13);
		v^=(v<<5);
	}
    return (v == 0) ? 0x6a27d958 : (v & 0x7FFFFFFF); // return at most a 31-bit seed
}

typedef pair<int,int> position;

inline position i2xy(int i, int mx, int my)
{
    assert((0 <= i) && (i < mx*my));
    return make_pair((i/my),(i%my));
}

inline int xy2i(int x, int y, int mx, int my) {
	assert(0 <= x && x < mx);
	assert(0 <= y && y < my);
	return x*my+y;
}

inline int xy2i(position xy, int mx, int my) {
	return xy2i(xy.first,xy.second,mx,my);
}

void Population::initialize(int nMaxX, int nMaxY, int nPollen, int nOvule, int nMarkers,float dSigmaP, float dSigmaS,  string si, string dist_name)
{
    m_nMaxX = nMaxX;
    m_nMaxY = nMaxY;
    m_nPollen = nPollen;
    m_nOvule = nOvule;
    m_nMarkers = nMarkers;
    m_nIndividuals = nMaxX * nMaxY;
    m_nSalleles = 0;
    m_nAlleles = 0;
    ostringstream out;
    m_dSigmaP = dSigmaP;
    m_dSigmaS = dSigmaS;
    pdisk.initialize((double)(2.0*m_dSigmaP));
    sdisk.initialize((double)(2.0*m_dSigmaS));
    if (dist_name != "disk"){
        dist.initialize(dist_name);
        pDisperse = &Population::pDisperseDist;
        sDisperse = &Population::sDisperseDist;
        out << "Dispersal distribution set to " << dist.getName() << ".\n";
    }
    else{
        pDisperse = &Population::pDisperseDisk;
        sDisperse = &Population::sDisperseDisk;
        out << "Dispersal distribution set to " << "Disk" << ".\n";
    }

    Individual::initialize(si);
    Individual::initDomRank(m_myrand,2*m_nIndividuals);

    // Initialize Population: each individual has unique allele
    for(int iii=0; iii<m_nIndividuals; iii++)
    {
        haplotype h0;
        haplotype h1;
        h0.push_back(m_nSalleles++);
        h1.push_back(m_nSalleles++);
        for(int jjj=0; jjj<nMarkers; jjj++)
        {
            h0.push_back(m_nAlleles++);
            h1.push_back(m_nAlleles++);
        }

        m_vPop1.emplace_back(iii,m_nOvule,m_nMarkers,h0,h1);
        m_vPop2.emplace_back(iii,m_nOvule,m_nMarkers);
    }


    out << "Self-Incompatibility set to " << Individual::getName() << ".\n";
    cout << out.str();
    pout << out.str();

}

void Population::param(double dSMut, double dMMut, double dDMut, unsigned int seed)
{
    //set Random seed
    if (seed==0) {
        seed = create_random_seed();
        pout << "Using Generated PRNG Seed: "<< seed << endl;
        cout << "Seed " << seed << endl;
    }
    m_myrand.seed(seed);
    m_pMut = -log((pow(1-dMMut, m_nMarkers))*(1-dSMut)*(1-dDMut));
    double totmut =  dDMut + dSMut + m_nMarkers*dMMut;
    m_pDMut = dDMut/totmut;
    m_pSMut = dSMut/totmut;
    setMutCount();
}

void Population::setMutCount() {

    m_nMutCount = floor(rand_exp(m_myrand, m_pMut));
}

void Population::mutCountDec() {
    if (--m_nMutCount > 0)
        return;
    setMutCount();
}


void Population::mutate(gamete &g)
{
    if (--m_nMutCount > 0)
        return;
    setMutCount();
    double rand = m_myrand.get_double52();
    if (rand < m_pDMut)
        {
        g.del = 1;
        return;
        }
    if ((rand < (m_pSMut+m_pDMut)) && (rand >= m_pDMut))
        {
        g.h.front() = m_nSalleles++;
        Individual::addDomRank(m_myrand.get_uint64());
        return;
        }
    if (rand >= (m_pSMut+m_pDMut))
        {
        g.h.at(1+m_myrand.get_uint(m_nMarkers)) = m_nAlleles++;
        return;
        }
}

void Population::evolve(int nBurnIn, int nGenerations, int nSample)
{
    //run burn-in period
    for(int ggg=0;ggg<nBurnIn;++ggg)
    {
        cout << "burn " << ggg << endl;
        for(int dad=0; dad<m_nIndividuals;dad++)
            pollenDispersal(dad);
        for(int mom=0; mom<m_nIndividuals;mom++)
            seedDispersal(mom);
        std::swap(m_vPop1,m_vPop2);
    }

    dout << "#Gen\tm\tibd\thibd\tko\tke\ts2\tthke\tNke\tNbke" << endl;
    for(int ggg=0;ggg<nGenerations;++ggg)
    {
        cout << "gen"<< ggg << endl;
        for(int dad=0; dad<m_nIndividuals;dad++)
            pollenDispersal(dad);
        for(int mom=0; mom<m_nIndividuals;mom++)
            seedDispersal(mom);
        if (ggg % nSample == 0)
            samplePop(ggg);
        std::swap(m_vPop1,m_vPop2);
    }
}


int Population::disperse(int x, int y, double sigma)
{
    double a = m_myrand.get_double52() * 2.0 * M_PI;
    double r = dist(m_myrand,sigma);
    double dX = static_cast<int>(floor(r*cos(a)+x+0.5));
    double dY = static_cast<int>(floor(r*sin(a)+y+0.5));
    if (dX >= 0 && dX < m_nMaxX && dY >= 0 && dY < m_nMaxY)
        return xy2i(dX,dY,m_nMaxX, m_nMaxY);
    return -1;
}

int Population::sDisperseDist(int x, int y)
{
    return disperse(x,y,m_dSigmaS);
}

int Population::pDisperseDist(int x, int y)
{
    return disperse(x,y,m_dSigmaP);
}

int Population::sDisperseDisk(int x, int y)
{
    return sdisk.disperse(x,y,m_myrand.get_uint64(),m_nMaxX,m_nMaxY);
}

int Population::pDisperseDisk(int x, int y)
{
    return pdisk.disperse(x,y,m_myrand.get_uint64(),m_nMaxX,m_nMaxY);
}


void Population::pollenDispersal(int dad)
{

    Individual dadHere = m_vPop1.at(dad);
    if(dadHere.weight() == 0)
        return;
    position xy = i2xy(dad,m_nMaxX,m_nMaxY);
    int nX = xy.first;
    int nY = xy.second;
    for (int p=0; p < m_nPollen; p++)
    {
        int nNewCell = (this->*pDisperse)(nX, nY);
        if (nNewCell == -1)
        {
            mutCountDec();
            continue;
        }
        Individual &momHere = m_vPop1[nNewCell];
        if (momHere.weight() == 0)
        {
            mutCountDec();
            continue;
        }
        gamete pollen = dadHere.makeGamete(m_myrand,m_nMarkers);
        mutate(pollen);
        if(!(momHere(pollen))) //check compatibility
            continue;
        unsigned int nPollenWeight = m_myrand.get_uint32();
        int ovule = m_myrand.get_uint(m_nOvule);
        if (nPollenWeight < momHere.ovuleWeight(ovule))
            continue;
        momHere.setOvuleWeight(nPollenWeight,ovule);
        momHere.setOvuleGenes(pollen,ovule);
    }
}

void Population::seedDispersal(int mom)
{
    Individual &momHere = m_vPop1[mom];
    if(momHere.weight() == 0)
    {
        momHere.clearOvuleWeight(m_nOvule);
        return;
    }
    momHere.setWeight(0);
    position xy = i2xy(mom,m_nMaxX,m_nMaxY);
    int nX = xy.first;
    int nY = xy.second;
    for (int seed=0; seed < m_nOvule; seed++)
    {
        if (momHere.ovuleWeight(seed)==0)
        {
            mutCountDec();
            continue;
        }
        momHere.setOvuleWeight(0,seed);
        int nNewCell = (this->*sDisperse)(nX, nY);
        if (nNewCell == -1)
        {
            mutCountDec();
            continue;
        }
        unsigned int nSeedWeight = m_myrand.get_uint32();
        if (nSeedWeight < m_vPop2[nNewCell].weight())
        {
            mutCountDec();
            continue;
        }
        gamete pollen = momHere.ovule(seed);
        gamete ovule = momHere.makeGamete(m_myrand,m_nMarkers);
        mutate(ovule);
        if (ovule.del==1 && pollen.del==1)
            {
                continue;
            }
        m_vPop2[nNewCell].newIndividual(pollen,ovule,nSeedWeight);
    }
}


double euclideanDist2(int i, int j, int mx, int my) {
	auto xy1 = i2xy(i,mx,my);
	auto xy2 = i2xy(j,mx,my);
	double dx = xy1.first - xy2.first;
	double dy = xy1.second - xy2.second;
	return (dx*dx+dy*dy);
}


struct popstats {
    typedef map<int,int> alleledb;
	alleledb num_allele;
	size_t num_homo;
	size_t num_ibd;
	size_t sum_dist2;

	popstats() : num_allele(), num_homo(0), num_ibd(0), sum_dist2(0)
		{ }
};

template<typename T>
inline T sq(const T& t) {
	return t*t;
}



void Population::samplePop(int gen)
{
    double M = 2.0*m_nIndividuals; //smaller sample??
    double s2 = 0.0;
    for(int m = 0; m < m_nMarkers+1; ++m)
    {
        popstats stats;
        for(int i=0; i<m_nIndividuals; i++)
        {
            Individual &I = m_vPop2[i];
            if(I.weight()== 0)
                continue;
            position xy = i2xy(i, m_nMaxX, m_nMaxY);
            int dX = xy.first;
            int dY = xy.second;
            stats.num_allele[I.dadGene(m)] += 1;
            stats.num_allele[I.momGene(m)] += 1;
            if(I.dadGene(m) == I.momGene(m))
                stats.num_homo += 1;
            if(I.dadGpar(m) == I.momGpar(m))
                stats.num_ibd += 1;
            if(m == 0)
            {
                stats.sum_dist2 += euclideanDist2(I.dadPos(),i,m_nMaxX, m_nMaxY);
                stats.sum_dist2 += euclideanDist2(I.momPos(),i,m_nMaxX, m_nMaxY);
            }
        }

        int dt = 0;
        foreach(popstats::alleledb::value_type &vv, stats.num_allele){
            dt += sq(vv.second);
        }

        if(m==0)
            s2 = 0.5*stats.sum_dist2/M;
        double ibd = stats.num_ibd/(1.0*m_nIndividuals);
        double hibd = stats.num_homo/(1.0*m_nIndividuals);
        double f = dt/sq(M);
        double Ke = 1.0/f;
        double theta_ke = Ke-1.0;
        double N_ke = 0.25*theta_ke/m_pMut;
        double Nb_ke = 4.0*M_PI*s2*N_ke/(m_nMaxX*m_nMaxY);
        double Ko = (double)stats.num_allele.size();
        string t = "\t";
        dout <<gen<<t<<m<<t<<ibd<<t<<hibd<<t<<Ko<<t<<Ke<<t<<s2<<t<<theta_ke<<t<<N_ke<<t<<Nb_ke<<endl;

    }
}







