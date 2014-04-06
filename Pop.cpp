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
    assert(0 <= i && i < mx*my);
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

void Population::initialize(int nMaxX, int nMaxY, int nPollen, int nOvule, int nMarkers, string si, string dist_name)
{
    m_nMaxX = nMaxX;
    m_nMaxY = nMaxY;
    m_nPollen = nPollen;
    m_nOvule = nOvule;
    m_nMarkers = nMarkers;
    m_nIndividuals = nMaxX * nMaxY;
    int m_nSalleles = 0;
    int m_nAlleles = 0;
    ostringstream out;
    dist.initialize(dist_name);
    Individual::initialize(si);
    Individual::initDomRank(m_myrand.get_uint64(),2*m_nIndividuals);
    cout << "I am here" << endl;

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
        h0.push_back(0);
        h1.push_back(0);
        m_vPop1.emplace_back(iii,m_nOvule,h0,h1);
        m_vPop1.emplace_back(iii,m_nOvule,m_nMarkers);
    }

    out << "Dispersal distribution set to " << dist.getName() << ".\n";
    out << "Self-Incompatibility set to " << Individual::getName() << ".\n";
    cout << out.str();
    mout << out.str() << "#" << endl;

}

void Population::param(float dSigmaP, float dSigmaS, double dSMut, double dMMut, double dDMut, unsigned int seed)
{
    //set Random seed
    if (seed==0) {
        seed = create_random_seed();
        mout << "Using Generated PRNG Seed: "<< seed << endl;
    }
    m_myrand.seed(seed);
    m_dSigmaP = dSigmaP;
    m_dSigmaS = dSigmaS;
    m_pMut = -log(1-(pow(1-dMMut, m_nMarkers)*(1-dSMut)*(1-dDMut)));
    double totmut =  dDMut + dSMut + m_nMarkers*dMMut;
    m_pDMut = dDMut/totmut;
    m_pSMut = dSMut/totmut;
    setMutCount();
}

void Population::setMutCount() {
    m_nMutCount = floor(rand_exp(m_myrand, m_pMut));
}

void Population::mutate(gamete &g)
{
    if (--m_nMutCount > 0)
        return;
    setMutCount();
    double rand = m_myrand.get_double52();
    if (rand < m_pDMut)
        g.del = 1;
    if (rand < m_pSMut)
        {
        g.h.front() = m_nSalleles++;
        Individual::addDomRank(m_myrand.get_uint64());
        }
    else
        g.h.at(1+m_myrand.get_uint(m_nMarkers)) = m_nAlleles++; //make sure this doesn't affect last marker
}

void Population::evolve(int nBurnIn, int nGenerations, int nSample)
{
    //run burn-in period
    for(int ggg=0;ggg<nBurnIn;++ggg)
    {
        for(int dad=0; dad<m_nIndividuals;dad++)
            pollenDispersal(dad);
        for(int mom=0; mom<m_nIndividuals;mom++)
            seedDispersal(mom);
        std::swap(m_vPop1,m_vPop2);
    }


    for(int ggg=0;ggg<nGenerations;++ggg)
    {
        for(int dad=0; dad<m_nIndividuals;dad++)
            pollenDispersal(dad);
        for(int mom=0; mom<m_nIndividuals;mom++)
            seedDispersal(mom);
        //if (ggg % nSample == 0)
            //samplePop(ggg);
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
        int nNewCell = disperse(nX,nY, m_dSigmaP);
        if (nNewCell == -1)
            continue;
        Individual &momHere = m_vPop1[nNewCell];
        if (momHere.weight() == 0)
            continue;
        dadHere.makeGamete(m_myrand, m_nMarkers);
        gamete pollen = dadHere.getGamete();
        mutate(pollen);
        if(!(momHere(dadHere))) //check compatibility
            continue;
        unsigned int nPollenWeight = m_myrand.get_uint32();
        int ovule = m_myrand.get_uint(m_nOvule);
        cout << ovule << "ovule" << endl;
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
        return;
    momHere.setWeight(0);
    position xy = i2xy(mom,m_nMaxX,m_nMaxY);
    int nX = xy.first;
    int nY = xy.second;
    for (int seed=0; seed < m_nOvule; seed++)
    {
        if (momHere.ovuleWeight(seed)==0)
            continue;
        momHere.setOvuleWeight(seed,0);
        int nNewCell = disperse(nX,nY, m_dSigmaS);
        if (nNewCell == -1)
            continue;
        unsigned int nSeedWeight = m_myrand.get_uint32();
        if (nSeedWeight < m_vPop2[nNewCell].weight())
            continue;
        gamete pollen = momHere.ovule(seed);
        momHere.makeGamete(m_myrand, m_nMarkers);
        gamete ovule = momHere.getGamete();
        mutate(ovule);
        if (ovule.del==1 && pollen.del==1)
            continue;
        m_vPop2[nNewCell].newIndividual(pollen,ovule);

    }
}


double minEuclideanDist2(int i, int j, int mx, int my) {
	auto xy1 = i2xy(i,mx,my);
	auto xy2 = i2xy(j,mx,my);
	double dx = abs(1.0*(xy1.first - xy2.first));
	double dy = abs(1.0*(xy1.second - xy2.second));
	dx = (dx < mx*0.5) ? dx : mx-dx;
	dy = (dy < my*0.5) ? dy : my-dy;
	return (dx*dx+dy*dy);
}

void Population::samplePop(int gen)
{
    vector<int> vIBD(1+m_nMaxY/2,0);
    vector<int> vN(1+m_nMaxY/2,0);
    typedef map<int,int> mapType;
    mapType alleleMap;
    int szSample = 0;
    double dSigma2 = 0.0;
    double ko = 0.0;
    double ke = 0.0;
/*
    int i0 = m_nTransPos * m_nMaxY;
    for(int i = i0; i < i0+m_nMaxY; ++i) {
    	Individual ind = m_vPop2[i];
        if(ind.weight() == 0)
    		continue;
    	szSample += 1;
    	alleleMap[ind.nAllele] += 1;
    	int p = ind.nParent_id;
    	dSigma2 += minEuclideanDist2(i,p,m_nMaxX,m_nMaxY);

    	for(int j=i; j < i0+m_nMaxY; ++j) {
    		if(m_vPop2[j].nWeight == 0)
    			continue;
    		int k = j-i;
    		k = (k <= m_nMaxY/2) ? k : m_nMaxY-k;
       		if(ind.nAllele == m_vPop2[j].nAllele) {
    			vIBD[k] += 1;
    		}
    		vN[k] += 1;
    	}
    }
    int df = 0;
    BOOST_FOREACH(mapType::value_type & vv, alleleMap){
            df += vv.second*vv.second;
        }
    ko = (double)alleleMap.size();
    double f = df/(double)(szSample*szSample);
    ke = 1.0/f;
    cout << "Ko: " << ko << " Ke: " << ke << endl;

    mout << "pIBD Gen " << gen << ": " << endl << "nIBD: \t";
    for(unsigned int k=0; k<vIBD.size();++k) {
        mout << vIBD[k] << ((k< vIBD.size()-1) ? "\t" : "\n");
    }
    mout << "Total: \t";
    for(unsigned int k=0; k< vN.size(); ++k) {
        mout << vN[k] << ((k < vN.size()-1 ? "\t" : "\n"));
    }
    mout << "sigma2: " << dSigma2/(2.0*szSample) << endl << endl;
*/
}


