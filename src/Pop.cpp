#include <fstream>
#include "Pop.h"

//urandom dev/urandom if it exists use it else use create random seed
using namespace std;

// random seed generator
inline unsigned int create_random_seed() {
	unsigned int v;
	//ifstream urandom("/dev/urandom", ios::in|ios::binary);
	//if(urandom.good()) {
		//urandom >> v;
	//} else {
	v = static_cast<unsigned int>(getpid());
	v += ((v << 15) + (v >> 3)) + 0x6ba658b3; // Spread 5-decimal PID over 32-bit number
	v^=(v<<17);
	v^=(v>>13);
	v^=(v<<5);
	v += static_cast<unsigned int>(time(NULL));
	v^=(v<<17);
	v^=(v>>13);
	v^=(v<<5);
	//}
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

void Population::initialize(int nMaxX, int nMaxY, string bound, int nPollen,\
 int nOvule, int nMarkers, int nDel, float dSigmaP, float dSigmaS,  string si, string dist_name, \
 float pp, float sp, bool fast, string mut_type, int nAlleles)
{
  m_nMaxX = nMaxX;
  m_nMaxY = nMaxY;
  m_sBound = bound;
  m_nPollen = nPollen;
  m_nOvule = nOvule;
  m_nMarkers = nMarkers;
  m_nDel = nDel;
  m_nIndividuals = nMaxX * nMaxY;
  m_nSalleles = nAlleles; //pop will be initialized wth nAlleles S-alleles
  m_nAlleles = 0;
  m_nMutAlleles = nAlleles;
  if(mut_type == "IAM")
    mutated_allele = &Population::mutate_iam;
  else if(mut_type == "KAM")
    mutated_allele = &Population::mutate_kam;
  ostringstream out;
  m_dSigmaP = dSigmaP;
  m_dSigmaS = dSigmaS;
  pDisp.initialize(dist_name, m_nMaxX, m_nMaxY, fast, m_sBound, m_dSigmaP, pp);
  sDisp.initialize(dist_name, m_nMaxX, m_nMaxY, fast, m_sBound, m_dSigmaS, sp);
  out << "Pollen dispersal distribution set to " << pDisp.getName() << ".\n";
  out << "Seed dispersal distribution set to " << sDisp.getName() << ".\n";
  out << "Landscape set to " << m_sBound << ".\n";

  Individual::initialize(si);
  Individual::initDomRank(m_myrand,2*m_nIndividuals);

 
  // Initialize Population: each individual has unique marker allele
  // Pop initialized with m_nMutAllele S-alleles
  for(int iii=0; iii<m_nIndividuals; iii++)
      {
      haplotype *h0 = new haplotype[nMarkers+1];
      haplotype *h1 = new haplotype[nMarkers+1];
      h0[0] = iii%m_nMutAlleles;
      h1[0] = (iii + 1) % m_nMutAlleles;
      //h0[0] = m_nSalleles++;
      //h1[0] = m_nSalleles++;
      for(int jjj=1; jjj<nMarkers+1; jjj++)
      {
          h0[jjj] = m_nAlleles++;
          h1[jjj] = m_nAlleles++;
      }
      
      m_vPop1.emplace_back(iii, m_nOvule, m_nMarkers, m_nDel, m_nMaxX, m_nMaxY, h0, h1);
      m_vPop2.emplace_back(iii, m_nOvule, m_nMarkers, m_nDel, m_nMaxX, m_nMaxY);
      m_nPopSample.push_back(iii);
      }

 
  out << "Self-Incompatibility set to " << Individual::getName() << ".\n";
  cout << out.str();
  pout << out.str();

}


void Population::param(double dSMut, double dMMut, double dDMut, double dPdel, unsigned int seed)
{
    //set Random seed
    if (seed==0) {
      seed = create_random_seed();
      pout << "Using Generated PRNG Seed: "<< seed << endl;
      cout << "Seed " << seed << endl;
    }
    m_myrand.seed(seed);
    m_pDel = dPdel;

    m_pMut = -log((pow(1-dMMut, m_nMarkers))*(1-dSMut)*pow(1-dDMut, m_nDel));
    double totmut =  m_nDel*dDMut + dSMut + m_nMarkers*dMMut;
    m_pDMut = (m_nDel * dDMut)/totmut;
    m_pMMut = (m_nMarkers * dMMut)/totmut;
    m_nLoci = 1;
    setMutCount();
}

void Population::setMutCount() {
    m_nMutCount = floor(rand_exp(m_myrand, m_pMut));
    //cout << "mut count: " << m_nMutCount << endl;
}

int Population::mutate_iam()
{
    return m_nSalleles ++;
}

int Population::mutate_kam()
{ // careful casting uint64_t to int, ok because n_alleles will typically not be large
    return (int)m_myrand.get_uint(m_nMutAlleles);
}

void Population::mutCountDec(){
    if(m_nMutCount-m_nLoci > 0){
      m_nMutCount -= m_nLoci;
      return;
    }
    setMutCount();
      return;
}

void Population::mutate(Individual &I)
{
    int count = m_nMutCount-m_nLoci;
    if (count > 0){
      m_nMutCount -= m_nLoci;
      return;
    }
    while(count< m_nLoci){
      bool chr = m_nMutCount == 1 ? 0 : 1;
      setMutCount();
      count += m_nMutCount;
      double ran = m_myrand.get_double52();
      if(ran < m_pDMut) // deleterious mutation occurs
      {
          if(chr)
            I.setMomDel(m_myrand.get_uint(m_nDel), 1);
          else
            I.setDadDel(m_myrand.get_uint(m_nDel), 1);
      }
      else if ((ran < (m_pMMut + m_pDMut)) && (ran >= m_pDMut))
      { // mutation in marker
          if(chr)
            I.setMomGene(1+m_myrand.get_uint(m_nMarkers),m_nAlleles++);
          else
            I.setDadGene(1+m_myrand.get_uint(m_nMarkers),m_nAlleles++);
      }
      else
      { //mutation in S-allele
          int allele = (this->*mutated_allele)();
          if(chr)
            I.setMomGene(0, allele);
          else
            I.setDadGene(0, allele);
          Individual::addDomRank(m_myrand.get_uint64());
      }
    }
}

void Population::evolve(int nBurnIn, int nGenerations, int nSample)
{
    //run burn-in period
    for(int ggg=0;ggg<nBurnIn;++ggg)
    {
        m_nLethal = 0;
        for(int dad=0; dad<m_nIndividuals;dad++)
            pollenDispersal(dad);
        m_nFecundity = 0;
        for(int mom=0; mom<m_nIndividuals;mom++)
            seedDispersal(mom);
        std::swap(m_vPop1,m_vPop2);
    }

    dout << "#Gen\tm\tibd\thibd\tko\tke\ts2\tthke\tNke\tNbke\tsterile\thoDel\the\thoDom\tdfreq\tN\tfecund" << endl;
    for(int ggg=0;ggg<nGenerations;++ggg)
    {
        m_nLethal = 0;
        for(int dad=0; dad<m_nIndividuals;dad++)
            pollenDispersal(dad);
        m_nFecundity = 0;
        for(int mom=0; mom<m_nIndividuals;mom++)
            seedDispersal(mom);
        if (ggg % nSample == 0)
            samplePop(ggg);

        std::swap(m_vPop1,m_vPop2);
      
    }
    
}

inline int wrap_around(int x, int w) {
    return ((x % w) + w) % w;
}


void Population::pollenDispersal(int dad)
{
    Individual &dadHere = m_vPop1.at(dad);
    if(dadHere.weight() == 0)
        return;
    for(int d=0; d<m_nDel; d++){
      if(dadHere.dadDel(d)!=0 && dadHere.momDel(d) != 0){
        if(m_pDel == 1.0 || m_myrand.get_double52() < m_pDel){
          m_nLethal ++;
          dadHere.setWeight(0);
          mutCountDec();
          return;
        }
      }
    }

    mutate(dadHere);
    int nX = dadHere.coordX();
    int nY = dadHere.coordY();
    ////position xy = i2xy(dad, m_nMaxX, m_nMaxY);
    ////int nX = xy.first;
    ////int nY = xy.second;

    for (int p=0; p < m_nPollen; p++)
    {
        int nNewCell = pDisp(m_myrand,nX, nY);
        if (nNewCell == -1)
            continue;
        Individual &momHere = m_vPop1[nNewCell];
        if (momHere.weight() == 0)
            continue;
        unsigned int nPollenWeight = m_myrand.get_uint32();
        //int ovule = m_myrand.get_uint(m_nOvule);
        if (nPollenWeight < momHere.getMinOvuleWeight()){
          continue;
        }
        gamete *pollen = dadHere.makeGamete(m_myrand,m_nMarkers,m_nDel);
        if(!(momHere(*pollen))){ //check compatibility
          delete pollen;
          continue;
        }
  

  // Moved conditional prior to make gamete for efficiency
  //unsigned int nPollenWeight = m_myrand.get_uint32();
        //int ovule = m_myrand.get_uint(m_nOvule);
        //if (nPollenWeight < momHere.ovuleWeight(ovule))
  //  {
  //    delete pollen;
        //    continue;
  //  }
        int ovule = momHere.getMinOvuleID();
        momHere.setOvuleWeight(nPollenWeight,ovule);
        momHere.setOvuleGenes(pollen,ovule);
  
    }
}

void Population::seedDispersal(int mom)
{
  Individual &momHere = m_vPop1[mom];
  if(momHere.weight() == 0) //No plant here
    {
      momHere.clearOvuleWeight(m_nOvule);  //clear out ovules for next generation
      return;
    }
  momHere.setWeight(0); //clear out weight to mark as completed

  int nX = momHere.coordX();
  int nY = momHere.coordY();
  //position xy = i2xy(mom,m_nMaxX,m_nMaxY);  
  //int nX = xy.first;
  //int nY = xy.second;
  for (int seed=0; seed < m_nOvule; seed++)
    {
      if (momHere.ovuleWeight(seed)==0) //No pollen landed in this ovule
          continue; //skip to next ovule
      m_nFecundity ++;
      momHere.setOvuleWeight(0, seed); //clear out weight 

      int nNewCell = sDisp(m_myrand,nX, nY); 
      if (nNewCell == -1) //reject if dispersal out of landscape
        continue;
      unsigned int nSeedWeight = m_myrand.get_uint32();
      if (nSeedWeight < m_vPop2[nNewCell].weight())
        continue;
      gamete *pollen = momHere.ovule(seed);
      gamete *ovule = momHere.makeGamete(m_myrand,m_nMarkers,m_nDel);
      // TODO: newIndividual should hold the reference to ovule and pollen instead
      // of making a copy.
      m_vPop2[nNewCell].newIndividual(pollen, ovule, nSeedWeight);
      delete ovule; 
    }
}


double euclideanDist2(int i, int j, int mx, int my) {
	auto xy1 = i2xy(i,mx,my);
	auto xy2 = i2xy(j,mx,my);
	double dx = xy1.first - xy2.first;
	double dy = xy1.second - xy2.second;
	return (dx*dx+dy*dy);
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


struct popstats {
  typedef map<int,int> alleledb; //map allele number to counts
	alleledb num_allele;
	size_t num_homo;
	size_t num_ibd;
	size_t sum_dist2;
  size_t num_del1; //number heterozygous for deleterious allele
  size_t num_del2; //number homozygous for deleterious allele
  size_t del_freq; //freq of deleterious allele

	popstats() : num_allele(), num_homo(0), num_ibd(0), sum_dist2(0), num_del1(0), num_del2(0), del_freq(0)
		{ }
};

template<typename T>
inline T sq(const T& t) {
	return t*t;
}




void Population::samplePop(int gen)
 {
     if (gen == 0)
         return;
     int sampleSz = 500; 
     double M = 2.0*sampleSz;
     double s2 = 0.0;
     int popCount = m_nIndividuals;
     for(int i=0; i<m_nIndividuals; i++){
        if(m_vPop2[i].weight() == 0)
          popCount -= 1;
     }

     random_shuffle(m_nPopSample.begin(),m_nPopSample.end());
     

     for(int m = 0; m < m_nMarkers+1; ++m)
     {
         int i = 0;
         int sample_count = 0;
         popstats stats;
         stats.del_freq = 0;
         while(i<(int)(m_nPopSample.size()) && sample_count < sampleSz)
         //for(int i=0; i<sampleSz; i++)
         {
             Individual &I = m_vPop2[m_nPopSample[i]];
             //cout << "m: " << m << "dad: " << I.dadGene(m) << "\n"
             //<< "mom: " << I.momGene(m) << "\n"
             //<< "dadPos: " << I.dadPos() << "\n"
             //<< "momPos: " << I.momPos() << "\n" 
             //<< "dadGpar: " << I.dadGpar(m) << "\n"
             //<< "momGpar: " << I.momGpar(m) << endl;
             if(I.weight()== 0){
                i++;
                continue;
             }

             int pos = I.position();
             //position xy = i2xy(i, m_nMaxX, m_nMaxY);
             //int dX = xy.first;
             //int dY = xy.second;
             stats.num_allele[I.dadGene(m)] += 1;
             stats.num_allele[I.momGene(m)] += 1;
             if(I.dadGene(m) == I.momGene(m))
                 stats.num_homo += 1;
             if(I.dadGpar(m) == I.momGpar(m)){
                 stats.num_ibd += 1;
               }
             if(m == 0)
             {
                 if(m_sBound=="rectangle"){
                   stats.sum_dist2 += euclideanDist2(I.dadPos(),pos,m_nMaxX, m_nMaxY);
                   stats.sum_dist2 += euclideanDist2(I.momPos(),pos,m_nMaxX, m_nMaxY);
                 }
                 else{
                   stats.sum_dist2 += minEuclideanDist2(I.dadPos(),pos,m_nMaxX, m_nMaxY);
                   stats.sum_dist2 += minEuclideanDist2(I.momPos(),pos,m_nMaxX, m_nMaxY);
                 }
                 if(I.dadDel(0)!=0)
                     stats.del_freq ++;
                 if(I.momDel(0)!=0)
                     stats.del_freq ++;
                 if(I.dadDel(0)!=0 && I.dadDel(0)==I.momDel(0)){
                     stats.num_del2++;
                 }
                 else if(I.dadDel(0)!=0 || I.momDel(0)!=0){
                     stats.num_del1++;
                 }
             }
             i++;
             sample_count ++;
         }
 
         int dt = 0;
         foreach(popstats::alleledb::value_type &vv, stats.num_allele){
             dt += sq(vv.second);
         }
 
         if(m==0)
             s2 = 0.5*stats.sum_dist2/M;
 
         double ibd = stats.num_ibd/(1.0*sampleSz);
         double hibd = stats.num_homo/(1.0*sampleSz);
         double f = dt/sq(M);
         double Ke = 1.0/f;
 
         double theta_ke = Ke-1.0;
         double N_ke = 0.25*theta_ke/m_pMut;
         double Nb_ke = 4.0*M_PI*s2*N_ke/(m_nMaxX*m_nMaxY);
         double Ko = (double)stats.num_allele.size();
         string t = "\t";
         int hoDom = sampleSz - stats.num_del1 - stats.num_del2;
 
  // TODO: Try using fopen/fwritef, the c routines are often faster than c++
  dout <<gen<<t<<m<<t<<ibd<<t<<hibd<<t<<Ko<<t<<Ke<<t<<s2<<t<<theta_ke<<t
  <<N_ke<<t<<Nb_ke<<t<<m_nLethal<<t<<stats.num_del2<<t<<stats.num_del1<<t
  <<hoDom<<t<<stats.del_freq/M<<t<<popCount<<t<<m_nFecundity<<endl;
 
     }
 }



