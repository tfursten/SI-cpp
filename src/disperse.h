#ifndef DISPERSE_INCLUDED
#define DISPERSE_INCLUDED

#include <boost/algorithm/string/predicate.hpp>

#include "xorshift64.h"
#include "rexp.h"
#include "rnormal.h"

template<class A, class B, int N>
int key_switch(A &ss, const B (&key)[N])
{
    using boost::algorithm::istarts_with;
    for(int iii=0; iii<N; ++iii)
    {
        if(istarts_with(key[iii],ss))
            return iii;
    }
    return (int) -1;

}

class Dispersal
{
protected:
    typedef double(Dispersal::*fptr)(xorshift64&, double);

public:
	template<class A>
	bool initialize(A &dist_name) {
		static const char name_keys[][16] = {
		    "exponential", "triangular", "normal", "rayleigh"
		};
		static const fptr dist_ops[] = {
            &Dispersal::dist_exponential,
            &Dispersal::dist_triangular,
            &Dispersal::dist_halfNormal,
            &Dispersal::dist_rayleigh
        };
        int pos = key_switch(dist_name, name_keys);
        if( pos == -1) {
        	std::cerr << "ERROR: Invalid dispersal distribution" << std::endl;
        	return false;
        }
        name = std::string(name_keys[pos]);
        op = dist_ops[pos];
        return true;
	}
	    
    inline std::string getName() { return name; }
    
    inline double operator()(xorshift64& rand, double sigma) {
    	return (this->*op)(rand,sigma);
    }
    
	double dist_exponential(xorshift64& rand, double sigma);
	double dist_triangular(xorshift64& rand, double sigma);
	double dist_halfNormal(xorshift64& rand, double sigma);
	double dist_rayleigh(xorshift64& rand, double sigma);
	
	Dispersal() {
		initialize("exponential");
	}
	
private:
    std::string name;
    fptr op;
};

#endif // DISPERSE_INCLUDED
