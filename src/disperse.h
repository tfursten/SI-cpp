#ifndef DISPERSE_INCLUDED
#define DISPERSE_INCLUDED

#include <boost/algorithm/string/predicate.hpp>

#include "xorshift64.h"
#include "rexp.h"
#include "rnormal.h"
#include "ring.h"
#include "ray.h"
#include "disk.h"

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

typedef pair<int,int> xyCoord;

class Dispersal
{
protected:
    typedef int(Dispersal::*fptr)(xorshift64&,int,int);
    int(Dispersal::*boundary)(int,int);


public:
	template<class A>
	bool initialize(A &dist_name, int x, int y, bool fast, std::string landscape, float p1, float p2=0) {
        bool ff = fast;
        if(landscape == "torus")
            boundary = &Dispersal::periodic;
        else
            boundary = &Dispersal::absorbing;
        maxX = x;
        maxY = y;

        if(ff){
            static const char name_keys[][16] = {
            "triangular", "rayleigh", "ring", "uniform"
            };
            static const fptr dist_ops[] = {
                &Dispersal::disc_triangular,
                &Dispersal::disc_rayleigh,
                &Dispersal::disc_ring,
                &Dispersal::disc_uniform
            };
            int pos = key_switch(dist_name, name_keys); 
            if(pos == -1) {
                std::cerr << "No fast version available" << std::endl;
                ff = false;
            }
            else{
            name = std::string(name_keys[pos]);
            op = dist_ops[pos];
            set_param(name,p1,p2);
            init_disc(name);
            }

        }
        if(!ff){
            static const char name_keys[][16] = {
            "exponential", "triangular", "normal", "rayleigh", "ring", "gamma"
            };
            static const fptr dist_ops[] = {
                &Dispersal::cont_exponential,
                &Dispersal::cont_triangular,
                &Dispersal::cont_halfNormal,
                &Dispersal::cont_rayleigh,
                &Dispersal::cont_ring,
                &Dispersal::cont_gamma
            };
            int pos = key_switch(dist_name, name_keys); 
            if( pos == -1) {
                std::cerr << "ERROR: Invalid dispersal distribution" << std::endl;
                return false;
            }
            name = std::string(name_keys[pos]);
            op = dist_ops[pos];
            set_param(name,p1,p2);

        }
        return true;
    }

	    
    inline std::string getName() { return name; }
    
    inline double operator()(xorshift64& rand, int x1, int y1) {
    	return (this->*op)(rand, x1, y1);
    }
    void set_param(std::string name, float p1, float p2);
	int cont_exponential(xorshift64& rand, int x1, int y1);
	int cont_triangular(xorshift64& rand, int x1, int y1);
	int cont_halfNormal(xorshift64& rand, int x1, int y1);
	int cont_rayleigh(xorshift64& rand, int x1, int y1);
    int cont_ring(xorshift64& rand, int x1, int y1);
    int cont_gamma(xorshift64& rand, int x1, int y1);
    int disc_triangular(xorshift64& rand, int x1, int y1);
	int disc_rayleigh(xorshift64& rand, int x1, int y1);
    int disc_ring(xorshift64& rand, int x1, int y1);
    //int disc_rectangular(xorshift64& rand, int x1, int y1);
    int disc_uniform(xorshift64& rand, int x1, int y1);

	
private:
    std::string name;
    int maxX, maxY;
    float param1, param2;
    void init_disc(std::string name);
    int disperse_cont(xorshift64& rand, int x, int y, double d);
    int disperse_disc(int x1, int y1, int x2, int y2);
    int absorbing(int x, int y);
    int periodic(int x, int y);
    fptr op;
    Disk disk;
    Ring ring;
    Ray ray1;
    Ray ray2;
    //Rect rect; not ready
    

};

#endif // DISPERSE_INCLUDED
