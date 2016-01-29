#include <iostream>


#include "disperse.h"

inline int xy2i(int x, int y, int mx, int my) {
    assert(0 <= x && x < mx);
    assert(0 <= y && y < my);
    return x*my+y;
}

inline int wrap_around(int x, int w) {
    return ((x%w)+w)%w;
}

void Dispersal::set_param(std::string name, float p1, float p2)
{
    assert(p1>0);
    if (name == "exponential")
        param1 = 1.0/p1;
    else if (name == "triangular")
        param1 = 2.0*p1;
    else if (name == "normal")
        param1 = p1 * sqrt(2.0);
    else if (name == "gamma"){
        //Using Marsaglia2000 method requires a>=1
        assert((p1 > 0) && (p2 > 0));
        param1 = p2; //alpha
        param2 = sqrt((param1*(1+param1))/(2.0*p1*p1)); //beta adjusted so second moment equals 2sigma^2
        assert(param2>0);
        if(param1<1){
            param5 = param1;
            param1 += 1.0;
        }
        else param5 = 1;
        param3 = param1-(1/3.0); //d
        param4 = (1.0/3.0)/sqrt(param3); //c   
    }
    else if (name == "pareto"){
        param1 = p2; //alpha
        param2 = sqrt((2*p1*p1*(param1-2))/param1); //xmax
        //param1 = p1;
        //param2 = p2;
    }
    else if (name == "lomax"){
        param1 = p2;
        param2 = sqrt((-2+p2)*(-1+p2)*p1*p1);
    }
    else if (name == "ring"){
        assert((p2>=0) && (p2<=1));
	param1 = sqrt((2.0*p1*p1)/(1.0-p2));
        param2 = p2; //center
    }
    else if (name == "rice"){
        param1 = p1*cos(M_PI/4.0);//x
        param2 = p1*sin(M_PI/4.0);//y
        param3 = sqrt((p1*p1)/2.0); //sigma
    }
    else if (name == "rayleigh"){
        assert(p2>=0);
        if(p2 != 0){
            param1 = p1*p2; //sigma_x
            param2 = sqrt(2*p1*p1-(param1*param1)); //sigma_y
        }
        else{ 
            param1 = p1;
            param2 = p1;
        }
    }
    else
        param1 = p1;

}

void Dispersal::init_disc(std::string name)
{
    if (name == "ring")
        ring.initialize(param1,param2);
    if (name == "triangular")
        disk.initialize(param1);
    if (name == "rayleigh"){
        ray1.initialize(param1,6);
        ray2.initialize(param2,6);
    }
}

int Dispersal::cont_exponential(xorshift64& rand, int x1, int y1)
{
    double d = rand_exp(rand, param1);
    return disperse_cont(rand, x1, y1, d);
}

int Dispersal::cont_triangular(xorshift64& rand, int x1, int y1)
{
    double d = param1*sqrt(rand.get_double52());
    return disperse_cont(rand, x1, y1, d);
}

int Dispersal::cont_halfNormal(xorshift64& rand, int x1, int y1)
{
    double d = rand_abs_normal(rand, 0.0, param1);
    return disperse_cont(rand, x1, y1, d);
}

int Dispersal::cont_rayleigh(xorshift64& rand, int x1, int y1)
{
    double dX = floor(rand_normal(rand,0.0,param1)+x1+0.5);
    double dY = floor(rand_normal(rand,0.0,param2)+y1+0.5);
    return (this->*boundary)(dX,dY);
}

int Dispersal::cont_rice(xorshift64& rand, int x1, int y1)
{
    double dX = floor(rand_normal(rand,param1,param3)+x1+0.5);
    double dY = floor(rand_normal(rand,param2,param3)+y1+0.5);
    return (this->*boundary)(dX,dY);
}

int Dispersal::cont_ring(xorshift64& rand, int x1, int y1)
{
    if((param2==0)||(rand.get_double52()>param2)){
        return disperse_cont(rand, x1, y1, param1);
    }
    else
        return (this->*boundary)(x1,y1);
}

int Dispersal::cont_gamma(xorshift64&rand, int x1, int y1)
{
    //Marsalgia and Tsang's Method
    double x, v, u, d;
    for(;;){
        do{
            x = rand_normal(rand,0.0,1);
            v = 1.0 + param4 * x;
        }
        while (v <= 0);
        v = v*v*v;
        u = rand.get_double52();
        if (u < 1 - 0.0331 * x * x * x * x) 
          break;
        if (log(u) < 0.5 * x * x + param3 * (1 - v + log(v)))
          break;
    }

    d = 1/param2 * param3 * v;
    if(param5 >= 1){
        return disperse_cont(rand, x1, y1, d);
    }
    else{
        u = rand.get_double52();
        d = d * pow(u,1.0/param5);
        return disperse_cont(rand, x1, y1, d);
    }
}

int Dispersal::cont_pareto(xorshift64& rand, int x1, int y1){
    
    //double d = param2/pow(1-rand.get_double53(),1/param1);
    double d = param2*exp(rand_exp(rand,param1)); //fasr
    return disperse_cont(rand,x1,y1,d);
}

int Dispersal::cont_lomax(xorshift64& rand, int x1, int y1){
    double d = param2*exp(rand_exp(rand,param1));
    return disperse_cont(rand,x1,y1,d-param2);
}


int Dispersal::disc_triangular(xorshift64& rand, int x1, int y1)
{
    xyCoord dXY = disk.disperse(rand.get_uint64());
    return disperse_disc(x1,y1,dXY.first,dXY.second);
}

int Dispersal::disc_rayleigh(xorshift64& rand, int x1, int y1)
{
    int dX = ray1.disperse(rand);
    int dY = ray2.disperse(rand);
    return disperse_disc(x1,y1,dX,dY);
}

int Dispersal::disc_ring(xorshift64& rand, int x1, int y1)
{
    xyCoord dXY = ring.disperse(rand.get_uint64());
    return disperse_disc(x1,y1,dXY.first,dXY.second);
}

int Dispersal::disc_uniform(xorshift64& rand, int x1, int y1)
{
    int dX = rand.get_uint64() % maxX;
    int dY = rand.get_uint64() % maxY;
    return (this->*boundary)(dX,dY);
}

int Dispersal::absorbing(int x, int y)
{
    if (x >= 0 && x < maxX && y >= 0 && y < maxY)
        return xy2i(x,y,maxX,maxY);
    return -1;
}

int Dispersal::periodic(int x, int y)
{
    int newX = wrap_around(x,maxX);
    int newY = wrap_around(y,maxY);
    return xy2i(newX,newY,maxX,maxY);
}


int Dispersal::disperse_cont(xorshift64& rand, int x, int y, double d)
{
    double a = rand.get_double53() * 2.0 * M_PI;  
    int newX = round(d*cos(a)+x);
    int newY = round(d*sin(a)+y);
    int i = (this->*boundary)(newX,newY);
    return i;
}

int Dispersal::disperse_disc(int x1, int y1, int x2, int y2)
{
    int newX = x1+x2;
    int newY = y1+y2;
    int i = (this->*boundary)(newX,newY);
    return i;
}

