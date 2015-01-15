#include <iostream>


#include "disperse.h"

inline int xy2i(int x, int y, int mx, int my) {
    assert(0 <= x && x < mx);
    assert(0 <= y && y < my);
    return x*my+y;
}

void Dispersal::set_param(std::string name, float p1, float p2)
{
    assert(p1>=0);
    if (name == "exponential")
        param1 = 1.0/p1;
    else if (name == "triangular")
        param1 = 2.0*p1;
    else if (name == "normal")
        param1 = p1 * sqrt(2.0);
    else if (name == "gamma"){
        assert((p1 != 0) && (p2 != 0));
        param1 = p1; //not working yet
        param2 = p2;
    }
    //else if (name == "rectangle") //not working yet
        //assert((p2>=0) && (p3>=0) && (p4>=0));
        //param1 = p1; //left
        //param2 = p2; //right
        //param3 = p3; //up
        //param4 = p4; //down
    else if (name == "ring"){
        param1 = p1;
        assert((p2>=0) && (p2<=1));
        param2 = p2; //center
    }
    else if (name == "rayleigh"){
        param1 = p1;
        assert(p2>=0);
        if(p2 != 0){
            param2 = p2;
        }
        else param2 = p1;
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
    //if (name == "rectangle")
        //rect.initialize(param1,param2,param3,param4);
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

int Dispersal::cont_ring(xorshift64& rand, int x1, int y1)
{
    if((param2==0)||(rand.get_double52()>param2)){
        return disperse_cont(rand, x1, y1, param1);
    }
    else
        return (this->*boundary)(x1,y1);
}

int Dispersal::cont_gamma(xorshift64& rand, int x1, int y1)
{
    double d = 0; //not finished
    return disperse_cont(rand, x1, y1, d);
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

//int Dispersal::disc_rectangular(xorshift& rand, int x1, int y1)
//{
    //xyCoord dXY = rect.disperse(rand.get_uint64());
    //return disperse_rect(x1,y1,dXY.first,dXY.second);
  //  return disperse_rect(x1,y1,x1,y1);
//}

int Dispersal::disc_uniform(xorshift64& rand, int x1, int y1)
{
    int dX = rand.get_uint64() % maxX;
    int dY = rand.get_uint64() % maxY;
    return (this->*boundary)(dX,dY);

}

inline int wrap_around(int x, int w) {
    return ((x%w)+w)%w;
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
    double a = rand.get_double52() * 2.0 * M_PI;
    double newX = floor(d*cos(a)+x+0.5);
    double newY = floor(d*sin(a)+y+0.5);
    int i = (this->*boundary)(newX,newY);
    return i;

}

int Dispersal::disperse_disc(int x1, int y1, int x2, int y2)
{
    double newX = x1+x2;
    double newY = y1+y2;
    int i = (this->*boundary)(newX,newY);
    return i;
}

