#include "disk.h"

inline int xy2i(int x, int y, int mx, int my){
assert(0 <= x && x < mx);
	assert(0 <= y && y < my);
	return x*my+y;
}
void Disk::initialize(double r)
{
    radius = r;
    rSq = r*r;
    halfRsq = rSq/2.0;
    getCellRange();
    vecDim = cellRange.size()-1;
    m_maxX = getMaxX();
    totalArea = M_PI*rSq;
    makeTables();

}

double Disk::integrate(double x1, double x2)
{
    return (0.5 * sqrt(rSq-x2*x2) * x2 + rSq*0.5* asin(x2/radius)) - (0.5 * sqrt(rSq-x1*x1) * x1 + rSq*0.5* asin(x1/radius));
}

void Disk::getCellRange()
{
    cellRange.push_back(0.0);
    for(int i=0; i < radius; i++)
    {
        double x = i + 0.5;
        if (x<radius)
            cellRange.push_back(x);
    }
    cellRange.push_back(radius);
}

double Disk::getMaxX(){
    double mX = pol2xy(M_PI/4.0).first;
    for (int i=vecDim; i>=1; i--){
        if((mX <= cellRange[i]) && (mX > cellRange[i-1]))
            return i;
    }
}

pair<double,double> Disk::pol2xy(double theta)
{
    double x = radius*cos(theta);
    double y = radius*sin(theta);
    return make_pair(x,y);
}

double Disk::circle(double &x)
{
    return (double)sqrt(radius*radius-x*x);

}

int Disk::getBin(double x){
    for (int i=0; i<vecDim; i++){
        if(x <= cellRange[i+1])
            return i;
    }
}

void Disk::areas(double x1, double x2, int i){
    double A = integrate(x1,x2);
    for (int y = 0; y<vecDim; y++){
        double y1 = cellRange[y];
        double y2 = cellRange[y+1];
        double a = (x2-x1)*(y2-y1);
        if(a<A){
            probMap[make_pair(i,y)] += a;
            A -= a;
        }
        else{
        probMap[make_pair(i,y)] += A;
        break;
        }
    }
}

void Disk::getAreas(int i){
    double x1 = cellRange[i];
    double x2 = cellRange[i+1];
    double fx1 = circle(x1);
    double fx2 = circle(x2);
    if(getBin(fx1) != getBin(fx2)){
        x2 = circle(cellRange[vecDim - getBin(fx2)]);
        areas(x1,x2,i);
        x1 = x2;
        x2 = cellRange[i+1];
        fx1 = circle(x1);
        fx2 = circle(x2);
        }
    areas(x1,x2,i);
}




void Disk::makeTables(){
    for(int i=0; i<m_maxX; i++)
        getAreas(i);
    for(int x=m_maxX; x<vecDim; x++){
        for(int y=0; y<vecDim; y++){
            if (probMap.count(make_pair(y,x)))
                probMap[make_pair(x,y)] = probMap[make_pair(y,x)];
        }
    }
    for(int x=0; x<vecDim; x++){
        for(int y=0; y<vecDim; y++){
            if (probMap.count(make_pair(x,y))){
                double prob = probMap[make_pair(x,y)]/totalArea;
                probMap[make_pair(x,y)] = prob;
                probMap[make_pair(x*-1,y*-1)] += prob;
                probMap[make_pair(x*-1,y)] += prob;
                probMap[make_pair(x,y*-1)] += prob;
            }
        }
    }
    makeVectors();
    makeAliasTable();
}

int Disk::disperse(int x, int y, uint64_t u, int maxX, int maxY){
    xyCoord xy = coordVec[xyTable(u)];
    x += xy.first;
    y += xy.second;
    if (x>=0 && x < maxX && y >= 0 && y < maxY)
        return xy2i(x,y,maxX,maxY);
    return -1;

}


void Disk::makeVectors(){
    for (map<pair<int,int>,double>::iterator it = probMap.begin(); it != probMap.end(); ++it){
        coordVec.push_back(it->first);
        probVec.push_back(it->second);
    }
}

void Disk::makeAliasTable(){
    xyTable.create(probVec.begin(),probVec.end());
}

void Disk::printTables(){
    for (map<pair<int,int>,double>::iterator it = probMap.begin(); it != probMap.end(); ++it)
            cout << "X: " << it->first.first<< " Y: " << it->first.second << " Prob: " << it->second << endl;
}

