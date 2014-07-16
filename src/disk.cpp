#include "disk.h"

void Disk::initialize(double r)
{
    radius = r;
    rSq = r*r;
    halfRsq = rSq/2.0;
    getCellRange();
    vecDim = cellRange.size()-1;
    maxX = getMaxX();
    totalArea = M_PI*rSq;
    cout << "radius: " << r << endl;
    cout << "dim: " << vecDim << endl;
    cout << "MaxX: " << maxX << endl;
    cout << "Total Area: " << totalArea << endl;

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
    cout << "X1,X2: " << x1 << " " << x2 << endl;
    cout << "A: " << A << endl;
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
    cout <<"Start:"<< x1 << " " << x2 << endl;
    if(getBin(fx1) != getBin(fx2)){
        cout << "INSIDE" << endl;
        cout << "bin1: " << getBin(fx1) << endl;
        cout << "bin2: " << getBin(fx2) << endl;
        cout << "vecDim: " << vecDim << endl;
        cout << cellRange[getBin(fx1)] << endl;
        cout << "CellRange: " ;
        for(int i=0; i<=vecDim; i++)
            cout << cellRange[i] << " ";
        cout << endl;
        x2 = circle(cellRange[getBin(fx1)]);
        cout << "X2 AFter: " << x2 << endl;
        areas(x1,x2,i);
        x1 = x2;
        x2 = cellRange[i+1];
        fx1 = circle(x1);
        fx2 = circle(x2);
        }
    areas(x1,x2,i);
}




void Disk::makeTables(){
    for(int i=0; i<maxX; i++)
        getAreas(i);
    for(int x=maxX; x<vecDim; x++){
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

xyCoord Disk::disperse(uint64_t u){
    return coordVec[xyTable(u)];
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

