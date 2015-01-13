#include "ring.h"


void Ring::initialize(double r, double c)
{
    radius = r;
    rSq = r*r;
    halfRsq = rSq/2.0;
    circum = 2*M_PI*r;
    pCenter = c*circum;
    tProb = pCenter + circum;
    getCellRange();
    vecDim = cellRange.size()-1;
    m_maxX = pol2xy(M_PI/4.0).first;
    getPoints();
    makeTables();

}

double Ring::arclen(double x1, double x2, double y1, double y2)
{
    double t1 = xy2Pol(x1,y1).second;
    double t2 = xy2Pol(x2,y2).second;
    return abs(t1-t2)*radius;
}

void Ring::getCellRange()
{
    cellRange.push_back(0.0);
    for(int i=0; i < radius; i++)
    {
        double x = i + 0.5;
        if (x<radius){
            cellRange.push_back(x);
        }
    }
    cellRange.push_back(radius);
}


pair<double,double> Ring::xy2Pol(double x, double y)
{
    double r = sqrt(x*x+y*y);
    double t = atan(double(y)/x);
    return make_pair(r,t);
}

pair<double,double> Ring::pol2xy(double theta)
{
    double x = radius*cos(theta);
    double y = radius*sin(theta);
    return make_pair(x,y);
}

double Ring::circle(double &x)
{
    return (double)sqrt(radius*radius-x*x);

}

void Ring::getPoints(){
    double x1 = cellRange[0];
    double y1 = circle(x1);
    int yyy = 1;
    points.push_back(make_pair(x1,y1));
    for(int i=1; i<=vecDim; i++){
        double x2 = cellRange[i];
        double y2 = circle(x2);
        if(getBin(y1) != getBin(y2)){
            x2 = circle(cellRange[vecDim-yyy]);
            y2 = cellRange[vecDim-yyy];
            yyy+=1;
        }
        if(x2>m_maxX){
            x2 = m_maxX;
            y2 = circle(m_maxX);
            points.push_back(make_pair(x2,y2));
            break;
        }
        points.push_back(make_pair(x2,y2));
        if(x2<cellRange[i]){
            x2 = cellRange[i];
            y2 = circle(x2);
            points.push_back(make_pair(x2,y2));
            }
        x1 = x2;
        y1 = y2;
    }
    cout << "pointsize: "<< points.size()<<endl;
    for(int i=0;i<points.size();i++)
        cout << points[i].first << " " << points[i].second << endl;
}



int Ring::getBin(double x){
    for (int i=0; i<vecDim; i++){
        if(x <= cellRange[i+1])
            return i;
    }
}


void Ring::length(double x1, double x2, double y1, double y2){
    double L = arclen(x1,x2,y1,y2);
    int x = getBin(x2);
    int y = getBin(y1);
    probMap[make_pair(x,y)]+= L;
    probMap[make_pair(y,x)]+= L;
}

void Ring::getLength(){
    for(int i=0; i<points.size()-1; i++){
        double x1 = points[i].first;
        double y1 = points[i].second;
        double x2 = points[i+1].first;
        double y2 = points[i+1].second;
        length(x1,x2,y1,y2);
    }
}


void Ring::makeTables(){
    getLength();
    for(int x=0; x<vecDim; x++){
        for(int y=0; y<vecDim; y++){
            if (probMap.count(make_pair(x,y))){
                double prob = probMap[make_pair(x,y)]/tProb;
                probMap[make_pair(x,y)] = prob;
                probMap[make_pair(x*-1,y*-1)] += prob;
                probMap[make_pair(x*-1,y)] += prob;
                probMap[make_pair(x,y*-1)] += prob;
            }
        }
    }
    if(pCenter)
        probMap[make_pair(0,0)] = pCenter/tProb;
    makeVectors();
    makeAliasTable();
}

xyCoord Ring::disperse(uint64_t u){
    return coordVec[xyTable(u)];
}

void Ring::makeVectors(){
    for (map<pair<int,int>,double>::iterator it = probMap.begin(); it != probMap.end(); ++it){
        coordVec.push_back(it->first);
        probVec.push_back(it->second);
    }
}

void Ring::makeAliasTable(){
    xyTable.create(probVec.begin(),probVec.end());
}

void Ring::printTables(){
    double sigma = 0.0;
    for (map<pair<int,int>,double>::iterator it = probMap.begin(); it != probMap.end(); ++it){
            double x = it->first.first;
            double y = it->first.second;
            double d = x*x+y*y;
            double p = it->second;
            sigma += p*d;
            cout << "X: " << it->first.first<< " Y: " << it->first.second << " Prob: " << it->second << endl;
            }
    cout << "sigma: " << sigma << endl;
}

