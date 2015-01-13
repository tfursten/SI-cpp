#ifndef RING_H_INCLUDED
#define RING_H_INCLUDED

#include <cmath>
#include <vector>
#include <map>
#include <iostream>
#include "aliastable.h"

using namespace std;
typedef pair<double,double> pts;
typedef pair<int,int> xyCoord;
typedef map<xyCoord, double> pMap;



class Ring
{
private:
    double radius;
    double rSq;
    double halfRsq;
    double pCenter;
    pMap probMap;
    vector<double> cellRange;
    double circum;
    double tProb;
    int vecDim;
    double m_maxX;
    vector<pts> points;
    vector<xyCoord> coordVec;
    vector<double> probVec;
    alias_table xyTable;
    double arclen(double x1, double x2, double y1, double y2);
    void getCellRange();
    pair<double,double> pol2xy(double theta);
    pair<double,double> xy2Pol(double x, double y);
    double circle(double &x);
    int getBin(double x);
    void getLength();
    void length(double x1, double x2, double y1, double y2);
    void makeVectors();
    void makeAliasTable();
    void getPoints();



public:
    void initialize(double r, double c);
    void makeTables();
    void printTables();
    xyCoord disperse(uint64_t u);

};

#endif // RING_H_INCLUDED
