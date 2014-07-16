#ifndef DISK_H_INCLUDED
#define DISK_H_INCLUDED

#include <math.h>
#include <vector>
#include <map>
#include <iostream>
#include "aliastable.h"

using namespace std;
typedef pair<int,int> xyCoord;
typedef map<xyCoord, double> pMap;



class Disk
{
private:
    double radius;
    double rSq;
    double halfRsq;
    pMap probMap;
    vector<double> cellRange;
    double m_maxX;
    double totalArea;
    int vecDim;
    vector<xyCoord> coordVec;
    vector<double> probVec;
    alias_table xyTable;

    double integrate(double x1, double x2);
    void getCellRange();
    double getMaxX();
    pair<double,double> pol2xy(double theta);
    double circle(double &x);
    int getBin(double x);
    void getAreas(int i);
    void areas(double x1, double x2, int i);
    void makeVectors();
    void makeAliasTable();



public:
    void initialize(double r);
    void makeTables();
    void printTables();
    int disperse(uint64_t u);

};

#endif // DISK_H_INCLUDED
