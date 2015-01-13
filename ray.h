#ifndef RAY_H_INCLUDED
#define RAY_H_INCLUDED

#include <cmath>
#include <vector>
#include <map>
#include <iostream>
#include "aliastable.h"
#include "xorshift64.h"


using namespace std;
typedef map<int, double> pMap;


class Ray
{
private:
	double sigma, range, xsigma;
    int tail, tcount;
    pMap probMap;
    vector<int> coordVec;
    vector<double> probVec;
    alias_table xyTable;
    double cdf(double x);
    double cdfRange(double x1, double x2);
    void getProb(double x1, double x2, int i);
    void makeTables();
    void makeVectors();
    void makeAliasTable();




public:
    void initialize(double s, double r);
    void printTables();
    int disperse(xorshift64& rand);
    void tailStats(int it);

};


#endif // RAY_H_INCLUDED
