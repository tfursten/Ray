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
    int rtail, ltail, tcount, vlen;
    pMap probMap;
    vector<int> coordVec;
    vector<double> probVec;
    alias_table xyTable;
    double cdf(double x);
    double getProb(double x1, double x2);
    void makeTables();
    void makeVectors();
    void makeAliasTable();
    int tail(xorshift64& rand);



public:
    void initialize(double s, double r);
    void printTables();
    int disperse(xorshift64& rand);
    void tailStats(int it);

};


#endif // RAY_H_INCLUDED
