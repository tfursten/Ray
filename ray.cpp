#include "ray.h"



void Ray::initialize(double s, double r)
{
    sigma = s;
    range = r;
    xsigma = range*sigma;
    tail = floor(xsigma)+1;
    tcount = 0;
    makeTables();
}


double Ray::cdf(double x){
    return 0.5*(1+erf(x/(sigma*sqrt(2.0))));
}

double Ray::cdfRange(double x1, double x2){
    assert(x2>x1);
    return cdf(x2)-cdf(x1);
}

void Ray::getProb(double x1, double x2, int i){
    double p = cdfRange(x1,x2);
    probMap[i]+= p;
    probMap[-i]+= p;
}

void Ray::makeTables(){
    double x1 = 0;
    double x2 = 0.5;
    int i=0;
    for(; x2<xsigma; i++){
        getProb(x1,x2,i);
        x1 = x2;
        x2 = x1+1.0;
    }
    getProb(x2-1.0,xsigma,i);
    double ptail = 1.0-cdf(xsigma);
    probMap[tail]=ptail*2.0;
    makeVectors();
    makeAliasTable();
}

int Ray::disperse(xorshift64& rand){
    uint64_t u = rand.get_uint64();
    int c = coordVec[xyTable(u)];
    if(c==tail){
        int64_t a = static_cast<int64_t>(u<<8);
        double aa = static_cast<double>(a);
        tcount += 1;
        double x,y;
        x = (xsigma + (-log(rand.get_double52())));
        x = floor(x+0.5);
        x = copysign(x,aa);
        return int(x);
    }
    return c;
}

void Ray::makeVectors(){
    for (map<int,double>::iterator it = probMap.begin(); it != probMap.end(); ++it){
        coordVec.push_back(it->first);
        probVec.push_back(it->second);
    }
}

void Ray::makeAliasTable(){
    xyTable.create(probVec.begin(),probVec.end());
}



void Ray::printTables(){
    for (map<int,double>::iterator it = probMap.begin(); it != probMap.end(); ++it){
            cout << "X: " << it->first << " Prob: " << it->second << endl;
            }
}

void Ray::tailStats(int it){
    cout << "Expected Tail: " << probMap[tail]*it*2 << endl;
    cout << "Current Tail: " << tcount << endl;
}

