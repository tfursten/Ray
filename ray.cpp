#include "ray.h"

template<typename T>
    inline static std::pair<T,int> round_up(T x) {
        T y = static_cast<T>(2);
        int k = 1;
        for(;y < x;y*=2,++k)
            /*noop*/;
        return std::make_pair(y,k);
    }

void Ray::initialize(double s, double r)
{
    sigma = s;
    range = r;
    xsigma = range*sigma;
    vlen = floor(xsigma);
    if(xsigma>vlen+0.5)
        vlen += 1;
    vlen = vlen*2 + 3;
    vlen = round_up(vlen).first - 1;
    coordVec.resize(vlen);
    probVec.resize(vlen,0.0);
    xsigma = (vlen-3)/2.0 + 0.5;
    tcount = 0;
    int tail = int(xsigma+0.5);
    ltail = -tail;
    rtail = tail;
    makeTables();
}


double Ray::cdf(double x){
    return 0.5*(1+erf(x/(sigma*sqrt(2.0))));
}

double Ray::getProb(double x1, double x2){
    assert(x2>x1);
    return cdf(x2)-cdf(x1);
}


void Ray::makeTables(){
    double x1 = 0;
    double x2 = 0.5;
    int i=rtail;
    int j=i;
    for(;j>0;i++,j--){
        double p = getProb(x1,x2);
        coordVec[j]=-floor(x2);
        coordVec[i]=floor(x2);
        probVec[i]+=p;
        probVec[j]+=p;
        x1 = x2;
        x2 = x1+1.0;
    }
    double ptail = 1.0-cdf(xsigma);
    coordVec[0] = ltail;
    coordVec[vlen-1] = rtail;
    probVec[0]=ptail;
    probVec[vlen-1]=ptail;
    makeAliasTable();
}

int Ray::tail(xorshift64& rand){
    double x;
    for(;;){
        x = (xsigma -log(rand.get_double52()))/double(xsigma);
        if(rand.get_double52()<exp(-0.5*(fabs(x)-xsigma)*(fabs(x)-xsigma)))
            break;
    }
    x = floor(x+0.5);
    return int(x);
}


int Ray::disperse(xorshift64& rand){
    uint64_t u = rand.get_uint64();
    int c = coordVec[xyTable(u)];
    if(c==ltail){
        tcount += 1;
        return -tail(rand);
    }
    if(c==rtail){
        tcount += 1;
        return tail(rand);
    }
    return c;
}


void Ray::makeAliasTable(){
    xyTable.create(probVec.begin(),probVec.end());
}



void Ray::printTables(){
    for (int i = 0; i<vlen; i++){
            cout << "X: " << coordVec[i] << " Prob: " << probVec[i] << endl;
            }
}

void Ray::tailStats(int it){
    cout << "Expected Tail: " << probVec[0]*2*it << endl;
    cout << "Current Tail: " << tcount << endl;
}

