#include <iostream>
#include <fstream>
#include "xorshift64.h"
#include "ray.h"
#include "rnormal.h"

using namespace std;
typedef pair<int,int> xyCoord;

inline int xy2i(int x, int y, int mx, int my) {
	return x*my+y;
}


int wrap_around(int x, int w) {
	return ((x % w) + w) % w;
}
int main()
{
    ofstream myfile, myfile2;
    myfile.open("test.txt");
    myfile.open("test2.txt");
    xorshift64 myrand;
    myrand.seed();
    Ray ray;
    //Ray ray2;
    double sigma = 3;
    ray.initialize(sigma,3);
    ray.printTables();
    //ray2.initialize(sigma+1,6);
    //ray2.printTables();
    int it = 10000000;
    map<int,int> dis;
    map<int,int> cont;
    myfile2 << "val\tcont\tdisc\n";
    for(int iii=0; iii<it; iii++)
    {
        int X = ray.disperse(myrand);
        dis[X] += 1;
        //int Y = ray2.disperse(myrand);
        int Y = floor(rand_normal(myrand, 0.0,2)+0.5);
        cont[Y] += 1;
        int dX = 25+X;
        int dY = 25+Y;
        int newX = wrap_around(static_cast<int>(dX), 100);
        int newY = wrap_around(static_cast<int>(dY), 100);
        myfile << newX << "\t" << newY << endl;

    }
    for(int i=-sigma*12; i<sigma*12; i++){
        double d = dis[i];
        double c = cont[i];
        myfile2 <<i<<"\t"<< c << "\t" << d << endl;
    }
    ray.tailStats(it);
    cout << myrand.get_count() << endl;
    myfile.close();
    myfile2.close();

    return 0;
}
