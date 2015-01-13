#include <iostream>
#include <fstream>
#include "xorshift64.h"
#include "ray.h"

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
    ofstream myfile;
    myfile.open("test.txt");
    xorshift64 myrand;
    myrand.seed(2348585);
    Ray rayl;
    Ray rays;
    //double sigma = 2.0;
    rayl.initialize(.25,3);
    rayl.printTables();
    rays.initialize(.25,3);
    rays.printTables();
    int it = 1000000;
    for(int iii=0; iii<it; iii++)
    {
        int X = rayl.disperse(myrand);
        int Y = rays.disperse(myrand);
        double dX = 25+X;
        double dY = 25+Y;
        int newX = wrap_around(static_cast<int>(dX), 100);
        int newY = wrap_around(static_cast<int>(dY), 100);
        myfile << newX << "\t" << newY << endl;
    }
    rayl.tailStats(it);
    rays.tailStats(it);
    myfile.close();

    return 0;
}
