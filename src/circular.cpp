/*
 * Generates synthetic data in two concentric, normally distributed clusters
 *
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/circular.cpp $
 */

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>

using namespace std;

void randCirc(double &x, double &y, double rc, double rc_range, double x1=0.0, double y1=0.0)
{
    double a= 2 * M_PI * ((rand()%100000) / 100000.0);
    double rc_range_up = rc_range * 100000;
    double r = sqrt(((rand()%(100000 - (int)rc_range_up)) / 100000.0) + (rc_range_up / 100000));
    x = (rc * r) * cos(a) + x1;
    y = (rc * r) * sin(a) + y1;
}

int main(int argc, char** argv)
{
    if(argc != 7)
    {
        cerr<<"ARGS: [outside radius] [outside tightness (0.0 to 1.0)] [outside count] [inside radius] [inside tightness (0.0 to 1.0)] [inside count]"<<endl;
        exit(1);
    }

    srand(time(NULL));
    double x, y;

    for(int i=0; i<atoi(argv[3]); ++i)
    {
        randCirc(x, y, atof(argv[1]), atof(argv[2]));
        cout<<x<<" "<<y<<endl;
    }

    for(int i=0; i<atoi(argv[6]); ++i)
    {
        randCirc(x, y, atof(argv[4]), atof(argv[5]));
        cout<<x<<" "<<y<<endl;
    }

    return 0;
}
