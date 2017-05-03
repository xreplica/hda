/*
 * Generates synthetic data in a number of circular groupings
 *
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/spots.cpp $
 */

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>

using namespace std;

void randCirc(double &x, double &y, double rc, double x1=0.0, double y1=0.0)
{
    double a= 2 * M_PI * ((rand()%100000) / 100000.0);
    double r;// = sqrt(rand()%100000 / 100000.0);

    bool valid = false;
    while(!valid)
    {
        double v1 = rand()%100000;// / 100000.0;
        bool use = true;;
        for(int i=0; i<10; ++i)
        {
            if(v1 > rand()%100000)
            {
                use=false;
                break;
            }
        }

        if(use)
        {
            r = sqrt(v1/100000.0);
            valid = true;
        }
    }

    x = (rc * r) * cos(a) + x1;
    y = (rc * r) * sin(a) + y1;
}

int main(int argc, char** argv)
{
    if(argc != 5)
    {
        cerr<<"ARGS: [area size] [number of spots] [outside radius] [outside count]"<<endl;
        exit(1);
    }

    srand(time(NULL));
    double x, y;

    for(int j=0; j<atoi(argv[2]); ++j)
    {
        double xCenter = rand()%atoi(argv[1]);
        double yCenter = rand()%atoi(argv[1]);
        for(int i=0; i<atoi(argv[4]); ++i)
        {
            randCirc(x, y, atof(argv[3]), xCenter, yCenter);
            cout<<x<<" "<<y<<endl;
        }
    }

    return 0;
}
