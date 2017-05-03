/*
 * Generates synthetic data in two mirrored and interwoven "C" shapes
 *
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/double_c.cpp $
 */

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>

using namespace std;

void randCircLeft(double &x, double &y, double rc, double rc_range, double x1=0.0, double y1=0.0)
{
    double a= 2 * M_PI * ((rand()%60000) / 100000.0 + 0.20);
    double rc_range_up = rc_range * 100000;
    double r = sqrt(((rand()%(100000 - (int)rc_range_up)) / 100000.0) + (rc_range_up / 100000));
    x = ((rc * r) * cos(a) + x1) * 1.1;
    y = ((rc * r) * sin(a) + y1) * 0.9;
}

int main(int argc, char** argv)
{
    double radius = 5.0;
    double tightness = 0.5;
    int count = 250;

    srand(time(NULL));
    double x, y;

    for(int i=0; i<count; ++i)
    {
        randCircLeft(x, y, radius, tightness, 0.0, 2.2);
        cout<<x<<" "<<y<<endl;
    }

    for(int i=0; i<count; ++i)
    {
        randCircLeft(x, y, radius, tightness, 0.0, -2.2);
        cout<<-x<<" "<<y<<endl;
    }

    return 0;
}
