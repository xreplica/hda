/*
 * CLI parameters handling implementation
 *
 * TODO:
 * -
 *
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/common/params.cpp $
 */

#ifndef __PARAMS_CPP__
#define __PARAMS_CPP__

#include <iostream>

#include "src/common/params.h"

using namespace std;

//CONSTRUCTORS

Params::Params() : vDimension(2), mergeMethod(1), mergeParameter(0.0), stopIterations(7), metric(EUCLIDEAN), displayClusters(false), displayIndices(false), displayU(false), noSplit(false), verbose(false)
{};

//FUNCTIONS

void Params::parse(int argc, const char* argv[])
{
    bool fProvided = false, dProvided=false;

    for(int i=1; i<argc; ++i)
    {
        if(argv[i][0] == '-')
        {
            switch(argv[i][1])
            {
            case 'f':
                if(argv[i+1][0] != '-')
                {
                    filename = argv[i+1];
                    ++i;
                    fProvided = true;
                }
                else
                    usage(argv[i][1], true);
                break;
            case 'd':
                if(argv[i+1][0] != '-')
                {
                    vDimension = atoi(argv[i+1]);
                    ++i;
                    dProvided = true;
                }
                else
                    usage(argv[i][1], true);
                break;
            case 'm':
                if(argv[i+1][0] != '-')
                {
                    mergeMethod = atoi(argv[i+1]);
                    ++i;
                }
                else
                    usage(argv[i][1], true);
                break;
            case 'l':
                if(argv[i+1][0] != '-')
                {
                    mergeParameter = atof(argv[i+1]);
                    ++i;
                }
                else
                    usage(argv[i][1], true);
                break;
            case 'D':
                if(argv[i+1][0] != '-')
                {
                    metric = (MetricType)atoi(argv[i+1]);
                    ++i;
                }
                else
                    usage(argv[i][1], true);
                break;
            case 's':
                if(argv[i+1][0] != '-')
                {
                    stopIterations = atoi(argv[i+1]);
                    ++i;
                }
                else
                    usage(argv[i][1], true);
                break;
            case 'u':
                displayU = true;
                break;
            case 'k':
                displayClusters = true;
                break;
            case 'i':
                displayIndices = true;
                break;
            case 'M':
                noSplit = true;
                break;
            case 'v':
                verbose = true;
                break;
            case 'h':
                usage(0);
                break;
            default:
                usage(argv[i][1]);
            }
        }
    }

    if(!(fProvided && dProvided))
        usage(0);

    if(mergeMethod == 1 && mergeParameter == 0.0)
        mergeParameter = 2.0;
    else if(mergeMethod == 2 && mergeParameter == 0.0)
        mergeParameter = 0.5;
    else if(mergeMethod == 3 && mergeParameter == 0.0)
        mergeParameter = 5.0;
    else if(mergeMethod == 4 && mergeParameter == 0.0)
        mergeParameter = 2.0;
    else if(mergeMethod == 5 && mergeParameter == 0.0)
        mergeParameter = 2.0;
}

void Params::usage(char invFlag, bool missingVal)
{
    if(missingVal)
        cout<<"Missing value for flag : "<<invFlag<<endl<<endl;
    else if(invFlag != 0)
        cout<<"Unknown flag: "<<invFlag<<endl<<endl;

    cout<<"USAGE: hda [options]"<<endl
        <<"  required flags:"<<endl
        <<"-f <input file name>"<<endl
        <<"-d <vector dimension>"<<endl
        <<"  optional flags:"<<endl
        <<"-m <merging method> 0=None, 1=NN(Nearest Neighbor), 2=K-means+(Chebychev), 3=x% within 5 sigma, 4=5%NN, 5=minimum nearest neighbor"<<endl
        <<"-l <merge parameter>"<<endl
        <<"-s <iterations for stop>"<<endl
        <<"-D <distance metric> 0=EUCLIDIAN, 2=(gaussian)KERNEL"<<endl
        <<"-M (no split/merge only, requires input be splitted"<<endl
        <<"-u (display partition matrix)"<<endl
        <<"-k (display clusters)"<<endl
        <<"-i (display vlidity indexes)"<<endl
        <<"-v (set verbose output)"<<endl
        <<"-h (display help)"<<endl<<endl;

    exit(1);
}

#endif
