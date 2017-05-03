/*
 * heuristic DIANA implementation
 *
 * TODO:
 * -
 *
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/cmds/ehda.cpp $
 */

#include <iostream>
#include <cstdlib>

#include "src/algorithms/ehdiana.h"
#include "src/core/validation.h"

using namespace std;

void usage(char, bool);
void parseParams(int, const char**, string &, unsigned int &, bool &, bool &, bool &, MetricType &);

int main(int argc, const char* argv[])
{
    string filename;
    unsigned int vDimension;
    MetricType metric = EUCLIDEAN;
    bool displayClusters = false;
    bool displayIndexes = false;
    bool displayU = false;

    parseParams(argc, argv, filename, vDimension, displayClusters, displayIndexes, displayU, metric);

    Point::setDimension(vDimension);
    Point::metric.setPreferredMetric(metric);

    EHDiana ehda(filename);

    try{
        ehda.initialize();
        ehda.run();
    }
    catch(matrixException &e)
    {
        cerr<<e.what()<<endl;
        exit(1);
    }

    if(displayU)
    {
        cout<<endl<<"Partition matrix:"<<endl;
        Matrix u = ehda.U().transposed();
        for(unsigned int i=0; i<u.rows();++i)
        {
            for(unsigned int j=0; j<u.cols();++j)
            {
                cout<<u[i][j]<<",";
            }
            cout<<endl;
        }
    }

    if(displayClusters)
    {
        ehda.printClusters();
    }

    if(displayIndexes)
    {
        cout<<"Xie-Beni index(-): "<<Validation::xieBeniIndex(ehda.clusterSet(), ehda.pointSet(), ehda.U(), 2)<<endl;
        cout<<"Fukuyama-Sugeno index(-): "<<Validation::fukuyamaSugenoIndex(ehda.clusterSet(), ehda.pointSet(), ehda.U(), 2)<<endl;
        cout<<"Kwon index(-): "<<Validation::kwonIndex(ehda.clusterSet(), ehda.pointSet(), ehda.U(), 2)<<endl;
        cout<<"CWB index*(-): "<<Validation::composeWithinBetweenIndex(ehda.clusterSet(), ehda.pointSet(), ehda.U(), 2)<<endl;
        cout<<"Zahid index*(+): "<<Validation::zahidIndex(ehda.clusterSet(), ehda.pointSet(), ehda.U(), 2)<<endl;
        cout<<"Fuzzy Hypervolume*(-): "<<Validation::fuzzyHypervolume(ehda.clusterSet(), ehda.pointSet(), ehda.U(), 2)<<endl;
        cout<<"PBM index*(+): "<<Validation::PBMIndex(ehda.clusterSet(), ehda.pointSet(), ehda.U())<<endl;
        cout<<"Silhouette index(+): "<<Validation::silhouette(ehda.clusterSet(), ehda.pointSet())<<endl;
    }

    return 0;
}

void usage(char invFlag = 0, bool missingVal = false)
{
    if(missingVal)
    {
        cout<<"Missing value for flag : "<<invFlag<<endl<<endl;
    }
    else if(invFlag != 0){
        cout<<"Unknown flag: "<<invFlag<<endl<<endl;
    }
    cout<<"USAGE: kmeansplus [options]"<<endl
        <<"  required flags:"<<endl
        <<"-f <input file name>"<<endl
        <<"-v <vector dimension>"<<endl
        <<"  optional flags:"<<endl
        <<"-d <distance metric> 0=EUCLIDEAN, 1=MAHALANOBIS, 2=(gaussian)KERNEL"<<endl
        <<"-u (display partition matrix)"<<endl
        <<"-k (display clusters)"<<endl
        <<"-i (display vlidity indexes)"<<endl
        <<"-h (display help)"<<endl<<endl;

    exit(1);
}

void parseParams(int argc, const char* argv[], string &filename, unsigned int &vDimension, bool &displayClusters, bool &displayIndexes, bool &displayU, MetricType &metric)
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
                {
                    usage(argv[i][1], true);
                }
                break;
            case 'v':
                if(argv[i+1][0] != '-')
                {
                    vDimension = atoi(argv[i+1]);
                    ++i;
                    dProvided = true;
                }
                else
                {
                    usage(argv[i][1], true);
                }
                break;
            case 'd':
                if(argv[i+1][0] != '-')
                {
                    metric = (MetricType)atoi(argv[i+1]);
                    ++i;
                }
                else
                {
                    usage(argv[i][1], true);
                }
                break;
            case 'u':
                displayU = true;
                break;
            case 'k':
                displayClusters = true;
                break;
            case 'i':
                displayIndexes = true;
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
    {
        usage(0);
    }
}
