/*
 * Fuzzy c-means+ implementation
 *
 * TODO:
 * -
 *
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/cmds/fcmp.cpp $
 */

#include <iostream>

#include "src/algorithms/fuzzycmeans.h"
#include "src/common/validation.h"

using namespace std;

void usage(char, bool);
void parseParams(int, const char**, string &, unsigned int &, unsigned int &, unsigned int &, bool &, bool &, bool &, MetricType &);

int main(int argc, const char* argv[])
{
    string filename;
    unsigned int vDimension;
    MetricType metric = EUCLIDEAN;
    bool displayClusters = false;
    bool displayIndexes = false;
    bool displayU = false;
    unsigned int fuzzyConst = 2;
    unsigned int numClust = 3;

    parseParams(argc, argv, filename, fuzzyConst, vDimension, numClust, displayClusters, displayIndexes, displayU, metric);

    Point::setDimension(vDimension);
    Point::metric.setPreferredMetric(metric);

    FuzzyCMeans fcm(filename, numClust, fuzzyConst);
    fcm.initialize();

    bool stable = false;
    while(!stable)
    {
        fcm.run();
        stable = ! fcm.fuzzyMerge();
        if(stable){
            stable = ! fcm.fuzzySplit();
        }
        //stable = true;
    }

    if(displayU)
    {
        cout<<endl<<"Partition matrix:"<<endl;
        Matrix u = fcm.U().transposed();
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
        cout<<"Clusters:"<<endl<<endl;
        fcm.printClusters();
    }

    if(displayIndexes)
    {
        cout<<"Xie-Beni index(-): "<<Validation::xieBeniIndex(fcm.clusterSet(), fcm.pointSet(), fcm.U(), fuzzyConst)<<endl;
        cout<<"Fukuyama-Sugeno index(-): "<<Validation::fukuyamaSugenoIndex(fcm.clusterSet(), fcm.pointSet(), fcm.U(), fuzzyConst)<<endl;
        cout<<"Kwon index(-): "<<Validation::kwonIndex(fcm.clusterSet(), fcm.pointSet(), fcm.U(), fuzzyConst)<<endl;
        cout<<"CWB index*(-): "<<Validation::composeWithinBetweenIndex(fcm.clusterSet(), fcm.pointSet(), fcm.U(), fuzzyConst)<<endl;
        cout<<"Zahid index*(+): "<<Validation::zahidIndex(fcm.clusterSet(), fcm.pointSet(), fcm.U(), fuzzyConst)<<endl;
        cout<<"Fuzzy Hypervolume*(-): "<<Validation::fuzzyHypervolume(fcm.clusterSet(), fcm.pointSet(), fcm.U(), fuzzyConst)<<endl;
        cout<<"PBM index*(+): "<<Validation::PBMIndex(fcm.clusterSet(), fcm.pointSet(), fcm.U())<<endl;
        cout<<"Silhouette index(+): "<<Validation::silhouette(fcm.clusterSet(), fcm.pointSet(), fcm.U())<<endl;
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
        <<"-m <fuzzification constant>"<<endl
        <<"-c <initial number of clusters>"<<endl
        <<"-u (display fuzzy partition)"<<endl
        <<"-k (display cluster center coordinates)"<<endl
        <<"-i (display vlidity indexes)"<<endl
        <<"-h (display help)"<<endl<<endl;

    exit(1);
}

void parseParams(int argc, const char* argv[], string &filename, unsigned int &fuzzyConst, unsigned int &vDimension, unsigned int &numClust, bool &displayClusters, bool &displayIndexes, bool &displayU, MetricType &metric)
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
            case 'm':
                if(argv[i+1][0] != '-')
                {
                    fuzzyConst = atoi(argv[i+1]);
                    ++i;
                }
                else
                {
                    usage(argv[i][1], true);
                }
                break;
            case 'c':
                if(argv[i+1][0] != '-')
                {
                    numClust = atoi(argv[i+1]);
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
                displayClusters = true;
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
