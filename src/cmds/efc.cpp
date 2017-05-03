/*
 * entropy-based fuzzy clustering implementation
 *
 * TODO:
 * -
 *
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/cmds/efc.cpp $
 */

#include <iostream>
#include <cstdlib>

#include "src/algorithms/fuzzyentropy.h"
#include "src/core/validation.h"

using namespace std;

void usage(char, bool);
void parseParams(int, const char**, string &, unsigned int &, bool &, bool &, bool &, MetricType &, double &, double &);

int main(int argc, const char* argv[])
{
    string filename;
    unsigned int vDimension;
    MetricType metric = EUCLIDEAN;
    bool displayClusters = false;
    bool displayIndexes = false;
    bool displayU = false;
    double beta = 0.7;
    double gamma = 5;

    parseParams(argc, argv, filename, vDimension, displayClusters, displayU, displayIndexes, metric, beta, gamma);

    Point::setDimension(vDimension);
    Point::metric.setPreferredMetric(metric);

    FuzzyEntropy efc(filename);

    efc.initialize(beta, gamma);
    efc.run();

    if(displayU)
    {
        cout<<endl<<"Partition matrix:"<<endl;
        Matrix u = efc.U().transposed();
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
        efc.printClusters();
    }

    if(displayIndexes)
    {
        cout<<"Xie-Beni index(-): "<<Validation::xieBeniIndex(efc.clusterSet(), efc.pointSet(), efc.U(), 2)<<endl;
        cout<<"Fukuyama-Sugeno index(-): "<<Validation::fukuyamaSugenoIndex(efc.clusterSet(), efc.pointSet(), efc.U(), 2)<<endl;
        cout<<"Kwon index(-): "<<Validation::kwonIndex(efc.clusterSet(), efc.pointSet(), efc.U(), 2)<<endl;
        cout<<"CWB index*(-): "<<Validation::composeWithinBetweenIndex(efc.clusterSet(), efc.pointSet(), efc.U(), 2)<<endl;
        cout<<"Zahid index*(+): "<<Validation::zahidIndex(efc.clusterSet(), efc.pointSet(), efc.U(), 2)<<endl;
        cout<<"Fuzzy Hypervolume*(-): "<<Validation::fuzzyHypervolume(efc.clusterSet(), efc.pointSet(), efc.U(), 2)<<endl;
        cout<<"PBM index*(+): "<<Validation::PBMIndex(efc.clusterSet(), efc.pointSet(), efc.U())<<endl;
        cout<<"Silhouette index(+): "<<Validation::silhouette(efc.clusterSet(), efc.pointSet())<<endl;
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
        <<"-b <beta threshold>"<<endl
        <<"-y <gamma threshold>"<<endl
        <<"  optional flags:"<<endl
        <<"-d <distance metric> 0=EUCLIDEAN, 1=MAHALANOBIS, 2=(gaussian)KERNEL"<<endl
        <<"-u (display partition matrix)"<<endl
        <<"-k (display clusters)"<<endl
        <<"-i (display vlidity indexes)"<<endl
        <<"-h (display help)"<<endl<<endl;

    exit(1);
}

void parseParams(int argc, const char* argv[], string &filename, unsigned int &vDimension, bool &displayClusters, bool &displayIndexes, bool &displayU, MetricType &metric, double &beta, double &gamma)
{
    bool fProvided = false, dProvided = false, bProvided = false, yProvided = false;

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
            case 'b':
                if(argv[i+1][0] != '-')
                {
                    beta = atof(argv[i+1]);
                    ++i;
                    bProvided = true;
                }
                else
                {
                    usage(argv[i][1], true);
                }
                break;
            case 'y':
                if(argv[i+1][0] != '-')
                {
                    gamma = atof(argv[i+1]);
                    ++i;
                    yProvided = true;
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

    if(!(fProvided && dProvided && bProvided && yProvided))
    {
        usage(0);
    }
}
