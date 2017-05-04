/*
 * Read in clustering result and validate
 *
 * TODO:
 * -
 *
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/cmds/validate.cpp $
 */

#include <iostream>
#include <cstdlib>
#include <fstream>

#include "src/common/validation.h"
#include "src/common/point.h"
#include "src/common/cluster.h"

using namespace std;

void usage(char, bool);
void parseParams(int, const char**, string &, unsigned int &, MetricType &);

int main(int argc, const char* argv[])
{
    string filename;
    unsigned int vDimension = 2;
    MetricType metric = EUCLIDEAN;
    vector<Cluster> _clusterSet;
    vector<Point> _pointSet;

    parseParams(argc, argv, filename, vDimension, metric);

    Point::setDimension(vDimension);
    Point::metric.setPreferredMetric(metric);

    //read file
    string line;
    ifstream inFile(filename.c_str(), ios::in);

    if(!inFile.is_open()){
        cout<<"Could not open file \""<<filename<<"\""<<endl;
        exit(1);
    }

    int groupNumber = -1;

    while(getline(inFile, line))
    {
        if(line.compare(0, 7, string("CLUSTER")) == 0)
        {
                ++groupNumber;
        }
        else if(line.length() > 0)
        {
            Point p = Point(line);
            p.setGroupNo(groupNumber);
            _pointSet.push_back(p);
        }
    }

    for(int i=0; i<=groupNumber; ++i)
    {
        Cluster c = Cluster();

        for(unsigned int j=0; j<_pointSet.size(); ++j)
        {
            if(_pointSet.at(j).groupNo() == i)
            {
                c.addPoint(&_pointSet.at(j));
            }
        }

        if(c.pointCount() > 0)
        {
            _clusterSet.push_back(c);
        }
    }

    inFile.close();

    Point::metric.setData(_pointSet);

    Matrix u(_pointSet.size(), _clusterSet.size());

    for(unsigned int i=0; i<_pointSet.size(); ++i)
    {
        u[i][_pointSet.at(i).groupNo()] = 1;
    }

    Matrix ut = u.transposed();

    for(unsigned int i=0; i<ut.rows(); ++i)
    {
        for(unsigned int j=0; j<ut.cols(); ++j)
        {
            cout<<ut[i][j]<<",";
        }

        cout<<endl;
    }

    cout<<"Xie-Beni index(-): "<<Validation::xieBeniIndex(_clusterSet, _pointSet, u, 2)<<endl;
    cout<<"Fukuyama-Sugeno index(-): "<<Validation::fukuyamaSugenoIndex(_clusterSet, _pointSet, u, 2)<<endl;
    cout<<"Kwon index(-): "<<Validation::kwonIndex(_clusterSet, _pointSet, u, 2)<<endl;
    cout<<"CWB index*(-): "<<Validation::composeWithinBetweenIndex(_clusterSet, _pointSet, u, 2)<<endl;
    cout<<"Zahid index*(+): "<<Validation::zahidIndex(_clusterSet, _pointSet, u, 2)<<endl;
    //cout<<"Fuzzy Hypervolume*(-): "<<Validation::fuzzyHypervolume(_clusterSet, _pointSet, u, 2)<<endl;
    cout<<"PBM index*(+): "<<Validation::PBMIndex(_clusterSet, _pointSet, u)<<endl;
    cout<<"Silhouette index(+): "<<Validation::silhouette(_clusterSet, _pointSet)<<endl;

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
        <<"-h (display help)"<<endl<<endl;

    exit(1);
}

void parseParams(int argc, const char* argv[], string &filename, unsigned int &vDimension, MetricType &metric)
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
