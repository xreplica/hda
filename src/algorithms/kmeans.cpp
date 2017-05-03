/*
 * K-means+ clustering implementation
 *
 * TODO:
 * - fix all points of all clusters fall to same point
 *
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/algorithms/kmeans.cpp $
 */

#ifndef __KMEANS_CPP__
#define __KMEANS_CPP__

#include <set>

#include "src/algorithms/kmeans.h"

//CONSTRUCTORS

KMeans::KMeans(string f, unsigned int d) : _dimension(d), _file(f), _finalClusterNum(0) {}

//FUNCTIONS

void KMeans::printClusters() const
{
    for (unsigned int i=0; i<_clusterSet.size(); ++i)
    {
        cout<<"The size of cluster [ "<<i<<" ] is: "<<_clusterSet.at(i).points().size()<<endl;
        cout<<_clusterSet.at(i).containerToString()<<endl;
    }
}

void KMeans::run(unsigned int num)
{
    _initialNum = num;
    initialize();
    assign();
}

void KMeans::initialize()
{
    string line;
    ifstream inFile(_file.c_str(), ios::in);
    set<unsigned int> clusterIndices;

    if(!inFile.is_open()){
        cout<<"Could not open file \""<<_file<<"\""<<endl;
        exit(1);
    }

    while(getline(inFile, line))
    {
        Point p = Point(line);
        _pointSet.push_back(p);
    }

    //randomly select initial centroid seeds
    srand((unsigned int)time(NULL));

    while(clusterIndices.size()<_initialNum)
    {
        int index = rand() % _pointSet.size();
        if(clusterIndices.find(index) == clusterIndices.end())
            clusterIndices.insert(index);
    }

    cout<<"initial clusters: "<<endl;
    for(set<unsigned int>::iterator i=clusterIndices.begin(); i!=clusterIndices.end(); ++i)
    {
        cout<<_pointSet.at(*i).toString()<<" "<<endl;
        Cluster c = Cluster(_pointSet.at(*i));
        _clusterSet.push_back(c);
    }
    cout<<endl;

    _numObjects = _pointSet.size();
    _numClusters = _clusterSet.size();

    inFile.close();

    Point::metric.setData(_pointSet);
}

void KMeans::assign()
{
    do
    {
        _stable = true;

        for(unsigned int i=0; i<_pointSet.size(); ++i)
        {
            double min = DBL_MAX;
            int clusterId = 0;
            double distance = 0; //distance from p to cluster
            for (unsigned int j=0; j<_clusterSet.size(); ++j)
            {
                distance = _pointSet.at(i).distance(_clusterSet.at(j));
                if (distance <= min)
                {
                    min = distance;
                    clusterId = j;
                }
            }

            _pointSet.at(i).setGroupNo(clusterId);

            _clusterSet.at(clusterId).addPoint(&_pointSet.at(i));
        }
        //remove the empty clusters
        for (unsigned int i=_clusterSet.size()-1; ; --i)
        {
            if (_clusterSet.at(i).isEmpty())
                _clusterSet.erase(_clusterSet.begin()+i);

            if(i == 0)
                break;
        }
        //check the stability of each centroid
        for (unsigned int j=0; j<_clusterSet.size(); ++j)
        {
            if (_clusterSet.at(j).isShifted())
            {
                _stable = false;
                j = _clusterSet.size();
                break;
            }
        }
        if (!_stable)
        {
            //update all clusters
            for (unsigned int j=0; j<_clusterSet.size(); ++j)
            {
                _clusterSet.at(j).update();
            }
        }
    }while (!_stable);

    //remove the empty clusters
    for (unsigned int i=_clusterSet.size()-1; ; --i)
    {
        if (_clusterSet.at(i).isEmpty())
            _clusterSet.erase(_clusterSet.begin()+i);

        if(i == 0)
            break;
    }
}

Matrix KMeans::U() const
{
    Matrix u(_numObjects, _clusterSet.size());

    for(unsigned int k=0; k<_clusterSet.size(); ++k)
    {
        for(unsigned int i=0; i<_clusterSet.at(k).pointCount(); ++i)
        {
            _clusterSet.at(k).point(i)->setGroupNo(k);
        }
    }

    for(unsigned int i=0; i<_numObjects; ++i)
    {
        u[i][_pointSet.at(i).groupNo()] = 1;
    }

    return u;
}

#endif
