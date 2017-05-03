/*
 * H-DIANA class implementation
 *
 * TODO:
 * -
 *
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/algorithms/fuzzyentropy.cpp $
 */

#ifndef __FUZZYENTROPY_CPP__
#define __FUZZYENTROPY_CPP__

#include <cmath>
#include <limits>

#include "src/algorithms/fuzzyentropy.h"

//CONSTRUCTORS

FuzzyEntropy::FuzzyEntropy() {}

FuzzyEntropy::FuzzyEntropy(string f) : _file(f) {}

//FUNCTIONS

void FuzzyEntropy::printClusters()
{
    for(unsigned int k=0; k<_clusterSet.size(); ++k)
    {
        cout<<"CLUSTER "<<k<<endl;
        for(unsigned int i=0; i<_clusterSet.at(k).pointCount(); ++i)
        {
            cout<<(*_clusterSet.at(k).point(i)).toString()<<endl;
        }
    }
}

void FuzzyEntropy::initialize(double b, double y)
{
    _beta = b;
    _gamma = y;
    double avgD = 0.0;
    string line;
    ifstream inFile(_file.c_str(), ios::in);

    if(!inFile.is_open()){
        cout<<"Could not open file \""<<_file<<"\""<<endl;
        exit(1);
    }

    while(getline(inFile, line))
    {
        Point p = Point(line);
        _pointSet.push_back(p);
    }

    inFile.close();

    Point::metric.setData(_pointSet);

    //Cluster c = Cluster();
    for(unsigned int i=0; i<_pointSet.size(); ++i)
    {
        _pointSet.at(i).setGroupNo(-1);
        //c.addPoint(&_pointSet.at(i));
    }
    //_clusterSet.push_back(c);

    _numObjects = _pointSet.size();
    //_numClusters = 1;

    _normalizedPointSet = normalize(_pointSet, 0.0, 1.0);

    for(unsigned int i=0; i<_normalizedPointSet.size(); ++i)
    {
        for(unsigned int j=i+1; j<_normalizedPointSet.size(); ++j)
        {
            avgD += _normalizedPointSet.at(i).distance(_normalizedPointSet.at(j));
        }
    }

    avgD /= _normalizedPointSet.size() * (_normalizedPointSet.size() - 1) / 2;

    _alpha = -log(0.5) / avgD;

    for(unsigned int i=0; i<_normalizedPointSet.size(); ++i)
    {
        double sum = 0.0;
        for(unsigned int j=0; j<_normalizedPointSet.size(); ++j)
        {
            if(i == j)
                continue;

            double Sij = similarity(_normalizedPointSet.at(i), _normalizedPointSet.at(j));
            if(Sij == 1)
            {
                continue;       //Sij == 1 wil cause division by zero ( log2(1) == 0 ), Ei == 0 when Sij == 1, so we skip it's calculation as it would simply add 0 to the sum
            }

            sum += (Sij * (log(Sij) / log(2))) + (1 - Sij) * (log(1 - Sij) / log(2));
        }

        _normalizedPointSet.at(i).entropy = -sum;
    }
}

void FuzzyEntropy::run()
{
    bool empty = false;

    while(!empty)
    {
        empty = true;
        double minEntropy = numeric_limits<double>::max( );
        Point *minPoint = NULL;

        for(unsigned int i=0; i<_normalizedPointSet.size(); ++i)
        {
            if(!_normalizedPointSet.at(i).flag())
                continue;

            empty = false;

            if(_normalizedPointSet.at(i).entropy < minEntropy)
            {
                minEntropy = _normalizedPointSet.at(i).entropy;
                minPoint = &_normalizedPointSet.at(i);
            }
        }

        if(empty)
            break;

        Cluster c;

        c.addPoint(minPoint);
        minPoint->unsetFlag();

        for(unsigned int i=0; i<_normalizedPointSet.size(); ++i)
        {
            if(similarity(*minPoint, _normalizedPointSet.at(i)) >= _beta)
            {
                c.addPoint(&_normalizedPointSet.at(i));
                _normalizedPointSet.at(i).unsetFlag();
            }
        }

        if(c.pointCount() >= _gamma)
            _clusterSet.push_back(c);
    }

    cout<<"Completed with "<<_clusterSet.size()<<" clusters"<<endl;
}

vector<Point> FuzzyEntropy::normalize(vector<Point> pointSet, double min, double max)
{
    vector<Point> normalizedSet;
    double mins[Point::dimension()], maxs[Point::dimension()];

    for(unsigned int j=0; j<Point::dimension(); ++j)
    {
        mins[j] = maxs[j] = pointSet.at(0).feature(j);
    }

    for(unsigned int i=1; i<pointSet.size(); ++i)
    {
        for(unsigned int j=0; j<Point::dimension(); ++j)
        {
            if(pointSet.at(i).feature(j) < mins[j])
                mins[j] = pointSet.at(i).feature(j);
            if(pointSet.at(i).feature(j) > maxs[j])
                maxs[j] = pointSet.at(i).feature(j);
        }
    }

    for(unsigned int i=1; i<pointSet.size(); ++i)
    {
        string pointVals = "";

        for(unsigned int j=0; j<Point::dimension(); ++j)
        {
            stringstream ss;
            ss << (pointSet.at(i).feature(j) - mins[j]) / (maxs[j] - mins[j]);
            pointVals = pointVals.append(ss.str()).append(",");
        }

        pointVals = pointVals.substr(0, pointVals.size()-1);

        normalizedSet.push_back(Point(pointVals));
    }

    return normalizedSet;
}

double FuzzyEntropy::similarity(Point pointA, Point pointB)
{
    return exp(-_alpha * pointA.distance(pointB));
}

/*
void FuzzyEntropy::splitCluster(Cluster &cluster)
{
    Cluster newCluster;
    double maxDistance = 0.0;
    unsigned int maxIndex = 0;
    bool change = true;

//cout<<cluster.containerToString()<<endl;

    if(cluster.getNumPoints() == 2)
    {
        newCluster.addPoint(cluster.getPoint(1));
        cluster.removePoint(1);
        _clusterSet.push_back(newCluster);
        return;
    }

    for(unsigned int i=0; i<cluster.getNumPoints(); ++i)
    {
        if(avgDistance(*cluster.getPoint(i), cluster) > maxDistance)
        {
            maxDistance = avgDistance((*cluster.getPoint(i)), cluster);
            maxIndex = i;
        }
    }

    newCluster.addPoint(cluster.getPoint(maxIndex));
    cluster.removePoint(maxIndex);

    while(change && cluster.getNumPoints() > 1)
    {
        change = false;
        for(unsigned int i=0; i<cluster.getNumPoints(); ++i)
        {
            if(DNewCluster(newCluster, *cluster.getPoint(i)) < DCluster(cluster, *cluster.getPoint(i)))
            {
                newCluster.addPoint(cluster.getPoint(i));
                cluster.removePoint(i);
                change = true;
                break;
            }
        }
    }

//cout<<endl<<"SPLIT"<<endl<<endl;
//cout<<cluster.containerToString()<<endl;;
//cout<<newCluster.containerToString()<<endl;

    _clusterSet.push_back(newCluster);
}
*/

Matrix FuzzyEntropy::U() const
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

