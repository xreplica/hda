/*
 * Cluster class implementation
 *
 * TODO:
 * -
 *
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/common/cluster.cpp $
 */

#ifndef __CLUSTER_CPP__
#define __CLUSTER_CPP__

#include "src/common/cluster.h"
#include "src/common/distance.h"

double Cluster::_outlier = 5;
double Cluster::_link = 0.5;
double Cluster::_threshold = 0.01;
unsigned int Cluster::_nextGroupNo = 0;

//CONSTRUCTORS

Cluster::Cluster() : _fpt(NULL), _longest(0), _population(0), _sigma(0)
{
    _groupNo = _nextGroupNo++;
}

Cluster::Cluster(string line) : _fpt(NULL), _longest(0), _population(0), _sigma(0)
{
    _groupNo = _nextGroupNo++;
    parseLine(line);
}

Cluster::Cluster(Point p) : _fpt(NULL), _longest(0), _population(0), _sigma(0)
{
    _groupNo = _nextGroupNo++;
    for(unsigned int i= 0; i<_dimension;++i)
        _featureSet.at(i) = p.feature(i);
}

Cluster::~Cluster() {}

//FUNCTIONS

//test if cluster has moved in any dimension by more than the threshold
bool Cluster::isShifted() const
{
    for(unsigned int i=0; i<_dimension; ++i)
    {
        if(_featureSet.at(i)-_featureSet.at(i) >= _threshold)
            return true;
    }

    return false;
}

//add a point to the cluster
void Cluster::addPoint(Point *point)
{
    double dist = distance(*point);

    for(unsigned int i=0; i<_dimension; ++i)
        _featureSet.at(i) = (_population * _featureSet.at(i) + (*point).feature(i)) / (_population + 1);

    if (dist > _longest)
    {
        _fpt = point;
        _longest = dist;
    }

    if(_points.size() == 0)
        point->avgClusterDist = 0;
    else
    {
        double avgDist = 0.0;
        for(unsigned long i=0; i<pointCount(); ++i)
        {
            double dist = point->distance(*_points.at(i));
            avgDist += (dist - avgDist) / (i + 1);

            _points.at(i)->avgClusterDist += (dist - _points.at(i)->avgClusterDist) / (_population);
        }
        point->avgClusterDist = avgDist;
    }

    _sigma = sqrt((_sigma * _sigma * _population + dist * dist) / (_population + 1));
    ++_population;
    _points.push_back(point);
}

void Cluster::removePoint(unsigned int index)
{
    double dist = distance(*_points.at(index));

    for(unsigned int i=0; i<_dimension; ++i)
        _featureSet.at(i) = (_population * _featureSet.at(i) - (*_points.at(index)).feature(i)) / (_population - 1);

    if(_longest == dist)
    {
        _longest = 0.0;

        for(unsigned long i=0; i<_points.size(); ++i)
        {
            if(i == index)
                continue;

            if(distance(*_points.at(i)) > _longest)
                _longest = distance(*_points.at(i));
        }
    }

    for(unsigned long i=0; i<pointCount(); ++i)
        _points.at(i)->avgClusterDist  -= (point(index)->distance(*point(i)) - _points.at(i)->avgClusterDist) / (_population - 2);

    _sigma = sqrt((_sigma * _sigma * _population - dist * dist) / (_population - 1));
    _population--;
    _points.erase(_points.begin() + index);
}

//update the cluster and reset dynamic parameters
void Cluster::update()
{
    for(unsigned int i=0; i<_dimension; ++i)
    {
        _featureSet.at(i) = _featureSet.at(i);
        _featureSet.at(i) = 0;
    }
}

//remove the farthest point from the cluster center
void Cluster::removeFpt()
{
    for(unsigned int i = 0; i<_dimension; ++i)
        _featureSet.at(i) = (_population * _featureSet.at(i) - (*_fpt).feature(i)) / (_population - 1);

    _sigma = sqrt((_sigma * _sigma * _population - _longest * _longest) / (_population - 1));
    _population--;

    for(unsigned int i=0; i<_points.size(); ++i)
    {
        Point* point = _points.at(i);
        if (_fpt == point)
            _points.erase(_points.begin()+(i));
    }

    _fpt = NULL;
}

//merge cluster with another given cluster
Cluster Cluster::mergeWith(const Cluster &cluster)
{
    Cluster newCluster = Cluster();
    double sigma1 = cluster.sigma();
    int p1 = _population;
    int p2 = cluster.population();

    for(unsigned int i=0; i<_dimension; ++i)
    {
        double temp = (cluster.feature(i)*p2 + _featureSet.at(i)*p1)/(p1+p2);
        _featureSet.at(i) = temp;
        newCluster.setFeature(i, temp);
    }

    for(unsigned int i=0; i<_points.size(); ++i)
        newCluster.addPoint(_points.at(i));

    for(unsigned int i=0; i<cluster.pointCount(); ++i)
        newCluster.addPoint(cluster.point(i));

    double newSigma = sqrt((_sigma*_sigma*p1 + sigma1*sigma1*p2)/(p1+p2));
    newCluster.setSigma(newSigma);
    newCluster.setPopulation(p1+p2);

    return newCluster;
}

//link cluster with another given cluster
void Cluster::linkWith(Cluster &cluster)
{
    int groupNo1 = _groupNo;
    int groupNo2 = cluster.groupNo();

    if (groupNo1 < groupNo2)
        cluster.setLinkNo(groupNo1);
    else if (groupNo2 < groupNo1)
        setLinkNo(groupNo2);

    _siblingSet.push_back(&cluster);
    cluster.addSibling(*this);
}

//set link number for cluster
void Cluster::setLinkNo(int no)
{
    if(_groupNo > no)
    {
        _groupNo = no;

        for(unsigned long i = 0; i<_siblingSet.size(); ++i)
            (*_siblingSet.at(i)).setGroupNo(no);
    }
}

//create string of points strings
string Cluster::containerToString() const
{
    if(!_points.empty()){
        Point c1 = *_points.at(0);
        string s = c1.toString() + "\n";

        for(unsigned long i=1; i<_points.size(); ++i)
        {
            Point c = *_points.at(i);
            s += c.toString() + "\n";
        }

        return s;
    }
    else
        return string("EMPTY");
}

#endif
