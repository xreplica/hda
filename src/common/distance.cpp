/*
 * Distance class implementation:
 *
 * TODO:
 * -
 *
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/common/distance.cpp $
 */

#ifndef __DISTANCE_CPP__
#define __DISTANCE_CPP__

//#define CACHE_DISTANCES

#include "src/common/distance.h"
#include "src/common/point.h"

//CONSTRUCTORS

Distance::Distance() : _dataset(vector<Point>()), _preferredMetric(EUCLIDEAN) {}

Distance::Distance(vector<Point> dataset) : _preferredMetric(EUCLIDEAN)
{
    setData(dataset);
}

Distance::Distance(MetricType metric) : _dataset(vector<Point>()), _preferredMetric(metric) {}

Distance::Distance(vector<Point> dataset, MetricType metric) : _preferredMetric(metric)
{
    setData(dataset);
}

Distance::~Distance() {}

//FUNCTIONS

void Distance::setData(const vector<Point> data)
{
    _dataset = data;

//  _covar = Matrix(Point::dimension());

    _avg = new double[(const int)Point::dimension()];
    for(unsigned int j=0; j<Point::dimension(); ++j)
        _avg[j] = 0.0;

    //calculate avg
    for(unsigned int i=0; i<_dataset.size(); ++i)
    {
        for(unsigned int j=0; j<Point::dimension(); ++j)
        {
            _avg[j] += _dataset.at(i).feature(j);
        }
    }
    for(unsigned int j=0; j<Point::dimension(); ++j)
    {
        _avg[j] /= _dataset.size();
    }

    //make a Point at avg
    stringstream ss;

    for(unsigned int i=0; i<Point::dimension(); ++i)
    {
        ss << _avg[i]<<",";
    }

    Point avgPoint = (ss.str().substr(0, ss.str().length() - 1));

    //calculate sigma
    _sigma = 0.0;

    for(unsigned int i=0; i<_dataset.size(); ++i)
    {
        _sigma += euclidianDistance(_dataset.at(i), avgPoint) * euclidianDistance(_dataset.at(i), avgPoint);
    }

    _sigma = sqrt(_sigma / _dataset.size());

    //set covar
/*  for(unsigned int i=0; i<Point::dimension(); ++i)
    {
        for(unsigned int j=0; j<Point::dimension(); ++j)
        {
            double sum = 0.0;

            for(unsigned int p=0; p<_dataset.size(); ++p)
            {
                sum += (_dataset.at(p).feature(i) - _avg[i]) * (_dataset.at(p).feature(j) - _avg[j]);
            }

            _covar[i][j] = sum / (_dataset.size() - 1);
        }
    }
*/
    //set matrices
/*  _covarInv = _covar.getInverse();
    _choleskyLower = _covarInv.getCholesky();
    _choleskyLowerInv = _choleskyLower.getInverse();
*/
#ifdef CACHE_DISTANCES
    _distanceCache.reset(_dataset.size(), _dataset.size()+1);

    for(unsigned int i=0; i<_dataset.size(); ++i)
    {
        for(unsigned int j=0; j<_dataset.size(); ++j)
        {
            double distance;

            if(_preferredMetric == EUCLIDEAN)
                distance = euclidianDistance(_dataset.at(i), _dataset.at(j));
//          else if(_preferredMetric == MAHALANOBIS)
//              distance = mahalanobisDistance(_dataset.at(i), _dataset.at(j));
            else if(_preferredMetric == KERNEL)
                distance = kernelDistance(_dataset.at(i), _dataset.at(j));
            else
                distance = 0;

            _distanceCache[_dataset.at(i).getPid()][_dataset.at(j).getPid()] = distance;
        }

        double distance;

        if(_preferredMetric == EUCLIDEAN)
            distance = euclidianDistance(_dataset.at(i));
//      else if(_preferredMetric == MAHALANOBIS)
//          distance = mahalanobisDistance(_dataset.at(i));
        else if(_preferredMetric == KERNEL)
            distance = kernelDistance(_dataset.at(i));
        else
            distance = 0;

        _distanceCache[_dataset.at(i).getPid()][_dataset.size()] = distance;
    }
#endif
}

double Distance::distance(const Point &a, const Point &b) const
{
#ifdef CACHE_DISTANCES
    if(a.getPid() < _dataset.size() && b.getPid() < _dataset.size())
        return _distanceCache[a.getPid()][b.getPid()];
#endif

    if(_preferredMetric == EUCLIDEAN)
        return euclidianDistance(a, b);
//  else if(_preferredMetric == MAHALANOBIS)
//      return mahalanobisDistance(a, b);
    else if(_preferredMetric == KERNEL)
        return kernelDistance(a, b);
    else
        return 0;
}

double Distance::distance(const Point &p) const
{
#ifndef _CACHE_DISTANCE

    if(_preferredMetric == EUCLIDEAN)
        return euclidianDistance(p);
//  else if(_preferredMetric == MAHALANOBIS)
//      return mahalanobisDistance(p);
    else if(_preferredMetric == KERNEL)
        return kernelDistance(p);
    else
        return 0;

#else

    return _distanceCache[p.getPid()][_dataset.size()];

#endif
}

double Distance::euclidianDistance(const Point &a, const Point &b) const
{
    double distance = 0.0;

    for (unsigned int j=0; j<Point::dimension(); ++j)
        distance += (a.feature(j) - b.feature(j)) * (a.feature(j) - b.feature(j));

    distance = sqrt(distance);

    return distance;
}

double Distance::euclidianDistance(const Point &a) const
{
    double distance = 0.0;

    for (unsigned int j=0; j<Point::dimension(); ++j)
        distance += (a.feature(j) - _avg[j]) * (a.feature(j) - _avg[j]);

    distance = sqrt(distance);

    return distance;
}

/*
double Distance::mahalanobisDistance(const Point &a, const Point &b) const
{
    if(_dataset.size() < 2)
    {
        cerr<<"insufficient data set, cannot calculate Mahalanobis distance"<<endl;
        exit(1);
    }

    double distance = 0.0;
    Matrix deltaV(1, Point::dimension());

    for(unsigned int j=0; j<Point::dimension(); ++j)
        deltaV[0][j] = a.feature(j) - b.feature(j);

    distance = (deltaV * _covarInv * deltaV.getTransposed())[0][0];

    return sqrt(distance);
}

double Distance::mahalanobisDistance(const Point &a) const
{
    if(_dataset.size() < 2)
    {
        cerr<<"insufficient data set, cannot calculate Mahalanobis distance"<<endl;
        exit(1);
    }

    double distance = 0.0;
    Matrix deltaV(1, Point::dimension());

    for(unsigned int j=0; j<Point::dimension(); ++j)
        deltaV[0][j] = a.feature(j) - _avg[j];

    distance = (deltaV.getTransposed() * _covarInv * deltaV)[0][0];

    return sqrt(distance);
}
*/

double Distance::gaussianKernelFcn(const Point &a, const Point &b) const
{
    return exp(-((euclidianDistance(a, b)) * (euclidianDistance(a, b))) / (_sigma * _sigma));
}

double Distance::kernelDistance(const Point &a, const Point &b) const
{
    return gaussianKernelFcn(a, a) + gaussianKernelFcn(b, b) - (2 * gaussianKernelFcn(a, b));
}

double Distance::kernelDistance(const Point &a) const
{
    stringstream ss;

    for(unsigned int i=0; i<Point::dimension(); ++i)
        ss << _avg[i]<<",";

    Point b = (ss.str().substr(0, ss.str().length() - 1));

    return gaussianKernelFcn(a, a) + gaussianKernelFcn(b, b) - (2 * gaussianKernelFcn(a, b));
}

#endif
