/*
 * Distance class header:
 *
 * TODO:
 * -
 *
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/common/distance.h $
 */

#ifndef __DISTANCE_H__
#define __DISTANCE_H__

#include <iostream>
#include <cstdlib>
#include <vector>
#include <sstream>

#include "src/common/matrix.h"

class Point;
class Cluster;

using namespace std;

enum MetricType {EUCLIDEAN = 0, /*MAHALANOBIS = 1, */ KERNEL = 2};

class Distance
{
private:

    vector<Point> _dataset;          //complete dataset
    double*       _avg;              //average/mean object (array of features)
    double        _sigma;            //standard deviation across all dataset
//  Matrix        _covar;            //covariance matrix
//  Matrix        _covarInv;         //inverse of _covar
//  Matrix        _choleskyLower;    //cholesky lower decomposition (of _covarInv)
//  Matrix        _choleskyLowerInv; //inverse of _choleskyLower
    MetricType    _preferredMetric;  //chosen distance metric

    Matrix        _distanceCache;    //matrix of cached distances between all objects in dataser

    //FUNCTIONS

    double gaussianKernelFcn(const Point &, const Point&) const;    //Gaussian kernel function

public:

    //CONSTRUCTORS

    Distance();                          //default constructor
    Distance(vector<Point>);             //constructor specifying dataset as vector of points
    Distance(MetricType);                //constructor specifying metric type
    Distance(vector<Point>, MetricType); //constructor specifying both dataset and metric type

    ~Distance();

    //ACCESSORS

    inline MetricType preferredMetric() const               { return _preferredMetric; }
    inline void       setPreferredMetric(MetricType metric) { _preferredMetric = metric; }
    void              setData(const vector<Point>);

    //FUNCTIONS

    double distance(const Point &,  const Point &) const;           //return distance between given points
    double distance(const Point &) const;                           //return distance between given point and mean
    double euclidianDistance(const Point &, const Point &) const;   //return euclidean distance between given points
    double euclidianDistance(const Point &) const;                  //return euclidean distance between given point and mean
//  double mahalanobisDistance(const Point &, const Point &) const; //return mahalanobis distance between given points
//  double mahalanobisDistance(const Point &) const;                //return mahalanobis distance between given point and mean
    double kernelDistance(const Point &, const Point &) const;      //return (gaussian) kernel distance between given points
    double kernelDistance(const Point &) const;                     //return (gaussian) kernel distance between given point and mean

    //hamming distance are too problematic to be useable in most cases
//  double hammingDistance(const Point &, const Point &) const;     //return hamming distance (for binary vectors) between given points
//  double hammingDistance(const Point &) const;                    //return hamming distance (for binary vectors) between given ^point and mean
};

#endif
