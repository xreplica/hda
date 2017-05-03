/*
 * H-DIANA class header
 *
 * TODO:
 * -
 *
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/algorithms/fuzzyentropy.h $
 */

#ifndef __FUZZYENTROPY_H__
#define __FUZZYENTROPY_H__

#include <cmath>
#include <queue>
#include <iostream>
#include <fstream>

#include "src/core/matrix.h"
#include "src/core/cluster.h"
#include "src/core/point.h"
#include "src/core/validation.h"

using namespace std;

class FuzzyEntropy
{
private:
    string          _file;
    unsigned int    _numClusters;
    unsigned int    _numObjects;
    vector<Cluster> _clusterSet;
    vector<Point>   _pointSet;
    vector<Point>   _normalizedPointSet;
//  Matrix          _u;
    double          _alpha;
    double          _beta;
    double          _gamma;

    vector<Point>   normalize(vector<Point>, double, double);
    double          similarity(Point, Point);

public:
    double beta;
    double gamma;

    //CONSTRUCTORS
    FuzzyEntropy();
    FuzzyEntropy(string);

    //ACCESSORS
    inline vector<Point>   pointSet() const   { return _pointSet; }
    inline vector<Cluster> clusterSet() const { return _clusterSet; }
    Matrix                 U() const;

    //FUNCTIONS
    void printClusters();
    void initialize(double, double);
    void run();
};

#endif

