/*
 * H-DIANA class header
 *
 * TODO:
 * -
 *
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/algorithms/hdiana.h $
 */

#ifndef __HDIANA_H__
#define __HDIANA_H__

#include <cmath>
#include <queue>
#include <iostream>
#include <fstream>

#include "src/common/matrix.h"
#include "src/common/cluster.h"
#include "src/common/point.h"
#include "src/common/validation.h"

using namespace std;

class HDiana
{
private:
    string _file;
    unsigned int    _mergeMethod;
    double          _mergeParameter;
    unsigned int    _numClusters;
    unsigned int    _numObjects;
    vector<Cluster> _clusterSet;
    vector<Point>   _pointSet;
    Matrix          _u;
    bool            _noSplit;
    bool            _verbose;
    unsigned short  _stopIterations;

    //FUNCTIONS
    double avgDistance(const Point&, const Cluster&);
    double avgNearestNeighbor(const Cluster&);
    double nearestPairDist(const Cluster&, const Cluster&);
    double clusterDiameter(const Cluster&);
    double clusterHistogramVariation(const Cluster&);
    void   splitCluster(Cluster&);
    double DNewCluster(Cluster&, Point&);
    double DCluster(Cluster&, Point&);

public:
    //CONSTRUCTORS
    HDiana(bool noSplit = false, unsigned short iterations = 7);
    HDiana(string, bool noSplit = false, unsigned short iterations = 7);

    //ACCESSORS/MUTATORS
    inline       vector<Cluster> clusterSet() const           { return _clusterSet; }
    inline       vector<Point> pointSet() const               { return _pointSet; }
    unsigned int mergeMethod() const                          { return _mergeMethod; }
    inline void  setMergeMethod(unsigned int m)               { _mergeMethod = m; }
    double       mergeParameter() const                       { return _mergeParameter; }
    void         setMergeParameter(double p)                  { _mergeParameter = p; }
    void         setVerbose(bool state=true)                  {  _verbose = state; }
    void         unsetVerbose()                               { _verbose = false; }
    void         setStopIterations(unsigned short iterations) { _stopIterations = iterations; }
    Matrix       U() const;

    //FUNCTIONS
    void printClusters();
    void initialize();
    void run();
};

#endif

