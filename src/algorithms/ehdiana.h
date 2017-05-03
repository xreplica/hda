/*
 * EH-DIANA class header
 * 
 * TODO:
 * - 
 * 
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/algorithms/ehdiana.h $
 */

#ifndef __EHDIANA_H__
#define __EHDIANA_H__

#include <cmath>
#include <queue>
#include <iostream>
#include <fstream>

#include "src/core/matrix.h"
#include "src/core/cluster.h"
#include "src/core/point.h"
#include "src/core/validation.h"

using namespace std;

class EHDiana
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

    //FUNCTIONS
    vector<Point> normalize(vector<Point>, double, double);
    double        similarity(Point, Point);
    double        similarity(Point, Point, double);
    double        avgDistance(Point, Cluster);
    double        EavgDistance(Point, Cluster);
    double        clusterDiameter(Cluster);
    double        EclusterDiameter(Cluster);
    void          splitCluster(Cluster&);
    void          EsplitCluster(Cluster&);
    double        DNewCluster(Cluster&, Point&);
    double        DCluster(Cluster&, Point&);
    double        EDNewCluster(Cluster&, Point&);
    double        EDCluster(Cluster&, Point&);
    bool          stop();
    Cluster*      toSplit();

public:
    //CONSTRUCTORS
    EHDiana();
    EHDiana(string);

    //ACCESSORS/MUTATORS    
    inline getClusterSet() const { return _clusterSet; }
    inline getPointSet() const   { return _pointSet; }
    Matrix getU() const;

    //FUNCTIONS
    void printClusters();
    void initialize();
    void run();
};

#endif

