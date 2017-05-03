/*
 * K-means+ class header:
 *
 * TODO:
 * -
 *
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/algorithms/kmeansplus.h $
 */

#ifndef __KMEANSPLUS_H__
#define __KMEANSPLUS_H__

#include <cfloat>
#include <ctime>

#include <fstream>

#include "src/core/cluster.h"
#include "src/core/matrix.h"

using namespace std;

class KMeansPlus
{
private:
    unsigned int    _dimension;
    string          _file;
    unsigned int    _initialNum;
    unsigned int    _finalClusterNum;
    unsigned int    _numClusters;
    unsigned int    _numObjects;
    vector<Cluster> _clusterSet;
    vector<Point>   _pointSet;
    int             _m1;              // count the # of interations of clustering
    int             _m2;              // count the # of interations of splitting
    int             _m2_m3;           // count the # of interations for stablizing clusters when splitting
    bool            _onSplitting;
    int*            _populationList;
    bool            _stable;

public:

    //CONSTRUCTORS
    KMeansPlus(string, unsigned int);

    //ACCESSORS
    inline vector<Cluster> clusterSet() const { return _clusterSet; }
    inline vector<Point>   pointSet() const   { return _pointSet; }
    Matrix                 U() const;

    //FUNCTIONS
    void   printClusters() const;
    void   run(unsigned int);
    void   initialize();
    void   assign();
    void   reAssign();
    Point* removeFarthestPoint();
    void   split();
    void   merge();
    void   link();
};


#endif
