/*
 * K-means+ class header:
 * 
 * TODO:
 * -
 * 
 * $Author: francois $
 * $Date: 2014-08-10 17:24:20 -0400 (Sun, 10 Aug 2014) $
 * $Revision: 36 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/algorithms/kmeans.h $
 */

#ifndef __KMEANS_H__
#define __KMEANS_H__

#include <cfloat>
#include <ctime>

#include <fstream>

#include "src/core/cluster.h"
#include "src/core/matrix.h"

using namespace std;

class KMeans
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
    int*            _populationList;
    bool            _stable;

public:
    
    //CONSTRUCTORS
    KMeans(string, unsigned int);

    //ACCESSORS/MUTATORS
    inline vector<Cluster> clusterSet() const { return _clusterSet; }
    inline vector<Point>   pointSet() const   { return _pointSet; }
    Matrix                 U() const;

    //FUNCTIONS
    void   printClusters() const;
    void   run(unsigned int);
    void   initialize();
    void   assign();    
    Point* removeFarthestPoint();
};


#endif
