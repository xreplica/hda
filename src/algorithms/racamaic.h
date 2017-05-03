/*
 * RACAMAIC class header
 *
 * TODO:
 * -
 *
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/algorithms/racamaic.h $
 */

#ifndef __RACAMAIC_H__
#define __RACAMAIC_H__

#include <cmath>
#include <queue>

#include "src/algorithms/fuzzycmeans.h"
#include "src/core/matrix.h"
#include "src/core/cluster.h"
#include "src/core/validation.h"

using namespace std;

class Racamaic
{
private:
    FuzzyCMeans     _fcm;
    string          _file;
    unsigned int    _numObjects;
    unsigned int    _fuzzyConstant;
    vector<Cluster> _clusterSet;
    vector<Point>   _pointSet;
    Matrix          _u;
    Matrix          _judgement;
    bool            reduceMatrix(Matrix &);

public:
    //CONSTRUCTORS
    Racamaic(string, unsigned int);

    //ACCESSORS/MUTATORS
    inline vector<Cluster> clusterSet() const {return _clusterSet; }
    inline vector<Point>   pointSet() const   { return _pointSet; }
    Matrix                 U() const;

    //FUNCTIONS
    void printClusters() const;
    void run();
};

#endif

