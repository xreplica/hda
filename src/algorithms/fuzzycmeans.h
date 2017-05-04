/*
 * Fuzzy c-means class header:
 * 
 * TODO:
 * - 
 * 
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/algorithms/fuzzycmeans.h $
 */

#ifndef __FUZZYCMEANS_H__
#define __FUZZYCMEANS_H__

#include <cfloat>
#include <ctime>

#include <fstream>
#include <iomanip>

#include "src/common/matrix.h"
#include "src/common/cluster.h"

using namespace std;

class FuzzyCMeans
{
private:
    string _file;
    static double   _variation_threshold;
    unsigned int    _numClusters;
    unsigned int    _numObjects;
    unsigned int    _fuzzyConstant;
    double          _epsilon;
    vector<Cluster> _clusterSet;
    vector<Point>   _pointSet;
    Matrix          _u;

public:

    //CONSTRUCTORS
    FuzzyCMeans();
    FuzzyCMeans(string, int, unsigned int);

    //ACCESSORS
    inline vector<Cluster> clusterSet() const             { return _clusterSet; }
    inline vector<Point>   pointSet() const               { return _pointSet; }
    inline Matrix          U() const                      { return _u; }
    inline void            setNumClusters(unsigned int c) { _numClusters = c; }

    //FUNCTIONS
    void   printU() const;
    void   printClusters() const;
    void   initialize(bool read = true);
    void   run();
    void   newRun(double threshold);
    double computeEpsilon(const Matrix oldU);
    void   updateU();
    void   updateC();
    void   newUpdateU();
    void   newUpdateC(double threshold);

    //functions belonging to fuzzy c-means+

    bool   fuzzyMerge();
    bool   fuzzySplit();
    void   removeCluster(unsigned int index);
    void   addCluster(Cluster);
};

#endif

