/*
 * Clustering validation class header:
 * Static class containing functions to calculate various validity indices. Functions are implemented for the folowing:
 * - LSE
 * - Xie & Beni validity index
 * - Fukuyama & Sugeno validity index
 * - Kwon validity index
 * - Compose Within and between (CWB) validity index
 * - Zahid's validity index
 * - Fuzzy Hypervolume
 *
 * TODO:
 * -
 *
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/common/validation.h $
 */

#ifndef __VALIDATION_H__
#define __VALIDATION_H__

#include <cfloat>
#include <algorithm>
#include <map>

#include "src/common/matrix.h"
#include "src/common/point.h"
#include "src/common/cluster.h"

using namespace std;

class Validation
{
private:
    //CONSTRUCTORS

    Validation(){}

    //FUNCTIONS

    static double  minInterClusterDistance(const vector<Cluster>&);
    static double  maxInterClusterDistance(const vector<Cluster>&);
    static double  avgInterClusterDistance(const vector<Cluster>&);
    static Cluster avgClusterCenter(const vector<Cluster>&);
    static Point   avgPoint(const vector<Point>&);
    static Point   patternVariance(const vector<Point>&);
    static double  scat(const vector<Cluster>&, const vector<Point>&, const Matrix&, unsigned int);
    static double  dist(const vector<Cluster>&, const vector<Point>&);
    static double  sc1(const vector<Cluster>&, const vector<Point>&, const Matrix&, unsigned int);
    static double  sc2(const vector<Cluster>&, const vector<Point>&, const Matrix&);
    static Matrix  fuzzyCovarianceMatrix(unsigned int, const Cluster&, unsigned int, const vector<Point>&, const Matrix&, unsigned int);
    static double  avgDistanceWithCluster(const Point&, const Cluster&);
public:
    static double  leastSquaredError(const vector<Cluster>&, const vector<Point>&, const Matrix&, unsigned int);
    static double  xieBeniIndex(const vector<Cluster>&, const vector<Point>&, const Matrix&, unsigned int);
    static double  fukuyamaSugenoIndex(const vector<Cluster>&, const vector<Point>&, const Matrix&, unsigned int);
    static double  kwonIndex(const vector<Cluster>&, const vector<Point>&, const Matrix&, unsigned int);
    static double  composeWithinBetweenIndex(const vector<Cluster>&, const vector<Point>&, const Matrix&, unsigned int);
    static double  zahidIndex(const vector<Cluster>&, const vector<Point>&, const Matrix&, unsigned int);
    static double  fuzzyHypervolume(const vector<Cluster>&, const vector<Point>&, const Matrix&, unsigned int);
    static double  randIndex(const vector<Cluster>&, const vector<Cluster>&, const vector<Point>&);
    static double  PBMIndex(const vector<Cluster>&, const vector<Point>&, const Matrix&);
    static double  silhouette(const vector<Cluster>&, const vector<Point>&);
    static double  silhouette(const vector<Cluster>&, const vector<Point>&, const Matrix&);
};

#endif
