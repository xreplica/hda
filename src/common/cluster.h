/*
 * Cluster class header:
 *
 * TODO:
 * -
 *
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/common/cluster.h $
 */

#ifndef __CLUSTER_H__
#define __CLUSTER_H__

#include <vector>

#include "src/common/point.h"
#include "src/common/matrix.h"

using namespace std;

class Cluster : public Point
{
protected:
    static double       _outlier;
    static double       _link;
    static double       _threshold;
    static unsigned int _nextGroupNo;
    Point*              _fpt;         //the farthest point from the center in the cluster
    double              _longest;     //distance of the farthest point from the cluster center
    int                 _population;  //number of data points assigned to this cluster
    double              _sigma;       //standard deviation
    double              _avgDist;     //average distance of points to centroid deviation
    vector<Cluster*>    _siblingSet;
    vector<Point*>      _points;

public:

    //CONSTRUCTORS
    Cluster();
    Cluster(string);
    Cluster(Point);

    ~Cluster();

    //ACCESSORS
    static inline void      setOutlier(double outlier)     { _outlier = outlier; }
    static inline void      setLink(double link)           { _link = link; }
    static inline void      setThreshold(double threshold) { _threshold = threshold; }
    inline Point*           fpt() const                    { return _fpt; }
    inline int              population() const             { return _population; }
    inline void             setPopulation(int population)  { _population = population; }
    inline double           sigma() const                  { return _sigma; }
    inline void             setSigma(double sigma)         { _sigma = sigma; }
    inline double           avgDist() const                { return _avgDist; }
    inline double           longest() const                { return _longest; }
    inline void             setLongest(double longest)     { _longest = longest; }
    inline vector<Point*>   points() const                 { return _points; }
    inline unsigned int     pointCount() const             { return _points.size(); }
    inline Point*           point(unsigned int i) const    { return _points.at(i); }
    inline vector<Cluster*> siblingSet() const             { return _siblingSet; }
    inline void             addSibling(Cluster &cluster)   { _siblingSet.push_back(&cluster); }

    //FUNCTIONS
    inline bool             isEmpty() const    { return _points.empty() ? true : false; }
    inline bool             hasOutlier() const { return _longest > _sigma * _outlier ? true : false; }
    bool                    isShifted() const;
    void                    addPoint(Point *);
    void                    removePoint(unsigned int);
    void                    update();
    void                    removeFpt();
    inline bool             canMerge(const Cluster &cluster)
    {
        double interDist = distance(cluster);
        double sigma1 = cluster.sigma();

        return interDist <= (_sigma + sigma1) * _link ? true : false;
    }
    Cluster                 mergeWith(const Cluster&);
    void                    linkWith(Cluster&);
    void                    setLinkNo(int);
    string                  containerToString() const;
};

#endif
