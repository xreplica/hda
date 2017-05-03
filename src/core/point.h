/*
 * Data point class header:
 *
 * TODO:
 * -
 *
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/core/point.h $
 */

#ifndef __POINT_H__
#define __POINT_H__

#include <cmath>
#include <cstdlib>

#include <iostream>
#include <set>
#include <sstream>
#include <vector>

#include "src/core/matrix.h"
#include "src/core/distance.h"

using namespace std;

class Point
{
protected:
    static unsigned int  _dimension;      //data vector dimension
    static unsigned long _nextPid;        //to generate unique IDs
    unsigned int         _pid;            //unqie ID of point
    int                  _groupNo;        //clustered group number (linked clusters for a single group)
    vector<double>       _featureSet;     //features of the data point (vector elements)
    string               _label;          //point label (or unused features)
    bool                 _flag;           //(temporary use)

public:
    static               Distance metric; //distance metric object
    double               entropy;         //temporary
    double               internalEntropy; //temporary
    double               avgClusterDist;

    //CONSTRUCTORS

    Point();
    Point(string);
    ~Point();

    //ACCESSORS

    static inline unsigned int dimension()                          { return _dimension; }
    static inline void         setDimension(unsigned int dimension) { _dimension = dimension; }
    inline double              feature(unsigned int i) const
    {
        double feature = 0;

        if (i>=_dimension)
            cout<<"***This feature index is out of range."<<endl;   //throw exception
        else
            feature = _featureSet.at(i);

        return feature;
    }
    inline void                setFeature(unsigned int i, double v)
    {
        if(i>=_dimension)
            cout<<"***This feature index is out of range."<<endl;   //throw exception
        else
            _featureSet.at(i)= v;
    }
    inline int                 groupNo() const         { return _groupNo; }
    inline void                setGroupNo(int groupNo) { _groupNo = groupNo; }
    inline vector<double>      featureSet() const      { return _featureSet; }
    inline bool                flag() const            { return _flag; }
    inline void                setFlag(bool val=true)  { _flag = val; }
    inline void                unsetFlag()             { _flag = false; }
    inline unsigned int        pid() const             { return _pid; }

    //FUNCTIONS

protected:
    void    parseLine(const string);     //parse input string (of features) into feature set

public:
    double  distance(const Point) const; //return distance between self and given point
    bool    isEqual(const Point) const;  //return true if given point matches self (all features have same value
    string  toString() const;            //return string representing self as feature vector (with label if present)
};

#endif
