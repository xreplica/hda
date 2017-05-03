/*
 * Data point class implementation
 * 
 * TODO:
 * - 
 * 
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/core/point.cpp $
 */

#ifndef __POINT_CPP__
#define __POINT_CPP__

#include "src/core/point.h"

unsigned int Point::_dimension = 2;             //2 dimensional objects by default
Distance Point::metric = Distance(EUCLIDEAN);   //euclidean distance by default
unsigned long Point::_nextPid = 0;              //nitialize poind IDs at 0

//CONSTRUCTORS

Point::Point() : _groupNo(-1),  _featureSet(vector<double>(Point::dimension(), 0.0)), _flag(true), avgClusterDist(0)            //initialize point with group no -1 (no group)
{
    _pid = _nextPid++;  //assign unique ID
}

Point::Point(string line) : _groupNo(-1), _flag(true), avgClusterDist(0)    //initialise point with input string and group no -1 (no group)
{
    parseLine(line);                                                        //parse input line to set features

    _pid = _nextPid++;                                                      //assign unique ID
}

Point::~Point() {}

//FUNCTIONS

void Point::parseLine(const string line)    //parse string formated point(vector)
{
    unsigned int currDimension = 0;         //track number of parsed dimensions
    stringstream ss;                        //string stream for parsed features
    stringstream label;                     //string stream for label
    set<char> delimiters;                   //set of delimiter characters

    char delims[]= {' ','\t',','};          //set (default) delimiters
    delimiters.insert(delims, delims+3);

    string::const_iterator it=line.begin();
    while(true)
    {
        char c;                 //current character

        if(it!=line.end()){     //if not eol
            c = *it;            //set c
            ++it;
        }
        else                                            //mirror: at eol
        {
            if(ss.str().length() > 0)                   //if ss is non-empty
            {
                if(currDimension >= Point::dimension()) //if parsed dimensions exceeds point dimensions (parsing label)
                {
                    label << " " << ss.str();           //append parsed stream to label
                    ss.str("");
                }
                else                                    //mirror: parsed dimension is within point dimensions
                {
                    if(currDimension < _featureSet.size())                                  //if featureSet has a position corresponding to the dimension
                        _featureSet.at(currDimension) = atof((char*)(ss.str().c_str()));    //insert parsed value at corresponding position
                    else                                                                    //mirror: featureSet's length is insufficient
                        _featureSet.push_back(atof((char*)(ss.str().c_str())));             //append parsed value at featureSet
                }
            }

            break;
        }

        if(delimiters.find(c) != delimiters.end())  //if current character is a delimiter
        {
            if(currDimension >= Point::dimension()) //if parsed dimension exceeds point dimensions (parsing label)
            {
                label << " " << ss.str();           //append parsed stream to label
                ss.str("");
            }
            else                                    //mirror: parsed dimensions is within point dimensions
            {
                if(currDimension < _featureSet.size())                                  //if featureSet has position corresponding to the dimensions
                {
                    _featureSet.at(currDimension) = atof((char*)(ss.str().c_str()));    //insert parsed value at corresponding position
                    ss.str("");
                    ++currDimension;
                }
                else                                //mirror: featureSet's length is insufficient
                {
                    _featureSet.push_back(atof((char*)(ss.str().c_str())));             //append parsed value to deatureSet
                    ss.str("");
                    ++currDimension;
                }
            }
        }
        else            //mirror: current character is not a delimiter
        {
            ss << c;    //append current character to parsed feature stream
        }
    }

    _label = label.str();
}

double Point::distance(const Point p) const
{
    return metric.distance(*this, p);
}

bool Point::isEqual(const Point p) const
{
    for (unsigned int i=0; i<_dimension; ++i)
    {
        if (_featureSet.at(i) != p.feature(i))
            return false;
    }

    return true;
}

string Point::toString() const
{
    stringstream oss;
    oss << _featureSet.at(0);
    for (unsigned int i=1; i<_featureSet.size(); ++i)
        oss << " " << _featureSet.at(i);

    oss << _label;

    return oss.str();
}

#endif
