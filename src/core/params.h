/*
 * CLI parameters handling class header
 *
 * TODO:
 * -
 *
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/core/params.h $
 */

#ifndef __PARAMS_H__
#define __PARAMS_H__

#include <string>

#include "src/core/distance.h"

using namespace std;

class Params
{
public:
    string         filename;
    unsigned int   vDimension;
    unsigned int   mergeMethod;
    double         mergeParameter;
    unsigned short stopIterations;
    MetricType     metric;
    bool           displayClusters;
    bool           displayIndices;
    bool           displayU;
    bool           noSplit;
    bool           verbose;

    //CONSTRUCTORS

    Params();

    //FUNCTIONS

    void parse(int argc, const char* argv[]);
    void usage(char invFlag = 0, bool missingVal = false);
};

#endif
