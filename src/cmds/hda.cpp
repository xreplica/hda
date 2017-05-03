/*
 * heuristic DIANA implementation
 *
 * TODO:
 * -
 *
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/cmds/hda.cpp $
 */

#include <iostream>
#include <cstdlib>

#include "src/algorithms/hdiana.h"
#include "src/core/validation.h"
#include "src/core/params.h"

using namespace std;

int main(int argc, const char* argv[])
{
    Params parameters;

    parameters.parse(argc, argv);

    Point::setDimension(parameters.vDimension);
    Point::metric.setPreferredMetric(parameters.metric);

    HDiana hda(parameters.filename, parameters.noSplit);

    hda.setMergeMethod(parameters.mergeMethod);
    hda.setMergeParameter(parameters.mergeParameter);
    hda.setStopIterations(parameters.stopIterations);
    hda.setVerbose(parameters.verbose);

    try
    {
        hda.initialize();
        hda.run();
    }
    catch(matrixException &e)
    {
        cerr<<e.what()<<endl;
        exit(1);
    }

    if(parameters.displayU)
    {
        cout<<endl<<"Partition matrix:"<<endl;
        Matrix u = hda.U().transposed();
        for(unsigned int i=0; i<u.rows();++i)
        {
            for(unsigned int j=0; j<u.cols();++j)
                cout<<u[i][j]<<",";
            cout<<endl;
        }
    }

    if(parameters.displayClusters)
        hda.printClusters();
    if(parameters.displayIndices)
    {
        cout<<"Xie-Beni index(-): "<<Validation::xieBeniIndex(hda.clusterSet(), hda.pointSet(), hda.U(), 2)<<endl;
        cout<<"Fukuyama-Sugeno index(-): "<<Validation::fukuyamaSugenoIndex(hda.clusterSet(), hda.pointSet(), hda.U(), 2)<<endl;
        cout<<"Kwon index(-): "<<Validation::kwonIndex(hda.clusterSet(), hda.pointSet(), hda.U(), 2)<<endl;
        cout<<"CWB index*(-): "<<Validation::composeWithinBetweenIndex(hda.clusterSet(), hda.pointSet(), hda.U(), 2)<<endl;
        cout<<"Zahid index*(+): "<<Validation::zahidIndex(hda.clusterSet(), hda.pointSet(), hda.U(), 2)<<endl;
        cout<<"PBM index*(+): "<<Validation::PBMIndex(hda.clusterSet(), hda.pointSet(), hda.U())<<endl;
        cout<<"Silhouette index(+): "<<Validation::silhouette(hda.clusterSet(), hda.pointSet())<<endl;
    }

    return 0;
}

