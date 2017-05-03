/*
 * RACAMAIC class implementation
 *
 * TODO:
 * -
 *
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/algorithms/racamaic.cpp $
 */

#ifndef __RACAMAIC_CPP__
#define __RACAMAIC_CPP__

#include "src/algorithms/racamaic.h"

//CONSTRUCTORS

Racamaic::Racamaic(string f, unsigned int fuzzy) : _file(f), _fuzzyConstant(fuzzy)
{
    _fcm = FuzzyCMeans(f, 2, fuzzy);
}

//FUNCTIONS

void Racamaic::printClusters() const
{
    for(unsigned int k=0; k<_clusterSet.size(); ++k)
    {
        cout<<"CLUSTER "<<k<<endl;
        for(unsigned int i=0; i<_clusterSet.at(k).pointCount(); ++i)
        {
            cout<<(*_clusterSet.at(k).point(i)).toString()<<endl;
        }
    }
}

bool Racamaic::reduceMatrix(Matrix &j)
{
    bool allZero = true;
    for(unsigned int i=0; i<j.rows(); ++i)
    {
        for(unsigned int i2=0; i2<j.cols(); ++i2)
        {
            if(j[i][i2] > 0)
            {
                j[i][i2]--;
                allZero = false;
            }
        }
    }

    return allZero;
}

void Racamaic::run()
{
    _fcm.initialize();
    _pointSet = _fcm.pointSet();
    _numObjects = _pointSet.size();
    _judgement = Matrix(_numObjects);
    unsigned int limit = (unsigned int)sqrt((double)_numObjects);

    for(unsigned int it=2; it<=limit; ++it)             //for each number of clusters
    {
        _fcm.setNumClusters(it);                            //set number of clusters
        _fcm.initialize(false);
        _fcm.run();

        _u = _fcm.U();                                      //retrieve relevant properties
        _pointSet = _fcm.pointSet();
        unsigned int numClusters = _fcm.clusterSet().size();
        for(unsigned int i=0; i<_numObjects; ++i)           //determine "discretized" memberships for all objects
        {
            unsigned int memberCluster = 0;
            double highestMembership = 0.0;

            for(unsigned int k=0; k<numClusters; ++k)
            {
                if(_u[i][k] > highestMembership)
                {
                    highestMembership = _u[i][k];
                    memberCluster = k;
                }
            }

            _pointSet.at(i).setGroupNo(memberCluster);
        }

        Matrix observation(_numObjects);                        //observation matrix
        for(unsigned int i=0; i<_numObjects; ++i)           //generate observation matrix
        {
            for(unsigned int i2=0; i2<_numObjects; ++i2)
            {
                if(_pointSet.at(i).groupNo() == _pointSet.at(i2).groupNo())
                    observation[i][i2] = 1.0;
            }
        }

        _judgement += observation;                          //add observation matrix to _judgement matrix (sum of all observation matrices)

        cout<<"Fuzzy c-means clustering for "<<it<<" clusters completed"<<endl;
    }

    //graph partitionning

    int mostStableGroupNumber = -1;
    int mostStableGroupCount = -1;
    Matrix mostStableGroupMatrix;
    int previousGroupNumber = -1;
    int currentGroupCount = 0;

    bool allZero = false;

    while(!allZero)                                         //while _judgement matrix is not all zeros
    {
        int numGroups = 0;

        queue<int> group;

        bool complete = false;

        for(unsigned int i=0; i<_numObjects; ++i)           //reset flag for all objects
        {
            _pointSet.at(i).setFlag();
        }

        while(!complete)                                    //while graph partitionning (via BFS) is not completed
        {
            int first = -1;

            for(unsigned int i=0; i<_numObjects; ++i)       //retrieve first valid object
            {
                if(_pointSet.at(i).flag())
                {
                    first = i;
                    group.push(i);
                    _pointSet.at(i).unsetFlag();
                    break;
                }
            }

            if(first == -1)                                 //if no object is valid, partitionning is complete
            {
                complete = true;
                break;
            }

            while(!group.empty())                           //BSF through objects using current _judgement matrix for adjacencies
            {
                for(unsigned int i=0; i<_numObjects; ++i)
                {
                    if(group.front() != (int)i && _pointSet.at(i).flag() && _judgement[group.front()][i] > 0)
                    {
                        group.push(i);
                        _pointSet.at(i).unsetFlag();
                    }
                }

                group.pop();
            }

            ++numGroups;
        }

        cout<<"Graph partitionning resulted in "<<numGroups<<" partitions"<<endl;

        if(previousGroupNumber == numGroups)
            ++currentGroupCount;
        else
        {
            currentGroupCount = 1;
            previousGroupNumber = numGroups;
        }

        if(currentGroupCount > mostStableGroupCount)
        {
            mostStableGroupNumber = numGroups;
            mostStableGroupCount = currentGroupCount;
            mostStableGroupMatrix = _judgement;
        }

        allZero = reduceMatrix(_judgement);                 //decrement all non-zero elements of _judgement matrix
    }

    cout<<mostStableGroupCount<<"/"<<limit<<" itterations of "<<mostStableGroupNumber<<" clusters"<<endl;

    _clusterSet.clear();
    queue<int> group;
    bool complete = false;

    for(unsigned int i=0; i<_numObjects; ++i)           //reset flag for all objects
    {
        _pointSet.at(i).setFlag();
    }

    while(!complete)                                    //while graph partitionning (via BFS) is not completed
    {
        Cluster cluster;
        int first = -1;

        for(unsigned int i=0; i<_numObjects; ++i)       //retrieve first valid object
        {
            if(_pointSet.at(i).flag())
            {
                first = i;
                group.push(i);
                _pointSet.at(i).unsetFlag();
                cluster.addPoint(&_pointSet.at(i));
                break;
            }
        }

        if(first == -1)                                 //if no object is valid, partitionning is complete
        {
            complete = true;
            break;
        }

        while(!group.empty())                           //BSF through objects using current _judgement matrix for adjacencies
        {
            for(unsigned int i=0; i<_numObjects; ++i)
            {
                if(group.front() != (int)i && _pointSet.at(i).flag() && mostStableGroupMatrix[group.front()][i] > 0)
                {
                    group.push(i);
                    _pointSet.at(i).unsetFlag();
                    cluster.addPoint(&_pointSet.at(i));
                }
            }

            group.pop();
        }

        _clusterSet.push_back(cluster);
    }
}

//ACCESSORS/MUTATORS

Matrix Racamaic::U() const
{
    Matrix u(_numObjects, _clusterSet.size());

    for(unsigned int k=0; k<_clusterSet.size(); ++k)
    {
        for(unsigned int i=0; i<_clusterSet.at(k).pointCount(); ++i)
        {
            _clusterSet.at(k).point(i)->setGroupNo(k);
        }
    }

    for(unsigned int i=0; i<_numObjects; ++i)
    {
        u[i][_pointSet.at(i).groupNo()] = 1;
    }

    return u;
}

#endif

