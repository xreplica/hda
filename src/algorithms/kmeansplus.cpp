/*
 * K-means+ clustering implementation
 *
 * TODO:
 * - fix all points of all clusters fall to same point
 *
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/algorithms/kmeansplus.cpp $
 */

#ifndef __KMEANSPLUS_CPP__
#define __KMEANSPLUS_CPP__

#include <set>

#include "src/algorithms/kmeansplus.h"

//CONSTRUCTORS

KMeansPlus::KMeansPlus(string f, unsigned int d) : _dimension(d), _file(f), _finalClusterNum(0), _m1(0), _m2(0), _m2_m3(0), _onSplitting(false) {}

//FUNCTIONS

void KMeansPlus::printClusters() const
{
    for (unsigned int i=0; i<_clusterSet.size(); ++i)
    {
        cout<<"The size of cluster [ "<<i<<" ] is: "<<_clusterSet.at(i).points().size()<<endl;
        cout<<_clusterSet.at(i).containerToString()<<endl;
    }
}

void KMeansPlus::run(unsigned int num)
{
    _initialNum = num;
    initialize();
    assign();
    split();
    //link();
    merge();
}

void KMeansPlus::initialize()
{
    string line;
    ifstream inFile(_file.c_str(), ios::in);
    set<unsigned int> clusterIndices;

    if(!inFile.is_open()){
        cout<<"Could not open file \""<<_file<<"\""<<endl;
        exit(1);
    }

    while(getline(inFile, line))
    {
        Point p = Point(line);
        _pointSet.push_back(p);
    }

    //randomly select initial centroid seeds
    srand((unsigned int)time(NULL));

    while(clusterIndices.size()<_initialNum)
    {
        int index = rand() % _pointSet.size();
        if(clusterIndices.find(index) == clusterIndices.end())
            clusterIndices.insert(index);
    }

    cout<<"initial clusters: "<<endl;
    for(set<unsigned int>::iterator i=clusterIndices.begin(); i!=clusterIndices.end(); ++i)
    {
        cout<<_pointSet.at(*i).toString()<<" "<<endl;
        Cluster c = Cluster(_pointSet.at(*i));
        _clusterSet.push_back(c);
    }
    cout<<endl;

    _numObjects = _pointSet.size();
    _numClusters = _clusterSet.size();

    inFile.close();

    Point::metric.setData(_pointSet);
}

void KMeansPlus::assign()
{
    do
    {
        ++_m1;
        _stable = true;

        for(unsigned int i=0; i<_pointSet.size(); ++i)
        {
            double min = DBL_MAX;
            int clusterId = 0;
            double distance = 0; //distance from p to cluster
            for (unsigned int j=0; j<_clusterSet.size(); ++j)
            {
                distance = _pointSet.at(i).distance(_clusterSet.at(j));
                if (distance <= min)
                {
                    min = distance;
                    clusterId = j;
                }
            }

            _pointSet.at(i).setGroupNo(clusterId);

            _clusterSet.at(clusterId).addPoint(&_pointSet.at(i));
        }
        //remove the empty clusters
        for (unsigned int i=_clusterSet.size()-1; ; --i)
        {
            if (_clusterSet.at(i).isEmpty())
                _clusterSet.erase(_clusterSet.begin()+i);

            if(i == 0)
                break;
        }
        //check the stability of each centroid
        for (unsigned int j=0; j<_clusterSet.size(); ++j)
        {
            if (_clusterSet.at(j).isShifted())
            {
                _stable = false;
                j = _clusterSet.size();
                break;
            }
        }
        if (!_stable)
        {
            //update all clusters
            for (unsigned int j=0; j<_clusterSet.size(); ++j)
            {
                _clusterSet.at(j).update();
            }
        }
    }while (!_stable);

    //remove the empty clusters
    for (unsigned int i=_clusterSet.size()-1; ; --i)
    {
        if (_clusterSet.at(i).isEmpty())
            _clusterSet.erase(_clusterSet.begin()+i);

        if(i == 0)
            break;
    }
}

void KMeansPlus::reAssign()
{
    for(unsigned int i=0; i<_clusterSet.size(); ++i)
    {
        _clusterSet.at(i).update();
    }

    assign();
}

Point* KMeansPlus::removeFarthestPoint()
{
    double maxDist = 0;
    int index=0;
    Point *pt = NULL;

    for(unsigned int i=0; i<_clusterSet.size(); ++i)
    {
        Point *temp = _clusterSet.at(i).fpt();
        double d = _clusterSet.at(i).longest();
        if(d>=maxDist)
        {
            maxDist = d;
            index = i;
            pt = temp;
        }
    }
    _clusterSet.at(index).setLongest(0);
    _clusterSet.at(index).removeFpt();

    return pt;
}

void KMeansPlus::split()
{
    unsigned int newClusterNum = 0; //number of clusters after splitting

    //flag of spltting: if no splitting, it is true
    bool _stable;

    do{
        _stable = true;

        for(unsigned int i=0; i<_clusterSet.size(); ++i)
        {
            Point *p = _clusterSet.at(i).fpt();
            if(_clusterSet.at(i).hasOutlier())
            {
                Cluster newC = Cluster((*p).toString());
                newC.addPoint(p);
                _clusterSet.at(i).removeFpt();
                _clusterSet.push_back(newC);
            }

            newClusterNum = _clusterSet.size();

            // if there is a new cluster created
            if(_clusterSet.size() < newClusterNum){
                ++_m2;
                reAssign(); // redo clustering
                _stable = false;
            }
        }
    }while(!_stable);

    cout<<"the number of iteration for splitting is: "<<_m2<<endl<<endl;
}

void KMeansPlus::merge()
{
    //check if two clusters can be merged
    for (unsigned int i = 0; i<_clusterSet.size(); ++i)
    {
        for (unsigned int j=i+1; j<_clusterSet.size(); ++j)
        {
            if(/*_clusterSet.at(j) != NULL && */_clusterSet.at(j).canMerge(_clusterSet.at(i)))
            {
                cout<<"Merging: "<<endl;

                Cluster newCluster = _clusterSet.at(i).mergeWith(_clusterSet.at(j));
                //set one merged cluster as newCluster;
                _clusterSet.at(i) = newCluster;//_clusterSet.setElementAt(newCluster, i);
                cout<<"Clusters [ "<<i<<", "<<j<<" ] are merged."<<endl;
                _clusterSet.erase(_clusterSet.begin()+j);
                if (j < _clusterSet.size())
                    j = j-1;

                cout<<"After merging, the population of newCluster [ "<<i<<" ] is:"<<newCluster.population()<<endl;
            }
        }
    }

    for (unsigned int i=0; i<_clusterSet.size(); ++i)
    {
        cout<<"Cluster [ "<<i<<" ] = "<<_clusterSet.at(i).population()<<endl;
    }
}

void KMeansPlus::link()
{

    for(unsigned int i=0; i<_clusterSet.size(); ++i)
    {
        _clusterSet.at(i).setGroupNo(i);//initally set group no. to each cluster
    }

    //check if two clusters can be linked
    for(unsigned int i=0; i<_clusterSet.size(); ++i)
    {
        //int grpNo1 = _clusterSet.at(i).groupNo();
        for(unsigned int j=i+1; j<_clusterSet.size(); ++j)
        {
            //int grpNo2 = _clusterSet.at(j).groupNo();
            if(_clusterSet.at(j).canMerge(_clusterSet.at(i)))
            {
                // link two close clusters using linkWith() module of Cluster class
                _clusterSet.at(i).linkWith(_clusterSet.at(j));
            }
        }
    }

    //record the population of each group, i.e., the merged cluster
    //delete _populationList;
    _populationList = new int[_clusterSet.size()];
    for(unsigned int i=0; i<_clusterSet.size(); ++i)
    {
        _populationList[i] = 0;
    }

    for(unsigned int j=0; j<_clusterSet.size(); ++j)
    {
        _populationList[_clusterSet.at(j).groupNo()] += _clusterSet.at(j).population();
    }

    //count the number of clusters after linking
    for(unsigned int j=0; j<_clusterSet.size(); ++j)
    {
        if(_populationList[j]!= 0)
        {
            cout<<"popuList["<<j<<"]="<<_populationList[j]<<endl;
            ++_finalClusterNum;
        }
    }
}

Matrix KMeansPlus::U() const
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
