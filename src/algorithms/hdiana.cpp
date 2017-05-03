/*
 * H-DIANA class implementation
 *
 * TODO:
 * -
 *
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/algorithms/hdiana.cpp $
 */

#ifndef __HDIANA_CPP__
#define __HDIANA_CPP__

#include <cmath>
#include <algorithm>

#include "src/algorithms/hdiana.h"

using namespace std;

//CONSTRUCTORS

HDiana::HDiana(bool noSplit, unsigned short iterations) : _noSplit(noSplit), _verbose(false), _stopIterations(iterations) {}

HDiana::HDiana(string f, bool noSplit, unsigned short iterations) : _file(f), _noSplit(noSplit), _verbose(false), _stopIterations (iterations) {}

//FUNCTIONS

void HDiana::printClusters()
{
    for(unsigned long k=0; k<_clusterSet.size(); ++k)
    {
        cout<<"CLUSTER "<<k<<endl;
        for(unsigned long i=0; i<_clusterSet.at(k).pointCount(); ++i)
            cout<<(*_clusterSet.at(k).point(i)).toString()<<endl;
    }
}

void HDiana::initialize()
{
    string line;

    unsigned long currClust = -1;

    ifstream inFile(_file.c_str(), ios::in);

    if(!inFile.is_open()){
        cout<<"Could not open file \""<<_file<<"\""<<endl;
        exit(1);
    }

    while(getline(inFile, line))
    {
        if(line.substr(0,7) == "CLUSTER")
            ++currClust;
        else
        {
            Point p = Point(line);
            p.setGroupNo(currClust);
            _pointSet.push_back(p);
        }
    }

    Point::metric.setData(_pointSet);

    if(_noSplit)
    {
        for(unsigned long k=0; k<currClust+1; ++k)
            _clusterSet.push_back(Cluster());
        for(unsigned long i=0; i<_pointSet.size(); ++i)
            _clusterSet.at(_pointSet.at(i).groupNo()).addPoint(&_pointSet.at(i));
    }
    else
    {
        Cluster c = Cluster();
        for(unsigned long i=0; i<_pointSet.size(); ++i)
            c.addPoint(&_pointSet.at(i));
        _clusterSet.push_back(c);
    }

    _numObjects = _pointSet.size();
    _numClusters = _clusterSet.size();

    if(_verbose)
    {
        time_t t = time(NULL);
        struct tm tm = *localtime(&t);
        cout<<"@"<<tm.tm_hour<<":"<<tm.tm_min<<":"<<tm.tm_sec<<" : "<<_numObjects<<" objects read into "<<_numClusters<<" clusters."<<endl<<endl;
    }

    inFile.close();
}

void HDiana::run()
{
    //splitting phase
    if(!_noSplit)
    {
        if(_verbose)
        {
            time_t t = time(NULL);
            struct tm tm = *localtime(&t);
            cout<<"@"<<tm.tm_hour<<":"<<tm.tm_min<<":"<<tm.tm_sec<<" : "<<"Splitting..."<<endl<<endl;
        }

        double minAvgCDist = DBL_MAX;
        vector<Cluster> optimalClusterSet = _clusterSet;

        for(unsigned long loop=_numObjects; loop>0; --loop)
        {
            double avgIntraClusterDist = 0.0;
            int maxIndex = -1;
            double maxDiameter = -1;//clusterDiameter(_clusterSet.at(0));

            for(unsigned long k=0; k<_clusterSet.size(); ++k)
            {
                if(clusterHistogramVariation(_clusterSet.at(k)) > maxDiameter && _clusterSet.at(k).pointCount() > 1)
                //if(clusterDiameter(_clusterSet.at(k)) > maxDiameter && _clusterSet.at(k).pointCount() > 1)
                {
                    maxIndex = k;
                    maxDiameter = clusterHistogramVariation(_clusterSet.at(k));
                    //maxDiameter = clusterDiameter(_clusterSet.at(k));
                }
            }

            if(maxIndex == -1)
            {
                cout<<"no more clusters are eligible for splitting"<<endl;
                break;
            }
            else if(_clusterSet.at(maxIndex).pointCount() < 2)
            {
                cout<<"over splitting: tried to split a single-object cluster"<<endl;//<<_clusterSet.size()<<" clusters."<<endl;
                break;
            }
            else
            {
                if(_verbose)
                {
                    time_t t = time(NULL);
                    struct tm tm = *localtime(&t);
                    cout<<"@"<<tm.tm_hour<<":"<<tm.tm_min<<":"<<tm.tm_sec<<" : "<<"Splitting cluster "<<maxIndex<<"."<<endl<<endl;
                }

                splitCluster(_clusterSet.at(maxIndex));
            }

            double avg = 0.0;
            for(unsigned long k=0; k<_clusterSet.size(); ++k)
            {
                double clusterAvg = 0.0;
                if(_clusterSet.at(k).pointCount() > 1)
                {
                    for(unsigned long i=0; i<_clusterSet.at(k).pointCount(); ++i)
                        clusterAvg += _clusterSet.at(k).point(i)->avgClusterDist;

                    clusterAvg /= _clusterSet.at(k).pointCount();
                    avg += clusterAvg;
                }
            }
            avg /= _clusterSet.size();
            avgIntraClusterDist = avg;

            if(avgIntraClusterDist < minAvgCDist)//if clustering is better
            {
                optimalClusterSet = _clusterSet;
                minAvgCDist = avgIntraClusterDist;
                loop = _stopIterations;

                if(_verbose)
                {
                    time_t t = time(NULL);
                    struct tm tm = *localtime(&t);
                    cout<<"@"<<tm.tm_hour<<":"<<tm.tm_min<<":"<<tm.tm_sec<<" : "<<"Current clustering is better than at previous iteration."<<endl<<endl;
                }
            }
            else if(_verbose)
            {
                time_t t = time(NULL);
                struct tm tm = *localtime(&t);
                cout<<"@"<<tm.tm_hour<<":"<<tm.tm_min<<":"<<tm.tm_sec<<" : "<<"Current clustering is worse than at previous iteration, "<<loop<<" iterations until stop."<<endl<<endl;
            }
        }

        _clusterSet = optimalClusterSet;
    }

    //merge phase
    if(_verbose && _mergeMethod > 0)
    {
        time_t t = time(NULL);
        struct tm tm = *localtime(&t);
        cout<<"@"<<tm.tm_hour<<":"<<tm.tm_min<<":"<<tm.tm_sec<<" : "<<"Merging..."<<endl;
    }

    //_mergeMethod == 0 implies no merging
    if(_mergeMethod == 1) //merge (NN)
    {
        bool done = false;

        while(!done)
        {
            bool merged = false;
            for(unsigned long k=0; k<_clusterSet.size(); ++k)
            {
                for(unsigned long m=k+1; m<_clusterSet.size(); ++m)
                {
                    bool mergeable = false;
                    double dist1 = avgNearestNeighbor(_clusterSet.at(k));
                    double dist2 = avgNearestNeighbor(_clusterSet.at(m));
                    double m_dist = 0.0;

                    if(dist1 == 0)
                        m_dist = dist2;
                    else if(dist2 == 0)
                        m_dist = dist1;
                    else if(dist1 > dist2)
                        m_dist = dist1;
                    else
                        m_dist = dist2;

                    for(unsigned long i=0; i<_clusterSet.at(k).pointCount(); ++i)
                    {
                        for(unsigned long j=0; j<_clusterSet.at(m).pointCount(); ++j)
                        {
                            if(_clusterSet.at(k).point(i)->distance(*_clusterSet.at(m).point(j)) < _mergeParameter*m_dist)
                            {
                                mergeable = true;
                                break;
                            }
                        }

                        if(mergeable)
                            break;
                    }

                    if(mergeable)
                    {
                        if(_verbose)
                        {
                            time_t t = time(NULL);
                            struct tm tm = *localtime(&t);
                            cout<<"@"<<tm.tm_hour<<":"<<tm.tm_min<<":"<<tm.tm_sec<<" : "<<"Merging clusters "<<k<<" and "<<m<<"."<<endl<<endl;
                        }

                        for(unsigned long i=0; i<_clusterSet.at(m).pointCount(); ++i)
                            _clusterSet.at(k).addPoint(_clusterSet.at(m).point(i));

                        _clusterSet.erase(_clusterSet.begin() + m);

                        merged = true;
                        break;
                    }
                }
                if(merged)
                    break;
            }
            if(merged)
                continue;

            done = true;
        }
    }
    else if(_mergeMethod == 2)  //alt merge (k-means+ based)
    {
        for(unsigned long k=0; k<_clusterSet.size(); ++k)
        {
            for(unsigned long m=k+1; m<_clusterSet.size(); ++m)
            {
                if(_clusterSet.at(k).distance(_clusterSet.at(m)) <= (_clusterSet.at(k).sigma() + _clusterSet.at(m).sigma()) * _mergeParameter)
                {
                    if(_verbose)
                    {
                        time_t t = time(NULL);
                        struct tm tm = *localtime(&t);
                        cout<<"@"<<tm.tm_hour<<":"<<tm.tm_min<<":"<<tm.tm_sec<<" : "<<"Linking clusters "<<k<<" and "<<m<<"."<<endl<<endl;
                    }

                    _clusterSet.at(k).linkWith(_clusterSet.at(m));
                }
            }
        }

        if(_verbose)
        {
            time_t t = time(NULL);
            struct tm tm = *localtime(&t);
            cout<<"@"<<tm.tm_hour<<":"<<tm.tm_min<<":"<<tm.tm_sec<<" : "<<"Merging clusters."<<endl<<endl;
        }

        for(unsigned long k=_clusterSet.size()-1; ; --k)
        {
            for(unsigned long m=k-1; ; --m)
            {
                if(_clusterSet.at(k).groupNo() == _clusterSet.at(m).groupNo())
                {
                    for(unsigned long i=0; i<_clusterSet.at(k).pointCount(); ++i)
                    {
                        _clusterSet.at(m).addPoint(_clusterSet.at(k).point(i));
                    }

                    _clusterSet.erase(_clusterSet.begin() + k);

                    if(k > 0)
                        --k;
                    else
                        break;
                }

                if(m == 0)
                    break;
            }

            if(k <= 1)
                break;
        }
    }
    else if(_mergeMethod == 3)  //alt merge (x% within y sigma)
    {
        for(unsigned long k=0; k<_clusterSet.size(); ++k)
        {
            unsigned long merge_threshold = (unsigned long)(_clusterSet.at(k).pointCount() * (_mergeParameter / 100.0)); //parameter for percentage, default 5%
            for(unsigned long m=k+1; m<_clusterSet.size(); ++m)
            {
                unsigned long m_count = 0;

                for(unsigned long i=0; i<_clusterSet.at(k).pointCount(); ++i)
                {
                    if(_clusterSet.at(k).point(i)->distance(_clusterSet.at(m)) < 5 * _clusterSet.at(m).sigma()) //parameter for #sigma, default 5 sigma
                    {
                        ++m_count;

                        if(m_count >= merge_threshold)
                        {
                            if(_verbose)
                            {
                                time_t t = time(NULL);
                                struct tm tm = *localtime(&t);
                                cout<<"@"<<tm.tm_hour<<":"<<tm.tm_min<<":"<<tm.tm_sec<<" : "<<"Linking clusters "<<k<<" and "<<m<<"."<<endl<<endl;
                            }

                            _clusterSet.at(k).linkWith(_clusterSet.at(m));
                            break;
                        }
                    }
                }
            }
        }

        if(_verbose)
        {
            time_t t = time(NULL);
            struct tm tm = *localtime(&t);
            cout<<"@"<<tm.tm_hour<<":"<<tm.tm_min<<":"<<tm.tm_sec<<" : "<<"Merging clusters."<<endl<<endl;
        }

        for(unsigned long k=_clusterSet.size()-1; ; --k)
        {
            for(unsigned long m=k-1; ; --m)
            {
                if(_clusterSet.at(k).groupNo() == _clusterSet.at(m).groupNo())
                {
                    for(unsigned long i=0; i<_clusterSet.at(k).pointCount(); ++i)
                    {
                        _clusterSet.at(m).addPoint(_clusterSet.at(k).point(i));
                    }

                    _clusterSet.erase(_clusterSet.begin() + k);

                    if(k > 0)
                        --k;
                    else
                        break;
                }

                if(m == 0)
                    break;
            }

            if(k <= 1)
                break;
        }
    }
    else if(_mergeMethod == 4)
    {
        bool done = false;

        while(!done)
        {
            bool merged = false;
            for(unsigned long k=0; k<_clusterSet.size(); ++k)
            {
                for(unsigned long m=k+1; m<_clusterSet.size(); ++m)
                {
                    bool mergeable = false;
                    unsigned int m_count = 0;
                    double dist1 = avgNearestNeighbor(_clusterSet.at(k));
                    double dist2 = avgNearestNeighbor(_clusterSet.at(m));
                    double m_dist = 0.0;

                    if(dist1 == 0 && dist2 == 0)
                    {
                        for(unsigned long i=0; i<=_numClusters; ++i)
                        {
                            m_dist += _clusterSet.at(i).sigma();
                        }

                        m_dist /= _numClusters;
                    }
                    else if(dist1 == 0)
                        m_dist = dist2;
                    else if(dist2 == 0)
                        m_dist = dist1;
                    else if(dist1 > dist2)
                        m_dist = dist1;
                    else
                        m_dist = dist2;

                    for(unsigned long i=0; i<_clusterSet.at(k).pointCount(); ++i)
                    {
                        for(unsigned long j=0; j<_clusterSet.at(m).pointCount(); ++j)
                        {
                            if(_clusterSet.at(k).point(i)->distance(*_clusterSet.at(m).point(j)) < _mergeParameter*m_dist)
                            {
                                //mergeable = true;
                                //break;
                                m_count++;
                            }
                        }

                        //if(mergeable)
                            //break;
                    }

                    if(_clusterSet.at(k).pointCount() < _clusterSet.at(m).pointCount())
                    {
                        if(m_count > _clusterSet.at(k).pointCount() * 0.05 || _clusterSet.at(k).pointCount() == 1)
                            mergeable = true;
                    }
                    else
                    {
                        if(m_count > _clusterSet.at(m).pointCount() * 0.05 || _clusterSet.at(m).pointCount() == 1)
                            mergeable = true;
                    }

                    if(mergeable)
                    {
                        if(_verbose)
                        {
                            time_t t = time(NULL);
                            struct tm tm = *localtime(&t);
                            cout<<"@"<<tm.tm_hour<<":"<<tm.tm_min<<":"<<tm.tm_sec<<" : "<<"Merging clusters "<<k<<" and "<<m<<"."<<endl<<endl;
                        }

                        for(unsigned long i=0; i<_clusterSet.at(m).pointCount(); ++i)
                            _clusterSet.at(k).addPoint(_clusterSet.at(m).point(i));

                        _clusterSet.erase(_clusterSet.begin() + m);

                        merged = true;
                        break;
                    }
                }
                if(merged)
                    break;
            }
            if(merged)
                continue;

            done = true;
        }
    }
    else if(_mergeMethod == 5) //merge (minNN)
    {
        bool done = false;

        while(!done)
        {
            bool merged = false;
            for(unsigned long k=0; k<_clusterSet.size(); ++k)
            {
                double min_dist = DBL_MAX;
                double min_index = -1;

                for(unsigned long m=k+1; m<_clusterSet.size(); ++m)
                {
                    double dist1 = avgNearestNeighbor(_clusterSet.at(k));
                    double dist2 = avgNearestNeighbor(_clusterSet.at(m));
                    double m_dist = 0.0;

                    if(dist1 == 0 && dist2 == 0)
                    {
                        for(unsigned long i=0; i<=_numClusters; ++i)
                        {
                            m_dist += _clusterSet.at(i).sigma();
                        }

                        m_dist /= _numClusters;
                    }
                    else if(dist1 == 0)
                        m_dist = dist2;
                    else if(dist2 == 0)
                        m_dist = dist1;
                    else if(dist1 > dist2)
                        m_dist = dist1;
                    else
                        m_dist = dist2;

                    for(unsigned long i=0; i<_clusterSet.at(k).pointCount(); ++i)
                    {
                        for(unsigned long j=0; j<_clusterSet.at(m).pointCount(); ++j)
                        {
                            double dist = _clusterSet.at(k).point(i)->distance(*_clusterSet.at(m).point(j));

                            if(dist < _mergeParameter*m_dist && dist < min_dist)
                            {
                                min_dist = _clusterSet.at(k).point(i)->distance(*_clusterSet.at(m).point(j));
                                min_index = m;
                            }
                        }
                    }
                }

                if(min_index >= 0)
                {
                    if(_verbose)
                    {
                        time_t t = time(NULL);
                        struct tm tm = *localtime(&t);
                        cout<<"@"<<tm.tm_hour<<":"<<tm.tm_min<<":"<<tm.tm_sec<<" : "<<"Merging clusters "<<k<<" and "<<min_index<<"."<<endl<<endl;
                    }

                    for(unsigned long i=0; i<_clusterSet.at(min_index).pointCount(); ++i)
                        _clusterSet.at(k).addPoint(_clusterSet.at(min_index).point(i));

                    _clusterSet.erase(_clusterSet.begin() + min_index);

                    merged = true;
                }
            }

            if(merged)
                continue;

            done = true;
        }
    }

    cout<<"finished with "<<_clusterSet.size()<<" clusters."<<endl;
}

double HDiana::avgDistance(const Point &point, const Cluster &cluster)
{
    double avgDist = 0.0;
    unsigned long t = 1;

    for(unsigned long i=0; i<cluster.pointCount(); ++i)
    {
        if(point.pid() == cluster.point(i)->pid())
            continue;
        avgDist += (point.distance(*cluster.point(i)) - avgDist) / t++;//avgDist += point.distance(*cluster.point(i));
    }

    return avgDist;//return avgDist / (cluster.pointCount() - 1);
}

double HDiana::avgNearestNeighbor(const Cluster &cluster)
{
    double avgDist = 0.0;

    for(unsigned long i=0; i<cluster.pointCount() - 1; ++i)
    {
        double closest = DBL_MAX;
        for(unsigned long j=i+1; j<cluster.pointCount(); ++j)
            if(cluster.point(i)->distance(*cluster.point(j)) < closest)
                closest = cluster.point(i)->distance(*cluster.point(j));
        avgDist += closest;
    }

    return avgDist / cluster.pointCount();
}

double HDiana::nearestPairDist(const Cluster &c1, const Cluster &c2)
{
    double nearest = DBL_MAX;
    for(unsigned long i=0; i<c1.pointCount(); ++i)
        for(unsigned long j=0; j<c2.pointCount(); ++j)
            if(c1.point(i)->distance(*c2.point(j)) < nearest)
                nearest = c1.point(i)->distance(*c2.point(j));

    return nearest;
}

double HDiana::clusterDiameter(const Cluster &cluster)
{
    double diameter = 0.0;

    for(unsigned long i=0; i<cluster.pointCount(); ++i)
    {
        for(unsigned long j=0; j<cluster.pointCount(); ++j)
        {
            if(i == j)
                continue;

            if((*cluster.point(i)).distance(*cluster.point(j)) > diameter)
                diameter = (*cluster.point(i)).distance(*cluster.point(j));
        }
    }

    return diameter;
}

double HDiana::clusterHistogramVariation(const Cluster &cluster)
{
    vector<double> max;
    vector<double> min;
    max.resize(Point::dimension());
    min.resize(Point::dimension());
    double diameter = 0.0;

    for(unsigned int j=0; j<Point::dimension(); ++j)
    {
        max.at(j) = cluster.points().at(0)->feature(j);
        min.at(j) = cluster.points().at(0)->feature(j);
    }

    for(unsigned int i=1; i<cluster.pointCount(); ++i)
    {
        for(unsigned int j=0; j<Point::dimension(); ++j)
        {
            if(cluster.points().at(i)->feature(j) > max.at(j))
                max.at(j) = cluster.points().at(i)->feature(j);
            if(cluster.points().at(i)->feature(j) < min.at(j))
                min.at(j) = cluster.points().at(i)->feature(j);
        }
    }

    for(unsigned int j=0; j<Point::dimension(); ++j)
        diameter += max.at(j) - min.at(j);

    return diameter;
}

void HDiana::splitCluster(Cluster &cluster)
{
    Cluster newCluster;
    double maxDistance = 0.0;
    unsigned int maxIndex = 0;
    bool change = true;

    if(cluster.pointCount() == 2)
    {
        newCluster.addPoint(cluster.point(1));
        cluster.removePoint(1);
        _clusterSet.push_back(newCluster);
        return;
    }

    for(unsigned long i=0; i<cluster.pointCount(); ++i)
    {
        if(cluster.point(i)->avgClusterDist > maxDistance)//if(avgDistance(*cluster.point(i), cluster) > maxDistance)
        {
            maxDistance = cluster.point(i)->avgClusterDist;
            maxIndex = i;
        }
    }

    newCluster.addPoint(cluster.point(maxIndex));
    cluster.removePoint(maxIndex);

    while(change && cluster.pointCount() > 1)
    {
        change = false;
        for(unsigned long i=0; i<cluster.pointCount(); ++i)
        {
            if(DNewCluster(newCluster, *cluster.point(i)) < DCluster(cluster, *cluster.point(i)))
            {
                newCluster.addPoint(cluster.point(i));
                cluster.removePoint(i);
                change = true;
                break;
            }
        }
    }

    _clusterSet.push_back(newCluster);
}

double HDiana::DNewCluster(Cluster &cluster, Point &point)
{
    double dist = 0.0;

    for(unsigned long i=0; i<cluster.pointCount(); ++i)
        dist += point.distance(*cluster.point(i));

    return dist / cluster.pointCount();
}

double HDiana::DCluster(Cluster &cluster, Point &point)
{
    return point.avgClusterDist;
}

//ACCESSORS/MUTATORS

Matrix HDiana::U() const
{
    Matrix u(_numObjects, _clusterSet.size());

    for(unsigned long k=0; k<_clusterSet.size(); ++k)
        for(unsigned long i=0; i<_clusterSet.at(k).pointCount(); ++i)
            _clusterSet.at(k).point(i)->setGroupNo(k);

    for(unsigned int i=0; i<_numObjects; ++i)
        u[i][_pointSet.at(i).groupNo()] = 1;

    return u;
}

#endif
