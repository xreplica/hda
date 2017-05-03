/*
 * EH-DIANA class implementation
 *
 * TODO:
 * -
 *
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/algorithms/ehdiana.cpp $
 */

#ifndef __EHDIANA_CPP__
#define __EHDIANA_CPP__

#include <cmath>
#include <cfloat>

#include "src/algorithms/ehdiana.h"

//CONSTRUCTORS

EHDiana::EHDiana() {}

EHDiana::EHDiana(string f) : _file(f) {}

//FUNCTIONS

void EHDiana::printClusters()
{
    for(unsigned int k=0; k<_clusterSet.size(); ++k)
    {
        cout<<"CLUSTER "<<k<<endl;
        for(unsigned int i=0; i<_clusterSet.at(k).getNumPoints(); ++i)
        {
            cout<<(*_clusterSet.at(k).getPoint(i)).toString()<<endl;
        }
    }
}

void EHDiana::initialize()
{
    double avgD = 0.0;
    string line;
    ifstream inFile(_file.c_str(), ios::in);

    if(!inFile.is_open()){
        cout<<"Could not open file \""<<_file<<"\""<<endl;
        exit(1);
    }

    while(getline(inFile, line))
    {
        Point p = Point(line);
        _pointSet.push_back(p);
    }

    inFile.close();

    _normalizedPointSet = _pointSet;//normalize(_pointSet, 0.0, 1.0);

    Point::metric.setData(_normalizedPointSet);

    Cluster c = Cluster("0,0,0,0");

    for(unsigned int i=0; i<_normalizedPointSet.size(); ++i)
    {
        _normalizedPointSet.at(i).setGroupNo(0);
        c.addPoint(&_normalizedPointSet.at(i));
        _normalizedPointSet.at(i).setGroupNo(0);//
    }
    _clusterSet.push_back(c);

    _numObjects = _normalizedPointSet.size();
    _numClusters = 1;

    for(unsigned int i=0; i<_normalizedPointSet.size(); ++i)
    {
        for(unsigned int j=i+1; j<_normalizedPointSet.size(); ++j)
        {
            avgD += _normalizedPointSet.at(i).distance(_normalizedPointSet.at(j));
        }
    }

    avgD /= _normalizedPointSet.size() * (_normalizedPointSet.size() - 1) / 2;

    _alpha = -log(0.5) / avgD;

    for(unsigned int i=0; i<_normalizedPointSet.size(); ++i)
    {
        double sum = 0.0;
        for(unsigned int j=0; j<_normalizedPointSet.size(); ++j)
        {
            if(i == j)
                continue;

            double Sij = similarity(_normalizedPointSet.at(i), _normalizedPointSet.at(j));
            if(Sij == 1)
            {
                continue;       //Sij == 1 wil cause division by zero ( log2(1) == 0 ), Ei == 0 when Sij == 1, so we skip it's calculation as it would simply add 0 to the sum
            }

            sum += (Sij * (log(Sij) / log(2))) + (1 - Sij) * (log(1 - Sij) / log(2));
        }

        _normalizedPointSet.at(i).entropy = -sum;
    }
}

void EHDiana::run()
{
    vector<Cluster> optimalClusterSet = _clusterSet;
    //double prevVariation = 0.0;
    //double currVariation = 0.0;
    double minEntropy = DBL_MAX;

    for(int loop=_numObjects; loop>0; --loop)
    {
        //double preSumEntropy = 0.0;
        double postSumEntropy = 0.0;

        unsigned int minIndex = 0;
        //TEST (invert min to max)
        double minDiameter = DBL_MAX;//clusterDiameter(_clusterSet.at(0)); //diameter here is similarity of the two most dissimilar objects
        //double minDiameter = 0.0;
        for(unsigned int k=0; k<_clusterSet.size(); ++k)
        {
            //TEST
            if(clusterDiameter(_clusterSet.at(k)) < minDiameter && _clusterSet.at(k).getNumPoints() > 1)
            {
                minIndex = k;
                minDiameter = clusterDiameter(_clusterSet.at(k));
            }
            /*if(EclusterDiameter(_clusterSet.at(k)) > minDiameter && _clusterSet.at(k).getNumPoints() > 1)
            {
                minIndex = k;
                minDiameter = EclusterDiameter(_clusterSet.at(k));
            }*/
        }
        //cout<<"index"<<minIndex<<endl;
        if(_clusterSet.at(minIndex).getNumPoints() < 2)
        {
            cout<<"over splitting... "<<_clusterSet.size()<<" clusters."<<endl;
            break;
        }

//-----------------------------------------------------BEGIN
        //get sum of the avg entropy of each cluster == sum(avgEntropy(cluster)) for each cluster
        /*
        for(unsigned int k=0; k<_clusterSet.size(); ++k)
        {
            double avgD = 0.0;
            double alpha = 0.0;

            for(unsigned int i=0; i<_clusterSet.at(k).getNumPoints(); ++i)
            {
                for(unsigned int j=i+1; j<_clusterSet.at(k).getNumPoints(); ++j)
                {
                    avgD += _clusterSet.at(k).getPoint(i)->distance(*_clusterSet.at(k).getPoint(j));
                }
            }

            avgD /= _clusterSet.at(k).getNumPoints() * (_clusterSet.at(k).getNumPoints() - 1) / 2;

            alpha = -log(0.5) / avgD;

            for(unsigned int i=0; i<_clusterSet.at(k).getNumPoints(); ++i)
            {
                double sum = 0.0;
                for(unsigned int j=0; j<_clusterSet.at(k).getNumPoints(); ++j)
                {
                    if(i == j)
                        continue;

                    double Sij = similarity(*_clusterSet.at(k).getPoint(i), *_clusterSet.at(k).getPoint(j), alpha);
                    if(Sij == 1)
                    {
                        continue;       //Sij == 1 wil cause division by zero ( log2(1) == 0 ), Ei == 0 when Sij == 1, so we skip it's calculation as it would simply add 0 to the sum
                    }

                    sum += (Sij * (log(Sij) / log(2))) + (1 - Sij) * (log(1 - Sij) / log(2));
                }

                _clusterSet.at(k).getPoint(i)->internalEntropy = -sum;
            }

            double avgE = 0.0;
            for(unsigned int i=0; i<_clusterSet.at(k).getNumPoints(); ++i)
            {
                avgE += _clusterSet.at(k).getPoint(i)->internalEntropy;
            }

            avgE /= _clusterSet.at(k).getNumPoints();

            preSumEntropy += avgE;
        }
        */

        //get avg of avg cluster entropy
        double avg = 0.0;
        for(unsigned int k=0; k<_clusterSet.size(); ++k)
        {
            double avgC = 0.0;

            for(unsigned int i=0; i<_clusterSet.at(k).getNumPoints(); ++i)
            {
                avgC += _clusterSet.at(k).getPoint(i)->entropy;
            }

            avgC /= _clusterSet.at(k).getNumPoints();
            avg += avgC;
        }
        avg /= _clusterSet.size();
        //preSumEntropy = avg;
//-----------------------------------------------------END

        if(_clusterSet.at(minIndex).getNumPoints() >= 2)
        {
            //TEST
            splitCluster(_clusterSet.at(minIndex));
            //EsplitCluster(_clusterSet.at(minIndex));
        }
        else
        {
            cout<<"tried to split a single-object cluster"<<endl;
            break;
        }

//-----------------------------------------------------BEGIN
        //get sum of the avg entropy of each cluster == sum(avgEntropy(cluster)) for each cluster
        /*
        for(unsigned int k=0; k<_clusterSet.size(); ++k)
        {
            double avgD = 0.0;

            for(unsigned int i=0; i<_clusterSet.at(k).getNumPoints(); ++i)
            {
                for(unsigned int j=i+1; j<_clusterSet.at(k).getNumPoints(); ++j)
                {
                    avgD += _clusterSet.at(k).getPoint(i)->distance(*_clusterSet.at(k).getPoint(j));
                }
            }

            avgD /= _clusterSet.at(k).getNumPoints() * (_clusterSet.at(k).getNumPoints() - 1) / 2;

            for(unsigned int i=0; i<_clusterSet.at(k).getNumPoints(); ++i)
            {
                double sum = 0.0;
                for(unsigned int j=0; j<_clusterSet.at(k).getNumPoints(); ++j)
                {
                    if(i == j)
                        continue;

                    double Sij = similarity(*_clusterSet.at(k).getPoint(i), *_clusterSet.at(k).getPoint(j));
                    if(Sij == 1)
                    {
                        continue;       //Sij == 1 wil cause division by zero ( log2(1) == 0 ), Ei == 0 when Sij == 1, so we skip it's calculation as it would simply add 0 to the sum
                    }

                    sum += (Sij * (log(Sij) / log(2))) + (1 - Sij) * (log(1 - Sij) / log(2));
                }

                _clusterSet.at(k).getPoint(i)->internalEntropy = -sum;
            }

            double avgE = 0.0;
            for(unsigned int i=0; i<_clusterSet.at(k).getNumPoints(); ++i)
            {
                avgE += _clusterSet.at(k).getPoint(i)->internalEntropy;
            }

            avgE /= _clusterSet.at(k).getNumPoints();

            postSumEntropy += avgE;
        }
        */

        //get avg of avg cluster entropy
        avg = 0.0;
        for(unsigned int k=0; k<_clusterSet.size(); ++k)
        {
            double avgC = 0.0;

            for(unsigned int i=0; i<_clusterSet.at(k).getNumPoints(); ++i)
            {
                avgC += _clusterSet.at(k).getPoint(i)->entropy;
            }

            avgC /= _clusterSet.at(k).getNumPoints();
            avg += avgC;
        }
        avg /= _clusterSet.size();
        postSumEntropy = avg;
//-----------------------------------------------------END

        /*currVariation = postSumEntropy - preSumEntropy;
        cout<<"prev:"<<prevVariation<<endl;
        cout<<"curr:"<<currVariation<<endl;
        if(currVariation < 0)//if clustering is better
        {
            cout<<"optimal = current at: "<<_clusterSet.size()<<" members"<<endl;
            optimalClusterSet = _clusterSet;
        }
        else if(prevVariation < 0)
        {
            loop = 7;
        }
        cout<<"loop: "<<loop<<endl;

        prevVariation = currVariation;*/

        //currVariation = postSumEntropy - preSumEntropy;
        //cout<<"prev:"<<prevVariation<<endl;
        //cout<<"curr:"<<currVariation<<endl;
        //cout<<"post:"<<postSumEntropy<<" min: "<<minEntropy<<endl;
        if(postSumEntropy < minEntropy)//if clustering is better
        {
            //cout<<"optimal = current at: "<<_clusterSet.size()<<" members"<<endl;
            optimalClusterSet = _clusterSet;
            minEntropy = postSumEntropy;
            loop = 7;
        }

        //prevVariation = currVariation;
    }

    _clusterSet = optimalClusterSet;
    cout<<"finished with "<<_clusterSet.size()<<" clusters."<<endl;
}

vector<Point> EHDiana::normalize(vector<Point> pointSet, double min, double max)
{
    vector<Point> normalizedSet;
    double mins[Point::getDimension()], maxs[Point::getDimension()];

    for(unsigned int j=0; j<Point::getDimension(); ++j)
    {
        mins[j] = maxs[j] = pointSet.at(0).getFeature(j);
    }

    for(unsigned int i=1; i<pointSet.size(); ++i)
    {
        for(unsigned int j=0; j<Point::getDimension(); ++j)
        {
            if(pointSet.at(i).getFeature(j) < mins[j])
                mins[j] = pointSet.at(i).getFeature(j);
            if(pointSet.at(i).getFeature(j) > maxs[j])
                maxs[j] = pointSet.at(i).getFeature(j);
        }
    }

    for(unsigned int i=0; i<pointSet.size(); ++i)
    {
        string pointVals = "";

        for(unsigned int j=0; j<Point::getDimension(); ++j)
        {
            stringstream ss;
            ss << (pointSet.at(i).getFeature(j) - mins[j]) / (maxs[j] - mins[j]);
            pointVals = pointVals.append(ss.str()).append(",");
        }

        pointVals = pointVals.substr(0, pointVals.size()-1);

        normalizedSet.push_back(Point(pointVals));
    }

    return normalizedSet;
}

double EHDiana::similarity(Point pointA, Point pointB)
{
    return exp(-_alpha * pointA.distance(pointB));
}

double EHDiana::similarity(Point pointA, Point pointB, double alpha)
{
    return exp(-alpha * pointA.distance(pointB));
}

double EHDiana::avgDistance(Point point, Cluster cluster)
{
    double avgDist = 0.0;

    for(unsigned int i=0; i<cluster.getNumPoints(); ++i)
    {
        avgDist += point.distance((*cluster.getPoint(i)));
    }

    return (1 / (cluster.getNumPoints() - 1)) * avgDist;
}

double EHDiana::EavgDistance(Point point, Cluster cluster)
{
    double avgDist = 0.0;

    for(unsigned int i=0; i<cluster.getNumPoints(); ++i)
    {
        avgDist += similarity(point, (*cluster.getPoint(i)));
    }

    return (1 / (cluster.getNumPoints() - 1)) * avgDist;
}

double EHDiana::clusterDiameter(Cluster cluster)
{
    if(cluster.getNumPoints() < 2)
        return 0.0;

    double diameter = DBL_MAX;

    for(unsigned int i=0; i<cluster.getNumPoints(); ++i)
    {
        for(unsigned int j=0; j<cluster.getNumPoints(); ++j)
        {
            if(i == j)
                continue;

            if(similarity(*cluster.getPoint(i), *cluster.getPoint(j)) < diameter)
                diameter = similarity(*cluster.getPoint(i), *cluster.getPoint(j));
        }
    }

    return diameter;
}

double EHDiana::EclusterDiameter(Cluster cluster)
{
    if(cluster.getNumPoints() < 2)
        return 0.0;

    double diameter = DBL_MAX;

    for(unsigned int i=0; i<cluster.getNumPoints(); ++i)
    {
        for(unsigned int j=0; j<cluster.getNumPoints(); ++j)
        {
            if(i == j)
                continue;

            if(similarity(*cluster.getPoint(i), *cluster.getPoint(j)) < diameter)
                diameter = similarity(*cluster.getPoint(i), *cluster.getPoint(j));
        }
    }

    return diameter;
}

void EHDiana::splitCluster(Cluster &cluster)
{
    Cluster newCluster;
    double maxDistance = 0.0;
    unsigned int maxIndex = 0;
    bool change = true;

    if(cluster.getNumPoints() == 2)
    {
        newCluster.addPoint(cluster.getPoint(1));
        cluster.getPoint(1)->setGroupNo(_clusterSet.size());//
        cluster.removePoint(1);
        _clusterSet.push_back(newCluster);
        return;
    }

    for(unsigned int i=0; i<cluster.getNumPoints(); ++i)
    {
        if(avgDistance(*cluster.getPoint(i), cluster) > maxDistance)
        {
            maxDistance = avgDistance((*cluster.getPoint(i)), cluster);
            maxIndex = i;
        }
    }

    newCluster.addPoint(cluster.getPoint(maxIndex));
    cluster.getPoint(maxIndex)->setGroupNo(_clusterSet.size());//
    cluster.removePoint(maxIndex);

    while(change && cluster.getNumPoints() > 1)
    {
        change = false;
        for(unsigned int i=0; i<cluster.getNumPoints(); ++i)
        {
            if(DNewCluster(newCluster, *cluster.getPoint(i)) < DCluster(cluster, *cluster.getPoint(i)))
            {
                newCluster.addPoint(cluster.getPoint(i));
                cluster.getPoint(i)->setGroupNo(_clusterSet.size());//
                cluster.removePoint(i);
                change = true;
                break;
            }
        }
    }

    _clusterSet.push_back(newCluster);
}

//TEST
void EHDiana::EsplitCluster(Cluster &cluster)
{
    double avgD = 0.0;
    Cluster newCluster;
    double maxDistance = 0.0;
    unsigned int maxIndex = 0;
    bool change = true;

    if(cluster.getNumPoints() == 2)
    {
        newCluster.addPoint(cluster.getPoint(1));
        cluster.getPoint(1)->setGroupNo(_clusterSet.size());//
        cluster.removePoint(1);
        _clusterSet.push_back(newCluster);
        return;
    }

    for(unsigned int i=0; i<cluster.getNumPoints(); ++i)
    {
        for(unsigned int j=i+1; j<cluster.getNumPoints(); ++j)
        {
            avgD += cluster.getPoint(i)->distance(*cluster.getPoint(j));
        }
    }

    avgD /= cluster.getNumPoints() * (cluster.getNumPoints() - 1) / 2;

    _alpha = -log(0.5) / avgD;

    for(unsigned int i=0; i<cluster.getNumPoints(); ++i)
    {
        double sum = 0.0;
        for(unsigned int j=0; j<cluster.getNumPoints(); ++j)
        {
            if(i == j)
                continue;

            double Sij = similarity(*cluster.getPoint(i), *cluster.getPoint(j));
            if(Sij == 1)
            {
                continue;       //Sij == 1 wil cause division by zero ( log2(1) == 0 ), Ei == 0 when Sij == 1, so we skip it's calculation as it would simply add 0 to the sum
            }

            sum += (Sij * (log(Sij) / log(2))) + (1 - Sij) * (log(1 - Sij) / log(2));
        }

        cluster.getPoint(i)->internalEntropy = -sum;
    }

    for(unsigned int i=0; i<cluster.getNumPoints(); ++i)
    {
        if(EavgDistance(*cluster.getPoint(i), cluster) > maxDistance)
        {
            maxDistance = EavgDistance((*cluster.getPoint(i)), cluster);
            maxIndex = i;
        }
    }

    newCluster.addPoint(cluster.getPoint(maxIndex));
    cluster.getPoint(maxIndex)->setGroupNo(_clusterSet.size());//
    cluster.removePoint(maxIndex);

    while(change && cluster.getNumPoints() > 1)
    {
        change = false;
        for(unsigned int i=0; i<cluster.getNumPoints(); ++i)
        {
            if(EDNewCluster(newCluster, *cluster.getPoint(i)) < EDCluster(cluster, *cluster.getPoint(i)))
            {
                newCluster.addPoint(cluster.getPoint(i));
                cluster.getPoint(i)->setGroupNo(_clusterSet.size());//
                cluster.removePoint(i);
                change = true;
                break;
            }
        }
    }

    _clusterSet.push_back(newCluster);
}

double EHDiana::DNewCluster(Cluster &cluster, Point &point)
{
    double dist = 0.0;

    for(unsigned int i=0; i<cluster.getNumPoints(); ++i)
    {
        dist += point.distance(*cluster.getPoint(i));
    }

    return dist / cluster.getNumPoints();
}

double EHDiana::DCluster(Cluster &cluster, Point &point)
{
    double dist = 0.0;

    for(unsigned int i=0; i<cluster.getNumPoints(); ++i)
    {
        dist += point.distance(*cluster.getPoint(i));
    }

    return dist / (cluster.getNumPoints() - 1);
}

//TEST
double EHDiana::EDNewCluster(Cluster &cluster, Point &point)
{
    double dist = 0.0;

    for(unsigned int i=0; i<cluster.getNumPoints(); ++i)
    {
        dist += similarity(point, (*cluster.getPoint(i)));
    }

    return dist / cluster.getNumPoints();
}

//TEST
double EHDiana::EDCluster(Cluster &cluster, Point &point)
{
    double dist = 0.0;

    for(unsigned int i=0; i<cluster.getNumPoints(); ++i)
    {
        dist += similarity(point, (*cluster.getPoint(i)));
    }

    return dist / (cluster.getNumPoints() - 1);
}

Matrix EHDiana::getU() const
{
    Matrix u(_numObjects, _clusterSet.size());

    for(unsigned int k=0; k<_clusterSet.size(); ++k)
    {
        for(unsigned int i=0; i<_clusterSet.at(k).getNumPoints(); ++i)
        {
            _clusterSet.at(k).getPoint(i)->setGroupNo(k);
        }
    }

    for(unsigned int i=0; i<_normalizedPointSet.size(); ++i)
    {
        u[i][_normalizedPointSet.at(i).getGroupNo()] = 1;
    }

    return u;
}

#endif

