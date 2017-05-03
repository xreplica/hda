/*
 * Fuzzy c-means class implementation
 *
 * TODO:
 * -
 *
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/algorithms/fuzzycmeans.cpp $
 */

#ifndef __FUZZYCMEANS_CPP__
#define __FUZZYCMEANS_CPP__

#include "src/algorithms/fuzzycmeans.h"

double FuzzyCMeans::_variation_threshold = 0.01;

//CONSTRUCTORS

FuzzyCMeans::FuzzyCMeans()
{
    srand((unsigned int)time(NULL));
}

FuzzyCMeans::FuzzyCMeans(string f, int c, unsigned int fuzzy) : _file(f), _numClusters(c), _fuzzyConstant(fuzzy)
{
    srand((unsigned int)time(NULL));
}

//FUNCTIONS

void FuzzyCMeans::printU() const
{
    for (unsigned int i=0; i<_numObjects; ++i)
    {
        for(unsigned int k=0; k<_numClusters; ++k)
        {
            cout<<fixed<<setprecision(8)<<_u[i][k]<<"\t";
        }
        cout<<endl;
    }
    cout<<endl;
}

void FuzzyCMeans::printClusters() const
{
    for (unsigned int k=0; k<_numClusters; ++k)
    {
        cout<<_clusterSet.at(k).toString()<<endl;
    }
    cout<<endl;
}

void FuzzyCMeans::initialize(bool read)
{
    if(read)
    {
        string line;
        //int lineCount = 0;
        ifstream inFile(_file.c_str(), ios::in);

        if(!inFile.is_open())
        {
            cout<<"Could not open file \""<<_file<<"\""<<endl;
            exit(1);
        }

        while(getline(inFile, line))
        {
            Point p = Point(line);
            _pointSet.push_back(p);
        }

        inFile.close();
    }

    _clusterSet.clear();

    for (unsigned int k=0; k<_numClusters; ++k)
    {
        _clusterSet.push_back(*(new Cluster()));
    }

    //_epsilon = DBL_MAX;
    _numObjects = _pointSet.size();
    //_numClusters = _clusterSet.size();

    _u.reset(_numObjects, _numClusters);
    //srand((unsigned int)time(NULL));

    for (unsigned int i=0; i<_numObjects; ++i)
    {
        double total = 0;

        for(unsigned int k=0; k<_numClusters-1; ++k)
        {
            _u[i][k] = ((double)rand() / RAND_MAX) * (1-total);
            total += _u[i][k];
        }
        _u[i][_numClusters-1] = (1-total);
    }

    Point::metric.setData(_pointSet);
}

void FuzzyCMeans::run()
{
    _epsilon = DBL_MAX;
    while(_epsilon >= _variation_threshold)
    {
        updateC();
        Matrix oldU = _u;

        updateU();
        _epsilon = computeEpsilon(oldU);
    }
}

void FuzzyCMeans::newRun(double threshold)
{
    _epsilon = DBL_MAX;
    while(_epsilon >= _variation_threshold)
    {
        newUpdateC(threshold);
        Matrix oldU = _u;
        newUpdateU();
        _epsilon = computeEpsilon(oldU);
    }
}

void FuzzyCMeans::updateC()
{
    for(unsigned int k=0; k<_numClusters; ++k)
    {
        double *v1 = new double[Point::dimension()];
        double s1 = 0.0;

        for(unsigned int j=0; j<Point::dimension(); ++j)
        {
            v1[j] = 0.0;
        }

        for(unsigned int i=0; i<_numObjects; ++i)
        {
            for(unsigned int j=0; j<Point::dimension(); ++j)
            {
                v1[j] += _pointSet.at(i).feature(j) * pow(_u[i][k], (double)_fuzzyConstant);
            }

            s1 += pow(_u[i][k], (double)_fuzzyConstant);
        }

        for(unsigned int j=0; j<Point::dimension(); ++j)
        {
            _clusterSet.at(k).setFeature(j, (v1[j] / s1));
        }

        //delete[] v1;
    }
}

void FuzzyCMeans::newUpdateC(double threshold)
{
    for(unsigned int k=0; k<_numClusters; ++k)
    {
        double *sum = new double[Point::dimension()];
        unsigned int count = 0;

        for(unsigned int j=0; j<Point::dimension(); ++j)
        {
            sum[j] = 0.0;
        }

        for(unsigned int i=0; i<_numObjects; ++i)
        {
            double membership = _u[i][k];

            if(membership > threshold)
            {
                bool greatestMember = true;
                for(unsigned int k2=0; k2<_numClusters; ++k2)
                {
                    if(_u[i][k2] > membership)
                    {
                        greatestMember = false;
                        break;
                    }
                }

                if(greatestMember)
                {
                    for(unsigned int j=0; j<Point::dimension(); ++j)
                    {
                        sum[j] += _pointSet.at(i).feature(j);// * pow(_u[i][k], (double)_fuzzyConstant);
                    }

                    ++count;
                }
            }
        }

        if(count != 0)
        /*{
            removeCluster(k);
        }
        else*/
        {
            for(unsigned int j=0; j<Point::dimension(); ++j)
            {
                _clusterSet.at(k).setFeature(j, (sum[j] / count));
            }
        }

        //delete[] sum;
    }
}

void FuzzyCMeans::updateU()
{
    for(unsigned int i=0; i<_numObjects; ++i)
    {
        //unsigned int maxOutIndex;
        //bool maxOut = false;

        for(unsigned int k=0; k<_numClusters; ++k)
        {
            double sum = 0;

            /*if(_pointSet.at(i).distance(_clusterSet.at(k)) == 0.0)
            {
                cout<<"ZERO DISTANCE ERROR"<<endl;
                maxOut = true;
                maxOutIndex = k;
                break;
            }*/

            for(unsigned int k2=0; k2<_numClusters; ++k2)
            {
                sum += pow(_pointSet.at(i).distance(_clusterSet.at(k)) / _pointSet.at(i).distance(_clusterSet.at(k2)), (2.0 / (_fuzzyConstant - 1)));
            }

            _u[i][k] = 1 / sum;
        }

        /*if(maxOut)
        {
            for(unsigned int k=0; k<_numClusters; ++k)
            {
                if(k == maxOutIndex)
                    _u[i][k] = 1.0 - (0.1);
                else
                    _u[i][k] = 0.1/(_numClusters-1);
            }
        }*/
    }
}

void FuzzyCMeans::newUpdateU()
{
    for(unsigned int i=0; i<_numObjects; ++i)
    {
        for(unsigned int k=0; k<_numClusters; ++k)
        {
            double maxDist = 0;

            for(unsigned int i2=0; i2<_numObjects; ++i2)
            {
                double dist = _pointSet.at(i2).distance(_clusterSet.at(k));
                if(dist > maxDist)
                    maxDist = dist;
            }

            _u[i][k] = (maxDist - _pointSet.at(i).distance(_clusterSet.at(k))) / maxDist;
        }
    }
}

double FuzzyCMeans::computeEpsilon(const Matrix oldU)
{
    double maxDiff = 0;

    for(unsigned int i=0; i<_u.rows(); ++i)
    {
        for(unsigned int k=0; k<_u.cols(); ++k)
        {
            if(abs(_u[i][k] - oldU[i][k]) > maxDiff)
            {
                maxDiff = abs(_u[i][k] - oldU[i][k]);
            }
        }
    }

    return maxDiff;
}

//===========================================================================================================================================================

bool FuzzyCMeans::fuzzyMerge()
{
cout<<"FuzzyMerge:"<<endl;
    double* sigmas = new double[_numClusters];

    for(unsigned int k=0; k<_numClusters; ++k)
    {
        double sigma = 0;
        for(unsigned int i=0; i<_numObjects; ++i)
        {
            sigma += pow(_clusterSet.at(k).distance(_pointSet.at(i)) * _u[i][k], 2.0);
        }
        sigma = sqrt(sigma/_numObjects);

        sigmas[k] = sigma;
    }

    for(unsigned int k1=0; k1<_numClusters; ++k1)
    {
        for(unsigned int k2=0; k2<_numClusters; ++k2)
        {
            if(k1 == k2)
                continue;

            if(_clusterSet.at(k1).distance(_clusterSet.at(k2)) < sqrt(2.0) * (sigmas[k1] + sigmas[k2]))
            {
                removeCluster(k2);

                //delete[] sigmas;

                return true;
            }
        }
    }

    //delete[] sigmas;

    return false;
}

bool FuzzyCMeans::fuzzySplit()
{
cout<<"FuzzySplit"<<endl;
    Point furthestPoint;
    double furthestDistance = 0.0;
    unsigned int furthestPointCluster = 0;

    double* sigmas = new double[_numClusters];

    for(unsigned int k=0; k<_numClusters; ++k)
    {
        double sigma = 0;
        for(unsigned int i=0; i<_numObjects; ++i)
        {
            unsigned int memberCluster = 0;
            double highestMembership = 0.0;

            sigma += pow(_clusterSet.at(k).distance(_pointSet.at(i)) * _u[i][k], 2.0);

            for(unsigned int k2=0; k2<_numClusters; ++k2)
            {
                if(_u[i][k2] > highestMembership)
                {
                    highestMembership = _u[i][k2];
                    memberCluster = k2;
                }
            }

            if(memberCluster == k)
            {
                if(_pointSet.at(i).distance(_clusterSet.at(k)) > furthestDistance)
                {
                    furthestPointCluster = k;
                    furthestPoint = _pointSet.at(i);
                    furthestDistance = _pointSet.at(i).distance(_clusterSet.at(k));
                }

            }

        }

        sigma = sqrt(sigma/_numObjects);
        sigmas[k] = sigma;
    }

    if(furthestDistance > 5*sigmas[furthestPointCluster])
    {
        addCluster(Cluster(furthestPoint));

        //delete[] sigmas;

        return true;
    }

    //delete[] sigmas;

    return false;
}

void FuzzyCMeans::removeCluster(unsigned int index)
{
    Matrix newU(_numObjects, _numClusters-1);
    unsigned int k2 = 0;

    for(unsigned int k=0; k<_numClusters; ++k)
    {
        if(k==index)
            continue;

        for(unsigned int i=0; i<_numObjects; ++i)
        {
            newU[i][k2] = _u[i][k] + ((1 / (_numClusters-1)) * _u[i][index]);
        }

        ++k2;
    }

    _u.swap(newU);
    _clusterSet.erase(_clusterSet.begin()+index);
    _numClusters--;
}

void FuzzyCMeans::addCluster(Cluster newCluster)
{
    Matrix newU(_numObjects, _numClusters+1);

    _u.swap(newU);
    _clusterSet.push_back(newCluster);
    ++_numClusters;
}

#endif
