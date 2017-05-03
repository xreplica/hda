/*
 * Cluster validation class implementation
 *
 * TODO:
 * -
 *
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/core/validation.cpp $
 */

#ifndef __VALIDATION_CPP__
#define __VALIDATION_CPP__

#include "src/core/validation.h"

double Validation::minInterClusterDistance(const vector<Cluster> &clusterSet){
    double minDist = DBL_MAX;

    for(unsigned int k1=0; k1<clusterSet.size(); ++k1)
    {
        for(unsigned int k2=0; k2<clusterSet.size(); ++k2)
        {
            if(k1 == k2)
                continue;

            double dist = clusterSet.at(k1).distance(clusterSet.at(k2));
            if(dist<minDist)
            {
                minDist = dist;
            }
        }
    }

    return minDist;
}

double Validation::maxInterClusterDistance(const vector<Cluster> &clusterSet){
    double maxDist = 0.0;

    for(unsigned int k1=0; k1<clusterSet.size(); ++k1)
    {
        for(unsigned int k2=0; k2<clusterSet.size(); ++k2)
        {
            if(k1 == k2)
                continue;

            double dist = clusterSet.at(k1).distance(clusterSet.at(k2));
            if(dist>maxDist)
            {
                maxDist = dist;
            }
        }
    }

    return maxDist;
}

double Validation::avgInterClusterDistance(const vector<Cluster> &clusterSet){
    double totDist = 0;
    //unsigned int num = 0;

    for(unsigned int k1=0; k1<clusterSet.size(); ++k1)
    {
        for(unsigned int k2=0; k2<clusterSet.size(); ++k2)
        {
            if(k1 == k2)
                continue;

            totDist  = clusterSet.at(k1).distance(clusterSet.at(k2));
            //++num;
        }
    }

    return totDist/((clusterSet.size() - 1) * (clusterSet.size() - 1));
}

Cluster Validation::avgClusterCenter(const vector<Cluster> &clusterSet){
    Cluster avg;

    for(unsigned int j=0; j<Point::dimension(); ++j){
        avg.setFeature(j, 0.0);
    }

    for(unsigned int i=0; i<clusterSet.size(); ++i)
    {
        for(unsigned int j=0; j<Point::dimension(); ++j){
            avg.setFeature(j, avg.feature(j) + clusterSet.at(i).feature(j));
        }
    }

    for(unsigned int j=0; j<Point::dimension(); ++j){
        avg.setFeature(j, avg.feature(j)/clusterSet.size());
    }

    return avg;
}

Point Validation::avgPoint(const vector<Point> &pointSet){
    Point avg;

    for(unsigned int j=0; j<Point::dimension(); ++j){
        avg.setFeature(j, 0.0);
    }

    for(unsigned int i=0; i<pointSet.size(); ++i)
    {
            for(unsigned int j=0; j<Point::dimension(); ++j){
                avg.setFeature(j, avg.feature(j) + pointSet.at(i).feature(j));
            }
    }

    for(unsigned int j=0; j<Point::dimension(); ++j){
        avg.setFeature(j, avg.feature(j)/pointSet.size());

    }

    return avg;
}

Point Validation::patternVariance(const vector<Point> &pointSet)
{
    Point avg = avgPoint(pointSet);
    Point variance;

    for(unsigned int j=0; j<Point::dimension(); ++j){
        avg.setFeature(j, 0.0);
    }

    for(unsigned int i=0; i<pointSet.size(); ++i)
    {
            for(unsigned int j=0; j<Point::dimension(); ++j){

                variance.setFeature(j, variance.feature(j) + pow(pointSet.at(i).feature(j) - avg.feature(j), 2.0));
            }
    }

    for(unsigned int j=0; j<Point::dimension(); ++j){
        avg.setFeature(j, avg.feature(j)/pointSet.size());
    }

    return variance;
}

double Validation::scat(const vector<Cluster> &clusterSet, const vector<Point> &pointSet, const Matrix &u, unsigned int fuzzyConstant)
{
    Point sigmaX = patternVariance(pointSet);
    double normSigmaX = 0;

    for(unsigned int j=0; j<Point::dimension(); ++j)
    {
        normSigmaX += sqrt(pow(sigmaX.feature(j), 2.0));
    }

    return leastSquaredError(clusterSet, pointSet, u, fuzzyConstant) / (clusterSet.size() * pointSet.size() * normSigmaX);
}

double Validation::dist(const vector<Cluster> &clusterSet, const vector<Point> &pointSet)
{
    double dMax = maxInterClusterDistance(clusterSet);
    double dMin = minInterClusterDistance(clusterSet);
    double s1 = 0.0;

    for(unsigned int k=0; k<clusterSet.size(); ++k)
    {
        double s2 = 0.0;

        for(unsigned int z=0; z<clusterSet.size(); ++z)
        {
            s2 += clusterSet.at(k).distance(clusterSet.at(z));
        }

        s1 += 1/s2;
    }

    return (dMax / dMin) * s1;;
}

double Validation::sc1(const vector<Cluster> &clusterSet, const vector<Point> &pointSet, const Matrix &u, unsigned int fuzzyConstant){
    Cluster avg = avgClusterCenter(clusterSet);

    double s1 = 0.0;
    for(unsigned int k=0; k<clusterSet.size(); ++k)
    {
        s1 += pow(clusterSet.at(k).distance(avg), 2.0);
    }

    double s2 = 0.0;
    for(unsigned int k=0; k<clusterSet.size(); ++k)
    {
        double s3 = 0.0, s4 = 0.0;

        for(unsigned int i=0; i<pointSet.size(); ++i)
        {
            s3 += pow(u[i][k], (double)fuzzyConstant) * pow(pointSet.at(i).distance(clusterSet.at(k)), 2.0);
        }

        for(unsigned int i=0; i<pointSet.size(); ++i)
        {
            s4 += u[i][k];
        }

        s2 += s3 / s4;
    }

    return (s1 / clusterSet.size()) / s2;
}

double Validation::sc2(const vector<Cluster> &clusterSet, const vector<Point> &pointSet, const Matrix &u){
    double s1 = 0.0;
    for(unsigned int k1=0; k1<clusterSet.size()-1; ++k1)
    {
        for(unsigned int k2=0; k2<clusterSet.size()-k1 ; ++k2)
        {
            double s2 = 0.0;
            for(unsigned int i=0; i<pointSet.size(); ++i)
            {
                s2 += pow(min(u[i][k1], u[i][k2]), 2.0);
            }

            double s3 = 0.0;
            for(unsigned int i=0; i<pointSet.size(); ++i)
            {
                s3 += min(u[i][k1], u[i][k2]);
            }

            s1 += s2 /s3;
        }
    }

    double s4 = 0.0;
    for(unsigned int i=0; i<pointSet.size(); ++i)
    {
        double max = 0.0;
        for(unsigned int k=0; k<clusterSet.size(); ++k)
        {
            if(u[i][k] > max)
                max = u[i][k];
        }

        s4 += max * max;
    }

    double s5 = 0.0;
    for(unsigned int i=0; i<pointSet.size(); ++i)
    {
        double max = 0.0;
        for(unsigned int k=0; k<clusterSet.size(); ++k)
        {
            if(u[i][k] > max)
                max = u[i][k];
        }

        s5 += max;
    }

    return s1 / (s4 / s5);
}

Matrix Validation::fuzzyCovarianceMatrix(unsigned int numClusters, const Cluster &cluster, unsigned int k, const vector<Point> &pointSet, const Matrix &u, unsigned int fuzzyConstant)
{
    Matrix m1(Point::dimension());
    double s1 = 0.0;
    Matrix mCluster(vector<double>(cluster.featureSet()));

    for(unsigned int i=0; i<pointSet.size(); ++i)
    {
        Matrix deltaM = Matrix(vector<double>(pointSet.at(i).featureSet())) - mCluster;

        m1 += (deltaM * deltaM.transposed()) * pow(u[i][k], (double)fuzzyConstant);
    }

    for(unsigned int i=0; i<pointSet.size(); ++i)
    {
        s1 += pow(u[i][k], (double)fuzzyConstant);
    }

    return m1 / s1;
}

double Validation::avgDistanceWithCluster(const Point &p, const Cluster &c)
{
    double dist = 0.0;

    for(unsigned int i=0; i<c.pointCount(); ++i)
    {
        dist += p.distance(*c.point(i));
    }

    return dist/c.pointCount();
}

/*double Validation::xieBeniIndex(vector<Cluster> clusterSet, unsigned int numClusters, vector<Point> pointSet, unsigned int numObjects, Matrix u, unsigned int fuzzyConstant){
    double s1 = 0;

    for(unsigned int k=0; k<numClusters; ++k)
    {
        for(unsigned int i=0; i<numObjects; ++i){
            s1 += pow(u[i][k], fuzzyConstant) * pow(pointSet.at(i).distance(clusterSet.at(k)), 2.0);
        }
    }

    return s1/(numObjects*pow(minInterClusterDistance(clusterSet, numClusters), 2.0));
}*/

double Validation::leastSquaredError(const vector<Cluster> &clusterSet, const vector<Point> &pointSet, const Matrix &u, unsigned int fuzzyConstant)
{
    double jm = 0;

    for(unsigned int i=0; i<pointSet.size(); ++i)
    {
        for(unsigned int k=0; k<clusterSet.size(); ++k){
            jm += pow(u[i][k], (double)fuzzyConstant) * pow(pointSet.at(i).distance(clusterSet.at(k)), 2.0);
        }
    }

    return jm;
}

double Validation::xieBeniIndex(const vector<Cluster> &clusterSet, const vector<Point> &pointSet, const Matrix &u, unsigned int fuzzyConstant)
{
    return leastSquaredError(clusterSet, pointSet, u, fuzzyConstant) / (pointSet.size() * pow(minInterClusterDistance(clusterSet), 2.0));
}

double Validation::fukuyamaSugenoIndex(const vector<Cluster> &clusterSet, const vector<Point> &pointSet, const Matrix &u, unsigned int fuzzyConstant)
{
    Cluster avg = avgClusterCenter(clusterSet);
    double km = 0;

    for(unsigned int i=0; i<pointSet.size(); ++i)
    {
        for(unsigned int k=0; k<clusterSet.size(); ++k){
            km += pow(u[i][k], (double)fuzzyConstant) * pow(pointSet.at(i).distance(avg), 2.0);
        }
    }

    return leastSquaredError(clusterSet, pointSet, u, fuzzyConstant) - km;
}

double Validation::kwonIndex(const vector<Cluster> &clusterSet, const vector<Point> &pointSet, const Matrix &u, unsigned int fuzzyConstant)
{
    Cluster avg = avgClusterCenter(clusterSet);
    double punish = 0;

    for(unsigned int k=0; k<clusterSet.size(); ++k)
    {
        punish += pow(clusterSet.at(k).distance(avg), 2.0);
    }

    return (leastSquaredError(clusterSet, pointSet, u, fuzzyConstant) + (punish / clusterSet.size())) / pow(minInterClusterDistance(clusterSet), 2.0);
}

double Validation::composeWithinBetweenIndex(const vector<Cluster> &clusterSet, const vector<Point> &pointSet, const Matrix &u, unsigned int fuzzyConstant)
{
    double alpha = 1.0; // = dist(cMax);

    return (alpha*scat(clusterSet, pointSet, u, fuzzyConstant))+dist(clusterSet, pointSet);
}

double Validation::zahidIndex(const vector<Cluster> &clusterSet, const vector<Point> &pointSet, const Matrix &u, unsigned int fuzzyConstant)
{
    return sc1(clusterSet, pointSet, u, fuzzyConstant) - sc2(clusterSet, pointSet, u);
}

double Validation::fuzzyHypervolume(const vector<Cluster> &clusterSet, const vector<Point> &pointSet, const Matrix &u, unsigned int fuzzyConstant)
{
    double fhv = 0.0;

    for(unsigned int k=0; k<clusterSet.size(); ++k)
    {
        fhv += sqrt(fuzzyCovarianceMatrix(clusterSet.size(), clusterSet.at(k), k, pointSet, u, fuzzyConstant).determinant());   //getDeterminant only valid for SPD Matrices
    }

    return fhv;
}

double Validation::randIndex(const vector<Cluster> &clusterSet1, const vector<Cluster> &clusterSet2, const vector<Point> &pointSet)
{
    int a=0, b=0, c=0, d=0;

    multimap<const Point*, unsigned int> x, y;

    for(unsigned int k=0; k<clusterSet1.size(); ++k)
    {
        for(unsigned int i=0; i<clusterSet1.at(k).pointCount(); ++i)
        {
            x.insert(pair<const Point*, int>(clusterSet1.at(k).point(i), k));
        }
    }

    for(unsigned int k=0; k<clusterSet2.size(); ++k)
    {
        for(unsigned int i=0; i<clusterSet2.at(k).pointCount(); ++i)
        {
            x.insert(pair<Point*, int>(clusterSet2.at(k).point(i), k));
        }
    }

    for(unsigned int i=0; i<pointSet.size()-1; ++i)
    {
        for(unsigned int j=i+1; j<pointSet.size(); ++j)
        {
            if(x.find(&pointSet.at(i)) == x.find(&pointSet.at(j)))      //same in X
            {
                if(y.find(&pointSet.at(i)) == y.find(&pointSet.at(j)))  //same in Y
                    ++a;
                else                                                    //different in Y
                    ++c;
            }
            else                                                        //different in X
            {
                if(y.find(&pointSet.at(i)) == y.find(&pointSet.at(j)))  //same in Y
                    ++d;
                else                                                    //different in Y
                    ++b;
            }
        }
    }

    return (a + b) / (a + b + c + d);
}

double Validation::PBMIndex(const vector<Cluster> &clusterSet, const vector<Point> &pointSet, const Matrix &u)
{
    double index = 0.0, E1 = 0.0, Ek = 0.0, Dk = 0.0;
    Cluster avg = avgClusterCenter(clusterSet);

    for(unsigned int i=0; i<clusterSet.size(); ++i)
    {
        for(unsigned int j=0; j<clusterSet.size(); ++j)
        {
            double dist = clusterSet.at(i).distance(clusterSet.at(j));

            if(dist > Dk)
                Dk = dist;
        }
    }

    for(unsigned int k=0; k<clusterSet.size(); ++k)
    {
        for(unsigned int i=0; i<pointSet.size(); ++i)
        {
            Ek += u[i][k] * clusterSet.at(k).distance(pointSet.at(i));
        }
    }

    for(unsigned int i=0; i<pointSet.size(); ++i)
    {
        E1 += avg.distance(pointSet.at(i));
    }

    index = ((1.0 / clusterSet.size()) * (E1 / Ek) * Dk);

    return index * index;
}

double Validation::silhouette(const vector<Cluster> &clusterSet, const vector<Point> &pointSet)
{
    double s = 0.0;

    for(unsigned int k=0; k<clusterSet.size(); ++k)
    {
        for(unsigned int i=0; i<clusterSet.at(k).pointCount(); ++i)
        {
            Point p = *clusterSet.at(k).point(i);
            double a = avgDistanceWithCluster(p, clusterSet.at(k));
            double b = DBL_MAX;

            for(unsigned int m=0; m<clusterSet.size(); ++m)
            {
                if(m != k)
                {
                    double avgDist = avgDistanceWithCluster(p, clusterSet.at(m));

                    if(avgDist < b)
                        b = avgDist;
                }
            }

            if(a<b)
                s += 1 - (a / b);
            else if(b<a)
                s += (b / a) - 1;
        }
    }

    return s / pointSet.size();
}

double Validation::silhouette(const vector<Cluster> &clusterSet, const vector<Point> &pointSet, const Matrix &u)
{
    vector<Cluster> newClusterSet;

    for(unsigned int k=0; k<clusterSet.size(); ++k)
    {
        newClusterSet.push_back(clusterSet.at(k));
    }

    for(unsigned int i=0; i<pointSet.size(); ++i)
    {
        int maxIndex = -1;
        double maxVal = 0.0;

        for(unsigned int k=0; k<newClusterSet.size(); ++k)
        {
            if(u[i][k] > maxVal)
            {
                maxVal = u[i][k];
                maxIndex = k;
            }
        }
        newClusterSet.at(maxIndex).addPoint((Point*)&pointSet.at(i));
    }

    return silhouette(newClusterSet, pointSet);
}

#endif
