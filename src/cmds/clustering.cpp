/*
 * Clustering guided command line tool
 *
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/cmds/clustering.cpp $
 */

 #include <iostream>
 #include <cstdlib>
 #include <string>
 #include <sstream>
 #include <fstream>

 using namespace std;

 int main()
 {
    int algorithm = -1;
    string filename = "";
    int vDimension = -1;
    int metric = -1;
    string s;
    stringstream ss("");

    #if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__) && !defined(__CYGWIN32__)
    ss<<"";
    #else
    ss<<"./";
    #endif

    cout<<"At any point enter \"0\" to quit"<<endl;

    //executable
    while(algorithm > 6 || algorithm < 0)
    {
        cout<<endl<<"Choose a clustering algorithm: "<<endl
            <<"1) k-means+ (y-means)"<<endl
            <<"2) fuzzy c-means"<<endl
            <<"3) racamaic"<<endl
            <<"4) entropy-based fuzzy clustering"<<endl
            <<"5) heuristic diana"<<endl
            <<"6) entropy-heuristic diana"<<endl;
        getline(cin, s);
        if(s.length() == 0)
            algorithm = -1;
        else
            algorithm = atoi(s.c_str());
        if(algorithm == 0)
            exit(0);
        if(algorithm > 6 || algorithm < 0)
            cout<<"Invalid selection."<<endl;
    }

    //-f
    while(filename.length() < 1)
    {
        cout<<endl<<"Enter the filename of the dataset: "<<endl;
        getline(cin, filename);
        if(filename == "0")
            exit(0);
        ifstream file(filename.c_str());
        if(!file.good())
        {
            cout<<"invalid file name"<<endl;
            filename = "";
        }
    }

    //-v
    while(vDimension < 0)
    {
        cout<<endl<<"Enter the data point dimension:"<<endl;
        getline(cin, s);
        if(s.length() == 0)
            vDimension = -1;
        else
            vDimension = atoi(s.c_str());
        if(vDimension == 0)
            exit(0);
        if(vDimension < 0)
            cout<<"Invalid input value."<<endl;
    }

    //-d
    while(metric > 3 || metric < 0)
    {
        cout<<endl<<"Choose distance metric: (default = euclidian)"<<endl
            <<"1) euclidian"<<endl
            <<"2) mahalanobis"<<endl
            <<"3) kernel (gaussian)"<<endl;
        getline(cin, s);
        if(s.length() == 0)
            metric = 1;
        else
            metric = atoi(s.c_str());
        if(metric == 0)
            exit(0);
        if(metric > 3 || metric < 0)
            cout<<"Invalid selection."<<endl;
    }

    switch(algorithm)
    {
    case 1:
        ss<<"kmp ";
        break;
    case 2 :
        ss<<"fcm ";
        break;
    case 3 :
        ss<<"rac ";
        break;
    case 4 :
        ss<<"efc ";
        break;
    case 5 :
        ss<<"hda ";
        break;
    case 6 :
        ss<<"ehda ";
        break;
    default:
        cerr<<"NO ALGORITHM SET!"<<endl;
        exit(1);
        break;
    }

    ss<<"-f "<<filename<<" -v "<<vDimension<<" -d "<<(metric - 1)<<" ";

    //-c
    if(algorithm < 3)
    {
        int numClust = -1;

        while(numClust < 0)
        {
            cout<<endl<<"Enter the initial number of clusters:"<<endl;
            getline(cin, s);
            numClust = atoi(s.c_str());
            if(numClust == 0)
                exit(0);
            if(numClust < 0)
                cout<<"Invalid value."<<endl;
        }

        ss<<"-c "<<numClust<<" ";
    }

    //-m
    if(algorithm < 4)
    {
        int fuzzyConst = -1;

        while(fuzzyConst < 0)
        {
            cout<<endl<<"Enter the fuzzyfication constant (default = 2)"<<endl;
            getline(cin, s);
            if(s.length() == 0)
                fuzzyConst = 2;
            else
                fuzzyConst = atoi(s.c_str());
            if(fuzzyConst == 0)
                exit(0);
            if(fuzzyConst < 0)
                cout<<"Invalid value."<<endl;
        }

        ss<<"-m "<<fuzzyConst<<" ";
    }

    //-b & -y
    if(algorithm == 4)
    {
        double beta = -1.0;
        int gamma = -1;

        //-b
        while(beta < 0.0 || beta > 1.0)
        {
            cout<<endl<<"Enter minimum similarity - beta (default = 0.7)"<<endl;
            getline(cin, s);
            if(s.length() == 0)
                beta = 0.7;
            else
                beta = atof(s.c_str());
            if(beta == 0.0)
                exit(0);
            if(beta < 0.0 || beta > 1.0)
                cout<<"Invalid value."<<endl;
        }

        //-y
        while(gamma < 0)
        {
            cout<<endl<<"Enter outlier threshold (minimum # of points in cluster)- gamma (default = 5)"<<endl;
            getline(cin, s);
            if(s.length() == 0)
                gamma = 5;
            else
                gamma = atoi(s.c_str());
            if(gamma == 0)
                exit(0);
            if(gamma < 0)
                cout<<"Invalid value."<<endl;
        }

        ss<<"-b "<<beta<<" -y "<<gamma<<" ";
    }

    //-u
    cout<<endl<<"Display fuzzy partition? (y/n)"<<endl;
    getline(cin, s);
    if(s.length() > 0 && (s.at(0) == 'y' || s.at(0) == 'Y'))
        ss<<"-u ";

    //-k
    cout<<endl<<"Display clusters? (y/n)"<<endl;
    getline(cin, s);
    if(s.length() > 0 && (s.at(0) == 'y' || s.at(0) == 'Y'))
        ss<<"-k ";

    //-i
    cout<<endl<<"Display validity indexes? (y/n)"<<endl;
    getline(cin, s);
    if(s.length() > 0 && (s.at(0) == 'y' || s.at(0) == 'Y'))
        ss<<"-i";

    cout<<endl<<"clustering..."<<endl<<endl;

    int result = system(ss.str().c_str());

    if(result != 0)
    {
        cout<<endl<<"To quickly run this clustering again, you can use the following command:"
            <<endl<<ss.str()<<endl<<endl;
    }
    else
    {
        cout<<"an error has occured"<<endl;
    }

    return 0;
 }

 //system(command)
