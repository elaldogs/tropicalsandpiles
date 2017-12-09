//============================================================================
// Name        : linearsandpile.cpp
// Author      : Aldo Guzmán-Sáenz, Nikita Kalinin
// Version     :
// Copyright   :
// Description : This program computes a linearized version of a sandpile, modeled
//				 with tropical curves
// to compile: g++ -std=c++11 -O3 linearsandpile.cpp -o linearsandpile
//============================================================================

#include <iostream>
#include <stack>
#include <vector>
#include <stdlib.h>
#include <fstream>
#include <map>
#include <set>
#include <exception>
#include <random>

using namespace std;

//Global variables and macros =================================================

#define CRITICAL 4								// Value at which points become unstable

int avalanchesize,volume,K;
map<pair<int, int>, int> current; 				// Map (dictionary) used to store the current (active) monomials as pairs and coefficients
map<pair<int, int>, set<int> > tocheck;         // map: monomial->unstable points contained in the part where this monomial is the minimal one
set<int> checkset;                              // indices of unstable points to check
vector< int> processed;                         // to estimate the size of the avalanche
vector<pair<int, int> > unstable;		        // Vector used to store the unstable points in the grid
int m, n;										// Sizes of the grid, m = # of rows; n = # of columns
pair<int, int> upper, lower, dexter, sinister;	// Current extreme monomials in each direction of the grid (dexter=right, sinister=left in latin)
int nunstable;									// Number of initial "unstable" cells in the grid
mt19937 mt;
uniform_int_distribution<int> dist2, dist1;
vector<pair<int, int> > curve;                  // tropical curve defined as the set where the min is attained twice
int curvesize;                                  // number of pixels in the curve
int touchboundary;                              // is -1 if the avalanche touched the boundary, 1 otherwise
int seed;

//=============================================================================

pair<int, int> operator+(						// Function to add pairs using the operator +
    const pair<int, int>& x,
    const pair<int, int>& y)
{
    return make_pair(x.first + y.first, x.second + y.second);
}

int ih(pair<int, int> index) 					// Index Helper function to convert from pairs of indices to a single index
{
    return index.first * m + index.second;
}

pair<int, int> ih(int index)					// Overload of ih to convert from an index to a pair of indices
{
    return make_pair(index / n, index % n);
}

int minimum(const vector<int>& collection)		// Returns the minimum value in collection
{
    int result = *collection.begin();
    for (vector<int>::const_iterator i = collection.begin(); i < collection.end();
            ++i)
    {
        if (*i < result)
        {
            result = *i;
        }
    }
    return result;
}

int coefficient(const pair<int, int>& element)  // initial coefficient of the monomial (element.first,element.second)
{
    int temp1[] =
    {
        element.first * 0 + element.second * 0,
        element.first * n + element.second * 0,
        element.first * 0 + element.second * m,
        element.first * n + element.second * m
    };
    vector<int> temp2(temp1, temp1 + sizeof(temp1) / sizeof(int));
    return -minimum(temp2);
}

vector<pair<int, int> > minimalmonomials(const pair<int, int>& cell) // The minimal polynomial at cell
{
    vector<int> temp1;
    map<pair<int, int>, int>::iterator i;
    vector<pair<int, int> > result;
    
    for (i = current.begin(); i != current.end(); ++i)
    {
        temp1.push_back(i->first.first * cell.first +
                        i->first.second * cell.second +
                        i->second);
    }
    int m = minimum(temp1);
    for (i = current.begin(); i != current.end(); ++i)
    {
        int val = i->first.first * cell.first +
                  i->first.second * cell.second +
                  i->second;
        if (val == m)
        {
            result.push_back(i->first);
        }
    }
    return result;
}

void add(pair<int, int> monomial)
{
    if (current.find(monomial) == current.end())
    {
        current[monomial] = coefficient(monomial);
    }
}
void operatorgp(const pair<int, int>& monomial, int pointnumber)
{
    bool flag = false;
    pair<int, int> temp1;
    vector<int> temp3;
    int temp2;
    if (monomial == upper)
    {
        upper = monomial + make_pair(0, 1);
        current[upper] = coefficient(upper);
        flag = true;
    }
    if (monomial == lower)
    {
        lower = monomial + make_pair(0, -1);
        current[lower] = coefficient(lower);
        flag = true;
    }
    if (monomial == sinister)
    {
        sinister = monomial + make_pair(-1, 0);
        current[sinister] = coefficient(sinister);
        flag = true;
    }
    if (monomial == dexter)
    {
        dexter = monomial + make_pair(1, 0);
        current[dexter] = coefficient(dexter);
        flag = true;
    }
    if (flag == true)
    {
        touchboundary = -1;
        for (int i = sinister.first; i != 0; ++i)
        {
            temp1 = make_pair(i,
                              int(-float(upper.second) / float(sinister.first) *
                                  float(i) +
                                  float(upper.second)));
            add(temp1);
            temp1 = make_pair(i,
                              int(-float(lower.second) / float(sinister.first) *
                                  float(i) +
                                  float(lower.second)));
            add(temp1);
        }
        for (int i = 0; i != dexter.first; ++i)
        {
            temp1 = make_pair(i,
                              int(-float(upper.second) / float(dexter.first) *
                                  float(i) +
                                  float(upper.second)));
            add(temp1);
            temp1 = make_pair(i,
                              int(-float(lower.second) / float(dexter.first) *
                                  float(i) +
                                  float(lower.second)));
            add(temp1);
        }
    }
    current.erase(monomial);
    temp1 = unstable[pointnumber];
    for (map<pair<int, int>, int>::iterator i = current.begin(); i != current.end(); ++i)
    {
        temp3.push_back(i->first.first * temp1.first +
                        i->first.second * temp1.second +
                        i->second);
    }
    temp2 = minimum(temp3);
    current[monomial] = temp2 -
                        monomial.first * temp1.first	-
                        monomial.second * temp1.second;
    
    vector<pair<int, int> > newmon = minimalmonomials(temp1);
    set<int> temp4;
    for(auto it : newmon)
    {
        if (tocheck.find(it)==tocheck.end())
        {
            temp4.insert(pointnumber);
            tocheck[it]=temp4;
        }
        else
        {
            tocheck[it].insert(pointnumber);
        }
    }
    
}
void init(int argc, char **argv)
{
    if (argc > 1)
    {
        if (argc == 5)
        {
            m = atoi(argv[1]);
            n = atoi(argv[2]);
            nunstable = atoi(argv[3]);
            mt19937 tempmt(atoi(argv[4]));
            seed = atoi(argv[4]);
            uniform_int_distribution<int> tempdist2(1,m-2), tempdist1(1,n-2);
            mt=tempmt;
            dist1=tempdist1;
            dist2=tempdist2;
        }
        else
        {
            cout << "Fatal error. Check number of parameters. Parameters should be: m,n,number of unstable points, seed";
            exit(-1);
        }
    }
    else
    {
        n = 1000;								// Default grid size
        m = n;									// By default, the grid is square
        nunstable = 900;                       // defaul the number of unstable points
        seed=2;                                 // default seed
        mt19937 tempmt(seed);
        uniform_int_distribution<int> tempdist2(1,m-2), tempdist1(1,n-2);
        mt=tempmt;
        dist1=tempdist1;
        dist2=tempdist2;
    }
    upper = make_pair(0, 1);
    lower = make_pair(0, -1);
    dexter = make_pair(1, 0);
    sinister = make_pair(-1, 0);
    current[make_pair(1, 0)] = coefficient(make_pair(1, 0));
    current[make_pair(-1, 0)] = coefficient(make_pair(-1, 0));
    current[make_pair(0, 1)] = coefficient(make_pair(0, 1));
    current[make_pair(0, -1)] = coefficient(make_pair(0, -1));
    current[make_pair(0, 0)] = coefficient(make_pair(0,0));
    unstable.resize(nunstable);
}

void pseudorelax()                          //Analogue of the relaxation function for sandpiles
{
    vector<pair<int, int> > monomial;
    int pointnumber; // index of the unstable point to relax
    for (int i=0;i< K+1 ; ++i)
    {
        processed.push_back(false);
    }
    while (!checkset.empty())
    {
        pointnumber = *(checkset.begin());
        checkset.erase(checkset.begin());   
        monomial = minimalmonomials(unstable[pointnumber]);
        if (monomial.size() == 1)
        {
            if (tocheck.find(monomial[0])!=tocheck.end())
            {
                for(auto it : tocheck[monomial[0]])
                {
                    checkset.insert(it);
                }
                tocheck[monomial[0]].clear();
            }
            operatorgp(monomial[0],pointnumber);
            ++volume;
            if (processed[pointnumber] == false)
            {
                ++avalanchesize;
                processed[pointnumber] = true;
            }
        }
        else
        {
            set<int> temp4;
            for(auto it : monomial)
            {
                if (tocheck.find(it)==tocheck.end())
                {
                    temp4.insert(pointnumber);
                    tocheck[it] = temp4;
                }
                else
                {
                    tocheck[it].insert(pointnumber);
                }
            }
        }
    }
    
}
void writeout()
{
    // Output of final state of the grid
    vector<pair<int, int> > temp1;
    int i = seed;
    std::string text = "./tsandpile/grid";
    //text += std::to_string(i);
    text += ".dat";
    for (int i = 0; i < n * m; ++i)
    {
        temp1 = minimalmonomials(ih(i));
        if (temp1.size() > 1)
        {
            curve.push_back(ih(i));
        }
    }
    curvesize = curve.size();

    string path(text);
    ofstream output(path.c_str(), ios::out | ofstream::binary);
    output.write(reinterpret_cast<const char *>(&n), sizeof(n));
    output.write(reinterpret_cast<const char *>(&nunstable), sizeof(nunstable));
    output.write(reinterpret_cast<const char *>(&curvesize),
                 sizeof(curvesize));
    for (vector<pair<int, int> >::iterator i = curve.begin(); i != curve.end();
         ++i)
    {
        output.write(reinterpret_cast<const char *>(&(i->first)),sizeof(i->first));
        output.write(reinterpret_cast<const char *>(&(i->second)),sizeof(i->second));
    }
    for (unsigned int i = 0; i < unstable.size(); ++i)
    {
        output.write(reinterpret_cast<const char *>(&unstable[i].first),
                     sizeof(unstable[i].first));
        output.write(reinterpret_cast<const char *>(&unstable[i].second),
                     sizeof(unstable[i].second));
    }
    output.close();
    
    // Output of map (i,j)->a_{i,j}
    text = "./tsandpile/active";
    //text += std::to_string(i);
    text += ".dat";
    path = text;
    int actlen = current.size();
    output.open(path.c_str(), ios::out | ofstream::binary);
    output.write(reinterpret_cast<const char *>(&actlen), sizeof(actlen));
    for (auto i = current.begin(); i != current.end(); ++i)
    {
        output.write(reinterpret_cast<const char *>(&(i->first.first)),sizeof(i->first.first));
        output.write(reinterpret_cast<const char *>(&(i->first.second)),sizeof(i->first.second));
        output.write(reinterpret_cast<const char *>(&(i->second)),sizeof(i->second));
    }
    output.close();
}

//============================================================================
// Parameters:
// m,n,number_of_added_points, seed
// -- m,n are the sides of the rectangular
// -- number_of_added_points,  number of initial unstable cells (at random positions)
// output:
// power_n_seed.txt -- sizes of the avalanches
// ...w.txt -- number of operations during avalanches
// ...a_00,a01,a10,a11,degree.txt -- files with such parameters of curves.
//============================================================================
int main(int argc, char **argv)
{
    init(argc,argv);
    
    string path("./tsandpile/power" + to_string(n) + "_" + to_string(nunstable) + "_" + to_string(seed) + ".txt");
    string pathw("./tsandpile/power" + to_string(n) + "_" + to_string(nunstable) + "_" + to_string(seed) + "w.txt");
    string patha00("./tsandpile/power" + to_string(n) + "_" + to_string(nunstable) + "_" + to_string(seed) + "a00.txt");
    string patha10("./tsandpile/power" + to_string(n) + "_" + to_string(nunstable) + "_" + to_string(seed) + "a10.txt");
    string patha01("./tsandpile/power" + to_string(n) + "_" + to_string(nunstable) + "_" + to_string(seed) + "a01.txt");
    string patha11("./tsandpile/power" + to_string(n) + "_" + to_string(nunstable) + "_" + to_string(seed) + "a11.txt");
    string pathdegree("./tsandpile/power" + to_string(n) + "_" + to_string(nunstable) + "_" + to_string(seed) + "degree.txt");
    
    ofstream output(path.c_str(), ios::out );
    ofstream outputw(pathw.c_str(), ios::out );
    ofstream outputa00(patha00.c_str(), ios::out );
    ofstream outputa10(patha10.c_str(), ios::out );
    ofstream outputa01(patha01.c_str(), ios::out );
    ofstream outputa11(patha11.c_str(), ios::out );
    ofstream outputdegree(pathdegree.c_str(), ios::out );
    
    K = 0;
    for (int i=0; i < nunstable; ++i)
    {
        touchboundary = 1;
        unstable[i].first = dist1(mt);
        unstable[i].second = dist2(mt);
        checkset.insert(K);
        ++K;
        avalanchesize = 0;
        volume = 0;
        pseudorelax();
        processed.clear();
        //cout<< i<<"\t"<<operationscount<<"\t"<<volume<<endl;
        output<<to_string(float(touchboundary)*float(avalanchesize)/float(K))+",";
        outputw<<to_string(float(touchboundary)*float(volume)/float(K))+",";
        outputa00<<to_string(current[make_pair(0, 0)])+",";
        outputa10<<to_string(current[make_pair(1, 0)])+",";
        outputa01<<to_string(current[make_pair(0, 1)])+",";
        outputa11<<to_string(current[make_pair(1, 1)])+",";
        outputdegree<<to_string(upper.second+dexter.first)+",";
    }
    // final check
    for (int i=0; i < nunstable; ++i)
    {
        vector<pair<int, int> > temp1;
        temp1 = minimalmonomials(unstable[i]);
        if (temp1.size() == 1)
        {
            std::cout << "did NOT stabilized!" << std::endl;
            std::cout <<i << std::endl;
        }
        
    }
    output.close();
    outputw.close();
    
    // produces a file with data with actual tropical curve to draw
    writeout();
    return 0;
}
