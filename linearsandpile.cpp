//============================================================================
// Name        : linearsandpile.cpp
// Author      : Aldo Guzmán-Sáenz
// Version     :
// Copyright   :
// Description : This program computes a linearized version of a sandpile, modeled
//				 with tropical curves
//============================================================================

#include <iostream>
#include <stack>
#include <vector>
#include <stdlib.h>
#include <fstream>
#include <map>
#include <exception>

using namespace std;

//Global variables and macros =================================================

#define CRITICAL 4								// Value at which points become unstable

map<pair<int, int>, int> current; 				// Map (dictionary) used to store the current (active) monomials as pairs and coefficients
vector<pair<int, int> > initialunstable;		// Vector used to store the initial unstable cells in the grid
int m, n;										// Sizes of the grid, m = # of rows; n = # of columns
pair<int, int> upper, lower, dexter, sinister;	// Current extreme monomials in each direction of the grid (dexter=right, sinister=left in latin)
int nunstable;									// Number of initial "unstable" cells in the grid

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

int coefficient(const pair<int, int>& element)  // Coefficient at (element.first,element.second)
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

void operatorgp(const pair<int, int>& cell, const pair<int, int>& monomial)
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
        for (int i = sinister.first; i != 0; ++i)
        {
            temp1 = make_pair(i,
                              int(-float(upper.second) / float(sinister.first) *
                                  float(i) +
                                  float(upper.second)));
            if (current.find(temp1) == current.end())
            {
                current[temp1] = coefficient(temp1);
            }
            temp1 = make_pair(i,
                              int(-float(lower.second) / float(sinister.first) *
                                  float(i) +
                                  float(lower.second)));
            if (current.find(temp1) == current.end())
            {
                current[temp1] = coefficient(temp1);
            }
        }
        for (int i = 0; i != dexter.first; ++i)
        {
            temp1 = make_pair(i,
                              int(-float(upper.second) / float(dexter.first) *
                                  float(i) +
                                  float(upper.second)));
            if (current.find(temp1) == current.end())
            {
                current[temp1] = coefficient(temp1);
            }
            temp1 = make_pair(i,
                              int(-float(lower.second) / float(dexter.first) *
                                  float(i) +
                                  float(lower.second)));
            if (current.find(temp1) == current.end())
            {
                current[temp1] = coefficient(temp1);
            }
        }
    }
    current.erase(monomial);
    for (map<pair<int, int>, int>::iterator i = current.begin(); i != current.end(); ++i)
    {
        temp3.push_back(i->first.first * cell.first +
                        i->first.second * cell.second +
                        i->second);
    }
    temp2 = minimum(temp3);
    current[monomial] = temp2 -
                        monomial.first * cell.first	-
                        monomial.second * cell.second;
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
            srand(atoi(argv[4]));
        }
        else
        {
            cout << "Fatal error. Check number of parameters.";
            exit(-1);
        }
    }
    else
    {
        n = 100;								// Default grid size
        m = n;									// By default, the grid is square
        nunstable = 250;
        srand(2);
    }
    //n=atoi(argv[1]);						// Read grid size from command line if available
    int temp2;
    vector<int> temp3;
    upper = make_pair(0, 1);
    lower = make_pair(0, -1);
    dexter = make_pair(1, 0);
    sinister = make_pair(-1, 0);
    current[make_pair(1, 0)] = coefficient(make_pair(1, 0));
    current[make_pair(-1, 0)] = coefficient(make_pair(-1, 0));
    current[make_pair(0, 1)] = coefficient(make_pair(0, 1));
    current[make_pair(0, -1)] = coefficient(make_pair(0, -1));
    initialunstable.resize(nunstable);
    for (unsigned int i = 0; i < initialunstable.size(); ++i)
    {
        initialunstable[i].first = rand() % (n-4)+2;
        initialunstable[i].second = rand() % (m-4)+2;
    }
    for (map<pair<int, int>, int>::iterator i = current.begin();
            i != current.end(); ++i)
    {
        temp3.push_back(
            i->first.first * initialunstable[0].first
            + i->first.second * initialunstable[0].second
            + i->second);
    }
    temp2 = minimum(temp3);
    current[make_pair(0, 0)] = temp2;
}

void pseudorelax()                          //Analogue of the relaxation function for sandpiles
{
    bool flag = true;
    vector<pair<int, int> > monomial;
    while (flag)
    {
        flag = false;
        for (unsigned int i = 0; i < initialunstable.size(); ++i)
        {
            monomial = minimalmonomials(initialunstable[i]);
            if (monomial.size() == 1)
            {
                operatorgp(initialunstable[i], monomial[0]);
                flag = true;
            }
        }
    }
}

//============================================================================
// TODO: Implement parameter parsing with custom names
// Parameters:
// --size N 			size of the grid
// --numberofgrains		number of initial unstable cells (at random positions)
//============================================================================
int main(int argc, char **argv)
{
    vector<pair<int, int> > temp1;
    vector<pair<int, int> > outputmonomials;
    int monomialcount;
    string path("./grid.dat");
    ofstream output(path.c_str(), ios::out | ofstream::binary);
    init(argc,argv);
    pseudorelax();
    output.write(reinterpret_cast<const char *>(&n), sizeof(n));
    output.write(reinterpret_cast<const char *>(&nunstable), sizeof(nunstable));
    for (int i = 0; i < n * m; ++i)
    {
        temp1 = minimalmonomials(ih(i));
        if (temp1.size() > 1)
        {
            outputmonomials.push_back(ih(i));
        }
    }
    monomialcount = outputmonomials.size();
    output.write(reinterpret_cast<const char *>(&monomialcount),
                 sizeof(monomialcount));
    for (vector<pair<int, int> >::iterator i = outputmonomials.begin(); i != outputmonomials.end();
            ++i)
    {
        output.write(reinterpret_cast<const char *>(&(i->first)),sizeof(i->first));
        output.write(reinterpret_cast<const char *>(&(i->second)),sizeof(i->second));
    }
    for (unsigned int i = 0; i < initialunstable.size(); ++i)
    {
        output.write(reinterpret_cast<const char *>(&initialunstable[i].first),
                     sizeof(initialunstable[i].first));
        output.write(reinterpret_cast<const char *>(&initialunstable[i].second),
                     sizeof(initialunstable[i].second));
    }
    return 0;
}
