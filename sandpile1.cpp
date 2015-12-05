//============================================================================
// Name        : sandpile1.cpp
// Author      : Aldo Guzmán-Sáenz
// Version     :
// Copyright   : 
// Description : This program simulates an abelian sandpile's evolution.
//============================================================================

#include <iostream>
#include <stack>
#include <vector>
#include <stdlib.h>
#include <fstream>


using namespace std;

#define CRITICAL 4								// Value at which points become unstable

stack< pair<int,int> > unstable1, unstable2; 	// Stacks used to store the current unstable cells in the grid
vector< pair<int,int> > initialunstable;		// Vector used to store the initial unstable cells in the grid
int n;											// Size of a size of the grid
vector<int> grid;								// Our sandpile, stored as an integer grid, the size is specified in main()

pair<int,int> operator+(						// Function to add pairs using the operator +
		const pair<int, int>& x,
		const pair<int, int>& y)
{
    return make_pair(x.first+y.first, x.second+y.second);
}

int ih(pair<int,int> index) 					// Index Helper function to convert from pairs of indices to a single index
{
	return index.first*n + index.second;
}

pair<int,int> ih(int index)						// Overload of ih to convert from an index to a pair of indices
{
	return make_pair(index/n,index%n);
}

vector< pair<int,int> > adjacent(const pair<int,int>& current) // Returns the positions of the 4 adjacent cells in our grid
{
	vector< pair<int,int> > result;
	result.push_back(current+make_pair(-1,0));
	result.push_back(current+make_pair(0,1));
	result.push_back(current+make_pair(1,0));
	result.push_back(current+make_pair(0,-1));
	return result;
}

bool issink(const pair<int, int>& index) 	// Determines if the cell at index is a sink (not affected by toppling)
{											// In our case, only cells at the boundary are sinks
	return (index.first>=n || index.first<0 || index.second>=n || index.second<0);
}

void relax()				// Main relaxation function (uses two stacks to keep track of unstable cells in our grid)
{
	vector< pair<int,int> > neighbors(CRITICAL);
	unsigned int i;
	pair<int,int> current;
	while(!unstable1.empty() || !unstable2.empty())
	{
		while(!unstable1.empty())
		{
			current=unstable1.top();
			unstable1.pop();
			if (grid[ih(current)] >= CRITICAL)
			{
				grid[ih(current)] -= CRITICAL;
				neighbors = adjacent(current);
				for(i=0;i<neighbors.size();i++)
					if (issink(neighbors[i])==false)
					{
						grid[ih(neighbors[i])]+=1;
						if (grid[ih(neighbors[i])] >= CRITICAL)
							unstable2.push(neighbors[i]);
					}
			}
		}
		while(!unstable2.empty())
		{
			current=unstable2.top();
			unstable2.pop();
			if (grid[ih(current)] >= CRITICAL)
			{
				grid[ih(current)] -= CRITICAL;
				neighbors = adjacent(current);
				for(i=0;i<neighbors.size();i++)
					if (issink(neighbors[i])==false)
					{
						grid[ih(neighbors[i])]+=1;
						if (grid[ih(neighbors[i])] >= CRITICAL)
							unstable1.push(neighbors[i]);
					}
			}
		}
	}
}

//============================================================================
// TODO: Implement parameter parsing
// Parameters:
// --size N 			size of the grid
// --numberofgrains		number of initial unstable cells (at random positions)
//============================================================================
int main(int argc, char **argv) {
	string path("./grid.dat");
	ofstream output(path.c_str(), ios::out | ofstream::binary);
	n=1000;											// Default grid size
	//n=atoi(argv[1]);								// Read grid size from command line if available
	int nunstable=10;
	grid.resize(n*n,3);								// Fill our grid with 3's
	srand(1);
	initialunstable.resize(nunstable);
	for(unsigned int i=0; i< initialunstable.size();i++)
	{
		initialunstable[i].first=rand()%n;
		initialunstable[i].second=rand()%n;
	}
	for(unsigned int i=0; i< initialunstable.size();i++)
	{
		grid[ih(initialunstable[i])]=4;
		unstable1.push(initialunstable[i]);
	}
	relax();
	output.write(reinterpret_cast<const char *>(&n), sizeof(n));
	output.write(reinterpret_cast<const char *>(&nunstable), sizeof(nunstable));
	for (unsigned int i=0;i<grid.size();i++)
	{
		output.write(reinterpret_cast<const char *>(&grid[i]),sizeof(grid[i]));
	}
	for (unsigned int i=0;i<initialunstable.size();i++)
	{
		output.write(reinterpret_cast<const char *>(&initialunstable[i].first),sizeof(initialunstable[i].first));
		output.write(reinterpret_cast<const char *>(&initialunstable[i].second),sizeof(initialunstable[i].second));
	}
	return 0;
}
