//============================================================================
// Name        : parallelsandpile.cpp
// Author      : Aldo Guzmán-Sáenz
// Version     :
// Copyright   :
// Description : This program simulates an abelian sandpile's evolution in a
//               parallel fashion
//============================================================================

//============================================================================
// Here is a drawing of a subgrid with its outer boundaries:
//          -------------------
//          ---+++++++++++++---
//          --+@############+--
//          --+#############+--
//          --+#############+--
//          --+#############+--
//          --+#############+--
//          ---+++++++++++++---
//          -------------------
//
// Where - are cells not in the subgrid, # are cells in the subgrid, and
// + are cells in the outer boundaries of the subgrid.
// @ is locationx,locationy (begins at (0,0))
//============================================================================


#include <iostream>
#include <stack>
#include <vector>
#include <stdlib.h>
#include <fstream>
#include <mpi.h>
#include <random>

using namespace std;

#define CRITICAL 4								// Value at which points become unstable
#define CRITICALMINUSONE 3
#define MASTERPROCESS 0

stack< pair<int,int> > unstable1, unstable2; 	// Stacks used to store the current unstable cells in the grid
vector< pair<int,int> > initialunstable;		// Vector used to store the initial unstable cells in the grid
vector< pair<int,int> > localinitialunstable;		// Vector used to store the initial unstable cells in each subgrid
vector<int> initialunstablesubgrids;
int partsx,partsy;                        // Number of parts to divide the grid in each axis, respectively
int n,m,nunstable;								// Size of the grid, m is size in x and n is size in y
vector<int> grid;								// Our sandpile, stored as an integer grid, the size is specified in main()
int stepx, stepy;
int world_rank;

mt19937 mt;
uniform_int_distribution<int> dist2, dist1;

vector<int> allneighborsleft;
vector<int> allneighborsright;
vector<int> allneighborstop;
vector<int> allneighborsbottom;

vector< vector<int> > allouterleft;
vector< vector<int> > allouterright;
vector< vector<int> > alloutertop;
vector< vector<int> > allouterbottom;


class subgrid                       // To greatly simplify calls for each thread, everything will be packed in a single object
{
    private:
        int locationx,locationy;
        int sizex,sizey;            // Size of our rectangular subgrid (all of them have the same sizes in this version)
    public:
        subgrid(){};                // Default constructor does nothing
        subgrid(int mi, int locationx, int locationy, int sizex, int sizey, int value);       // Constructor with options
        int getsizex() const;             // Getters for sizes
        int getsizey() const;
        pair<int,int> getlocation();
        int getlocationx() const;
        int getlocationy() const;
        stack< pair<int,int> > unstable1, unstable2;
        vector<int> actual;         // Actual values of cells in our subgrid
        vector<int> outerleft;      // Values of the outer boundary of our subgrid (divided in four parts for ease of use)
        vector<int> outerright;     // All of these should have sizes equal to the size of the corresponding side of our subgrid
        vector<int> outertop;       // The direction of iteration is left->right or top->bottom
        vector<int> outerbottom;
        int myid;                   // id of the subgrid's creator
        int neighborleft;           // Process id's of the neighbors of each subgrid, -1 means it doesn't have a neighbor in that
        int neighborright;          // direction
        int neighbortop;
        int neighborbottom;
        int isboundary(const pair<int,int>& p) const;
        int& operator ()  (const pair<int,int>& ncolrow) ;     // Some syntactic sugar to get/set the values of cells our subgrid
        int& operator ()  (const int ncol,const int nrow) ;
        subgrid& operator= (const subgrid& rhs);
};

subgrid::subgrid(int mi, int lx, int ly, int sx, int sy, int value)
{
    myid = mi;
    locationx=lx;
    locationy=ly;
    sizex=sx;
    sizey=sy;
    actual.resize(sx*sy,value);
    outerleft.resize(sy,0);
    outerright.resize(sy,0);
    outertop.resize(sx,0);
    outerbottom.resize(sx,0);
    if (lx/stepx < partsx - 1)
    {
        neighborright = myid + 1;
    }
    else
    {
        neighborright = -1;
    }
    if (lx/stepx >= 1)
    {
        neighborleft = myid - 1;
    }
    else
    {
        neighborleft = -1;
    }
    if (ly/stepy < partsy - 1)
    {
        neighborbottom = myid + partsx;
    }
    else
    {
        neighborbottom = -1;
    }
    if (ly/stepy >= 1)
    {
        neighbortop = myid - partsx;
    }
    else
    {
        neighbortop = -1;
    }
}

subgrid& subgrid::operator= (const subgrid&rhs)
{
        myid = rhs.myid;
        locationx=rhs.locationx;
        locationy=rhs.locationy;
        sizex=rhs.sizex;
        sizey=rhs.sizey;            // Size of our rectangular subgrid
        unstable1=rhs.unstable1;
        unstable2=rhs.unstable2;
        actual=rhs.actual;
        outerleft=rhs.outerleft;
        outerright=rhs.outerright;
        outertop=rhs.outertop;
        outerbottom=rhs.outerbottom;
        neighborbottom=rhs.neighborbottom;
        neighborleft=rhs.neighborleft;
        neighborright=rhs.neighborright;
        neighbortop=rhs.neighbortop;
        return *this;
}

int subgrid::isboundary(const pair <int,int>& p) const // Returns 0 if (x,y) belongs to outertop, 1 outerright, 2 outerbottom,
{                                                       // 3 outerleft, else it returns -1; this is in LOCAL coordinates of the
    int x=p.first;                                      // subgrid
    int y=p.second;
    if (y==-1 && x<sizex && 0<=x ) // Outertop
    {
        return 0;
    }
    if (x==sizex && y<sizey && 0<=y ) // Outerright
    {
        return 1;
    }
    if (y==sizey && x<sizex && 0<=x  ) // Outerbottom
    {
        return 2;
    }
    if (x==-1 && y<sizey && 0<=y ) // Outerleft
    {
        return 3;
    }
    return -1;
}

pair<int,int> subgrid::getlocation()
{
    return make_pair(locationx,locationy);
}

int& subgrid::operator() (const pair<int,int>& ncolrow)
{
    return actual[ncolrow.first + ncolrow.second*sizex];
}

int& subgrid::operator() (const int ncol,const int nrow)
{
    return actual[ncol + nrow*sizex];
}


int subgrid::getlocationx() const
{
    return locationx;
}

int subgrid::getlocationy() const
{
    return locationy;
}

int subgrid::getsizex() const
{
    return sizex;
}

int subgrid::getsizey() const
{
    return sizey;
}


int whichsubgrid(const pair<int, int> & cell ) // subgrids are numbered by putting all rows of subgrids in a single row and counting
{
    for(int i=0; i<partsy;i++)
    {
        for(int j=0; j<partsx; j++)
        {
            if(cell.first< stepx*(j + 1) && stepx*j<=cell.first && cell.second< stepy*(i + 1) && stepy*i<=cell.second)
            {
                return i*partsx + j;
            }
        }
    }
    return -1;
}

int nonzerocount(vector<int> x )   // How many nonzero elements we have in x
{
    int result=0;
    for(unsigned int i=0;i<x.size();++i)
    {
        if(x[i]!=0)
        {
            result++;
        }
    }
    return result;
}

pair<int,int> operator+(						// Function to add pairs using the operator +
		const pair<int, int>& x,
		const pair<int, int>& y)
{
    return make_pair(x.first+y.first, x.second+y.second);
}

int ih(pair<int,int> index) 					// Index Helper function to convert from pairs of indices to a single index
{
	return index.first*m + index.second;
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

bool issink(const subgrid& s, const pair<int, int>& index) 	// Determines if the cell at index is a sink (not affected by toppling)
{											// In our case, only cells at the boundary are sinks
	return (index.first+s.getlocationx()>=m || index.first+s.getlocationx()<0 || index.second+s.getlocationy()>=n || index.second+s.getlocationy()<0);
}

void checkcriticals(subgrid& s)
{
    for(int i=0; i<s.getsizey();i++)
    {
        for(int j=0; j<s.getsizex();j++)
        {
            if (s(make_pair(j,i))>=CRITICAL)
                s.unstable1.push(make_pair(j,i));
        }
    }
}

void relax(subgrid& s)				// Main relaxation function (uses two stacks to keep track of unstable cells in our grid)
{
	vector< pair<int,int> > neighbors(CRITICAL);
	unsigned int i;
	pair<int,int> current;
	fill(s.outerbottom.begin(),s.outerbottom.end(),0);
	fill(s.outertop.begin(),s.outertop.end(),0);
	fill(s.outerleft.begin(),s.outerleft.end(),0);
	fill(s.outerright.begin(),s.outerright.end(),0);
	while(!s.unstable1.empty() || !s.unstable2.empty())
	{
		while(!s.unstable1.empty())
		{
			current=s.unstable1.top();
			while (s(current) >= CRITICAL)
			{
				s(current) -= CRITICAL;
				neighbors = adjacent(current);
				for(i=0;i<neighbors.size();++i)
					if (issink(s,neighbors[i])==false)
					{
                        switch( s.isboundary(neighbors[i]) )
                        {
                            case 0:
                                s.outertop[neighbors[i].first]+=1;
                                break;
                            case 1:
                                s.outerright[neighbors[i].second]+=1;
                                break;
                            case 2:
                                s.outerbottom[neighbors[i].first]+=1;
                                break;
                            case 3:
                                s.outerleft[neighbors[i].second]+=1;
                                break;
                            case -1:
                                s(neighbors[i])+=1;
                                if (s(neighbors[i]) >= CRITICAL)
                                    s.unstable2.push(neighbors[i]);
                                break;
                        }
					}
			}
			s.unstable1.pop();
		}
		while(!s.unstable2.empty())
		{
			current=s.unstable2.top();
			while (s(current) >= CRITICAL)
			{
				s(current) -= CRITICAL;
				neighbors = adjacent(current);
				for(i=0;i<neighbors.size();++i)
					if (issink(s,neighbors[i])==false)
					{
					switch( s.isboundary(neighbors[i]) )
                        {
                            case 0:
                                s.outertop[neighbors[i].first]+=1;
                                break;
                            case 1:
                                s.outerright[neighbors[i].second]+=1;
                                break;
                            case 2:
                                s.outerbottom[neighbors[i].first]+=1;
                                break;
                            case 3:
                                s.outerleft[neighbors[i].second]+=1;
                                break;
                            case -1:
                                s(neighbors[i])+=1;
                                if (s(neighbors[i]) >= CRITICAL)
                                    s.unstable1.push(neighbors[i]);
                                break;
                        }
					}
			}
			s.unstable2.pop();
		}
	}
}

void sanitycheck(int argc, char **argv)
{
    int world_size;
    if (argc>1)
    {
        if (argc==7)
        {
            m=atoi(argv[1]);
            n=atoi(argv[2]);
            nunstable=atoi(argv[3]);
            mt19937 tempmt(atoi(argv[4]));
            uniform_int_distribution<int> tempdist2(1,m-2), tempdist1(1,n-2);
            mt=tempmt;
            dist1=tempdist1;
            dist2=tempdist2;
            partsx=atoi(argv[5]);
            partsy=atoi(argv[6]);
        }
        else
        {
            cout<<"Fatal error. Check number of arguments."<<endl;
            exit(-1);
        }
    }
    else
    {
        m=100;      // Default grid size. m is size in x and n is size in y
        n=100;
        nunstable=150;  // Default number of critical cells
        mt19937 tempmt(2);
        uniform_int_distribution<int> tempdist2(1,m-2), tempdist1(1,n-2);
        mt=tempmt;
        dist1=tempdist1;
        dist2=tempdist2;
        partsx=1;    // Default number of parts (no divisions)
        partsy=1;
    }
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    if(m%partsx!=0 || n%partsy!=0)
    {
        cout<<"Fatal error. Length of side of grid not divisible by number of parts."<<endl;
        exit(-1);
    }
    if(world_size!=partsx*partsy)
    {
        cout<<"Fatal error. Not enough processes for all parts, or not enough parts for all processes. "<<endl;
        cout<<"Check that number of processes is equal to partsx*partsy."<<endl;
        exit(-1);
    }
    stepx=m/partsx;
    stepy=n/partsy;
}

void init(subgrid& s)
{
	initialunstable.resize(nunstable);
	initialunstablesubgrids.resize(nunstable);
	for(unsigned int i=0; i < initialunstable.size();++i)
	{
		initialunstable[i].first=dist1(mt);
		initialunstable[i].second=dist2(mt);
	}
    for(unsigned int i=0; i < initialunstable.size();++i)
    {
        initialunstablesubgrids[i]=whichsubgrid(initialunstable[i]);
    }
    for (unsigned int j=0; j < initialunstable.size(); ++j)
    {
        if (0 == initialunstablesubgrids[j])
        {
            s(initialunstable[j])=CRITICAL;
        }
    }
    allneighborsbottom.resize(partsx*partsy);
    allneighborstop.resize(partsx*partsy);
    allneighborsleft.resize(partsx*partsy);
    allneighborsright.resize(partsx*partsy);
    for( int i=0; i < partsx*partsy; ++i)
    {
        if (i%partsx < partsx - 1)
        {
            allneighborsright[i] = i + 1;
        }
        else
        {
            allneighborsright[i] = -1;
        }
        if (i%partsx >= 1)
        {
            allneighborsleft[i] = i - 1;
        }
        else
        {
            allneighborsleft[i] = -1;
        }
        if (i/partsx < partsy - 1)
        {
            allneighborsbottom[i] = i + partsx;
        }
        else
        {
            allneighborsbottom[i] = -1;
        }
        if (i/partsx >= 1)
        {
            allneighborstop[i] = i - partsx;
        }
        else
        {
            allneighborstop[i] = -1;
        }
    }
    allouterbottom.resize(partsx*partsy);
    alloutertop.resize(partsx*partsy);
    allouterleft.resize(partsx*partsy);
    allouterright.resize(partsx*partsy);
    for( int i=0; i< partsx*partsy; ++i)
    {
        allouterbottom[i].resize(stepx);
        alloutertop[i].resize(stepx);
        allouterleft[i].resize(stepy);
        allouterright[i].resize(stepy);
    }
}

void addsubgridtototal(subgrid& total, const vector<int>& sub, const int lx, const int ly)
{
    for(int i=0;i<stepy;++i)
    {
        for(int j=0;j<stepx;++j)
        {
            total(ly+i,lx+j)+=sub[i*stepx + j];
        }
    }
}


void writeout(const subgrid& s)
{
    string path("./grid.dat");
	ofstream output(path.c_str(), ios::out | ofstream::binary);
    output.write(reinterpret_cast<const char *>(&n), sizeof(n));
	output.write(reinterpret_cast<const char *>(&nunstable), sizeof(nunstable));
	for (unsigned int i=0;i<(s.actual).size();++i)
	{
		output.write(reinterpret_cast<const char *>(&s.actual[i]),sizeof(s.actual[i]));
	}
	for (unsigned int i=0;i<initialunstable.size();++i)
	{
		output.write(reinterpret_cast<const char *>(&initialunstable[i].first),sizeof(initialunstable[i].first));
		output.write(reinterpret_cast<const char *>(&initialunstable[i].second),sizeof(initialunstable[i].second));
	}
}

void receiveallouters()
{
    for (int i=1;i<partsx*partsy;++i)   // This is the only part where the sizes being equal for all subgrids matter
    {
        if (allneighborsbottom[i]!=-1)
        {
            MPI_Recv(&(allouterbottom[i].front()),allouterbottom[i].size(),MPI_INT,i,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        if (allneighborstop[i]!=-1)
        {
            MPI_Recv(&(alloutertop[i].front()),alloutertop[i].size(),MPI_INT,i,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        if (allneighborsleft[i]!=-1)
        {
            MPI_Recv(&(allouterleft[i].front()),allouterleft[i].size(),MPI_INT,i,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        if (allneighborsright[i]!=-1)
        {
            MPI_Recv(&(allouterright[i].front()),allouterright[i].size(),MPI_INT,i,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
    }
}

void sendallouterstoadd()
{
    for (int i=1;i<partsx*partsy;++i)   // This is the only part where the sizes being equal for all subgrids matter
    {
        if (allneighborsbottom[i]!=-1)
        {
            MPI_Send(&(alloutertop[allneighborsbottom[i]].front()),alloutertop[allneighborsbottom[i]].size(),MPI_INT,i,0,MPI_COMM_WORLD);
        }
        if (allneighborstop[i]!=-1)
        {
            MPI_Send(&(allouterbottom[allneighborstop[i]].front()),allouterbottom[allneighborstop[i]].size(),MPI_INT,i,0,MPI_COMM_WORLD);
        }
        if (allneighborsleft[i]!=-1)
        {
            MPI_Send(&(allouterright[allneighborsleft[i]].front()),allouterright[allneighborsleft[i]].size(),MPI_INT,i,0,MPI_COMM_WORLD);
        }
        if (allneighborsright[i]!=-1)
        {
            MPI_Send(&(allouterleft[allneighborsright[i]].front()),allouterleft[allneighborsright[i]].size(),MPI_INT,i,0,MPI_COMM_WORLD);
        }
    }
}


void receiveouterstoaddfrommaster(subgrid &s)
{
    vector<int> tempouterbottom(s.outerbottom.size(),0);
    vector<int> tempoutertop(s.outertop.size(),0);
    vector<int> tempouterleft(s.outerleft.size(),0);
    vector<int> tempouterright(s.outerright.size(),0);
    if (s.neighborbottom!=-1)
    {
        MPI_Recv(&(tempouterbottom.front()),tempouterbottom.size(),MPI_INT,MASTERPROCESS,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        for(unsigned int i=0; i< tempouterbottom.size();++i)
        {
            s(i,stepy-1)+=tempouterbottom[i];
        }
    }
    if (s.neighbortop!=-1)
    {
        MPI_Recv(&(tempoutertop.front()),tempoutertop.size(),MPI_INT,MASTERPROCESS,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        for(unsigned int i=0; i< tempoutertop.size();++i)
        {
            s(i,0)+=tempoutertop[i];
        }
    }
    if (s.neighborleft!=-1)
    {
        MPI_Recv(&(tempouterleft.front()),tempouterleft.size(),MPI_INT,MASTERPROCESS,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        for(unsigned int i=0; i< tempouterleft.size();++i)
        {
            s(0,i)+=tempouterleft[i];
        }
    }
    if (s.neighborright!=-1)
    {
        MPI_Recv(&(tempouterright.front()),tempouterright.size(),MPI_INT,MASTERPROCESS,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        for(unsigned int i=0; i< tempouterright.size();++i)
        {
            s(stepx-1,i)+=tempouterright[i];
        }
    }
}

void sendouterstomaster(subgrid& s) // The master process should receive these in the same order!
{
    if (s.neighborbottom!=-1)
    {
        MPI_Send(&(s.outerbottom.front()),s.outerbottom.size(),MPI_INT,MASTERPROCESS,0,MPI_COMM_WORLD);
    }
    if (s.neighbortop!=-1)
    {
        MPI_Send(&(s.outertop.front()),s.outertop.size(),MPI_INT,MASTERPROCESS,0,MPI_COMM_WORLD);
    }
    if (s.neighborleft!=-1)
    {
        MPI_Send(&(s.outerleft.front()),s.outerleft.size(),MPI_INT,MASTERPROCESS,0,MPI_COMM_WORLD);
    }
    if (s.neighborright!=-1)
    {
        MPI_Send(&(s.outerright.front()),s.outerright.size(),MPI_INT,MASTERPROCESS,0,MPI_COMM_WORLD);
    }
}


void addoutersinmaster(subgrid& s)
{
    if (s.neighborbottom!=-1)
    {
        vector<int> tempouterbottom=alloutertop[s.neighborbottom];
        for(unsigned int i=0; i< tempouterbottom.size();++i)
        {
            s(i,stepy-1)+=tempouterbottom[i];
        }
    }
    if (s.neighborright!=-1)
    {
        vector<int> tempouterright=allouterleft[s.neighborright];
        for(unsigned int i=0; i< tempouterright.size();++i)
        {
            s(stepx-1,i)+=tempouterright[i];
        }
    }

}

void waitformaster()
{
    int a;
    cout<<world_rank<<": "<<"Waiting for master"<<endl;
    MPI_Recv(&a,1,MPI_INT,MASTERPROCESS,-1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    cout<<world_rank<<": "<<"Master is done"<<endl;
}

void masterdone()
{
    int a=1;
    for(int i=1;i<partsx*partsy;i++)
    {
        MPI_Send(&a,1,MPI_INT,i,-1,MPI_COMM_WORLD);
    }

}

void debug_messages(int messagenumber, int showmessages)
{
    if (showmessages==1)
    {
        switch(messagenumber)
        {
            case 1: cout<<"MASTER: Done with initializations"<<endl;break;
            case 2: cout<<world_rank<<": Done receiving initial unstables."<<endl; break;
            case 3: cout<<world_rank<<": Done with initializations."<<endl; break;
            case 4: cout<<world_rank<<": Checking criticals."<<endl;break;
            case 5: cout<<world_rank<<": Relaxing..."<<endl;break;
            case 6: cout<<world_rank<<": Relaxed."<<endl; break;
            case 7: cout<<"MASTER: Receiving pending counts."<<endl; break;
            case 8: cout<<"MASTER: Done receiving pending counts."<<endl; break;
            case 9: cout<<"MASTER: Done adding outers."<<endl; break;
            case 10: cout<<"MASTER: Done receiving all outers."<<endl; break;
            case 11: cout<<"MASTER: Done sending all outers to add to slaves."<<endl; break;
            case 12: cout<<world_rank<<": Receiving outers to add."<<endl;break;
            case 13: cout<<world_rank<<": Done receiving outers to add."<<endl;break;
            case 14: cout<<world_rank<<": Sending outers to master."<<endl;break;
            case 15: cout<<world_rank<<": Done sending outers to master."<<endl;break;
        }
    }
}

void debug_messages(int messagenumber, int showmessages, int process)
{
    if (showmessages==1)
    {
        switch(messagenumber)
        {
            case 0: cout<<"MASTER: Started sending initial unstables to "<<process<<"."<<endl; break;
        }
    }
}

//============================================================================
// TODO: Implement parameter parsing
//
//============================================================================
int main(int argc, char **argv) {
    int pendingcount;
    int numberinitialunstable;
    int accum = 1;
    int debugging=0; // Verbosity is OFF by default
    int donemessage=0;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    sanitycheck(argc,argv);
    subgrid s(        world_rank,
                                (world_rank%partsx) * stepx,
                                int(world_rank/partsx) * stepy,
                                stepx,
                                stepy,
                                CRITICALMINUSONE);
    if(world_rank == 0)
    {
        init(s);
        for (int i=1; i<partsx*partsy; ++i)
        {
            vector <int> tempunstablex;
            vector <int> tempunstabley;
            for (unsigned int j=0; j < initialunstable.size(); ++j)
            {
                if (i==initialunstablesubgrids[j])
                {
                    tempunstablex.push_back(initialunstable[j].first);
                    tempunstabley.push_back(initialunstable[j].second);
                }
            }
            numberinitialunstable = tempunstablex.size();
            debug_messages(0,debugging,i);
            MPI_Send(&numberinitialunstable,1,MPI_INT,i,1,MPI_COMM_WORLD);
            MPI_Send(&(tempunstablex.front()),numberinitialunstable,MPI_INT,i,2,MPI_COMM_WORLD);
            MPI_Send(&(tempunstabley.front()),numberinitialunstable,MPI_INT,i,3,MPI_COMM_WORLD);
        }
        debug_messages(1,debugging);

    }
    else
    {
        MPI_Recv(&numberinitialunstable,1,MPI_INT,MASTERPROCESS,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        vector <int> tempunstablex(numberinitialunstable,0);
        vector <int> tempunstabley(numberinitialunstable,0);
        // We send &(vector.front()) instead of vector because we are using vectors, not arrays!
        MPI_Recv(&(tempunstablex.front()),numberinitialunstable,MPI_INT,MASTERPROCESS,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Recv(&(tempunstabley.front()),numberinitialunstable,MPI_INT,MASTERPROCESS,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        debug_messages(2,debugging);
        for (int i=0; i<numberinitialunstable; i++)
        {
            s(tempunstablex[i]-s.getlocationx(),tempunstabley[i]-s.getlocationy())=CRITICAL;
        }
        debug_messages(3,debugging);
    }
    if(world_rank==0)
    {
        MPI_Bcast(&donemessage,1,MPI_INT,MASTERPROCESS,MPI_COMM_WORLD);
    }
    else
    {
        MPI_Bcast(&donemessage,1,MPI_INT,MASTERPROCESS,MPI_COMM_WORLD);
    }
    while(accum!=0)
    {
        accum=0;
        debug_messages(4,debugging);
        checkcriticals(s);
        debug_messages(5,debugging);
        relax(s);
        debug_messages(6,debugging);
        if (world_rank==0)
        {
            int mypendingcount=nonzerocount(s.outerbottom)
                            +nonzerocount(s.outertop)
                            +nonzerocount(s.outerleft)
                            +nonzerocount(s.outerright);
            accum += mypendingcount;
            allouterbottom[MASTERPROCESS]=s.outerbottom;
            alloutertop[MASTERPROCESS]=s.outertop;
            allouterleft[MASTERPROCESS]=s.outerleft;
            allouterright[MASTERPROCESS]=s.outerright;
            debug_messages(7,debugging);
            for (int i=1;i<partsx*partsy;++i)
            {
                MPI_Recv(&pendingcount,1,MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                accum += pendingcount;
            }
            debug_messages(8,debugging);
            addoutersinmaster(s);
            debug_messages(9,debugging);
            receiveallouters();
            debug_messages(10,debugging);
            sendallouterstoadd();
            debug_messages(11,debugging);
            for (int i=1;i<partsx*partsy;++i)
            {
                MPI_Send(&accum,1,MPI_INT,i,0,MPI_COMM_WORLD);
            }
        }
        else
        {
            pendingcount =   nonzerocount(s.outerbottom)
                            +nonzerocount(s.outertop)
                            +nonzerocount(s.outerleft)
                            +nonzerocount(s.outerright) ;
            MPI_Send(&pendingcount,1,MPI_INT,MASTERPROCESS,0,MPI_COMM_WORLD);
            debug_messages(14,debugging);
            sendouterstomaster(s);
            debug_messages(15,debugging);
            debug_messages(12,debugging);
            receiveouterstoaddfrommaster(s);
            debug_messages(13,debugging);
            MPI_Recv(&accum,1,MPI_INT,MASTERPROCESS,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
    }
    if (world_rank==0)
    {
        subgrid total(  MASTERPROCESS,
                        0,
                        0,
                        m,
                        n,
                        0);
        vector<int> tempactual(stepx*stepy,0);
        int templx,temply;
        addsubgridtototal(total,s.actual,s.getlocationx(),s.getlocationy());
        for (int i=1;i<partsx*partsy;++i)
        {
            MPI_Recv(&(tempactual.front()),tempactual.size(),MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            MPI_Recv(&templx,1,MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            MPI_Recv(&temply,1,MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            addsubgridtototal(total,tempactual,templx,temply);
        }
        writeout(total);
    }
    else
    {
        int templx=s.getlocationx();
        int temply=s.getlocationy();
        MPI_Send(&(s.actual.front()),s.actual.size(),MPI_INT,MASTERPROCESS,0,MPI_COMM_WORLD);
        MPI_Send(&templx,1,MPI_INT,MASTERPROCESS,0,MPI_COMM_WORLD);
        MPI_Send(&temply,1,MPI_INT,MASTERPROCESS,0,MPI_COMM_WORLD);
    }
    MPI_Finalize();
    return 0;
}
