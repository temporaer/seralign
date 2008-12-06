#include <signal.h>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <factory/factory.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/foreach.hpp>
#include <configuration.hpp>
#include "random_pat_adjmat_gen.hpp"
#include <nana.h>
#undef A

using namespace std;
namespace ublas = boost::numeric::ublas;


#define SQR(X) ((X)*(X))

RandomPatAdjmatGen::RandomPatAdjmatGen()
	:mSize(20),
	 mPatSize(6),
	 mEdgeInclusionProb(0.5),
	 mRunningID(0)
{
}


void RandomPatAdjmatGen::configure()
{
	AdjMatGen::configure(); 
	this->setSize(gCfg().getInt("rand_pat_adjmat_gen.size"));
	this->setPatSize(gCfg().getInt("rand_pat_adjmat_gen.pattern_size"));
	this->setEdgeInclusionProb(gCfg().getFloat("rand_pat_adjmat_gen.edge_inclusion_prob"));
	this->setGeneratePattern(gCfg().getBool("rand_pat_adjmat_gen.gen_pat"));
}

void RandomPatAdjmatGen::fillNodeArrayCirclePattern()
{
	int n = mSize;
	int m = mPatSize;
	if(n<m)
		throw logic_error("Pattern-Size is smaller than matrix size.");
	mNodes.clear();

	const float radius = 15.f;
	const float maxX = 100.f;
	const float maxY = 100.f;
	RNode center;
	center.x = radius + (maxX-2.0f*radius) * drand48();
	center.y = radius + (maxY-2.0f*radius) * drand48();

	// Put Circle
	if(mGeneratePattern){
		for(int i=0;i<m;i++){
			RNode node;
			float ang = i*2.0*3.14/m;
			node.x = center.x + radius * sin(ang);
			node.y = center.y + radius * cos(ang);
			node.part_of_pattern = true;
			mNodes.push_back(node);
		}
	}
	else{
		// assign part_of_pattern property to random nodes
		for(int i=0;i<m;i++){
			RNode node;
			node.x = maxX * drand48();
			node.y = maxY * drand48();
			node.part_of_pattern = true;
			mNodes.push_back(node);
		}
	}
	// Put arbitrary nodes.
	for(int i=m;i<n;){
		RNode node;
		node.x = maxX * drand48();
		node.y = maxY * drand48();
		node.part_of_pattern = false;
		if(SQR(node.x-center.x)+SQR(node.y-center.y) < SQR(radius))
			continue;
		i++;
		mNodes.push_back(node);
	}
}
void RandomPatAdjmatGen::triangulateNodes(std::vector<Edge>& edges)
{
	const char* outfile = "/tmp/.rpag_nodes";
	const char*  infile = "/tmp/.rpag_triangles";
	ofstream os(outfile);
	for(unsigned int i=0;i<mNodes.size();i++){
		os<<i<<"\t"<<mNodes[i].x<<"\t"<<mNodes[i].y<<endl;
	}
	os.close();

	// execute triangulator
	char cmd[255];
	sprintf(cmd,"%s %s > %s 2> /dev/null",TRIANGULATION_COMMAND,outfile,infile);
	LG(isVerbose(),"Exec Cmd: %s...",cmd);
	int res = system(cmd);
	LG(isVerbose(),"done.\n");
	if(res == -1)
		throw runtime_error(std::string("RandomPatAdjmatGen could not execute triangulator."));
	if(WIFSIGNALED(res) &&
			(WTERMSIG(res) == SIGINT || WTERMSIG(res) == SIGQUIT)){
		cerr << "Got interrupt, calling exit." << endl;
		exit(0);
	}

	// read triangulator output
	ifstream is(infile);
	while(!is.eof()){
		int a=-1,b=-1;
		is >> a >> b;
		if(a<0 || b<0)
			break;
		Edge e;
		e.first = a; e.second = b;
		edges.push_back(e);
	}
}


ProbAdjPerm RandomPatAdjmatGen::operator()()
{
	// shuffle nodes
	fillNodeArrayCirclePattern();
	random_shuffle(mNodes.begin(),mNodes.end());

	// get edges
	vector<Edge> edges;
	triangulateNodes(edges);

	boost::shared_ptr<AdjMat::AdjMatT> A_ptr;
	A_ptr.reset(new ublas::matrix<double>(mSize,mSize));
	AdjMat::AdjMatT& A = *A_ptr;

	A = ublas::zero_matrix<double>(mSize, mSize);

	BOOST_FOREACH(Edge& e, edges){
		bool incl = false;
		bool bothInPattern = mNodes[e.first].part_of_pattern && mNodes[e.second].part_of_pattern;
		incl |= bothInPattern;
		incl |= drand48()<mEdgeInclusionProb; 
		if(!incl)
			continue;
		double d = 1;
		d *= sqrt( SQR(mNodes[e.first].x-mNodes[e.second].x)
			+      SQR(mNodes[e.first].y-mNodes[e.second].y));
		//d *=      (mNodes[e.first ].part_of_pattern || mNodes[e.second].part_of_pattern)? 2 : 1;
		A(e.first,e.second) = d;
		A(e.second,e.first) = d;
	}

	ProbAdjPerm pap;
	pap.setAdjMat(A_ptr);

	mRunningID++;
	return pap;
}

RandomPatAdjmatGen::~RandomPatAdjmatGen()
{
}
bool RandomPatAdjmatGen::hasNext()
{
	return true;
}

string RandomPatAdjmatGen::getPlainDescription(int ser_idx, const Serialization&s)
{
	int idx = s[ser_idx];
	return mNodes[idx].part_of_pattern ? string("X") : string("_");
}
string RandomPatAdjmatGen::getPrologDescription(int ser_idx, const Serialization&s)
{
	int idx = s[ser_idx];
	stringstream str;
	str	<<"xpos("<<(int)(mNodes[idx].x)<<");"
		<<"ypos("<<(int)(mNodes[idx].y)<<")";
	return str.str();
}

string RandomPatAdjmatGen::getGraphID()
{
	stringstream str;
	str<<"randpatgraph_"<<mSize<<"_"<<mPatSize<<"_"<<mRunningID;
	return str.str();
}
std::string RandomPatAdjmatGen::getGraphVizNodeAttribs(int idx)
{
	stringstream str;

	str<<",";
	if(mNodes[idx].part_of_pattern) str<<"color=red";
	else                            str<<"color=black";

	str<<",pos=\"" << 10*(int)mNodes[idx].x <<"," << 10*(int)mNodes[idx].y <<  "\"";
	return str.str();
}


namespace{
	registerInFactory<AdjMatGen, RandomPatAdjmatGen> registerBase("RandomPatAdjmatGen");
}


