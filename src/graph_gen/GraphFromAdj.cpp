// 2008-12-19 Hannes Schulz <mail at hannes-schulz dot de>

#include "GraphFromAdj.hpp"
#include <boost/graph/graph_traits.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <ProbData.hpp>
#include <nana.h>
#undef A
#define V(X) #X << "=" << (X) << " "

namespace bgl = boost;
using namespace std;


/**********************************************************
 *          GraphFromAdj Implementation                *
 **********************************************************/
struct GraphFromAdj::Impl{
	typedef double cost;
	typedef bgl::adjacency_matrix<bgl::undirectedS,bgl::no_property,bgl::property<bgl::edge_weight_t, cost> > Graph;
	typedef bgl::property_map<Graph, bgl::edge_weight_t>::type WeightMap;
	typedef Graph::vertex_descriptor vertex_descriptor;
	typedef Graph::edge_descriptor edge_descriptor;
	typedef Graph::vertex_iterator vertex_iterator;

	Impl(const ProbAdjPerm&, unsigned int idx1, unsigned int idx2);
	~Impl();
	double  getDist(int i);
	Graph               mG;
	WeightMap           mWeightMap;
	double              mTotalDistance;
	std::vector<pair<vertex_descriptor,vertex_descriptor> > mPivots; 
	std::vector<double> mDistanceMapA;
	std::vector<double> mDistanceMapB;
};

double
GraphFromAdj::Impl::getDist(int idx){
	vertex_descriptor v = bgl::vertex(idx, mG);
	if(v == mPivots[0].first)  return 0.0;
	if(v == mPivots[0].second) return mTotalDistance;
	
	double a = mDistanceMapA[v];
	double b = mDistanceMapB[v];
#ifndef NDEBUG
	if(a+b<mTotalDistance)
		cerr << "Warning: Triangle Inequality does not hold!"<<endl;
#endif
	double tmp = (mTotalDistance*mTotalDistance - a*a - b*b) / (-2*a*b);
	double gamma;
	if(0);
	else if(tmp<-1) gamma = M_PI;
	else if(tmp> 1) gamma = 0;
	else            gamma = acos( tmp );
	double s_gamma = sin(gamma);
	tmp = b/mTotalDistance*s_gamma;
	double x = (tmp>1)?0:(a * sqrt(1-tmp*tmp));
	return x;
}
// Impl constructor
GraphFromAdj::Impl::Impl(const ProbAdjPerm& pap, unsigned int idx1, unsigned int idx2)
	:
		 mG(pap.getAdjMat()->size1())
		,mDistanceMapA(pap.getAdjMat()->size1())
		,mDistanceMapB(pap.getAdjMat()->size1())
{
	AdjMat::AdjMatT& A = *pap.getAdjMat(); 
	unsigned int n=A.size1();
	mWeightMap = get(bgl::edge_weight,mG);
	for(unsigned int i=0;i<n;i++)
		for(unsigned int j=0;j<n;j++){
			edge_descriptor e; bool inserted;
			boost::tie(e, inserted) = bgl::add_edge(i,j,mG);
			if(A(i,j)<0.0000001) mWeightMap[e] = INT_MAX;
			else                 mWeightMap[e] = A(i,j);
		}
	vertex_descriptor s = vertex(idx1, mG);
	bgl::dijkstra_shortest_paths(mG, s, bgl::distance_map(&mDistanceMapA[0]));
	vertex_descriptor t = vertex(idx2, mG);
	bgl::dijkstra_shortest_paths(mG, t, bgl::distance_map(&mDistanceMapB[0]));
	I(mDistanceMapB.size() == mDistanceMapA.size());
	mPivots.push_back(make_pair(s,t)); 
	mTotalDistance = mDistanceMapA[t];
#ifndef NDEBUG
	if(fabs(mTotalDistance - mDistanceMapB[s]) > 0.000001){
		cerr << "Warning: Graph Distance Metric assymetric!"<< endl;
	}
#endif
}

// Impl destructor
GraphFromAdj::Impl::~Impl(){
}


/**********************************************************
 *          GraphFromAdj Interface                     *
 **********************************************************/
GraphFromAdj::GraphFromAdj(const ProbAdjPerm& pap, unsigned int idx1, unsigned int idx2)
	:mImpl(new Impl(pap,idx1,idx2))
{
}

GraphFromAdj::~GraphFromAdj()
{
  // cleanup
}
double GraphFromAdj::getDist(int i)
{
	return mImpl->getDist(i);
}

