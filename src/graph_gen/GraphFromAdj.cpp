// 2008-12-19 Hannes Schulz <mail at hannes-schulz dot de>

#include "GraphFromAdj.hpp"
#include <limits.h>
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
	typedef std::vector<double> Map;
	typedef bgl::adjacency_matrix<bgl::undirectedS,bgl::no_property,bgl::property<bgl::edge_weight_t, cost> > Graph;
	typedef bgl::property_map<Graph, bgl::edge_weight_t>::type WeightMap;
	typedef Graph::vertex_descriptor vertex_descriptor;
	typedef Graph::edge_descriptor edge_descriptor;
	typedef Graph::vertex_iterator vertex_iterator;

	struct DistAndHeight{
		vertex_descriptor A;
		vertex_descriptor B;
		Map distA;
		Map distB;
		Map height;
		double totalDist;
		DistAndHeight(){}
		DistAndHeight(int n) :A(0),B(n),distA(n),distB(n),height(n),totalDist(1)
		{
		}
	};

	Impl(const ProbAdjPerm&, unsigned int idx1, unsigned int idx2, int dim=3);
	Impl(const ProbAdjPerm&, int dim=3);
	~Impl();
	inline double  getDist(unsigned int k, const int& i, const int& j){
		vertex_descriptor f = bgl::vertex(i, mG);
		vertex_descriptor t = bgl::vertex(j, mG);
		return getDist(k,f,t);
	}
	double  getDist(unsigned int k, const vertex_descriptor& i, const vertex_descriptor& j);
	inline double  getProjection(unsigned int k, const int& i){
		vertex_descriptor f = bgl::vertex(i, mG);
		return getProjection(k,f);
	}
	double  getProjection(unsigned int k, const vertex_descriptor& i);
	double  heightFoot(double a, double b, double g);
	void    setA(unsigned int k, int i);
	void    setB(unsigned int k, int i);
	unsigned int     getFarthestFrom(unsigned int k, const Map&);
	unsigned int     getFarthestFromA(unsigned int k);
	unsigned int     getFarthestFromB(unsigned int k);
	Graph            mG;
	WeightMap        mWeightMap;
	std::vector<DistAndHeight> mDistAndHeights;
};

double GraphFromAdj::Impl::getProjection(unsigned int k, const vertex_descriptor& i){
	I(mDistAndHeights.size()>k); 
	DistAndHeight& dh = mDistAndHeights[k];
	double dai = getDist(k,dh.A,i);
	double dbi = getDist(k,dh.B,i);
	double dab = getDist(k,dh.A,dh.B);
	double xi  = heightFoot(dai, dbi, dab);
	return xi;
}
double GraphFromAdj::Impl::getDist(unsigned int k, const vertex_descriptor& i, const vertex_descriptor& j){
	I(mDistAndHeights.size()>k); 
	if(i==j) return 0.0;
	if(k==0){
		DistAndHeight& dh = mDistAndHeights[k];
		if(dh.A==i) return dh.distA[j];
		if(dh.B==i) return dh.distB[j]; 
		I(false); // user makes sure we never get here *cough*
	}
	DistAndHeight& dh = mDistAndHeights[k];
	double td;
	if(0);
	else if(dh.A==i) td = dh.distA[j]; 
	else if(dh.B==i) td = dh.distB[j];
	else{
		I(false);
	}
	double dist2 = td*td;
	for(unsigned int dim=0;dim<k;dim++){
		const double xi  = getProjection(dim,i);
		const double xj  = getProjection(dim,j);
		const double dxij = xi -xj;
		dist2 -= dxij*dxij;
	}
#ifndef NDEBUG
	if(dist2<0) cerr << "Warning: getDist: dist2<0"<<endl;
#endif
	return sqrt(dist2);
}

double 
GraphFromAdj::Impl::heightFoot(double a, double b, double g){
#ifndef NDEBUG
	if(a!=a) cerr << "Warning: NaN found in distA!"<<endl;
	if(b!=b) cerr << "Warning: NaN found in distB!"<<endl;
	if(a+b<g-1E10)
		cerr << "Warning: Triangle Inequality does not hold: "<< (g-a-b)<<endl;
	if(g!=g) cerr << "Warning: Base is NaN!"<<endl;
	if(g==0) cerr << "Warning: Base is 0!"<<endl;
#endif
	//cout << "heightFoot: " << V(a)<<V(b)<<V(g)<<endl;
	if(b==0) return a;
	if(a==0) return g-b;
#if 0
	double cosgamma = (g*g-a*a-b*b)/(2*a*b), gamma;
	if(0);
	else if(cosgamma < -1.0){ cout<<V(cosgamma)<<endl;gamma = 0;   }
	else if(cosgamma >  1.0){ cout<<V(cosgamma)<<endl;gamma = M_PI;}
	else gamma = acos(cosgamma);
	double h = a*b*sin(gamma)/g;
	return sqrt(a*a - h*h);
#else
	return (a*a + g*g - b*b)/(2*g);
#endif
}

unsigned int
GraphFromAdj::Impl::getFarthestFrom(unsigned int k, const Map& m){
	return distance(m.begin(),max_element(m.begin(),m.end()));
}

unsigned int
GraphFromAdj::Impl::getFarthestFromA(unsigned int k){
	I(mDistAndHeights.size() > k);
	return getFarthestFrom(k, mDistAndHeights[k].distA);
}
unsigned int
GraphFromAdj::Impl::getFarthestFromB(unsigned int k){
	I(mDistAndHeights.size() > k);
	return getFarthestFrom(k, mDistAndHeights[k].distB);
}
void
GraphFromAdj::Impl::setA(unsigned int k, int i){ 
	I(mDistAndHeights.size() > k);
	vertex_descriptor s = bgl::vertex(i, mG);
	DistAndHeight& dh = mDistAndHeights[k];
	I(dh.distA.size()>0);
	dh.A = s;
	bgl::dijkstra_shortest_paths(mG, s, bgl::distance_map(&dh.distA[0]));
	dh.totalDist = dh.distA[dh.B];
}
void
GraphFromAdj::Impl::setB(unsigned int k, int i){
	I(mDistAndHeights.size() > k);
	vertex_descriptor s = bgl::vertex(i, mG);
	DistAndHeight& dh = mDistAndHeights[k];
	I(dh.distB.size()>0);
	dh.B = s;
	bgl::dijkstra_shortest_paths(mG, s, bgl::distance_map(&dh.distB[0]));
	dh.totalDist = dh.distB[dh.A];
}


// Impl constructor
GraphFromAdj::Impl::Impl(const ProbAdjPerm& pap, unsigned int idx1, unsigned int idx2, int dim)
	:mG(pap.getAdjMat()->size1())
	,mDistAndHeights(dim, DistAndHeight(pap.getAdjMat()->size1()))
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
	setA(0,idx1);
	setB(0,idx2);
}

// Impl constructor
GraphFromAdj::Impl::Impl(const ProbAdjPerm& pap, int dim)
	: mG(pap.getAdjMat()->size1())
	 ,mDistAndHeights(dim,DistAndHeight(pap.getAdjMat()->size1()))
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
}

// Impl destructor
GraphFromAdj::Impl::~Impl(){
}


/**********************************************************
 *          GraphFromAdj Interface                     *
 **********************************************************/
GraphFromAdj::GraphFromAdj(const ProbAdjPerm& pap, unsigned int idx1, unsigned int idx2, int dim)
	:mImpl(new Impl(pap,idx1,idx2,dim))
{
}
GraphFromAdj::GraphFromAdj(const ProbAdjPerm& pap, int dim)
	:mImpl(new Impl(pap,dim))
{
}

GraphFromAdj::~GraphFromAdj()
{
  // cleanup
}
double GraphFromAdj::getProjection(unsigned int k, int i)
{
	return mImpl->getProjection(k, i);
}
double GraphFromAdj::getDist(unsigned int k, int i, int j)
{
	return mImpl->getDist(k, i, j);
}
double GraphFromAdj::getTotalDist(unsigned int k){
	return mImpl->mDistAndHeights[k].totalDist;
}
void   GraphFromAdj::setA(unsigned int k, int i)
{ mImpl->setA(k, i);}
void   GraphFromAdj::setB(unsigned int k, int i)
{ mImpl->setB(k, i);}
unsigned int    GraphFromAdj::getFarthestFromA(unsigned int k)
{ return mImpl->getFarthestFromA(k); }
unsigned int    GraphFromAdj::getFarthestFromB(unsigned int k)
{ return mImpl->getFarthestFromB(k); }


