#define BOOST_TEST_DYN_LINK 
#define BOOST_TEST_MODULE graph_from_adj

#include <numeric>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/shared_ptr.hpp>
#include <DegreeSort.hpp>
#include <GraphFromAdj.hpp>
using namespace std;
using namespace boost;
namespace ublas = boost::numeric::ublas;

struct Fixture{
	const int n;
	AdjMat::AdjMatT adj;
	ProbAdjPerm pap;
	
	Fixture()
		:  n(4)
		 , adj(n,n)
	{
		adj(0,0) = 0;
		adj(1,1) = 0;
		adj(2,2) = 0;
		adj(3,3) = 0;
		adj(1,0) = adj(0,1) = 10;
		adj(2,0) = adj(0,2) =  6;
		adj(2,1) = adj(1,2) =  5;
		adj(3,0) = adj(0,3) =  5;
		adj(3,1) = adj(1,3) =  6;
		adj(3,2) = adj(2,3) =  7.21110255092798;
		pap.setAdjMat(boost::shared_ptr<AdjMat::AdjMatT>(new AdjMat::AdjMatT(adj)));
	}
	~Fixture(){
	}
};

BOOST_FIXTURE_TEST_SUITE( suite, Fixture )



BOOST_AUTO_TEST_CASE( testFastmapping )
{
	int from, to;
	from = 0;
	to   = 1;
	GraphFromAdj g(pap,from,to);
	int i = g.getFarthestFromA(0); // level 0: find node far away from `from'=0
	BOOST_CHECK_EQUAL(i,1);        // level 0: should be node 1
	g.setA(1,2);                 // level 1: set root A=2
	i = g.getFarthestFromA(1);   // level 1: find node far away from A=2
	BOOST_CHECK_EQUAL(i,3);      // level 1: should be node 3
}
BOOST_AUTO_TEST_CASE( testProjection )
{
	int from,to;
	if(1){
		from = 0;
		to   = 1;
		GraphFromAdj g(pap,from,to);
		cout << "========= projecting on line " << from << " to " << to<<endl;
		BOOST_CHECK_LT(fabs(g.getProjection(0,0)), 1E-9);
		BOOST_CHECK_CLOSE(g.getProjection(0,1),10.0,0.1);
		BOOST_CHECK_CLOSE(g.getProjection(0,2), 5.55,0.1);
		BOOST_CHECK_CLOSE(g.getProjection(0,3), 4.45,0.1);

		cout << "========= projecting on line 2 to 3"<<endl; 
		g.setA(1,2);
		g.setB(1,3);

		BOOST_CHECK_LT(fabs(g.getProjection(1,2)), 1E-9);
		BOOST_CHECK_CLOSE(g.getProjection(1,0),3.563,0.1);
		BOOST_CHECK_CLOSE(g.getProjection(1,1),3.563,0.1);
		BOOST_CHECK_CLOSE(g.getProjection(1,3),7.127,0.1);

		//for(int i=0;i<n;i++){
			//cout << i << "  " << g.getProjection(1,i)<<endl;
		//}
	}
	if(1){
		from = 3;
		to   = 2;
		GraphFromAdj g(pap,from,to);
		cout << "======= projecting on line " << from << " to " << to<<endl;
		BOOST_CHECK_CLOSE(g.getProjection(0,0), 2.843, 0.1);
		BOOST_CHECK_CLOSE(g.getProjection(0,1), 4.368, 0.1);
		BOOST_CHECK_CLOSE(g.getProjection(0,2), 7.2111, 0.1);
		BOOST_CHECK_LT(fabs(g.getProjection(0,3)), 1E-9);
	}
}



BOOST_AUTO_TEST_SUITE_END()

