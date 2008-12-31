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
		:  n(6)
		 , adj(n,n)
	{
		for(int i=0;i<n;i++)
			for(int j=i;j<n;j++)
			{
				adj(i,j) = adj(j,i) = 0.01+0.01* (int)(100*drand48());
			}
		pap.setAdjMat(boost::shared_ptr<AdjMat::AdjMatT>(new AdjMat::AdjMatT(adj)));
	}
	~Fixture(){
	}
};

BOOST_FIXTURE_TEST_SUITE( suite, Fixture )



BOOST_AUTO_TEST_CASE( testIncreasingDegree )
{
	DegreeSort ds;
	ds.sort( pap );
	
	for(int s=0;s<n;s++){
		GraphFromAdj g(pap,s,n-s-1);
		BOOST_CHECK_CLOSE(g.getDist(0,s),0.0, 0.01);
		for(int i=0;i<n;i++){
			if(i==s)
				continue;
			BOOST_CHECK_GT(g.getDist(0,i),0.0099);
		}
	}
}



BOOST_AUTO_TEST_SUITE_END()

