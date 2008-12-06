#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE origidx

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <configuration.hpp>
#include <jumbled_adjmat_gen.hpp>
#include <boost/program_options.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <DegreeSort.hpp>
#include <matlab_io.hpp>
using namespace boost::program_options;
namespace ublas = boost::numeric::ublas;
using boost::any_cast;

using namespace std;
#define N 8

struct Fixture{
	AdjMat::AdjMatT adjmat; 
	ProbAdjPerm pap;
	
	Fixture()
		:adjmat(N,N)
	{
		for(int i=0;i<N;i++)
			for(int j=i;j<N;j++){
				switch(i){
					case 0:  adjmat(i,j) = adjmat(j,i) = 3; break;
					case 1:  adjmat(i,j) = adjmat(j,i) = 2; break;
					case 3:  adjmat(i,j) = adjmat(j,i) = 1; break;
					default: adjmat(i,j) = adjmat(j,i) = 0;
				}
			}
		boost::shared_ptr<AdjMat::AdjMatT> ptr(new AdjMat::AdjMatT(adjmat));
		pap.setAdjMat(ptr);
	}
	~Fixture(){
	}
};

BOOST_FIXTURE_TEST_SUITE( suite, Fixture )

BOOST_AUTO_TEST_CASE( testOrigIdx )
{
	matlab_matrix_out(cout,"sortedmat",*pap.getAdjMat());
	JumbledAdjMatGen jumb(pap);
	for(int iter=0;iter<10;iter++){
		ProbAdjPerm mypap = jumb();
		DegreeSort ds;
		ds.sort(mypap);
		matlab_matrix_out(cout,"sortedmat",*mypap.getAdjMat());
		BOOST_CHECK_EQUAL(mypap.getOriginalIndex(N-1), 0);
		BOOST_CHECK_EQUAL(mypap.getOriginalIndex(N-2), 1);
		BOOST_CHECK_EQUAL(mypap.getOriginalIndex(N-3), 3);
	}
}

BOOST_AUTO_TEST_SUITE_END()

