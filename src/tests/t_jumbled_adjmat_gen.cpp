#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE 

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
			for(int j=i;j<N;j++)
				adjmat(i,j) = adjmat(j,i) = drand48()>0.5 ? 1 : 0;
		boost::shared_ptr<AdjMat::AdjMatT> ptr(new AdjMat::AdjMatT(adjmat));
		pap.setAdjMat(ptr);
	}
	~Fixture(){
	}
};

BOOST_FIXTURE_TEST_SUITE( suite, Fixture )

BOOST_AUTO_TEST_CASE( testIdentical )
{
	JumbledAdjMatGen jumble(pap);
	for(int tries = 0; tries < 10; tries ++ ){
		ProbAdjPerm        jpap  =  jumble();
		AdjMat::AdjMatT&   jadj  = *jpap.getAdjMat();
		PermMat::PermMatT& jperm = *jpap.getPermMat();
		AdjMat::AdjMatT    tmp   = prod(jperm, jadj);
		tmp = prod(tmp, trans(jperm));
		for(int i=0; i<N;i++)
			for(int j=0; j<N; j++)
				BOOST_CHECK_CLOSE( tmp(i,j), adjmat(i,j), 0.01 ) ;
	}
}

BOOST_AUTO_TEST_CASE( testIdenticalOldPerm )
{
	DegreeSort ds;
	ds.sort(pap);
	JumbledAdjMatGen jumble(pap);
	for(int tries = 0; tries < 10; tries ++ ){
		ProbAdjPerm        jpap  =  jumble();
		AdjMat::AdjMatT&   jadj  = *jpap.getAdjMat();
		PermMat::PermMatT& jperm = *jpap.getPermMat();
		AdjMat::AdjMatT    tmp   = prod(jperm, jadj);
		tmp = prod(tmp, trans(jperm));
		for(int i=0; i<N;i++)
			for(int j=0; j<N; j++)
				BOOST_CHECK_CLOSE( tmp(i,j), adjmat(i,j), 0.01 ) ;
	}
}

BOOST_AUTO_TEST_CASE( testIdenticalOldPerm2 )
{
	DegreeSort ds;
	ds.sort(pap);
	JumbledAdjMatGen jumble(pap);
	for(int tries = 0; tries < 10; tries ++ ){
		ProbAdjPerm        jpap  =  jumble();
		AdjMat::AdjMatT&   jadj  = *jpap.getAdjMat();
		PermMat::PermMatT& jperm = *jpap.getPermMat();
		AdjMat::AdjMatT    tmp   = prod(jperm, jadj);
		tmp = prod(tmp, trans(jperm));
		for(int i=0; i<N;i++)
			for(int j=0; j<N; j++)
				BOOST_CHECK_CLOSE( tmp(i,j), adjmat(i,j), 0.01 ) ;
	}
}

BOOST_AUTO_TEST_SUITE_END()

