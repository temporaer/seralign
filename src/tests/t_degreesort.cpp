#define BOOST_TEST_DYN_LINK 
#define BOOST_TEST_MODULE degreesort

#include <numeric>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/shared_ptr.hpp>
#include <DegreeSort.hpp>
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
				adj(i,j) = adj(j,i) = 0.01* (int)(100*drand48());
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

	AdjMat::AdjMatT&   sadj  = *pap.getAdjMat();
	PermMat::PermMatT& sperm = *pap.getPermMat();
	AdjMat::AdjMatT     tmp  = prod(sperm, sadj);

	double lasts = -1;
	for(int i=0;i<n;i++){
		ublas::vector<double> col = ublas::column(sadj, i);
		double s = accumulate(col.begin(),col.end(),0.0);
		BOOST_CHECK_LT(lasts,s);
		lasts = s;
	}
	
}

BOOST_AUTO_TEST_CASE( testNewPermCorrect )
{
	DegreeSort ds;
	ds.sort( pap );

	AdjMat::AdjMatT&   sadj  = *pap.getAdjMat();
	PermMat::PermMatT& sperm = *pap.getPermMat();
	AdjMat::AdjMatT     tmp  = prod(sperm, sadj);
	tmp = prod(tmp, trans(sperm));

	for(int i=0;i<n;i++)
		for(int j=0; j<n; j++)
			BOOST_CHECK_CLOSE(tmp(i,j),adj(i,j),0.01);
	
}

BOOST_AUTO_TEST_CASE( testOldPermCorrect )
{
	pap.setPermMat(
			boost::shared_ptr<PermMat::PermMatT>(
				new PermMat::PermMatT(ublas::identity_matrix<int>(n))));
	DegreeSort ds;
	ds.sort( pap );

	AdjMat::AdjMatT&   sadj  = *pap.getAdjMat();
	PermMat::PermMatT& sperm = *pap.getPermMat();
	AdjMat::AdjMatT     tmp  = prod(sperm, sadj);
	tmp = prod(tmp, trans(sperm));

	for(int i=0;i<n;i++)
		for(int j=0; j<n; j++)
			BOOST_CHECK_CLOSE(tmp(i,j),adj(i,j),0.01);
	
}

BOOST_AUTO_TEST_CASE( testOldPermCorrect2 )
{
	pap.setPermMat(
			boost::shared_ptr<PermMat::PermMatT>(
				new PermMat::PermMatT(ublas::identity_matrix<int>(n))));
	DegreeSort ds;
	ds.sort( pap );
	ds.sort( pap ); // second sorting, should not change anything

	AdjMat::AdjMatT&   sadj  = *pap.getAdjMat();
	PermMat::PermMatT& sperm = *pap.getPermMat();
	AdjMat::AdjMatT     tmp  = prod(sperm, sadj);
	tmp = prod(tmp, trans(sperm));

	for(int i=0;i<n;i++)
		for(int j=0; j<n; j++)
			BOOST_CHECK_CLOSE(tmp(i,j),adj(i,j),0.01);
	
}

BOOST_AUTO_TEST_SUITE_END()

