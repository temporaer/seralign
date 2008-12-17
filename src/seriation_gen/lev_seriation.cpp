// 2008-11-17 Hannes Schulz <mail at hannes-schulz dot de>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/bindings/lapack/lapack.hpp>
#include <factory/factory.h>
#include <matlab_io.hpp>
#include "lev_seriation.hpp"
#include <nana.h>
#undef C

namespace lapack = boost::numeric::bindings::lapack;
namespace ublas = boost::numeric::ublas;
using namespace std;

/**********************************************************
 *          LEVSeriationGen Implementation                *
 **********************************************************/
struct LEVSeriationGen::Impl{
	Impl();
	~Impl();
	Serialization operator()(const ProbAdjPerm& pap);
	Serialization readout_plain(ublas::vector<double>& x,const AdjMat::AdjMatT& adj);
};

// Impl constructor
LEVSeriationGen::Impl::Impl(){
}

// Impl destructor
LEVSeriationGen::Impl::~Impl(){
}

Serialization LEVSeriationGen::Impl::operator()(const ProbAdjPerm& pap)
{
	ProbAdjLapPerm pap2(pap);
	LG(isVerbose(),"Generating Laplacian\n");
	pap2.calculateLaplacian();
	Laplacian::LaplacianT& L = *pap2.getLaplacian(); 
	int n = L.size1();
	Serialization ret(n);

	LG(isVerbose(),"Determine eigenvalues of laplacian\n");
	ublas::matrix<double,ublas::column_major> Eigv(L);
	ublas::vector<double> lambda(n);
	ublas::vector<double>::iterator best_lambda;
	int best_lambda_idx=0;

	lapack::syev( 'V', 'L', Eigv, lambda, lapack::minimal_workspace() );
	//best_lambda     = min_element(lambda.begin(),lambda.end());   // lambdas are sorted!
	//best_lambda_idx = std::distance(lambda.begin(),best_lambda);  // lambdas are sorted!
	best_lambda_idx   = 1;                                          // 1==lambda_2==fiedler vector
	ublas::vector<double> leading = ublas::column(Eigv,best_lambda_idx);

	LG(isVerbose(),"eigenvalue-readout\n");
	//L("Best Lambda: %3.3f\n",lambda[best_lambda_idx]);
	//matlab_vector_out(cout,"Lambdas",lambda);
	//matlab_vector_out(cout,"Leading",leading);
	return readout_plain(leading,*pap.getAdjMat());
}

#   define BEST_ELEM(X) max_element(X.begin(),X.end())
Serialization LEVSeriationGen::Impl::readout_plain(ublas::vector<double>& x,const AdjMat::AdjMatT& adj)
{
	unsigned int n = x.size();

	// tricky: make sure x > 0 at all times.
	x += ublas::scalar_vector<double>(n, 1 - (*min_element(x.begin(),x.end())));

	Serialization ret(n);
	std::vector<bool> done(n,false);

	// find highest component of x
	ublas::vector<double>::iterator it = BEST_ELEM(x);
	int idx = std::distance(x.begin(),it);

	LG(isVerbose(),"Determine Actual Path through Graph.\n");
	for(unsigned int i=0;i<n;i++){
		ret[i] = idx;
		done[idx] = true;
		*it = 0.0; 
		it = BEST_ELEM(x);
		idx = std::distance(x.begin(),it);
	}
	return ret;
}

/**********************************************************
 *          LEVSeriationGen Interface                     *
 **********************************************************/
LEVSeriationGen::LEVSeriationGen()
	:mImpl(new Impl())
{
}

LEVSeriationGen::~LEVSeriationGen()
{
  // cleanup
}
Serialization LEVSeriationGen::operator()(const ProbAdjPerm& pap)
{
	return (*mImpl)(pap);
}


namespace{ registerInFactory<SerGenAdj, LEVSeriationGen> registerBase("LEVSeriationGen"); }
