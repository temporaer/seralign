#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <sdp_prob.hpp>
#include <sdp_seriation_gen.hpp>
#include <sdp_seriation_prob_gen.hpp>
#include <sdp_wrapper.hpp>
#include <factory/factory.h>
#include <stats.hpp>
#include <configuration.hpp>
#include <boost/numeric/bindings/lapack/lapack.hpp>
#include <cholesky.hpp> // from third_party
#include <arg_max.hpp>
#include <boost/lambda/lambda.hpp>
#include <matlab_io.hpp>
#include <cmath>
#include <nana.h>
#undef C
#undef A

#undef L_DEFAULT_GUARD
#define L_DEFAULT_GUARD isVerbose()

namespace l = boost::lambda;
namespace lapack = boost::numeric::bindings::lapack;
using namespace boost::numeric;
using namespace std;

double my_sqrt(double d){return sqrt(d);}

// Implementation of SeriationGen

struct SDPSeriationGen::Impl{
    std::auto_ptr<SDPWrapper> mSDPWrapper;
    typedef AdjMat::AdjMatT AdjMatT;
    Serialization operator()(const ProbAdjPerm& pap);
    void setSDPWrapper(std::auto_ptr<SDPWrapper> w);
    ~Impl();

	Serialization readout_connected(boost::numeric::ublas::vector<double>& x,const AdjMat::AdjMatT& adj);
	Serialization readout_plain(boost::numeric::ublas::vector<double>& x,const AdjMat::AdjMatT& adj);
};
SDPSeriationGen::Impl::~Impl(){
}

void SDPSeriationGen::configure()
{
	SerGenAdj::configure();
	string sdp_wrapper_name          = gCfg().getString("ser-gen.sdp-wrapper");
	auto_ptr<SDPWrapper> sdp_wrapper = genericFactory<SDPWrapper>::instance().create(sdp_wrapper_name);
	if(!sdp_wrapper.get())
		throw logic_error(string("Supplied SDPWrapper `") + sdp_wrapper_name + "' does not exist");
	sdp_wrapper->configure();
	setSDPWrapper(sdp_wrapper);
}


Serialization SDPSeriationGen::Impl::operator()(const ProbAdjPerm& pap)
{
	L("Generating Seriation using SDP.\n");
	const AdjMatT& adj = *pap.getAdjMat();
	unsigned int n = adj.size1();
	I(adj.size1() == adj.size2());
	I(mSDPWrapper.get() != NULL);

	L("Generate the SDP-Problem\n");
	SDPProb prob;
	SDPSeriationProbGen sdp_probgen(pap.getAdjMat());
	sdp_probgen(prob);

	L("Solve the SDP-Problem using an SDP Wrapper\n");
	SDPWrapper::AnswerT _X(n,n);
	noalias(_X) = (*mSDPWrapper)(prob);
	ublas::symmetric_adaptor<SDPWrapper::AnswerT, ublas::lower> X(_X);
	//SDPWrapper::AnswerT& X=_X;

#ifndef NDEBUG
	I(ublas::is_symmetric(X));
	//L("Make sure we got the right thing: tr(EX)=1\n");
	using ublas::range;
	using ublas::prod;
	for(unsigned int i=0;i<prob.F.size();i++){
		ublas::matrix<double> R = prod(prob.F[i],X);
		ublas::matrix_vector_range<SDPWrapper::AnswerT> diag(R, range (0,n), range (0,n));
		double trace = ublas::sum(diag);
		//L("Trace should be: %2.3f, Trace is: %2.3f\n",prob.b[i],trace);
		I(fabs(trace-prob.b[i])<0.001);
	}
#endif

// TODO: Determine Method for rounding
	int Y_METHOD = gCfg().getInt("ser-gen.sdp-rounding-method");
	SDPWrapper::AnswerT V (X);
	if (Y_METHOD == 1 || Y_METHOD == 3){
		// cholesky-decomposition
		L("Decompose X = V*V'\n");
		if(!ulapack::chol_checked_inplace(V, true))
			throw runtime_error("Could not cholesky-decompose matrix, seems not to be positive-semidefinite.");
	}
	ublas::matrix<double,ublas::column_major> Eigv(X);
	ublas::vector<double> lambda(n);
	ublas::vector<double>::iterator max_lambda;
	int max_lambda_idx=0;
	if (Y_METHOD == 2 || Y_METHOD == 3){
		// eigen-decomposition
		L("Decompose X = Q * D * Q'\n");
		lapack::syev( 'V', 'L', Eigv, lambda, lapack::minimal_workspace() );
		max_lambda = max_element(lambda.begin(),lambda.end());
		max_lambda_idx = std::distance(lambda.begin(),max_lambda);
		ublas::vector<double> lambda_sqrt(n);
		std::transform(lambda.begin(),lambda.end(),lambda_sqrt.begin(),my_sqrt);
	}

	if(Y_METHOD == 3) // FIXME: ??? war 1
		V = prod(trans(V),Eigv);

    // What do we want to maximize?
	// maximize symmetric version of C
	ublas::symmetric_adaptor<SDPProb::MatT, ublas::upper> B(prob.C);
	// maximize original C
	//SDPProb::MatT& B = prob.C;

	L("Rank reduction using random hyperplanes\n");
	ublas::vector<double> best_y(n);
	double best_y_val = -1E6;
	ublas::vector<double> r(n);
	ublas::vector<double> y(n);
	int tries = gCfg().getInt("ser-gen.sdp-rand-plane-tries");
	if(Y_METHOD==2)
		tries = 1;
	ExactDescriptiveStatistics yvalstats("y_val ");
	for(int iter=0;iter<tries;iter++){
		// generate unit-vector
		generate(r.begin(), r.end(), drand48);
		r -= ublas::scalar_vector<double>(n,0.5);

		// normieren ist unnoetig, da hinterher y auf laenge 1 gebracht wird
		//r /= ublas::norm_2(r); 

		if(Y_METHOD == 1){
			//use combination of cholesky-vectors
			noalias(y) = prod(ublas::trans(V),r);
		}else if(Y_METHOD == 2){
			// use combination of eigen-vectors
			//noalias(y) = prod(Eigv,r);
			noalias(y) = sqrt(lambda(max_lambda_idx)) * ublas::column(Eigv,max_lambda_idx);
		}else if(Y_METHOD == 3){
			// use technique from Nemirovski, Roos, Terlaky
			for(unsigned int i=0;i<r.size();i++)
				r(i) = r(i) > 0 ? 1 : -1;
			noalias(y) =  prod(V,r);
		}
		y /= ublas::norm_2(y);

		double y_val = inner_prod(y,prod(B,y));
		yvalstats += y_val;

		if(y_val> best_y_val){
		    // the solver always _maximizes_
			best_y_val = y_val;
			best_y     = y;
		}
	}
	L("best_y_val = %2.10f\n",best_y_val);
	if(isVerbose()) cout << yvalstats<<endl;

	ublas::vector<double> x = best_y;

	//return readout_connected(x,adj);
	return readout_plain(x,adj);
}

#   define BEST_ELEM(X) max_element(X.begin(),X.end())
Serialization SDPSeriationGen::Impl::readout_plain(ublas::vector<double>& x,const AdjMat::AdjMatT& adj)
{
	unsigned int n = x.size();

	// tricky: make sure x > 0 at all times.
	x += ublas::scalar_vector<double>(n, 1 - (*min_element(x.begin(),x.end())));

	Serialization::RankT ranks(n);
	std::vector<bool> done(n,false);

	// find highest component of x
	ublas::vector<double>::iterator it = BEST_ELEM(x);
	int idx = std::distance(x.begin(),it);

	L("Determine Actual Path through Graph.\n");
	for(unsigned int i=0;i<n;i++){
		ranks[i] = idx;
		done[idx] = true;
		*it = 0.0; 
		it = BEST_ELEM(x);
		idx = std::distance(x.begin(),it);
	}
	return Serialization(ranks);
}

Serialization SDPSeriationGen::Impl::readout_connected(ublas::vector<double>& x,const AdjMat::AdjMatT& adj)
{
	unsigned int n = x.size();

	// tricky: make sure x > 0 at all times.
	x += ublas::scalar_vector<double>(n, 1 - (*min_element(x.begin(),x.end())));

	Serialization::RankT ranks(n);
	std::vector<bool> done(n,false);

	// find highest component of x
	ublas::vector<double>::iterator it = BEST_ELEM(x);
	int idx = std::distance(x.begin(),it);


	L("Determine Actual Path through Graph.\n");
	for(unsigned int i=0;i<n;i++){
		// mark as visited
		ranks[i] = idx;
		done[idx] = true;

		// make sure we do not visit again
		*it = 0.0; 

		// [m,i] = max(x.*A(:,i));
		ublas::vector<double> adjcol = ublas::column(adj,idx);
		for(unsigned int j=0;j<adjcol.size();j++)
			if(adjcol[j]>0.00001) adjcol[j]=1;
		ublas::vector<double> tmp = ublas::element_prod(x,adjcol);
		it = BEST_ELEM(tmp);

		if( *it < 0.000000001 && i<n-1)
		{
			// if *it small, then either x[it] visited or adj(old_idx,idx) not connected
			// --> we reached a dead end, find next best start point
			it  = BEST_ELEM(x);
			idx = std::distance(x.begin(),it);
		}else{
			idx = std::distance(tmp.begin(),it);
			// point it in x, not tmp:
			it  = x.begin();
			it += idx;
		}


	}
	return Serialization(ranks);
}


void SDPSeriationGen::Impl::setSDPWrapper(std::auto_ptr<SDPWrapper> w){
  mSDPWrapper = w;
}

// Wrap SeriationGen::Impl
SDPSeriationGen::SDPSeriationGen()
:mImpl(new Impl)
{
}


Serialization SDPSeriationGen::operator()(const ProbAdjPerm& pap)
{
    return (*mImpl)(pap);
}

void SDPSeriationGen::setSDPWrapper(std::auto_ptr<SDPWrapper> w){
    mImpl->setSDPWrapper(w);
}

SDPSeriationGen::~SDPSeriationGen()
{
}

namespace{ registerInFactory<SerGenAdj, SDPSeriationGen> registerBase("SDPSeriationGen"); }
