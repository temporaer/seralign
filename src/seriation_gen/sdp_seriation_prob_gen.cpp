#include <sdp_seriation_prob_gen.hpp>
#include <sdp_prob.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <nana.h>
#undef C

namespace ublas = boost::numeric::ublas;
using namespace std;

SDPSeriationProbGen::SDPSeriationProbGen(boost::shared_ptr<AdjMatT>a)
	: mAdjMat(a)
{
}

void SDPSeriationProbGen::operator()(SDPProb& prob)
{
	using ublas::identity_matrix;
	typedef AdjMatT Matrix ;
	
	Matrix& adj = *mAdjMat;

	I(adj.size1() == adj.size2());

	int n = adj.size1();

	// C = F_0
	// negative, since we want to _minimize_ instead of maximize as in normal SDP
	prob.C = adj;

	// F_1
	identity_matrix<double> Id(n);
	prob.F.push_back(Id);

	// b_1
	prob.b = ublas::vector<double>(1);
	prob.b(0) = 1;
}

