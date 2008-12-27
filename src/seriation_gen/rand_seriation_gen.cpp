#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <factory/factory.h>
#include "rand_seriation_gen.hpp"

using namespace std;

// 2008-11-22 Hannes Schulz <mail at hannes-schulz dot de>

RandSerGen::RandSerGen()
{
}

RandSerGen::~RandSerGen()
{
  // cleanup
}

Serialization RandSerGen::operator()(const ProbAdjPerm& pap)
{
	int n = pap.getAdjMat()->size1();
	Serialization::RankT ranks(n);
	for(int i=0;i<n;i++)
		ranks[i]=i;
	random_shuffle(ranks.begin(),ranks.end());
	return Serialization(ranks);
}

void RandSerGen::configure()
{
}

namespace{ registerInFactory<SerGenAdj, RandSerGen> registerBase("RandSeriationGen"); }
