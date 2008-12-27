#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <factory/factory.h>
#include "plain_seriation_gen.hpp"

using namespace std;

// 2008-11-22 Hannes Schulz <mail at hannes-schulz dot de>

PlainSerGen::PlainSerGen()
{
}

PlainSerGen::~PlainSerGen()
{
  // cleanup
}

Serialization PlainSerGen::operator()(const ProbAdjPerm& pap)
{
	int n = pap.getAdjMat()->size1();
	Serialization::RankT ranks(n);
	for(int i=0;i<n;i++)
		ranks[i]=i;
	return Serialization(ranks);
}

void PlainSerGen::configure()
{
}

namespace{ registerInFactory<SerGenAdj, PlainSerGen> registerBase("PlainSeriationGen"); }
