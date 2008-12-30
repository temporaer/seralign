#include "postproc.hpp"
#include <factory/factory.h>

using namespace std;

// 2008-11-21 Hannes Schulz <mail at hannes-schulz dot de>

PostProc::PostProc()
{
}

PostProc::~PostProc()
{
  // cleanup
}

void PostProc::atStart()
{
}

void PostProc::atEnd()
{
}

void PostProc::atSeriation(AdjMatGen& gen, Serialization& ser, ProbAdjPerm& prob)
{
}
void PostProc::atSeriation(AdjMatGen& gen, std::vector<boost::numeric::ublas::vector<double> >& cloud, ProbAdjPerm& prob)
{
}


namespace{ registerInFactory<PostProc, PostProc> registerBase("PostProc"); }
