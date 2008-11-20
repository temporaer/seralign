#include <stdexcept>
#include <factory/factory.h>
#include "adjmat_gen.hpp"
using namespace std;

void AdjMatGen::configure()
{
}

ProbAdjPerm AdjMatGen::operator()()
{
	throw logic_error("Called AdjMatGen() w/o subclassing");
}

AdjMatGen::~AdjMatGen()
{
}
bool AdjMatGen::hasNext()
{
	return true;
}

string AdjMatGen::getPrologDescription(int idx)
{
	throw runtime_error("AdjMatGen::getPrologDescription() not implemented for subclass"); 
}

string AdjMatGen::getGraphID()
{
	throw runtime_error("AdjMatGen::getGraphID() not implemented for subclass"); 
}


namespace{
	registerInFactory<AdjMatGen, AdjMatGen> registerBase("AdjMatGen");
}


