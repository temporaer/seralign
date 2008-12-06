#include <stdexcept>
#include <factory/factory.h>
#include "adjmat_gen.hpp"
#include <Serialization.hpp>
#include <configuration.hpp>
using namespace std;

void AdjMatGen::configure()
{
	setVerbose(gCfg().getBool("adjmat_gen.verbose"));
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

string AdjMatGen::getPlainDescription(int ser_idx, const Serialization&)
{
	throw runtime_error("AdjMatGen::getPlainDescription() not implemented for subclass"); 
}
string AdjMatGen::getPrologDescription(int ser_idx, const Serialization&)
{
	throw runtime_error("AdjMatGen::getPrologDescription() not implemented for subclass"); 
}

string AdjMatGen::getGraphID()
{
	throw runtime_error("AdjMatGen::getGraphID() not implemented for subclass"); 
}

std::string AdjMatGen::getGraphVizNodeAttribs(int idx)
{
	return string("");
}

bool AdjMatGen::mVerbose;

namespace{
	registerInFactory<AdjMatGen, AdjMatGen> registerBase("AdjMatGen");
}


