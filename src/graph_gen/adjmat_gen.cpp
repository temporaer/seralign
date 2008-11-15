#include <stdexcept>
#include <factory/factory.h>
#include "adjmat_gen.hpp"
using namespace std;

void AdjMatGen::configure()
{
}

boost::shared_ptr<AdjMatGen::AdjMatT> AdjMatGen::operator()()
{
	throw logic_error("Called AdjMatGen() w/o subclassing");
}

AdjMatGen::~AdjMatGen()
{
}

namespace{
	registerInFactory<AdjMatGen, AdjMatGen> registerBase("AdjMatGen");
}
