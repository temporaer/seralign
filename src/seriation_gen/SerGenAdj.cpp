#include <stdexcept>
#include "SerGenAdj.hpp"
#include <factory/factory.h>
#include <boost/numeric/ublas/vector.hpp>

using namespace std;

// 2008-11-15 Hannes Schulz <mail at hannes-schulz dot de>

SerGenAdj::SerGenAdj()
{
}

SerGenAdj::~SerGenAdj()
{
  // cleanup
}

Serialization SerGenAdj::operator()(const ProbAdjPerm& pap)
{
	throw logic_error("SerGenAdj::operator() called w/o subclassing!");
}

void SerGenAdj::configure()
{
}

namespace{ registerInFactory<SerGenAdj, SerGenAdj> registerBase("SerGenAdj"); }
