#include <stdexcept>
#include <factory/factory.h>
#include "graph_embedder.hpp"

using namespace std;

// 2009-01-17 Hannes Schulz <mail at hannes-schulz dot de>

GraphEmbedder::GraphEmbedder()
{
}
void GraphEmbedder::configure()
{
}

GraphEmbedder::cloud_type 
GraphEmbedder::operator()(const ProbAdjPerm& pap, int dim)
{
	throw runtime_error("GraphEmbedder operator() not overloaded by subclass!");
}


GraphEmbedder::~GraphEmbedder()
{
  // cleanup
}

namespace{ registerInFactory<GraphEmbedder, GraphEmbedder> registerBase("GraphEmbedder"); }
