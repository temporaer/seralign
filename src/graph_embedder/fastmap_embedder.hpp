#ifndef __FASTMAPGRAPHEMBEDDER_HPP__
#define __FASTMAPGRAPHEMBEDDER_HPP__

#include <boost/tuple/tuple.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <Serialization.hpp>
#include <configuration.hpp>
#include "graph_embedder.hpp"

class FastmapGraphEmbedder : public GraphEmbedder
{
  public:
    FastmapGraphEmbedder();
	virtual cloud_type operator()(const ProbAdjPerm& pap, int dim);
	virtual void configure();
    virtual ~FastmapGraphEmbedder();

  private:
	int               mFastmapTries;

	// how good this is for mapping (=longest distance from node i)
	boost::tuple<float,unsigned int,unsigned int>  
		getFastmapQuality(GraphFromAdj&,int k, int seed);

};

#endif /* #ifndef __FASTMAPGRAPHEMBEDDER_HPP__ */
