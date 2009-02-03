#ifndef __HEATKERNELGRAPHEMBEDDER_HPP__
#define __HEATKERNELGRAPHEMBEDDER_HPP__
#include <graph_embedder.hpp>

class HeatkernelGraphEmbedder : public GraphEmbedder
{
  public:
    HeatkernelGraphEmbedder();
	virtual cloud_type operator()(const ProbAdjPerm& pap, unsigned int dim);
	virtual void configure();
    virtual ~HeatkernelGraphEmbedder();

  private:
	double mT;

};

#endif /* #ifndef __HEATKERNELGRAPHEMBEDDER_HPP__ */
