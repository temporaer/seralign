#ifndef __SPECTRALGRAPHEMBEDDER_HPP__
#define __SPECTRALGRAPHEMBEDDER_HPP__
#include <graph_embedder.hpp>

class SpectralGraphEmbedder : public GraphEmbedder
{
  public:
    SpectralGraphEmbedder();
	virtual cloud_type operator()(const ProbAdjPerm& pap, int dim);
	virtual void configure();
    virtual ~SpectralGraphEmbedder();

  private:

};

#endif /* #ifndef __SPECTRALGRAPHEMBEDDER_HPP__ */
