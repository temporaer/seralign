#ifndef __FEATUREGRAPHEMBEDDER_HPP__
#define __FEATUREGRAPHEMBEDDER_HPP__
#include <vector>
#include <fstream>
#include <boost/numeric/ublas/matrix.hpp>
#include <adjmat_gen.hpp>
#include <graph_embedder.hpp>

class FeatureGraphEmbedder : public GraphEmbedder
{
  public:
    FeatureGraphEmbedder();
	virtual cloud_type operator()(const ProbAdjPerm& pap, unsigned int dim);
	virtual void configure();
    virtual ~FeatureGraphEmbedder();
	void    setAdjMatGen(boost::shared_ptr<AdjMatGen> g, unsigned int dim);

  private:
	boost::shared_ptr<AdjMatGen> mGen;
	std::vector<boost::numeric::ublas::matrix<double> > mMats;
	std::ofstream mOS;
};

#endif /* #ifndef __FEATUREGRAPHEMBEDDER_HPP__ */
