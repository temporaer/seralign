#ifndef __GRAPHEMBEDDER_HPP__
#define __GRAPHEMBEDDER_HPP__
#include <ProbData.hpp>
#include <boost/numeric/ublas/fwd.hpp>
#include <vector>

class GraphEmbedder
{
  public:
	  typedef boost::numeric::ublas::vector<double> point_type;
	  typedef std::vector<point_type>      cloud_type;

    GraphEmbedder();
	virtual cloud_type operator()(const ProbAdjPerm& pap, unsigned int dim);
    virtual void configure();
    virtual ~GraphEmbedder();

  private:

};

#endif /* #ifndef __GRAPHEMBEDDER_HPP__ */
