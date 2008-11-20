#ifndef __LEVSERIATION_HPP__
#define __LEVSERIATION_HPP__

#include <boost/shared_ptr.hpp>
#include <SerGenAdj.hpp>

/**
 * @brief LEVSeriation Serialize using leading eigenvector of Laplacian
 * @author Hannes Schulz <mail at hannes-schulz dot de>
 * @date 2008-11-17
 */

class LEVSeriationGen : public SerGenAdj
{
  public:

    /**
     * Default constructor
     */
    LEVSeriationGen();

    /**
     * Destructor
     */
    virtual ~LEVSeriationGen();

	/**
	 * Serialize
	 */
	virtual Serialization operator()(const ProbAdjPerm& pap);

  private:

  	struct Impl;
	boost::shared_ptr<Impl> mImpl;

};

#endif /* #ifndef __LEVSERIATION_HPP__ */
