#ifndef __PLAINSERGEN_HPP__
#define __PLAINSERGEN_HPP__

#include "SerGenAdj.hpp"

/**
 * @brief PlainSerGen generates a seriation in the same order as the adjacency matrix
 * @author Hannes Schulz <mail at hannes-schulz dot de>
 * @date 2008-11-22
 */

class PlainSerGen : public SerGenAdj
{
  public:

    /**
     * Default constructor
     */
    PlainSerGen();

    /**
     * Destructor
     */
    virtual ~PlainSerGen();

	/**
	 * Serialize
	 */
	virtual Serialization operator()(const ProbAdjPerm& pap);

	virtual void configure();

  private:

};

#endif /* #ifndef __PLAINSERGEN_HPP__ */
