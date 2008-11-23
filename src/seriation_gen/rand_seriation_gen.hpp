#ifndef __RANDSERGEN_HPP__
#define __RANDSERGEN_HPP__

#include "SerGenAdj.hpp"

/**
 * @brief RandSerGen generates a random seriation
 * @author Hannes Schulz <mail at hannes-schulz dot de>
 * @date 2008-11-22
 * <+long description+>
 */

class RandSerGen : public SerGenAdj
{
  public:

    /**
     * Default constructor
     */
    RandSerGen();

    /**
     * Destructor
     */
    virtual ~RandSerGen();

	/**
	 * Serialize
	 */
	virtual Serialization operator()(const ProbAdjPerm& pap);

	virtual void configure();

  private:

};

#endif /* #ifndef __RANDSERGEN_HPP__ */
