#ifndef __GDISTSERGEN_HPP__
#define __GDISTSERGEN_HPP__

#include "SerGenAdj.hpp"

/**
 * @brief GDistSeriationGen generates a seriation based on the distances of
 * nodes in the graph from a start node n. The start node is determined by
 * another serialization method, which can be configured.
 * @author Hannes Schulz <mail at hannes-schulz dot de>
 * @date 2008-12-20
 */

class GDistSeriationGen : public SerGenAdj
{
  public:

    /**
     * Default constructor
     */
    GDistSeriationGen();

    /**
     * Destructor
     */
    virtual ~GDistSeriationGen();

	/**
	 * Serialize
	 */
	virtual Serialization operator()(const ProbAdjPerm& pap);

	virtual void configure();

  private:
	std::string mSeriationGenName;
	bool        mNormalize;

};

#endif /* #ifndef __GDISTSERGEN_HPP__ */
