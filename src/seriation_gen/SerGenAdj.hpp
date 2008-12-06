#ifndef __SERGENADJ_HPP__
#define __SERGENADJ_HPP__

/**
 * @brief SerGenAdj base class for serializers based on ProbAdjPerm
 * @author Hannes Schulz <mail at hannes-schulz dot de>
 * @date 2008-11-15
 */

#include <ProbData.hpp>
#include <Serialization.hpp>

class SerGenAdj
{
  public:

    /**
     * Default constructor
     */
    SerGenAdj();

    /**
     * Destructor
     */
    virtual ~SerGenAdj();

	/**
	 * Serialize
	 */
	virtual Serialization operator()(const ProbAdjPerm& pap);

	virtual void configure();

	/// whether to be verbose
	static inline void setVerbose(bool b){mVerbose=b;}
	/// whether to be verbose
	static inline bool isVerbose()       {return mVerbose;}

  private:
	static bool mVerbose;

};

#endif /* #ifndef __SERGENADJ_HPP__ */
