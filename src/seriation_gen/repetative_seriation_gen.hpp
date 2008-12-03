#ifndef __REPETATIVESERGEN_HPP__
#define __REPETATIVESERGEN_HPP__

#include <string>
#include "SerGenAdj.hpp"

/**
 * @brief RepetativeSerGen Generates a Seriation by repetatively calling
 *                         another seriation on a jumbled adjacency matrix.
 * @author Hannes Schulz <mail at hannes-schulz dot de>
 * @date 2008-12-02
 * <+long description+>
 */

class RepetativeSerGen: public SerGenAdj
{
  public:

    /**
     * Default constructor
     */
    RepetativeSerGen();

    /**
     * Destructor
     */
    virtual ~RepetativeSerGen();

	/**
	 * Serialize
	 */
	virtual Serialization operator()(const ProbAdjPerm& pap);

	virtual void configure();

	inline void setSerializer(std::string s){mSerializer=s;}
	inline void setRepetitions(int i){mRepetitions=i;}

  private:

	std::string mSerializer;
	int         mRepetitions;

};

#endif /* #ifndef __REPETATIVESERGEN_HPP__ */
