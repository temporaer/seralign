#ifndef __LOCALIGN_DISTMAT_HPP__
#define __LOCALIGN_DISTMAT_HPP__

#include <iostream>
#include <fstream>
#include <vector>
#include <seqan/align.h>
#include <boost/numeric/ublas/fwd.hpp>

#include "postproc.hpp"

class AdjMatGen;
class ProbAdjPerm;

/**
 * @brief LocAlignDistmat aligns all serializations and outputs a distance matrix
 * @author Hannes Schulz <mail at hannes-schulz dot de>
 * @date 2008-11-21
 */

class LocAlignDistmat : public PostProc
{
  public:
	typedef boost::numeric::ublas::matrix<double> DistMatT;

    /**
     * Default constructor
     */
    LocAlignDistmat();

	virtual void atStart();
	virtual void atSeriation(AdjMatGen& gen, Serialization& ser, ProbAdjPerm& prob);
	virtual void atEnd();


    /**
     * Destructor
     */
    virtual ~LocAlignDistmat();

  private:

	struct Node{
		std::string s;
		inline bool operator== (const Node& n) const {return s == n.s;}
		inline operator unsigned char () const {return s[0];}
		Node(const std::string& t);
		Node();
	};

	std::vector<seqan::String<Node> > mStrings;
	std::vector<std::string>          mGraphIDs;
};

#endif /* #ifndef __LOCALIGN_DISTMAT_HPP__ */

