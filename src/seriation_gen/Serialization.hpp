#ifndef __SERIALIZATION_HPP__
#define __SERIALIZATION_HPP__

/**
 * @brief Serialization A Serialization
 * @author Hannes Schulz <mail at hannes-schulz dot de>
 * @date 2008-11-15
 */
#include <boost/numeric/ublas/vector.hpp>

class Serialization
{
  public:
	  typedef boost::numeric::ublas::vector<int>   RankT;
	  typedef boost::numeric::ublas::vector<float> PosT;
  private:
	  RankT mRanks;
	  PosT  mPositions;
  public:

	  Serialization(int n):mRanks(n),mPositions(n){}
	  Serialization(const RankT&r, const PosT&p):mRanks(r),mPositions(p){}
	  Serialization(const RankT&r):mRanks(r),mPositions(r.size()){}

	  inline void setRanks    (const RankT& r){ mRanks = r; }
	  inline void setPositions(const PosT& p) { mPositions = p; }

	  inline const RankT& getRanks    ()const{ return mRanks; }
	  inline const PosT&  getPositions()const{ return mPositions; }

	  virtual ~Serialization();

};

#endif /* #ifndef __SERIALIZATION_HPP__ */
