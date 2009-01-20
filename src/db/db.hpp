#ifndef __DB_HPP__
#define __DB_HPP__
#include <stdexcept>
#include <ProbData.hpp>
#include <vector>
#include <boost/numeric/ublas/fwd.hpp>
class GraphDB{
	public:
		typedef boost::numeric::ublas::vector<double> point_type;
		typedef std::vector<point_type>               TCloud;
		virtual void add(const ProbAdjPerm& pap, const TCloud& cloud){};
		virtual boost::numeric::ublas::vector<int>    getFeatures(int i)
		{
			throw std::runtime_error("getFeatures() not overloaded by your DB");
		}
		virtual boost::numeric::ublas::vector<int>    getFeatures(const TCloud& cloud)
		{
			throw std::runtime_error("getFeatures() not overloaded by your DB");
		}
		virtual void init(int dim){};
		virtual void finish(){};
		inline unsigned int size() { return mPaps.size(); }
		inline const ProbAdjPerm& getPap(unsigned int i) { return mPaps[i]; }
		virtual void configure(){};
	protected:
		std::vector<ProbAdjPerm> mPaps;
};
#endif /* __DB_HPP__ */
