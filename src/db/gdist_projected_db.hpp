#ifndef __GDIST_PROJECTED_DB_HPP__
#define __GDIST_PROJECTED_DB_HPP__
#include <ProbData.hpp>
#include <icp.hpp>
#include <boost/tuple/tuple.hpp>

#define DIM 3
class GraphFromAdj;
class GDistProjectedDB{
	public: // types
		typedef boost::numeric::ublas::vector<double> point_type;
		typedef util::ICP<DIM,point_type>             TICP;
		typedef TICP                                  value_type;
		typedef std::vector<TICP>                     TICPVec;
		typedef TICPVec::iterator                     iterator;
	public:
		inline iterator begin(){return mDB.begin();}
		inline iterator end()  {return mDB.end();}
		void add(const std::string& id, const ProbAdjPerm& pap);
		void configure();
	private:
		TICPVec           mDB;
		std::string       mSeriationGenName;
		int               mICPMaxIters;
		float             mICPLambda;
		int               mFastmapTries;
		std::vector<std::string> mIDs;

		// how good this is for mapping (=longest distance from node i)
		boost::tuple<float,unsigned int,unsigned int>  
			getFastmapQuality(GraphFromAdj&, int seed);
};
#endif /* __GDIST_PROJECTED_DB_HPP__ */