#ifndef __GDIST_PROJECTED_DB_HPP__
#define __GDIST_PROJECTED_DB_HPP__
#include <ProbData.hpp>
#include <icp.hpp>
#include <boost/tuple/tuple.hpp>
#include "db.hpp"
#include "DBPoint.hpp"

#define DIM 3
class GraphFromAdj;
class GDistProjectedDB : GraphDB{
	public: // types
		//typedef boost::numeric::ublas::vector<double> point_type;
		typedef double TPrecision;
		typedef DBPoint<TPrecision,DIM>               point_type;
		typedef util::ICP<DIM,point_type,TPrecision>  TICP;
		typedef TICP                                  value_type;
		typedef std::vector<TICP>                     TICPVec;
		typedef TICPVec::iterator                     iterator;
		typedef std::vector<point_type>               TCloud;
	public:
		inline iterator begin(){return mDB.begin();}
		inline iterator end()  {return mDB.end();}
		inline const TICP& back()const{ return mDB.back();}
		inline       TICP& back()     { return mDB.back();}
		inline const TICP& operator[](unsigned int i)const{ return mDB[i];}
		inline       TICP& operator[](unsigned int i){ return mDB[i];}
		TCloud add(const std::string& id, const ProbAdjPerm& pap);
		TCloud addSpectral(const std::string& id, const ProbAdjPerm& pap);
		TCloud addCompleteSpectral(const std::string& id, const ProbAdjPerm& pap);
		void configure();
		inline const std::vector<std::string>& getIDs(){return mIDs;}
	private:
		TICPVec           mDB;
		std::string       mSeriationGenName;
		int               mICPMaxIters;
		float             mICPLambda;
		int               mFastmapTries;
		std::vector<std::string> mIDs;

		// how good this is for mapping (=longest distance from node i)
		boost::tuple<float,unsigned int,unsigned int>  
			getFastmapQuality(GraphFromAdj&,int k, int seed);
};
#endif /* __GDIST_PROJECTED_DB_HPP__ */
