#ifndef __BUILD_DB_HPP__
#define __BUILD_DB_HPP__

#include <memory>
#include "action.hpp"
#include <db.hpp>

class AdjMatGen;
class PostProc;
class BuildDB:public Action{
	public:
		virtual void operator()();
		virtual ~BuildDB();
	private:
		std::auto_ptr<GraphDB> mDB;
		//void printDistMatrix(int cnt);
		//int knn_classify(int cnt, int id, AdjMatGen&, PostProc&, int k);
		//double match(GDistProjectedDB::TICP&,GDistProjectedDB::TICP&);
		void spatialAnalysis(int cnt, AdjMatGen&);
		void spatialAnalysisCloud(int cnt, AdjMatGen&);
};

#endif /* __BUILD_DB_HPP__ */
