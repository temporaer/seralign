#ifndef __BUILD_DB_HPP__
#define __BUILD_DB_HPP__

#include "action.hpp"
#include <gdist_projected_db.hpp>

class AdjMatGen;
class PostProc;
class BuildDB:public Action{
	public:
		virtual void operator()();
		virtual ~BuildDB();
	private:
		GDistProjectedDB mDB;
		void printDistMatrix(int cnt);
		void printClosest(int cnt, int id, AdjMatGen&, PostProc&);
};

#endif /* __BUILD_DB_HPP__ */
