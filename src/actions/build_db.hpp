#ifndef __BUILD_DB_HPP__
#define __BUILD_DB_HPP__

#include "action.hpp"
#include <gdist_projected_db.hpp>

class BuildDB:public Action{
	public:
		virtual void operator()();
		virtual ~BuildDB();
	private:
		GDistProjectedDB mDB;
};

#endif /* __BUILD_DB_HPP__ */
