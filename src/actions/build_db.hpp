#ifndef __BUILD_DB_HPP__
#define __BUILD_DB_HPP__

#include "action.hpp"

class BuildDB:public Action{
	virtual void operator()();
	virtual ~BuildDB();
};

#endif /* __BUILD_DB_HPP__ */
