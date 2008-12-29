#ifndef __SERIALIZE_HPP__
#define __SERIALIZE_HPP__

#include "action.hpp"

class Serialize:public Action{
	virtual void operator()();
	virtual ~Serialize();
};

#endif /* __SERIALIZE_HPP__ */
