#ifndef __SDP_SIMPLE_HPP__
#define __SDP_SIMPLE_HPP__

#include "action.hpp"

class Serialize:public Action{
	virtual void operator()();
	virtual ~Serialize();
};

#endif /* __SDP_SIMPLE_HPP__ */
