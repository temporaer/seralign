#ifndef __FASTMAP_DISTANCEMAT_HPP__
#define __FASTMAP_DISTANCEMAT_HPP__

#include "action.hpp"

class FastmapDistancemat:public Action{
	virtual void operator()();
	virtual ~FastmapDistancemat();
};

#endif /* __FASTMAP_DISTANCEMAT_HPP__ */
