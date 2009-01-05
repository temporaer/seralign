#ifndef __DBPOINT_HPP__
#define __DBPOINT_HPP__
#include <boost/numeric/ublas/vector.hpp>
#include <boost/math_fwd.hpp>

namespace boost{
	namespace math{
		template<class T> T conj(T);
	}
}

template <class T, int DIM=3>
struct DBPoint {
	typedef boost::numeric::ublas::vector<T, boost::numeric::ublas::bounded_array<T,DIM> > pos_type;
	typedef T value_type;
	pos_type pos;
	int id;
	DBPoint() : pos(DIM,0){}
	inline T operator[](unsigned int i)const{return pos[i];}
	inline operator const pos_type& ()const {return pos;}
	struct Translator;

	template<class _Quat>
	struct Rotator;
};

template <class T, int DIM>
struct DBPoint<T,DIM>::Translator
{	
	typedef DBPoint<T,DIM> TPoint;
	typedef typename TPoint::pos_type TVec;
	inline TPoint operator()(const TPoint& orig, const TVec& t){
		TPoint res = orig;
		res.pos   += t;
		return res;
	} 
};
template <class T, int DIM>
template <class _Quat>
struct DBPoint<T,DIM>::Rotator
{ 
	typedef DBPoint<T,DIM> TPoint;
	typedef _Quat          TQuat;
	//inline TPoint operator()(const TPoint& orig, const TMat& R){return prod(R, orig);} 
	inline TPoint operator()(const TPoint& orig, const TQuat& q){
		TQuat a=q * TQuat(0,orig[0],orig[1],orig[2]) * boost::math::conj(q);
		TPoint res = orig;
		res.pos[0] = a.R_component_2();
		res.pos[1] = a.R_component_3();
		res.pos[2] = a.R_component_4();
		return res;
	} 
};

#endif /* __DBPOINT_HPP__ */
