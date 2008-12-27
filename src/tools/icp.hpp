#ifndef __ICP_HPP__
#define __ICP_HPP__
#include <kdtree++/kdtree.hpp>
#include <arg_max.hpp>
#include <numeric>
#include <functional>
#include "boost/tuple/tuple.hpp"
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/bind.hpp>
#include <boost/math/quaternion.hpp>
#include <boost/numeric/bindings/lapack/lapack.hpp>
#include <quad_helper.hpp>


/**
 * @brief ICP The iterative Closest Point Algorithm
 * @author Hannes Schulz <mail at hannes-schulz dot de>
 * @date 2008-12-26
 */

namespace util{

#define ICP_TMPL    template<int PtDim, class T, class AccessorT>
#define ICP_FQCLASS ICP<PtDim,T,AccessorT>

template<int i,class T> 
struct TupleNComp
{ inline bool operator()(const T& a, const T& b){
		return a.template get<i>() < b.template get<i>();
	}
};
template<int i,class T> 
struct TupleN
{ inline bool operator()(const T& a, const T& b){
		return a.template get<i>() + b.template get<i>();
	}
};

ICP_TMPL
class ICP
{
	public:
		typedef float                                               TPrecision;
		typedef T                                                   TPoint;
		typedef KDTree::KDTree<PtDim,TPoint, AccessorT>             TTree;
		typedef boost::numeric::ublas::vector<TPrecision>           TVec; 
		typedef boost::numeric::ublas::matrix<TPrecision, boost::numeric::ublas::column_major>           TMat; 
		typedef boost::math::quaternion<TPrecision>                 TQuat; 
	public:
		struct VectorTranslator
		{ inline TPoint operator()(const TPoint& orig, const TVec& t){return orig+t;} };
		struct VectorRotator
		{ inline TPoint operator()(const TPoint& orig, const TMat& R){return prod(R, orig);} };

		ICP(){ }

		template<class Iter>
		void registerModel(Iter begin, Iter end);
		template<class Iter, class Translator>
		void registerModel(Iter begin, Iter end, Translator t);

		template<class Iter>
		void match(Iter begin, Iter end);
		template<class Iter, class Translator>
		void match(Iter begin, Iter end, Translator t);

		inline TVec getTrans(){return mDeterminedTranslation;}
		inline TMat getRot(){return mDeterminedRotation;}

	private:
		TTree     mModelTree;
		TVec      mModelCentroid;

		TMat      mDeterminedRotation;
		TVec      mDeterminedTranslation;

		template<class Iter>
		boost::tuple<TVec,int>   getCentroid(Iter begin, Iter end);
		TQuat getQuat(const TMat&);
};

ICP_TMPL
template<class Iter>
boost::tuple<typename ICP_FQCLASS::TVec,int>   
ICP_FQCLASS::getCentroid(Iter begin, Iter end){
	Iter it = begin;

	// determine mean of data
	TVec v(PtDim,0);
	int cnt=0;
	while(it!=end, cnt++)
		v += (TVec)(*it++);
	v /= (TPrecision) cnt;
	return boost::make_tuple(v,cnt);
}

/**
 * Register the model points
 */
ICP_TMPL
template<class Iter>
void ICP_FQCLASS::registerModel(Iter begin, Iter end){ 
	registerModel(begin,end,VectorTranslator());
}
/**
 * Register the model points
 */
ICP_TMPL
template<class Iter, class Translator>
void ICP_FQCLASS::registerModel(Iter begin, Iter end, Translator translate){ 
	int cnt;

	// determine mean of data
	boost::tie(mModelCentroid, cnt) = getCentroid(begin,end);

	// put mean-less (^^) points in tree
	mModelTree = TTree(); 
	Iter it = begin;
	while(it!=end)
		mModelTree.insert( translate(*it++,-mModelCentroid) );
}

ICP_TMPL
typename ICP_FQCLASS::TQuat
ICP_FQCLASS::getQuat(const TMat& M){
	// as in Horn 1987
	TPrecision Sxx = M(0,0), Syy = M(1,1), Szz = M(2,2),
			   Syx = M(1,0), Szx = M(2,0), Szy = M(2,1),
			   Sxy = M(0,1), Sxz = M(0,2), Syz = M(1,2);
	TMat N(4,4);
	// diagonal
	N(0,0) = Sxx+Syy+Szz;
	N(1,1) = Sxx-Syy-Szz;
	N(2,2) =-Sxx+Syy-Szz;
	N(3,3) =-Sxx-Syy+Szz;
	N(1,0) = Syx-Szx;        // 1st column
	N(2,0) = Szx-Sxz;
	N(3,0) = Sxy-Syx;
	N(0,1) = Syz-Szy;        // 2nd column
	N(2,1) = Sxy+Syz;
	N(3,1) = Szx+Sxz;
	N(0,2) = Szx-Sxz;        // 3rd column
	N(1,2) = Sxy+Syx;
	N(3,2) = Syz+Szy;
	N(0,3) = Sxy-Syx;        // 4th column
	N(1,3) = Szx+Sxz;
	N(2,3) = Syz+Szy;

	TMat Eigv(4,4);
	TVec lambda(4);
	boost::numeric::bindings::lapack::syev( 'V', 'L', Eigv, lambda, boost::numeric::bindings::lapack::minimal_workspace() );
	typename TVec::iterator best_lambda     = min_element(lambda.begin(),lambda.end());   
	int best_lambda_idx = std::distance(lambda.begin(),best_lambda);  
	TVec leading = boost::numeric::ublas::column(Eigv,best_lambda_idx);
	return TQuat(leading[0], leading[1], leading[2], leading[3]);
}

/**
 * Match points to the model
 */
ICP_TMPL
template<class Iter>
void ICP_FQCLASS::match(Iter begin, Iter end){
	match(begin, end, VectorTranslator());
}
ICP_TMPL
template<class Iter, class Translator>
void ICP_FQCLASS::match(Iter begin, Iter end, Translator translate){
	int nIter = 10;
	typedef boost::tuple<TPrecision, Iter,typename TTree::iterator>  TMatch; 
	typedef std::vector<TMatch>                                      TMatchList;
	TMatchList matchlist;

	// determine centroid of query
	TVec queryCentroid; int p;
	boost::tie(queryCentroid, p) = getCentroid(begin,end);
	queryCentroid = -queryCentroid;
	std::for_each(begin,end,boost::bind<void>(translate,::_1,queryCentroid));

	matchlist.resize(p); 
	typename TMatchList::iterator oit = matchlist.begin();
	TPrecision old_err = INT_MAX;
	TMat R(PtDim,PtDim);
	TVec t(PtDim);

	while(true){
		// step 1: find best match
		Iter it = begin;
		pair<typename TTree::iterator, TPrecision> res; 
		while(it!=end) {
			res = mModelTree.find_nearest(*it);
			*oit++ = boost::tie(res.second, it,res.first);
			it++;
		}
		// TODO: we want squared distance -- what is in find_nearest?

		// step 2: sort according to dist. 
		sort(matchlist.begin(), matchlist.end(),TupleNComp<0,TMatch>());

		// step 3: stop if conditions apply
		if(nIter==0) break;
		int   Npo = p/2;
		typename TMatchList::iterator mit, mend = matchlist.begin()+Npo;
		TPrecision sts = std::accumulate(matchlist.begin(),mend,0,TupleN<0,TMatch>());
		TPrecision   e = sts/nIter;
		if(e < 0.0001) break;
		if(fabs(e-old_err)<0.0001) break;
		old_err = e;

		// step 4: compute optimal motion (Horn, 1987)
		// convention: left: Query;  right: Model
		// 4a: Determine Rotation Matrix
		TMat M(boost::numeric::ublas::zero_matrix<TPrecision>(PtDim, PtDim));
		mit = matchlist.begin();
		while(mit!=mend){
			M += boost::numeric::ublas::outer_prod((TVec)(*mit->template get<1>()), (TVec)(*mit->template get<2>()));
			mit++;
		}
		TQuat q(getQuat(M));
		R = quaternion_to_R3_rotation<TPrecision>(q);

		// TODO: Determine Scale factor (necessary ???)
		// 4b: Determine Translation.
		mit = matchlist.begin();
		VectorRotator rotate;
		t = boost::numeric::ublas::scalar_vector<TPrecision>(PtDim,0);
		while(mit!=matchlist.end()) {
			rotate(*mit->template get<1>(), R);
			t += (TVec)(*mit->template get<2>()) - (TVec)(*mit->template get<1>());
			mit++;
		}
		nIter--;
	}
	mDeterminedRotation = R;
	mDeterminedTranslation = t;
}

}  // namespace util

#endif /* #ifndef __ICP_HPP__ */
