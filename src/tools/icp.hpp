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
#include <boost/lambda/lambda.hpp>
#include <boost/math/quaternion.hpp>
#include <boost/numeric/bindings/lapack/lapack.hpp>
#include <quad_helper.hpp>

#define V(X) #X << "=" << (X) <<" "

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
{ inline double operator()(double d, const T& b){
		return d + b.template get<i>();
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
		template<class Iter, class Translator, class Rotator>
		void match(Iter begin, Iter end, Translator t, Rotator r);

		inline TVec  getTrans(){return mDeterminedTranslation;}
		inline TQuat getRot(){return mDeterminedRotation;}

	public:
		TTree     mModelTree;
		TVec      mModelCentroid;

		TQuat     mDeterminedRotation;
		TVec      mDeterminedTranslation;

		template<class Iter>
		inline boost::tuple<TVec,int>   getCentroid(Iter begin, Iter end);

		template<class Iter, class OIt>
		inline void determineBestMatches(Iter begin, Iter end, OIt out);

		template<class Iter>
		inline TQuat getQuat(Iter begin, Iter end);
};

ICP_TMPL
template<class Iter, class OIt>
void 
ICP_FQCLASS::determineBestMatches(Iter begin, Iter end, OIt oit){
	std::pair<typename TTree::iterator, TPrecision> res; 
	for(Iter it=begin; it!=end; it++) {
		res = mModelTree.find_nearest(*it);
		//cout << "Match: " << (*it) << "  --  " << (*res.first)<< " dist="<<res.second<<endl;
		*oit++ = boost::make_tuple(res.second, it,res.first);  // TODO: we want squared distance -- what is in find_nearest?
	}
}

ICP_TMPL
template<class Iter>
boost::tuple<typename ICP_FQCLASS::TVec,int>   
ICP_FQCLASS::getCentroid(Iter begin, Iter end){
	Iter it = begin;

	// determine mean of data
	TVec v(PtDim,0);
	int cnt=0;
	for(cnt=0; it!=end; cnt++, it++)
		v += (TVec)(*it);
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
	for(Iter it = begin; it!=end; it++)
		mModelTree.insert( translate(*it,-mModelCentroid) );
}

ICP_TMPL
template<class Iter>
typename ICP_FQCLASS::TQuat
ICP_FQCLASS::getQuat(Iter begin, Iter end){
	TMat M(boost::numeric::ublas::zero_matrix<TPrecision>(PtDim, PtDim));
	for(Iter mit = begin; mit!=end; mit++){
		M += boost::numeric::ublas::outer_prod((TVec)(*mit->template get<1>()), (TVec)(*mit->template get<2>()));
	}
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
	N(1,0) = Syz-Szy;        // 1st column
	N(2,0) = Szx-Sxz;
	N(3,0) = Sxy-Syx;
	N(0,1) = Syz-Szy;        // 2nd column
	N(2,1) = Sxy+Syx;
	N(3,1) = Szx+Sxz;
	N(0,2) = Szx-Sxz;        // 3rd column
	N(1,2) = Sxy+Syx;
	N(3,2) = Syz+Szy;
	N(0,3) = Sxy-Syx;        // 4th column
	N(1,3) = Szx+Sxz;
	N(2,3) = Syz+Szy;

	TVec lambda(4);
	boost::numeric::bindings::lapack::syev( 'V', 'L', N, lambda, boost::numeric::bindings::lapack::minimal_workspace() );
	typename TVec::iterator best_lambda     = max_element(lambda.begin(),lambda.end());   
	int best_lambda_idx = std::distance(lambda.begin(),best_lambda);  
	TVec leading = boost::numeric::ublas::column(N,best_lambda_idx);
	return TQuat(leading[0], leading[1], leading[2], leading[3]);
}

/**
 * Match points to the model
 */
ICP_TMPL
template<class Iter>
void ICP_FQCLASS::match(Iter begin, Iter end){
	match(begin, end, VectorTranslator(), VectorRotator());
}
ICP_TMPL
template<class Iter, class Translator, class Rotator>
void ICP_FQCLASS::match(Iter begin, Iter end, Translator translate, Rotator rotate){
	int nIter = 10000;
	mDeterminedRotation    = TQuat(1,0,0,0);
	mDeterminedTranslation = TVec(PtDim);
	typedef boost::tuple<TPrecision, Iter,typename TTree::iterator>  TMatch; 
	typedef std::vector<TMatch>                                      TMatchList;
	TMatchList matchlist;

	// determine centroid of query
	TVec queryCentroid; int p;
	boost::tie(queryCentroid, p) = getCentroid(begin,end);
	mDeterminedTranslation = -queryCentroid;
	std::transform(begin,end,begin,boost::bind<TPoint>(translate,::_1,mDeterminedTranslation));

	matchlist.resize(p); 
	typename TMatchList::iterator oit = matchlist.begin(), mit;
	TPrecision old_err = INT_MAX;
	TMat R(boost::numeric::ublas::identity_matrix<TPrecision>(PtDim,PtDim));
	TVec t(boost::numeric::ublas::scalar_vector<TPrecision>(PtDim, 0));

	while(true){
		// step 1: find best match
		determineBestMatches(begin,end,matchlist.begin());

		// step 2: sort according to dist. 
		sort(matchlist.begin(), matchlist.end(),TupleNComp<0,TMatch>());

		// step 3: stop if conditions apply
		if(nIter==0)               { cout << "break: nIter 0"<<endl;                           break; }

		float alpha = 0.4, beta = 1.0;

		int Npo;
		const float lambda = 2;
		float e;
		TupleN<0,TMatch> acc;
		for(int iter=0;iter<10;iter++){
			float x1    = alpha + 0.38197 * (beta - alpha);
			Npo = alpha * p + 0.5f;
			TPrecision e_alpha   = e =std::accumulate(matchlist.begin(),matchlist.end()+Npo,0.0,acc)/Npo;
			TPrecision phi_alpha = e_alpha/pow(alpha,1+lambda);
			Npo = x1    * p + 0.5f;
			TPrecision e_x1    = e =std::accumulate(matchlist.begin(),matchlist.end()+Npo,0.0,acc)/Npo;
			TPrecision phi_x1  = e_x1/pow(x1,1+lambda);
			if(phi_alpha < phi_x1) {
				beta = x1;
				continue;
			}
			float x2    = beta  - 0.38197 * (beta - alpha);
			Npo = x2    * p + 0.5f;
			TPrecision e_x2    = e =std::accumulate(matchlist.begin(),matchlist.end()+Npo,0.0,acc)/Npo;
			TPrecision phi_x2  = e_x2/pow(x2,1+lambda);
			if(phi_x1 < phi_x2) {
				alpha = x1; 
				beta = x2;
				continue;
			}
			alpha = x2;
		}
		cout << V(Npo) <<endl;

		if(fabs(e) < 0.0001)       { cout << "break: error small: "   <<e<<endl;               break; }
		if(fabs(e-old_err)<0.0001) { cout << "break: err diff small: "<<fabs(e-old_err)<<endl; break; }
		old_err = e;

		// step 4: compute optimal motion (Horn, 1987)
		// - convention: left: Query;  right: Model
		// - 4a: Determine Rotation Matrix
		TQuat q(getQuat(matchlist.begin(), matchlist.end()));
		mDeterminedRotation = mDeterminedRotation * q;
		R = quaternion_to_R3_rotation<TPrecision>(q);

		// TODO: Determine Scale factor (necessary ???)
		// - 4b: Determine Translation.
		t = boost::numeric::ublas::scalar_vector<TPrecision>(PtDim,0);
		for(mit = matchlist.begin(); mit!=matchlist.end(); mit++) {
			t += (TVec)(*mit->template get<2>()) - (TVec)(rotate(*mit->template get<1>(), R));
		}
		t/=(TPrecision)matchlist.size();
		mDeterminedTranslation += t;

		// step 5: Transform all points according to R, t
		std::transform(begin, end, begin, boost::bind<TPoint>(rotate, ::_1, R));
		std::transform(begin, end, begin, boost::bind<TPoint>(translate, ::_1, t));
		nIter--;
	}
}

}  // namespace util

#endif /* #ifndef __ICP_HPP__ */
