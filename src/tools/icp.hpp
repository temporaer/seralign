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

//#define V(X) #X << "=" << (X) <<" "

/**
 * @brief ICP The iterative Closest Point Algorithm
 * @author Hannes Schulz <mail at hannes-schulz dot de>
 * @date 2008-12-26
 */

namespace util{

#define ICP_TMPL    template<int PtDim, class T, class AccessorT>
#define ICP_FQCLASS ICP<PtDim,T,AccessorT>


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
		/// Moves a TPoint by adding a TVec to it
		struct VectorTranslator
		{ inline TPoint operator()(const TPoint& orig, const TVec& t){return orig+t;} };

		/// Rotates a TPoint by multiplying it with a quaternion q
		struct VectorRotator
		{ //inline TPoint operator()(const TPoint& orig, const TMat& R){return prod(R, orig);} 
		  inline TPoint operator()(const TPoint& orig, const TQuat& q){
			  TQuat a=q * TQuat(0,orig[0],orig[1],orig[2]) * boost::math::conj(q);
			  TPoint res(3);
			  res[0] = a.R_component_2();
			  res[1] = a.R_component_3();
			  res[2] = a.R_component_4();
			  return res;
		  } 
		};
		inline TVec rotateVector(const TVec& orig, const TQuat& q){
			TQuat a=q * TQuat(0,orig[0],orig[1],orig[2]) * boost::math::conj(q);
			TVec res(3);
			res[0] = a.R_component_2();
			res[1] = a.R_component_3();
			res[2] = a.R_component_4();
			return res;
		} 

		/// compare two tuples w.r.t. i-th component
		template<int i,class Tuple> struct TupleNComp
		{ inline bool operator()(const Tuple& a, const Tuple& b)
			{ return a.template get<i>() < b.template get<i>(); }
		};

		/// accumulate i-th component of a tuple
		template<int i,class Tuple> struct TupleNAcc
		{ inline double operator()(double d, const Tuple& b)
			{ return d + b.template get<i>(); } };

		/// Construct an ICP object
		ICP(float lambda=2.f, int overlapiter=10, float trim_alpha=0.4, float trim_beta=1.0) 
			: mVerbose(true), mLambda(lambda), mOverlapIterations(overlapiter), mTrimAlpha(trim_alpha), mTrimBeta(trim_beta)
		{ }

		/// Register a model. 
		/// Assumes that the TPoints can be moved by adding a TVec to them.
		template<class Iter>
		void registerModel(Iter begin, Iter end);

		/// Register a model. 
		/// Assumes that the TPoints can be moved using the supplied Translator t.
		template<class Iter, class Translator>
		void registerModel(Iter begin, Iter end, Translator t);

		/// Match something to the model.
		/// Assumes that the TPoints can be moved by adding a TVec to them.
		/// Assumes that the TPoints can be rotated by multiplication with a rotation matrix.
		template<class Iter>
		void match(Iter begin, Iter end);

		/// Match something to the model.
		/// Assumes that the TPoints can be moved using supplied Translator.
		/// Assumes that the TPoints can be rotated using supplied Rotator.
		template<class Iter, class Translator, class Rotator>
		void match(Iter begin, Iter end, Translator t, Rotator r);

		/// Returns the determined translation.
		inline TVec  getTrans(){return mDeterminedTranslation;}

		/// Returns the determined rotation as a quaternion.
		inline TQuat getRot(){return mDeterminedRotation;}

	public:
		const bool mVerbose;
		TTree      mModelTree;             ///< the registered model in a kD-Tree
		TVec       mModelCentroid;         ///< the model centroid
                                           
		TQuat      mDeterminedRotation;    ///< the determined rotation
		TVec       mDeterminedTranslation; ///< the determined translation

		const float  mLambda;              ///< parameter of phi-function in Chetverikov (Robust Euclidean alignment...)
		const int    mOverlapIterations;   ///< how many iterations of golden section search for optimal overlap
		const float  mTrimAlpha;           ///< trimming parameter: lower bound of interval to search for minimum of phi(xi)
		const float  mTrimBeta;            ///< trimming parameter: upper bound of interval to search for minimum of phi(xi)

		/// determine the centroid of the points between begin and end. 
		/// Assumes that they are convertible to TVec.
		template<class Iter>
		inline boost::tuple<TVec,int>   getCentroid(Iter begin, Iter end);

		/// determine the best match of the points in begin, end to model.
		/// The result is saved in the output-iterator out.
		template<class Iter, class OIt>
		inline void determineBestMatches(Iter begin, Iter end, OIt out);

		/// determine the quaternion for the matches between begin and end.
		template<class Iter>
		inline TQuat getQuat(Iter begin, Iter end);

		/// determine number of points to use and resulting error
		/// @param begin start of match-vector
		/// @param     p total number of points in match-vector
		/// @param   acc an accumulator to count errors in overlap section
		template<class Iter, class Acc>
		boost::tuple<TPrecision,int>
		getTrimmingParams(Iter begin, int p, Acc acc);
};

/*-----------------------------------------------------------------------------
 *  Implementation of class ICP
 *-----------------------------------------------------------------------------*/
ICP_TMPL
template<class Iter, class Accu>
boost::tuple<typename ICP_FQCLASS::TPrecision,int>
ICP_FQCLASS::getTrimmingParams(Iter begin, int p, Accu acc){
	TPrecision alpha = mTrimAlpha, beta=mTrimBeta, e;
	const TPrecision goldenSection = 0.38197;
	int n;

	for(int iter=0;iter<mOverlapIterations;iter++){
		TPrecision x1    = alpha + goldenSection * (beta - alpha);
		n = alpha * p + 0.5f;
		TPrecision e_alpha   = e =std::accumulate(begin,begin+n,0.0,acc)/n;
		TPrecision phi_alpha = e_alpha/pow(alpha,1+mLambda);
		n = x1    * p + 0.5f;
		TPrecision e_x1      = e =std::accumulate(begin,begin+n,0.0,acc)/n;
		TPrecision phi_x1    = e_x1/pow(x1,1+mLambda);
		if(phi_alpha < phi_x1) {
			beta = x1;
			continue;
		}
		TPrecision x2    = beta  - goldenSection * (beta - alpha);
		n = x2    * p + 0.5f;
		TPrecision e_x2      = e =std::accumulate(begin,begin+n,0.0,acc)/n;
		TPrecision phi_x2    = e_x2/pow(x2,1+mLambda);
		if(phi_x1 < phi_x2) {
			alpha = x1; 
			beta = x2;
			continue;
		}
		alpha = x2;
	}
	return boost::tie(e,n);
}

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
	mModelTree.optimize(); 
}

ICP_TMPL
template<class Iter>
typename ICP_FQCLASS::TQuat
ICP_FQCLASS::getQuat(Iter begin, Iter end){
	// as in Horn (1987)
	TMat M(boost::numeric::ublas::zero_matrix<TPrecision>(PtDim, PtDim));
	for(Iter mit = begin; mit!=end; mit++){
		M += boost::numeric::ublas::outer_prod((TVec)(*mit->template get<1>()), (TVec)(*mit->template get<2>()));
	}
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
	typedef boost::tuple<TPrecision, Iter,typename TTree::iterator>  TMatch; 
	typedef std::vector<TMatch>                                      TMatchList;
	TMatchList matchlist;

	// determine centroid of query
	TVec queryCentroid; int p;
	boost::tie(queryCentroid, p) = getCentroid(begin,end);
	mDeterminedTranslation = TVec(PtDim,0);
	std::transform(begin,end,begin,boost::bind<TPoint>(translate,::_1,-queryCentroid));

	matchlist.resize(p); 
	typename TMatchList::iterator oit = matchlist.begin(), mit, mend;
	TPrecision old_err = INT_MAX;
	TMat R(PtDim,PtDim);
	TVec t(PtDim);

	while(true){
		// step 1: find best match
		determineBestMatches(begin,end,matchlist.begin());

		// step 2: sort according to dist. 
		sort(matchlist.begin(), matchlist.end(),TupleNComp<0,TMatch>());

		// step 3: stop if conditions apply
		if(nIter==0)               { if(mVerbose) cout << "ICP break: nIter 0"<<endl;                           break; }

		// trimming: determine which overlap parameters to use
		int Npo; TPrecision e;
		boost::tie(e,Npo) = getTrimmingParams(matchlist.begin(),p, TupleNAcc<0,TMatch>()); 
		if(mVerbose) cout << "ICP Trimmed error ("<<Npo<<"/"<<p <<" is "<< e<<endl;
		mend = matchlist.begin() + Npo;

		if(fabs(e) < 0.0001)       { if(mVerbose) cout << "ICP break: error small: "   <<e<<endl;               break; }
		if(fabs(e-old_err)<0.0001) { if(mVerbose) cout << "ICP break: err diff small: "<<fabs(e-old_err)<<endl; break; }
		old_err = e;

		// step 4: compute optimal motion (Horn, 1987)
		// - convention: left: Query;  right: Model
		// - 4a: Determine Rotation Matrix from Npo matches
		TQuat q(getQuat(matchlist.begin(), mend));
		mDeterminedRotation = q * mDeterminedRotation ;

		// TODO: Determine Scale factor (necessary ???)
		// - 4b: Determine Translation from Npo matches.
		t = boost::numeric::ublas::scalar_vector<TPrecision>(PtDim,0);
		for(mit = matchlist.begin(); mit!=mend; mit++) {
			t += (TVec)(*mit->template get<2>()) - (TVec)(rotate(*mit->template get<1>(), q));
		}
		t/=(TPrecision)matchlist.size();
		mDeterminedTranslation += t;

		// step 5: Transform _all_ points according to R, t
		std::transform(begin, end, begin, boost::bind<TPoint>(rotate, ::_1, q));
		std::transform(begin, end, begin, boost::bind<TPoint>(translate, ::_1, t));
		nIter--;
	}
	mDeterminedTranslation += -rotateVector(queryCentroid,mDeterminedRotation) +mModelCentroid;
}
#undef ICP_TMPL    
#undef ICP_FQCLASS 

}  // namespace util

#endif /* #ifndef __ICP_HPP__ */
