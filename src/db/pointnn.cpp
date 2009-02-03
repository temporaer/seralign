#include <stdexcept>
#include <fstream>
#include <parallel/algorithm>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/ref.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <configuration.hpp>
#include <kdtree++/kdtree.hpp>
#include <adjmat_gen.hpp>
#include <progressbar.hpp>
#include "pointnn.hpp"

#define PARALLELIZE 0
#if PARALLELIZE
#  include <tbb/task_scheduler_init.h>
#  include <tbb/parallel_for.h>
#  include <tbb/blocked_range.h>
#  include <tbb/parallel_do.h>
#endif


#include <nana.h>

using namespace std;
namespace ublas = boost::numeric::ublas;
namespace ll    = boost::lambda;

// 2009-01-17 Hannes Schulz <mail at hannes-schulz dot de>

#define DIM 3
template<class _Prec, int _Dim>
struct Point{
	typedef Point<_Prec, _Dim>         self_t;
	typedef _Prec                      value_type;
	typedef ublas::vector<value_type>  pos_t;
	typedef unsigned char              klass_t;
	static const int  sDim = _Dim;
	Point(){ 
		throw runtime_error("Point() called");
	}
	Point(const pos_t& pos, const boost::any& ref, const unsigned int& nodeid)
		:mPos(pos),mRef(ref),mNodeID(nodeid)
	{
	}
	inline value_type operator[](unsigned int i)const{ return mPos[i]; }
	inline value_type dist (const self_t& p) const { return norm_2(mPos-p.mPos); }
	inline value_type dist2(const self_t& p) const { 
		const pos_t difference = mPos - p.mPos;
		return inner_prod(difference, difference);
	}
	pos_t          mPos;
	boost::any     mRef;
	unsigned int   mNodeID;
};

struct PointNNDB::Impl{
	typedef Point<float,DIM>             point_t;
	typedef KDTree::KDTree<DIM,point_t>  tree_t;
	typedef tree_t::iterator             tree_iterator;
	
	Impl() :mOS("/tmp/classrate")
	{
	}
	void add(const ProbAdjPerm& pap, const TCloud& cloud);
	inline void finish(){ mTree.optimize(); }
	inline void clear(){ mTree.clear(); }
	void evaluate(AdjMatGen&, int);
	int classify(const TCloud& cloud, AdjMatGen& gen, int k);

	tree_t mTree;
	ublas::vector<int>    mFeatureDir;
	ublas::vector<double> mFeatureMomentum;
	ofstream              mOS;
};

inline double my_fabs(const double& d){return fabs(d);}

template <class _Tree, class _Point>
struct NNHelper{
	typedef _Point point_t;
	typedef _Tree  tree_t;
	NNHelper(const tree_t& t, AdjMatGen& gen, unsigned int k, ProgressBar& pb):
		mTree(t),mGen(gen), mK(k), in_class_cnt(0), btw_class_cnt(0), class0cnt(0), class1cnt(0), mPb(pb)
		,btw_class_feat(mGen.getNumFeatures(),0)
		, in_class_feat(mGen.getNumFeatures(),0)
		, mRange(0.06)
	{
	}
	const _Tree& mTree;
	mutable AdjMatGen&   mGen; // make sure to call de-facto const methods _only_
	unsigned int mK;
	mutable int in_class_cnt, btw_class_cnt, class0cnt, class1cnt;
	mutable ProgressBar& mPb;  // this one is not so important...
	mutable AdjMatGen::feature_t btw_class_feat, in_class_feat;
	mutable float mRange;
	void operator()(const point_t& p)const{
		mPb.inc();
		vector<point_t> neigh;
		vector<point_t> k_nearest(mK+1,p);
		float range = mRange;

		// fill vector with k nearest neighbours
		while(true){
			unsigned int cnt = mTree.count_within_range(p,range);
			if(cnt>=mK+1) break;
			range *= 1.1;
		}
		mTree.find_within_range(p,range,back_inserter(neigh));
		if(neigh.size()>mK+5)
			range *= 0.7f;

		mRange = 0.9f*mRange + 0.1f*range;
		//cout <<"Query: " << p.mPos <<"Range: "<< mRange << " --> " <<neigh.size()<<endl;
		if(neigh.size() > 1.5f*mK)
		{
			int i=0;
			vector<float> dists(neigh.size());
			vector<int>   idx(neigh.size());
			transform(neigh.begin(), neigh.end(), dists.begin(), ll::bind(&point_t::dist2,p,ll::_1));
			generate(idx.begin(),idx.end(),ll::var(i)++);
			partial_sort(idx.begin(), idx.begin()+mK+1, idx.end(), 
					ll::var(dists)[ll::_1] < ll::var(dists)[ll::_2]);
			transform(idx.begin(), idx.begin()+mK+1, k_nearest.begin(), ll::var(neigh)[ll::_1]);
		}else
		{
			partial_sort_copy(neigh.begin(), neigh.end(), k_nearest.begin(), k_nearest.end(), ll::bind(&point_t::dist2,p,ll::_1) < ll::bind(&point_t::dist2,p,ll::_2));
		}

		// determine features of neighbours, dependent on between-class and in-class
		typename point_t::klass_t pklass = mGen.getClassID(p.mRef);
		if(pklass == 0) __sync_fetch_and_add(&class0cnt,1);
		else            __sync_fetch_and_add(&class1cnt,1);

		typename vector<point_t>::iterator it = k_nearest.begin(); 
		it++; // NOTE: ignore ourselves. Not necessary for new queries!
		for(;it!=k_nearest.end();it++){
			const point_t& q = *it;
			typename point_t::klass_t qklass = mGen.getClassID(q.mRef);
			bool in_class = qklass == pklass;
			//if(in_class){
			if(in_class){ 
				__sync_fetch_and_add(& in_class_cnt,1);
				in_class_feat += mGen.getFeatures(q.mNodeID,q.mRef);
			}
			else{
				__sync_fetch_and_add(&btw_class_cnt,1);
				btw_class_feat += mGen.getFeatures(q.mNodeID,q.mRef);
			}
		}
	}
};

template <class _Tree, class _Point>
struct NNClassifier{
	typedef _Point point_t;
	typedef _Tree  tree_t;
	NNClassifier(const tree_t& t, AdjMatGen& gen, unsigned int k):
		mTree(t),mGen(gen), mK(k), class0cnt(0), class1cnt(0), 
		mRange(0.06)
	{
	}
	const _Tree& mTree;
	mutable AdjMatGen&   mGen; // make sure to call de-facto const methods _only_
	unsigned int mK;
	mutable float class0cnt, class1cnt;
	mutable float mRange;
	void operator()(const point_t& p)const{
		vector<point_t> neigh;
		vector<point_t> k_nearest(mK,p);
		float range = mRange;
		if(fabs(accumulate(p.mPos.begin(), p.mPos.end(),0.0))<1E-6)
			return;

		// fill vector with k nearest neighbours
		while(true){
			unsigned int cnt = mTree.count_within_range(p,range);
			if(cnt>=mK+1) break;
			range *= 1.1;
		}
		mTree.find_within_range(p,range,back_inserter(neigh));
		if(neigh.size()>mK+5)
			range *= 0.7f;

		mRange = 0.9f*mRange + 0.1f*range;
		if(neigh.size() > 1.5f*mK)
		{
			int i=0;
			vector<float> dists(neigh.size());
			vector<int>   idx(neigh.size());
			transform(neigh.begin(), neigh.end(), dists.begin(), ll::bind(&point_t::dist2,p,ll::_1));
			generate(idx.begin(),idx.end(),ll::var(i)++);
			partial_sort(idx.begin(), idx.begin()+mK, idx.end(), 
					ll::var(dists)[ll::_1] < ll::var(dists)[ll::_2]);
			transform(idx.begin(), idx.begin()+mK, k_nearest.begin(), ll::var(neigh)[ll::_1]);
		}else
		{
			partial_sort_copy(neigh.begin(), neigh.end(), k_nearest.begin(), k_nearest.end(), ll::bind(&point_t::dist2,p,ll::_1) < ll::bind(&point_t::dist2,p,ll::_2));
		}

		double klass0weight = 0, klass1weight = 0;
		for(unsigned int i=0;i<k_nearest.size();i++){
			const point_t& q = k_nearest[i];
			typename point_t::klass_t qklass = mGen.getClassID(q.mRef);
			double d = p.dist2(q);
			//double d = 1.0;
			if(qklass == 0) klass0weight += d>0?1.0/d:1E6;
			else            klass1weight += d>0?1.0/d:1E6;
		}
		//cout<<boost::format("klasses: %6.6f   vs   %6.6f\n")%klass0weight%klass1weight;
		if(klass1weight < klass0weight)
			class0cnt += 1.0;
		else
			class1cnt += 1.0;
	}
};

void PointNNDB::Impl::evaluate(AdjMatGen& gen, int k){
	if(mFeatureDir.size()==0){
		mFeatureDir      = ublas::vector<int>  (gen.getNumFeatures(),1);
		mFeatureMomentum = ublas::vector<double>(gen.getNumFeatures(),0.01);
	}
	ProgressBar pb(mTree.size(), "Determine NN", 50);
	NNHelper<tree_t,point_t> nnh(mTree, gen, k,pb);
#if PARALLELIZE
	tbb::task_scheduler_init init;
	tbb::parallel_do(mTree.begin(),mTree.end(),nnh);
#else
	BOOST_FOREACH(const point_t& p, mTree){nnh(p);};
#endif
	//double c_rat = nnh.class0cnt/nnh.class1cnt;
	AdjMatGen::feature_t in_c = nnh.in_class_feat/nnh.in_class_cnt, 
						btw_c = nnh.btw_class_feat/nnh.btw_class_cnt;
	cout <<endl;
	cout << boost::format("Found %d/%d in-/btw-class links, %d/%d points (%3.3f cmp %3.3f)\n")%nnh.in_class_cnt%nnh.btw_class_cnt%nnh.class0cnt%nnh.class1cnt%((float)nnh.in_class_cnt/nnh.btw_class_cnt)%((float)nnh.class0cnt/nnh.class1cnt);
	cout << " InClassFeat: " << in_c<<endl;
	cout << "BtwClassFeat: " << btw_c<<endl;
	ublas::vector<double> fw = gen.getFeatureWeights();
	for(unsigned int i=0;i<fw.size();++i){
		// RPROP
		int dir    = 0;
		if(0);
		else if(in_c(i)-btw_c(i)>0) dir = +1; 
		else if(in_c(i)-btw_c(i)<0) dir = -1;
		int sgn = dir*mFeatureDir(i);
		if(0);
		else if(sgn>0)
			mFeatureMomentum(i) = min(1.2*mFeatureMomentum(i),50.0);
		else if(sgn<0)
			mFeatureMomentum(i) = max(0.5*mFeatureMomentum(i),1E-6);

		fw(i) += mFeatureMomentum(i) * sgn;
		mFeatureDir(i) = dir;
	}
	cout << "    Features: " << fw<<endl;
	cout << "    Moments : " << mFeatureMomentum<<endl;
	cout << "    Dirs    : " << mFeatureDir<<endl;
	
	gen.setFeatureWeights(fw);
	mOS << (float)nnh.in_class_cnt/nnh.btw_class_cnt;
	mOS << "\t";
	copy(fw.begin(),fw.end(),ostream_iterator<double>(mOS,"\t"));
	mOS << endl;
}
void PointNNDB::Impl::add(const ProbAdjPerm& pap, const TCloud& cloud){
	const boost::any& ref = pap.getBackground();

	for(unsigned int i=0;i<cloud.size();i++){
		I(cloud[i].size() == DIM);
		if(fabs(accumulate(cloud[i].begin(), cloud[i].end(),0.0)) < 1E-6) 
			continue;
		mTree.insert(point_t(cloud[i],ref,i));
	}
}

int PointNNDB::Impl::classify(const TCloud& cloud, AdjMatGen& gen, int k){
	NNClassifier<tree_t,point_t> nnc(mTree, gen, k);
	for(unsigned int i=0; i<cloud.size(); i++){
		point_t p(cloud[i], boost::any(), i);
		nnc(p);
	};
	//cout << boost::format("Classifying as %d   (%d vs. %d)\n")%(nnc.class1cnt>nnc.class0cnt)%nnc.class0cnt%nnc.class1cnt;
	return nnc.class0cnt>nnc.class1cnt ? 0 : 1;
}


// -------------------------------------------------
//                   Wrappers
// -------------------------------------------------

PointNNDB::PointNNDB()
	:mImpl(new Impl())
{
}
int PointNNDB::classify(const TCloud& cloud, AdjMatGen& gen, int k)
{
	return mImpl->classify(cloud, gen, k);
}


void PointNNDB::findNearest(const TCloud& cloud)
{
	throw runtime_error("PointNNDB::findNearest(cloud) not implemented.");
}


void PointNNDB::add(const ProbAdjPerm& pap, const TCloud& cloud)
{
	mPaps.push_back(pap);
	mImpl->add(pap, cloud);
}


void PointNNDB::init(int dim)
{
	I(dim==DIM);
}


PointNNDB::~PointNNDB()
{
}

void PointNNDB::configure()
{
	GraphDB::configure();
}

void PointNNDB::finish(){
	mImpl->finish();
}

void PointNNDB::clear(){
	mPaps.clear();
	mImpl->clear();
}

void PointNNDB::evaluate(AdjMatGen&gen, int k){
	mImpl->evaluate(gen, k);
}
