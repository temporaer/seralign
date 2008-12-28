#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE icptest

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <icp.hpp>
#include <quat_helper.hpp>
using namespace util;
using namespace boost::numeric;
using namespace std;

#define V(X) #X << "=" << (X) <<" "
struct Fixture{
	
	Fixture(){
	}
	~Fixture(){
	}
};

typedef ublas::vector<double> Point;

#define DIM 3

struct F{
	typedef double result_type;
	inline result_type operator()(const Point& p, int i)const {return p[i];}
};


BOOST_FIXTURE_TEST_SUITE( suite, Fixture )

BOOST_AUTO_TEST_CASE( testTrans )
{
	ICP<DIM,Point,F> icp;
	vector<Point> v, w;
	const int S=9;
	for(int i=0;i<S;i++){
		Point p(DIM), q(DIM);
		p[0] = i; p[1] = i*i; p[2] = 3*i;
		v.push_back(p);
		w.push_back(p);
	}
	for(int i=0;i<S; i++){
		v[i] += ublas::scalar_vector<double>(3,0.1);
		//swap(v[i][1],v[i][2]);
		v[i][0] += 0.01*drand48();
		v[i][1] += 0.01*drand48();
		v[i][2] += 0.01*drand48();
	}
	icp.registerModel(v.begin(), v.end());
	icp.match(w.begin(),w.end());
	ICP<DIM,Point,F>::TQuat q = icp.getRot();
	BOOST_CHECK_CLOSE(q.R_component_1(),1.0f,1.0);
	BOOST_CHECK_CLOSE(q.R_component_2(),0.0f,1.0);
	BOOST_CHECK_CLOSE(q.R_component_3(),0.0f,1.0);
	BOOST_CHECK_CLOSE(q.R_component_4(),0.0f,1.0);
}
BOOST_AUTO_TEST_CASE( testRot )
{
	ICP<DIM,Point,F> icp;
	vector<Point> model, query_orig, query;
	const int S=10;
	for(int i=0;i<S;i++){
		Point p(DIM), q(DIM);
		p[0] = i; p[1] = 2*i; p[2] = 3*i;
		model.push_back(p);
	}
	copy(model.begin(),model.end(),back_inserter(query_orig));
	copy(model.begin(),model.end(),back_inserter(query));

	// modify model a bit
	for(int i=0;i<S; i++) {
		model[i] += ublas::scalar_vector<double>(3,0.1);
		swap(model[i][1],model[i][2]);
	}

	// register model and match query
	icp.registerModel(model.begin(), model.end());
	icp.match(query.begin(),query.end());

	cout << "trans: " <<icp.getTrans() <<endl;
	cout << "rotmat:" << quaternion_to_R3_rotation(icp.getRot()) <<endl;
	cout << "rot:   " << icp.getRot() <<endl;

	// test transformation properties
	ICP<DIM,Point,F>::TQuat q = icp.getRot();
	ICP<DIM,Point,F>::TVec  v = icp.getTrans();
	cout << V(v)<<endl;
	for(int i=0; i<S; i++){
		Point x = query_orig[i];
		Point y = model[i];
		ICP<DIM,Point,F>::TQuat qx(0,x[0],x[1],x[2]);
		ICP<DIM,Point,F>::TQuat qy(0,y[0],y[1],y[2]);
		ICP<DIM,Point,F>::TQuat qv(0,v[0],v[1],v[2]);
		ICP<DIM,Point,F>::TQuat qxtrans = q * qx * boost::math::conj(q) + qv;
		BOOST_CHECK_CLOSE(qxtrans.R_component_2(),qy.R_component_2(),1.0);
		BOOST_CHECK_CLOSE(qxtrans.R_component_3(),qy.R_component_3(),1.0);
		BOOST_CHECK_CLOSE(qxtrans.R_component_4(),qy.R_component_4(),1.0);
	}

}







BOOST_AUTO_TEST_SUITE_END()

