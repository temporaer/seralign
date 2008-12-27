#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE icptest

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <icp.hpp>
using namespace util;
using namespace boost::numeric;
using namespace std;

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

BOOST_AUTO_TEST_CASE( test1 )
{
	ICP<DIM,Point,F> icp;
	vector<Point> v, w;
	for(int i=0;i<10;i++){
		Point p(DIM), q(DIM);
		p[0] = i; p[1] = i*i; p[2] = i*i*i*i;
		v.push_back(p);
		q[0] = i; q[1] = i*i; q[2] = i*i*i*i;
		w.push_back(p);
	}
	for(int i=0;i<10; i++){
		//v[i] += ublas::scalar_vector<double>(3,100+0.3*i+0.01*i*i);
		swap(v[i][1],v[i][2]);
		//v[i][2] += 50*drand48();
	}
	icp.registerModel(v.begin(), v.end());
	icp.match(w.begin(),w.end());
	cout << "trans: " <<icp.getTrans() <<endl;
	//cout << "rot:   " << quaternion_to_R3_rotation(icp.getRot()) <<endl;
	cout << "rot:   " << icp.getRot() <<endl;
	ICP<DIM,Point,F>::TQuat q = icp.getRot();
	BOOST_CHECK_CLOSE(q.R_component_1(),0.13f,1.0);
	BOOST_CHECK_CLOSE(q.R_component_2(),0.13f,1.0);
	BOOST_CHECK_CLOSE(q.R_component_3(),-0.69f,1.0);
	BOOST_CHECK_CLOSE(q.R_component_4(),-0.69f,1.0);
}







BOOST_AUTO_TEST_SUITE_END()

