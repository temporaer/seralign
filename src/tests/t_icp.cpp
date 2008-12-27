#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE icptest

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
	inline result_type operator()(Point p, int i)const {return p[i];}
};


BOOST_FIXTURE_TEST_SUITE( suite, Fixture )

BOOST_AUTO_TEST_CASE( test1 )
{
	ICP<DIM,Point,F> icp;
	vector<Point> v, w;
	for(int i=0;i<10;i++){
		Point p(DIM);
		p[0] = p[1] = p[2] = i;
		v.push_back(p);
	}
	w=v;
	icp.registerModel(v.begin(), v.end());
	icp.match(w.begin(),w.end());
	cout << "trans: " <<icp.getTrans() <<endl;
	cout << "rot:   " <<icp.getRot() <<endl;
}







BOOST_AUTO_TEST_SUITE_END()

