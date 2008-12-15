#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE fastmaptest

#include <vector>
#include <fstream>
#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>
#include <maart/utility/fastmap.h>

using namespace std;
using namespace utility;
struct S{
	vector<float> feat;
	float x, y, z;
	S(){}
	S(float f, float g, float h, float fact){
		x = f + fact*drand48();
		y = g + fact*drand48();
		z = h + fact*drand48();
	}
	operator vector<float> () {
		vector<float> v(3);
		v[0] = x;
		v[1] = y;
		v[2] = z;
		return v;
	}
};
#define SQR(x) fabs((x)*(x))
class D : public abstract_distance_functor<vector<float>, float>{
	virtual bool is_symmetric()const {return true;}
	virtual bool zero_means_same() const {return true;}
	virtual score_type operator() (const feature_type &i, const feature_type &j) const{
		float sum = 0.f;
		for(unsigned int f=0;f<i.size();f++)
			sum += SQR(i[f] - j[f]);
		return sqrt(sum);
	}
};

struct Fixture{
	fastmap<vector<float> > fm;
	D d;
	
	Fixture(){
		fm.set_distance_function(&d);
	}
	~Fixture(){
	}
};


BOOST_FIXTURE_TEST_SUITE( suite, Fixture )

BOOST_AUTO_TEST_CASE( test1 )
{
	vector<S> v;
	for(int i=0; i<10;i++) v.push_back(S(1,2,2, 0.2));
	for(int i=0; i<10;i++) v.push_back(S(2,2,1, 0.2));
	for(int i=0; i<10;i++) v.push_back(S(2,1,2, 0.2));
	for(int i=0; i<10;i++) v.push_back(S(1,1,2, 0.2));
	for(int i=0; i<10;i++) v.push_back(S(1,2,1, 0.2));
	for(int i=0; i<10;i++) v.push_back(S(2,1,1, 0.2));
	vector<vector<float> > map;
	for(int i=0;i<v.size();i++) map.push_back(v[i]);

	fm.make_map(2, map);
	map = fm.get_map();
	int i=0;
	BOOST_FOREACH(vector<float>& v, map){
		if(i!=0 && (i % 10 == 0))
			cout << endl;
		cout << v[0] << "\t" << v[1]<<endl;
		i++;
	}
}


BOOST_AUTO_TEST_SUITE_END()

