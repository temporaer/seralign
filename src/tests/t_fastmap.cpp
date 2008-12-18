#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE fastmaptest

#include <exception>
#include <vector>
#include <fstream>
#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>
#include <maart/utility/fastmap.h>

using namespace std;
using namespace utility;
struct S{
	vector<float> feat;
	float w, x, y, z;
	S(){}
	S(float e, float f, float g, float h, float fact){
		w = e + fact*drand48();
		x = f + fact*drand48();
		y = g + fact*drand48();
		z = h + fact*drand48();
	}
	unsigned int size()const{return 4;}
	float operator [](unsigned int i)const{
		switch(i){
			case 0:return w;
			case 1:return x;
			case 2:return y;
			case 3:return z;
		}
		throw exception();
	}
};
#define SQR(x) fabs((x)*(x))
class D : public abstract_distance_functor<S, float>{
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
	fastmap<S> fm;
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
	for(int i=0; i<10;i++) v.push_back(S(1,2,2,0, 0.2));
	for(int i=0; i<10;i++) v.push_back(S(2,2,1,0, 0.2));
	for(int i=0; i<10;i++) v.push_back(S(2,1,2,0, 0.2));
	for(int i=0; i<10;i++) v.push_back(S(1,1,2,0, 0.2));
	for(int i=0; i<10;i++) v.push_back(S(1,2,1,0, 0.2));
	for(int i=0; i<10;i++) v.push_back(S(2,1,1,0, 0.2));

	for(int i=0; i<10;i++) v.push_back(S(1,2,2,1, 0.2));
	for(int i=0; i<10;i++) v.push_back(S(2,2,1,1, 0.2));
	for(int i=0; i<10;i++) v.push_back(S(2,1,2,1, 0.2));
	for(int i=0; i<10;i++) v.push_back(S(1,1,2,1, 0.2));
	for(int i=0; i<10;i++) v.push_back(S(1,2,1,1, 0.2));
	for(int i=0; i<10;i++) v.push_back(S(2,1,1,1, 0.2));

	fm.make_map(3, v);

	vector<vector<float> > map = fm.get_map();
	int i=0;
	ofstream os(".dat");
	BOOST_FOREACH(vector<float>& v, map){
		if(i!=0 && (i % 10 == 0))
			os << endl;
		copy(v.begin(),v.end(),ostream_iterator<float>(os,"\t"));
		os << i - (i%10);
		os <<endl;
		i++;
	}
	cout << fm.evaluate_stress(v,2);
}


BOOST_AUTO_TEST_SUITE_END()

