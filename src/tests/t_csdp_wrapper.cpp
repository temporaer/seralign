#define BOOST_TEST_DYN_LINK 
#define BOOST_TEST_MODULE csdp_wrapper_test

#include <iostream>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <factory/factory.h>
#include <sdp_wrapper.hpp>
#include <sdp_prob.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <memory>
using namespace std;
using namespace boost::numeric;

struct Fixture{
	typedef ublas::matrix<double> Mat;
	typedef ublas::vector<double> Vec;
	Mat mC;
	Mat mA1;
	Mat mA2;
	Vec mb;
	SDPWrapper::AnswerT mAnswer;
	SDPProb mProb;
	
	Fixture()
	{
		double C[7][7] = {
			{2,1,0,0,0,0,0},
			{1,2,0,0,0,0,0},
			{0,0,3,0,1,0,0},
			{0,0,0,2,0,0,0},
			{0,0,1,0,3,0,0},
			{0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0}};
		double A1[7][7] = {
			{3,1,0,0,0,0,0},
			{1,3,0,0,0,0,0},
			{0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0},
			{0,0,0,0,0,1,0},
			{0,0,0,0,0,0,0}};
		double A2[7][7] = {
			{0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0},
			{0,0,3,0,1,0,0},
			{0,0,0,4,0,0,0},
			{0,0,1,0,5,0,0},
			{0,0,0,0,0,0,0},
			{0,0,0,0,0,0,1}};

		double X[7][7] = {
			{0.125, 0.125, 0, 0, 0, 0, 0},
			{0.125, 0.125, 0, 0, 0, 0, 0},
			{0, 0, 0.667, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0}};


		mC  = Mat(7,7);
		mA1 = Mat(7,7);
		mA2 = Mat(7,7);
		mAnswer = SDPWrapper::AnswerT(7,7);
		for(int i=0;i<7;i++)
			for(int j=0;j<7;j++){
				mC(i,j) = C[i][j];
				mA1(i,j) = A1[i][j];
				mA2(i,j) = A2[i][j];
				mAnswer(i,j) = X[i][j];
			}
		mb = Vec(2);
		mb(0) = 1;
		mb(1) = 2;

		mProb.b = mb;
		mProb.C = mC;
		mProb.F.push_back(mA1);
		mProb.F.push_back(mA2);

	}
	~Fixture(){
	}
};

BOOST_FIXTURE_TEST_SUITE( suite, Fixture )

BOOST_AUTO_TEST_CASE( testSDPLR )
{
	auto_ptr<SDPWrapper> wrap = genericFactory<SDPWrapper>::instance().create("CSDPWrapper");
	SDPWrapper* solv = (SDPWrapper*) wrap.get();
	SDPWrapper::AnswerT ret = (*solv)(mProb);
	BOOST_REQUIRE_EQUAL(ret.size1(),mAnswer.size1());
	BOOST_REQUIRE_EQUAL(ret.size2(),mAnswer.size2());
	for(int i=0;i<ret.size1(); i++){
		for(int j=i;j<ret.size2(); j++){
			BOOST_CHECK_SMALL(ret(i,j)-mAnswer(i,j),0.001);
		}
	}
}


BOOST_AUTO_TEST_SUITE_END()

