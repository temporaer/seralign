#define BOOST_TEST_DYN_LINK 
#define BOOST_TEST_MODULE sdp_ser_gen_test

#include <boost/shared_ptr.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <sdp_wrapper.hpp>
#include <sdp_prob.hpp>
#include <sdp_seriation_gen.hpp>
#include <sdp_seriation_prob_gen.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <factory/factory.h>
#include <configuration.hpp>
using namespace std;
using namespace boost;

struct Fixture{
	const int n;
	shared_ptr<SDPSeriationGen> sdpsergen;
	shared_ptr<SDPSeriationProbGen::AdjMatT> am;
	
	Fixture()
		:  n(6)
		 , am(new SDPSeriationProbGen::AdjMatT(n,n))
	{
		for(int i=0;i<n;i++)
			for(int j=0;j<n;j++)
			{
				(*am)(i,j) = ((i+j)%2>0)?1.0:1.0;
			}
		sdpsergen.reset(new SDPSeriationGen());
	}
	~Fixture(){
	}
};

void testSerGen(const char* sdpwrapper, SDPSeriationGen* sdpsergen, boost::shared_ptr<AdjMat::AdjMatT> adj){
	auto_ptr<SDPWrapper> wrap = genericFactory<SDPWrapper>::instance().create(sdpwrapper);
	wrap->configure();
	sdpsergen->setSDPWrapper(wrap);
	ProbAdjPerm pap;
	pap.setAdjMat(adj);
	(*sdpsergen)(pap);
}

BOOST_FIXTURE_TEST_SUITE( suite, Fixture )

#ifdef HAVE_SDPA
BOOST_AUTO_TEST_CASE( testSDPA )
{
	gCfg().parsecfg(0,NULL);
	testSerGen("SDPAWrapper",sdpsergen.get(),am);
}
#endif
#ifdef HAVE_DSDP
BOOST_AUTO_TEST_CASE( testDSDP )
{
	gCfg().parsecfg(0,NULL);
	testSerGen("DSDPWrapper",sdpsergen.get(),am);
}
#endif
#ifdef xHAVE_SDPLR
BOOST_AUTO_TEST_CASE( testSDPLR )
{
	gCfg().parsecfg(0,NULL);
	testSerGen("SDPLRWrapper",sdpsergen.get(),am);
}
#endif


BOOST_AUTO_TEST_SUITE_END()

