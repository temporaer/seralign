#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE configtest

#include <boost/test/unit_test.hpp>
#include <configuration.hpp>
#include <boost/program_options.hpp>
using namespace boost::program_options;
using boost::any_cast;

using namespace std;

struct Fixture{
	Configuration cfg;
	
	Fixture(){
	}
	~Fixture(){
	}
};

BOOST_FIXTURE_TEST_SUITE( suite, Fixture )

BOOST_AUTO_TEST_CASE( testVerboseT )
{
	char * argv[] = {"prog","action", "-v" };
	cfg.parsecfg(3,argv);
	bool verbose = any_cast<bool>(cfg.getAny("verbose"));
	BOOST_CHECK_EQUAL( true, verbose );
}

BOOST_AUTO_TEST_CASE( testVerboseF )
{
	char * argv[] = {"prog","action", "-q" };
	cfg.parsecfg(3,argv);
	bool verbose = any_cast<bool>(cfg.getAny("verbose"));
	BOOST_CHECK_EQUAL( false, verbose );
}

BOOST_AUTO_TEST_CASE( testIncompat )
{
	char * argv[] = {"prog","action", "-q","-v" };
	BOOST_CHECK_THROW( cfg.parsecfg(4,argv), logic_error ) ;
}

BOOST_AUTO_TEST_CASE( testDependent )
{
	char * argv[] = {"prog","action"};
	cfg.dependent_options("action","quiet");
	BOOST_CHECK_THROW( cfg.parsecfg(2,argv), logic_error ) ;
}

BOOST_AUTO_TEST_CASE( testAction )
{
	char * argv[] = {"prog","action", "-v" };
	cfg.parsecfg(3,argv);
	string action = any_cast<string>(cfg.getAny("action"));
	BOOST_CHECK_EQUAL( string("action"), action );
}
BOOST_AUTO_TEST_CASE( testAction2 )
{
	char * argv[] = {"prog","action", "-v" };
	cfg.parsecfg(3,argv);
	string action = cfg.getString("action");
	BOOST_CHECK_EQUAL( string("action"), action );
}

class Module{
	public:
		options_description od(){
			options_description o("Module Blubb");
			o.add_options()
				("modopt,m", value<string>(), "module blubb opt")
				;
			return o;
		};
};
BOOST_AUTO_TEST_CASE( testExtension )
{
	char * argv[] = {"prog","action", "-v", "-m", "holla"};
	Module m;
	cfg.addModuleOptions(m.od());
	cfg.parsecfg(5,argv);
}


BOOST_AUTO_TEST_SUITE_END()

