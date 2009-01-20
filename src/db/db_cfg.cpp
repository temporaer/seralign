#include <configuration.hpp>
#include <boost/program_options.hpp>
#include <factory/factory.h>

using namespace boost::program_options;
using namespace std;

class DBCfg{
	public:
		DBCfg();
};

DBCfg::DBCfg(){
	options_description od("========== DB Options ==========");

	options_description gdpdb("  GDistProjectedDB Options");
	gdpdb.add_options()
		("GDistProjectedDB.seriation_gen", value<string>()->default_value("LEVSeriationGen"), "Seriation Generator")
		("GDistProjectedDB.maxiter", value<int>()->default_value(100), "ICP max iterations")
		("GDistProjectedDB.fastmap_tries", value<int>()->default_value(5), "How hard to look for longest dist")
		("GDistProjectedDB.lambda", value<float>()->default_value(10), "ICP lambda")
		;
	options_description anndb("  ANNDB Options");
	anndb.add_options()
		("ANNDB.bucket_size", value<int>()->default_value(10), "KD-Tree Bucket Size")
		;
	od.add(gdpdb);
	od.add(anndb);
	gCfg().addModuleOptions(od);
}

namespace {
	DBCfg _dbcfg;
}
