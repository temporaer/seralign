#include <configuration.hpp>
#include <boost/program_options.hpp>
#include <factory/factory.h>

using namespace boost::program_options;
using namespace std;

class ActionCfg{
	public:
		ActionCfg();
};

ActionCfg::ActionCfg(){
	options_description od("========== Action Options ==========");

	options_description serialize("  Serialize Options");
	serialize.add_options()
		("serialize.want_degree_sorting",    value<bool>()->default_value(true), "whether to sort by degree before calling serializer")
		("serialize.adjmat_gen",    value<string>()->default_value("RandomAdjMatGen"), "Adjacency Matrix Generator")
		("serialize.seriation_gen", value<string>()->default_value("SDPSeriationGen"), "Seriation Generator")
		("serialize.max_num", value<int>()->default_value(0), "Maximum number of graphs to serialize")
		;
	options_description fastmap("  Fastmap Options");
	fastmap.add_options()
		("fastmap.matrices",    value<vector<string> >(), "list of csv files to convert")
		("fastmap.classlabels", value<string>()->default_value("../../data/mutagenesis/188/is_active.txt"), "tsv files with graphname and label on each line")
		("fastmap.filename_replace_what", value<string>()->default_value(".csv"), "what to replace in the filename")
		("fastmap.filename_replace_with", value<string>()->default_value(".map"), "what to replace 'what' with")
		;
	options_description builddb("  BuildDB Options");
	builddb.add_options()
		("BuildDB.knn_k", value<int>()->default_value(3), "k for knn")
		("BuildDB.query_id,p", value<int>()->default_value(0), "ID of object to query")
		("BuildDB.evalmode", value<int>()->default_value(1), "0printDistMat 1spatAnalysis 2knn")
		("BuildDB.embed_meth", value<int>()->default_value(0), "0spectral 1fastmap")
		;
	od.add(serialize);
	od.add(fastmap);
	od.add(builddb);
	gCfg().addModuleOptions(od);
}

namespace {
	ActionCfg _actioncfg;
}
