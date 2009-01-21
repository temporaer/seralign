#include <sstream>
#include <fstream>
#include <configuration.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include <boost/foreach.hpp>
#include <factory/factory.h>
#include <matlab_io.hpp>
#include "sdf_adjmat_gen.hpp"
//#include "jumbled_adjmat_gen.hpp"
#include <nana.h>

using namespace std;
using namespace boost;
namespace ublas = boost::numeric::ublas;

void SDFAdjmatGen::configure()
{
	AdjMatGen::configure(); 
	string files = gCfg().getString("SDFAdjmatGen.files");
	mFixedSize   = gCfg().getInt("SDFAdjmatGen.fixed_size");
	setInputFiles(files);
}

SDFAdjmatGen::SDFAdjmatGen()
	:mFixedSize(0)
{
}
bool SDFAdjmatGen::hasNext()
{
	return true;
}
SDFAdjmatGen::~SDFAdjmatGen()
{
}
std::string SDFAdjmatGen::getGraphID(const string& ref)
{
	stringstream str;
	return str.str();
}
std::string SDFAdjmatGen::getGraphVizNodeAttribs(int idx, const string& ref)
{
	I(ref == "");
	return string("");
}
std::string SDFAdjmatGen::getPlainDescription(int ser_idx, const Serialization& ser, const string& ref)
{
	I(ref == "");
	stringstream str;
	return str.str();
}
int SDFAdjmatGen::getClassID(const string& ref)
{
	return 0;
}

ProbAdjPerm SDFAdjmatGen::operator()()
{


	ProbAdjPerm pap;
	//pap.setAdjMat(adj);
	pap.setId(getGraphID(""));

	return pap;
}

void SDFAdjmatGen::setInputFiles(const std::string& inputfileline)
{
	string line = boost::trim_copy(inputfileline);

	vector<string> file_class_pairs;
	boost::split( file_class_pairs, line, boost::is_any_of(";") );

	BOOST_FOREACH(const string& pair, file_class_pairs){
		vector<string> tmp;
		boost::split( tmp, pair, boost::is_any_of(":") );
		if(!tmp.size()==2)
		{
			cerr << "Warning: SDFAdjmatGen: InputFormatError: "<<pair<< " in " << line<< " should have format filename:classid"<<endl;
			continue;
		}
		FileDescriptor fd;
		fd.name = tmp[0];
		fd.classid = boost::lexical_cast<int>(tmp[1]);
		mInputFiles.push_back(fd);
	}
}

void SDFAdjmatGen::readInputFiles()
{
	BOOST_FOREACH(const FileDescriptor& fd, mInputFiles){
		ifstream is(fd.name.c_str());
		readMolekule(is);
	}
}

void SDFAdjmatGen::readMolekule(std::istream& is)
{
	string line, id;
	getline(is,line); // contains filename
	getline(is,line); // contains id of openbabel converter
	getline(is,line); // is empty
	getline(is,line); // space separated: numberAtoms, numberBonds, unknown stuff....
	boost::trim(line);
	vector<string> strvec;
	boost::split( strvec, line, boost::is_any_of("	") );
	if(!strvec.size()==11){
		cerr << "Warning: SDFAdjmatGen: InputFileFormat Error: input description line `"<<line<<"' has unexpected format."<<endl;
		return;
	}

	// create a matrix with numberAtoms size
	int numberAtoms = boost::lexical_cast<int>(strvec[0]);
	int numberBonds = boost::lexical_cast<int>(strvec[1]);
	int n = (mFixedSize==0) ? numberAtoms : mFixedSize;
	shared_ptr<AdjMat::AdjMatT> adj(new AdjMat::AdjMatT(n, n));

	// read all atoms
	for(int i=0;i<numberAtoms; i++){
		getline(is, line);
		boost::trim(line);
		boost::split( strvec, line, boost::is_any_of("	") );
		if(strvec.size() != 9){
			cerr << "Warning: SDFAdjmatGen: InputFileFormat Error: atom description line `"<<line<<"' has unexpected format."<<endl;
			return;
		}
		string atomType = strvec[3]; // C or N or...
	}
	
	// read all bonds
	for(int i=0;i<numberBonds; i++){
		getline(is, line);
		boost::trim(line);
		boost::split( strvec, line, boost::is_any_of("	") );
		if(strvec.size() != 6){
			cerr << "Warning: SDFAdjmatGen: InputFileFormat Error: bond description line `"<<line<<"' has unexpected format."<<endl;
			return;
		}
		int atomFrom = boost::lexical_cast<int>(strvec[0])-1;
		int atomTo   = boost::lexical_cast<int>(strvec[1])-1;
		int bondType = boost::lexical_cast<int>(strvec[2]);
	}

	// read foot
	getline(is, line);
	static const boost::regex e("END");
	if(!regex_match(line, e))
	{
		cerr << "Warning: SDFAdjmatGen: InputFileFormat Error: end-line `"<<line<<"' has unexpected format."<<endl;
		return;
	}

}


namespace{ registerInFactory<AdjMatGen, SDFAdjmatGen> registerBase("SDFAdjmatGen"); }
