#ifndef __SDF_ADJMAT_HPP__
#define __SDF_ADJMAT_HPP__

#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <adjmat_gen.hpp>
#include <vector>
#include <map>

class SDFAdjmatGen: public AdjMatGen{
	private:
		int          mFixedSize;
	public:
		virtual void configure();
		SDFAdjmatGen();
		virtual ProbAdjPerm operator()();
		virtual ~SDFAdjmatGen();
		virtual bool hasNext();
		virtual std::string getPlainDescription(int ser_idx, const Serialization&, const std::string&);
		virtual std::string getGraphID(const std::string&);
		virtual std::string getGraphVizNodeAttribs(int idx, const std::string&); //< should start with a comma!
		virtual int getClassID(const std::string&);          ///< which class the graph belongs to
		
	private:

		struct FileDescriptor{
			std::string name;
			int         classid;
		};
		std::vector<FileDescriptor> mInputFiles; 

		struct Descriptor{
			boost::shared_ptr<AdjMat::AdjMatT> mA_ptr;
			int                                mClassID;
		};
		std::map<std::string, Descriptor> mDescriptors;

		void setInputFiles(const std::string& s);
		void readInputFiles();
		void readMolekule(std::istream& is);
};

#endif /* __SDF_ADJMAT_HPP__ */

