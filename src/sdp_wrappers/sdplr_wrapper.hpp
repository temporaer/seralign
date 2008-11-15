#ifndef __SDPLR_WRAPPER_HPP__
#define __SDPLR_WRAPPER_HPP__
#include <sdp_wrapper.hpp>
//! Wraps the SDPLR solver to give it an interface specified in SDPWrapper
class SDPLRWrapper: public SDPWrapper
{
    public:
        //! solve an @see SDPProb SDP-problem
        virtual AnswerT operator()(const SDPProb&);
        virtual ~SDPLRWrapper();
		inline void setParamFile(const std::string& s){mParamFile=s;}
	private:
		//! The parameter file
		std::string mParamFile;

        /*! Write an input file for sdpa to read.
         *  \param p   the problem to print
         *  \param in  the file to write to */
        void writeSDPLRInputFile(const SDPProb& p, const char* in);

        /*! Read the output of sdpa.
         *  \param out the name of the sdpa-output file to read
         *  \param ret the vector where we should record the answer */
        bool readSDPLROutputFile(const char* out, AnswerT& ret);

        /*! Run the sdpa-binary.
         *  \param in where the problem was written to
         *  \param out where the solution shall be written to
         *  \param par the parameter file for sdpa (see its doc) */
        void runSDPLR(const char* in, const char* out,const char* par);
};
#endif

