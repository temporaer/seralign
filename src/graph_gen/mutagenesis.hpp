#ifndef __MUTAGENESIS_HPP__
#define __MUTAGENESIS_HPP__

#include <string>
#include <boost/shared_ptr.hpp>
#include <adjmat_gen.hpp>

/**
 * @brief Mutagenesis generate a graph for a molecule in the mutagenesis dataset.
 * @author Hannes Schulz <mail at hannes-schulz dot de>
 * @date 2008-11-17
 */

class Mutagenesis :public AdjMatGen
{
  public:

    /**
     * Default constructor
     */
    Mutagenesis();

    /**
     * Destructor
     */
    virtual ~Mutagenesis();

	virtual void configure();
	virtual ProbAdjPerm operator()();
	virtual bool hasNext();
	virtual std::string getPrologDescription(int ser_idx,const Serialization&,const boost::any&ref);
	virtual std::string getPlainDescription(int ser_idx,const Serialization&,const boost::any&ref);
	virtual std::string getGraphID(const boost::any&ref);
	virtual boost::shared_ptr<AdjMat::AdjMatT> getAdjMat(const boost::any&ref);
	virtual int         getClassID(const boost::any&ref);

	void setInputFilename(std::string&);
	void open();

  private:

  	struct Impl;
	boost::shared_ptr<Impl> mImpl;

};

#endif /* #ifndef __MUTAGENESIS_HPP__ */
