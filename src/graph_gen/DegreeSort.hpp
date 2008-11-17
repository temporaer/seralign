#ifndef __DEGREESORT_HPP__
#define __DEGREESORT_HPP__

/**
 * @brief DegreeSort sort members of adjacency matrix by degree
 * @author Hannes Schulz <mail at hannes-schulz dot de>
 * @date 2008-11-17
 */


#include	<ProbData.hpp>

class DegreeSort
{
  public:

    /**
     * Default constructor
     */
    DegreeSort();

	/**
	 * Sort a the adjacency matrix of a pap by degree.
	 * Records the changes in associated Permutation matrix.
	 */
	void sort(ProbAdjPerm& pap);

    /**
     * Destructor
     */
    virtual ~DegreeSort();

  private:

};

#endif /* #ifndef __DEGREESORT_HPP__ */
