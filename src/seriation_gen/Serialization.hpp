#ifndef __SERIALIZATION_HPP__
#define __SERIALIZATION_HPP__

/**
 * @brief Serialization A Serialization
 * @author Hannes Schulz <mail at hannes-schulz dot de>
 * @date 2008-11-15
 */

#if 1
# include <vector>
 typedef std::vector<int> Serialization;
#else
class Serialization
{
  public:

    /**
     * Default constructor
     */
    Serialization();

    /**
     * Destructor
     */
    virtual ~Serialization();

  private:

};
#endif

#endif /* #ifndef __SERIALIZATION_HPP__ */
