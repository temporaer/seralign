/**
 *  @file utility/distance.h
 *
 *  Defines distance functions.
 *
 *  Copyright 2005 Steven Blackburn, http://www.beeka.org
 *
 *  This file is part of the Music and Audio Retreival Tools (MaART).
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
#ifndef maart_distanceH
#define maart_distanceH

#include <iostream>
#include <cmath>
#include <string>
#include <vector>


/**
 * A function to make a vector from a standard C array of features.
 * This is particularly useful for test data.
 */
template <class T>
const std::vector<T> to_vector(const T array[], size_t n)
{
    std::vector<T> v(n);
    for (unsigned int i=0; i<n; i++)
      v[i] = array[i];
}


/**
 * Function to stream a vector.
 */
template <class T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &v)
{
    os << '(';

    if (v.size() > 0)
    {
      os << v[0];
      for (size_t n=1; n<v.size(); n++)
        os << ',' << v[n];
    }
    os << ')';

    return os;
}


//
// Comparison functions
//




/**
 * Implements basic difference function.
 * Returns a normalised difference (mean difference).
 * Zero means that there is no difference. The more positive the number,
 * the more difference there is between the two histograms.
 */
template <class T>
float basic_distance(const std::vector<T> &i, const std::vector<T> &j)
{
    float diff = 0;

    for (size_t n=0; n<i.size(); n++)
      diff += fabs(i[n] - j[n]);

    return diff / i.size();
}

/**
 * Sum of Squares of the differences.
 * This is has most of the benefits of the Euclidean distance but without
 * the expense of calculating the square root. If comparing two distances,
 * rather than iterpreting them then this is a good candidate. The square
 * root of this function can be taken to arrive at the Euclidean distance,
 * if you wish to maintain the same sort of 'scale' as the input vectors.
 *
 * Zero means that there is no difference. The more positive the number,
 * the more difference there is between the two histograms.
 */
template <class T>
float sum_of_squares_distance(const std::vector<T> &i, const std::vector<T> &j)
{
    float diff = 0;

    for (size_t n=0; n<i.size(); n++)
      diff += (i[n] - j[n]) * (i[n] - j[n]);

    return diff;
}


/**
 * Euclidean distance of the differences.
 * In practice, the distinguishing ability of this is about the same as the
 * basic average. It is a more expensive function to computate.
 *
 * Zero means that there is no difference. The more positive the number,
 * the more difference there is between the two histograms.
 */
template <class T>
float euclidean_distance(const std::vector<T> &i, const std::vector<T> &j)
{
    float diff = 0;

    for (size_t n=0; n<i.size(); n++)
      diff += (i[n] - j[n]) * (i[n] - j[n]);

    return sqrt(diff);
}


// TODO: Make this function a template
float cosine_distance(const std::vector<float> &i, const std::vector<float> &j);


//
// Comparison functors
//

/**
 * The first template parameter defines the type of feature.
 *
 * A second template parameter defines the type returned by the distance function (i.e. the score).
 */
template <class featureT, class scoreT = float>
class abstract_distance_functor
{
  public:
    typedef featureT feature_type; ///< The type being compared
    typedef scoreT score_type; ///< The value returned by the comparison


    virtual ~abstract_distance_functor() {};

    virtual score_type operator() (const feature_type &i, const feature_type &j) const = 0;

    /**
     * @return True if the function is symmetric, i.e. fn(i, j) == fn(j, i).
     */
    virtual bool is_symmetric() const = 0;

    /**
     * @return True if the result of the operator() is zero if the two
     *         histograms are the same. This implies that this is a
     *         distance measure.
     */
    virtual bool zero_means_same() const = 0;
};


/**
 * Defines a functor which compares vectors.
 * This is currently intended purely to tidy the code - there might not
 * be a need for this level of abstraction.
 */
template <class featureT, class scoreT = float>
class abstract_vector_distance_functor
  : public abstract_distance_functor<std::vector<featureT>, scoreT>
{};



/// Define the signature of the histogram distance function
typedef float (*distance_fn_t)(const std::vector<float> &i, const std::vector<float> &j);


/**
 * This class can be used to wrap a functor around a standard comparison
 * function.
 */
template <distance_fn_t fn, bool symmetric, bool zero_is_same>
class concrete_distance_functor
  : public abstract_vector_distance_functor<float, float>
{
  public:
    concrete_distance_functor() {}

    virtual score_type operator() (const feature_type &i, const feature_type &j) const
      { return fn(i, j); }

    virtual bool is_symmetric() const { return symmetric; }
    virtual bool zero_means_same() const { return zero_is_same; }
};


const concrete_distance_functor<basic_distance, true, true> basic_distance_functor;
const concrete_distance_functor<sum_of_squares_distance, true, true> sum_of_squares_distance_functor;
const concrete_distance_functor<euclidean_distance, true, true> euclidean_distance_functor;
const concrete_distance_functor<cosine_distance, true, false> cosine_distance_functor;


// Utility function
const abstract_vector_distance_functor<float> * make_distance_functor(const std::string &fn);

#endif // End of sentry
