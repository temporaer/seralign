/**
 *  @file utility/distance.cpp
 *
 *  Implements distance functions.
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

#include <cctype>
#include <cmath>

#include "distance.h"

using namespace std;


/**
 * Utility function to convert a string into a function pointer to a
 * histogram comparison function.
 *
 * @param fn The name of the function.
 *           Can be: "basic"; "euclidean"; "cosine".
 *
 * @return A pointer to the function, or NULL if the name is not
 *         recognised.
 */
const abstract_vector_distance_functor<float> * make_distance_functor(const std::string &fn)
{
    string lower;
    lower.resize(fn.size());
    for (size_t n=0; n<lower.size(); n++)
      lower[n] = tolower(fn[n]);

    if (lower == "basic")
      return &basic_distance_functor;
    else if (lower == "sumofsquares")
      return &sum_of_squares_distance_functor;
    else if (lower == "euclidean")
      return &euclidean_distance_functor;
    else if (lower == "cosine")
      return &cosine_distance_functor;
    else
      return NULL;
}




#include <iostream>
using namespace std;
// Cosine distance.
// Divides the dot product by the magnitudes multiplied together.
// This does not cope well with discrete vectors which share no common
// signals. E.g. A = [0, 64, 0]; B = [0, 0, 64]; The dot product is zero
// and so the answer will be zero, regardless of the maginitude.
// Special case is when one of the vectors contains all zeros, as this
// makes the divisor also zero.
float cosine_distance(const std::vector<float> &i, const std::vector<float> &j)
{
    float dot_product = 0.0;
    float magnitude_i = 0.0;
    float magnitude_j = 0.0;

    for (size_t n=0; n<i.size(); n++)
    {
      dot_product += i[n] * j[n];
      magnitude_i += i[n] * i[n];
      magnitude_j += j[n] * j[n];
    }

    magnitude_i = sqrt(magnitude_i);
    magnitude_j = sqrt(magnitude_j);
//cout << dot_product << " / (" << magnitude_i << " * " << magnitude_j << ")" << endl;
    if (magnitude_i != 0 && magnitude_j != 0)
      return static_cast<float>( 1.0 - fabs(dot_product / (magnitude_i * magnitude_j)) );
    else if (magnitude_i == 0 && magnitude_j == 0)
      return 0.0; // They must be the same if the magnitude of both vectors is zero
    else
      // Would return NAN as one of the vectors is empty.
      // As far as cosine distance is concerned, they couldn't be
      // more different.
      return 1.0;
}

