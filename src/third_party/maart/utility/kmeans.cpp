/**
 *  @file kmeans.cpp
 *
 *  Implements the k-means clustering algorithm.
 *
 *  Copyright 2006 Steven Blackburn, http://www.beeka.org
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

#include <iostream>

#include "utility/kmeans.h"

namespace utility {

using namespace std;


/**
 * Implements a progress reporting function which outputs to cout.
 *
 * @param iteration The iteration this progress comes from
 *
 * @param points The total number of data points
 *
 * @param stable The current number of data points which are stable (i.e. haven't moved this iteration)
 */
void display_cluster_progress(size_t iteration, size_t points, size_t stable)
{
    const float fidelity = 10.0; // 1.0 will update on whole numbers, 10.0 on tenths...
    static float last_displayed = 1000.0;

    const float progress = static_cast<float>((100.0 * stable) / points);

    // Only print out changes greater than 1/fidelity
    if (static_cast<int>((progress - last_displayed) * fidelity) != 0)
    {
//      cout.precision(3);
      cout << "Iteration " << iteration << ": fit quality = " << progress << "%       \r";
      cout.flush();
      last_displayed = progress;
    }
}


}; //  End of the utility namespace
