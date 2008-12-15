/**
 *  @file kmeans.h
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

#ifndef maart_kmeans_h
#define maart_kmeans_h

#include <cassert>
#include <cfloat>
#include <cmath>
#include <vector>

#include "utility/distance.h"
#include "utility/matrix.h"

#include <iostream>

namespace utility {

/**
 * Define the signature of the progress reporting function.
 * The reporting function must not expect the total to stay the
 * same. Although in many situations it will remain constant,
 * some functionality might dynamically calculate the total.
 * This could cause the
 *
 * @param iteration The iteration this progress comes from
 *
 * @param points The total number of data points
 *
 * @param stable The current number of data points which are stable (i.e. haven't moved this iteration)
 */
typedef void (*cluster_progress_fn_t)(size_t iteration, size_t points, size_t stable);


void display_cluster_progress(size_t iteration, size_t points, size_t stable);


/**
 * Implements the k-means clustering algorithm.
 *
 * @todo Implement a different starting algorithm, which finds disperate initial clusters.
 * @todo Change the datapoints to be a vector, rather than a matrix (saves spending time calling to_vector).
 * @todo Once changed to input vectors, consider using iterators, so any container could be used.
 */
template <class featureT, class scoreT>
void kmeans(std::vector< std::vector<featureT> > *cluster_centers, std::vector<size_t> *cluster_assignment, const matrix<featureT> &records, size_t k,
    const abstract_distance_functor<std::vector<featureT>, scoreT> &distance_fn,
    cluster_progress_fn_t progress_fn = NULL)
{
    const size_t R = records.rows();
    const size_t d = records.columns();

    // Read something saying that the average case converges in 'much less 
    // than the number of data points', so using double the number of points 
    // as a limit on the iterations
    const size_t max_iterations = R * 2;

    assert(k > 0); // Check we have been asked to calculate something

    // Clear the cluster assignment record (setting to k => not assigned)
    cluster_assignment->resize(R);
    for (size_t n=0; n<R; ++n)
      (*cluster_assignment)[n] = k;

    //
    // Determine random starting points
    //
    std::vector< std::vector<featureT> > centers(k);
    std::vector<size_t> center_populations(k);
    for (unsigned int c=0; c < k; ++c)
    {
      size_t record;

      do {
        double r = (double)rand() / (double)(RAND_MAX);
        record = (size_t) (r * R);
      } while ((*cluster_assignment)[record] != k);

      centers[c] = to_vector(records.row(record));
      (*cluster_assignment)[record] = c;
      center_populations[c] = 1;
    }


    //
    // Now iterate:
    //
    unsigned long iteration = 0;
    size_t assignment_changed;

    do
    {
      assignment_changed = 0;
      iteration++;

      // Assign each row of the matrix to the cluster they are closest to
      for (size_t n=0; n<R; ++n)
      {
        // Find closest cluster
        size_t best_cluster = 0;
        scoreT best_score = distance_fn(to_vector(records.row(n)), centers[0]);
        for (unsigned int c=1; c < k; ++c)
        {
          const scoreT this_score = distance_fn(to_vector(records.row(n)), centers[c]);
          if (this_score < best_score)
          {
            best_score = this_score;
            best_cluster = c;
          }
        }

        // Assign to closest cluster
        const size_t current_cluster = (*cluster_assignment)[n];
        if (current_cluster != best_cluster)
        {
          if (current_cluster == k || center_populations[current_cluster] > 1)  // Ensure the cluster is never empty...
          {
            if (current_cluster != k)
              center_populations[current_cluster]--;
            (*cluster_assignment)[n] = best_cluster;
            center_populations[best_cluster]++;
            assignment_changed++;
//std::cout << n << ": switching from " << current_cluster << " to " << best_cluster << std::endl;
//std::cout << '.';
          }
        }

      }

      // Cope with instances of empty clusters.
      // Currently the code prevents empty clusters by preventing the last point leaving it
      // but the very next row might be added.
      // Solutions include:
      //   * picking a point furthest away from all the non-empty centres (after calculating new centroids)
      //   * pick a point at random
      //   * note the point which was last closest to the centroid (assuming centroids don't move too much)
//std::cout << std::endl;  
      // Calculate the centroid of the cluster (catching the special case of an empty cluster)
      for (unsigned int c=0; c < k; ++c)
      {
        // Initialise variables for mean
        for (unsigned int i=0; i<d; ++i)
          centers[c][i] = 0.0;
  
        for (size_t n=0; n<R; ++n)
        {
          const std::vector<featureT> &record = to_vector(records.row(n));
          for (unsigned int i=0; i<d; ++i)
            if ((*cluster_assignment)[n] == c)
              centers[c][i] += record[i];
        }

        for (unsigned int i=0; i<d; ++i)
          centers[c][i] /= center_populations[c];

//std::cout << "Population of cluster " << c << " is " << center_populations[c];
//        std::cout << ", mean [0] is " << centers[c][0] << std::endl;
      }

      if (progress_fn != NULL)
        progress_fn(iteration, R, R - assignment_changed);

      // Stop if no reassignments take place (or a limit on the iterations is reached)
    } while (assignment_changed > 0 && iteration < max_iterations);

    
    (*cluster_centers) = centers;
}


/**
 * Calculates the Residual Sum of Squares of the cluster assignments.
 * This gives a measure of how well the clusters represent the dataset.
 */
template <class featureT, class scoreT>
scoreT calculate_rss(const std::vector<std::vector<featureT> > &cluster_centers, 
    const std::vector<size_t> &cluster_assignment, const matrix<featureT> &records,
    const abstract_distance_functor<std::vector<featureT>, scoreT> &distance_fn)
{
    scoreT distortion = 0;

    const size_t R = records.rows();
    for (size_t n=0; n<R; ++n)   
    {
      const size_t this_cluster = cluster_assignment[n];
      const scoreT this_score = distance_fn(to_vector(records.row(n)), cluster_centers[this_cluster]);
      distortion += this_score * this_score;
    }
    
    return distortion;
}


/**
 * Determine the best number for k, when using k-means cluster, using the Schwarz criterion.
 *
 * @param max_k The maximum value to try for k
 *
 * @return The best value for k.
 */
template <class featureT, class scoreT>
size_t determine_best_kmeans_by_sc(std::vector< std::vector<featureT> > *cluster_centers, 
                                   std::vector<size_t> *cluster_assignment, 
                                   const matrix<featureT> &records, 
                                   size_t max_k, 
                                   const abstract_distance_functor<std::vector<featureT>, scoreT> &distance_fn,
                                   cluster_progress_fn_t progress_fn = NULL)
{
    const size_t R = records.rows();
    const size_t d = records.columns();

    size_t best_k = 0;
    double best_rss = DBL_MAX;
    double best_sic = DBL_MAX;
    double best_sc = DBL_MAX;

    std::vector< std::vector<featureT> > best_cluster_centers;
    std::vector<size_t> best_cluster_assignment;

    bool found_best = false;
    for (size_t this_k=2; !found_best && this_k<=max_k; this_k++)
    {
      kmeans(cluster_centers, cluster_assignment, records, this_k, distance_fn, progress_fn);

      const double rss = calculate_rss(*cluster_centers, *cluster_assignment, records, distance_fn);
      const double sic = std::log((double)R) * d * this_k  +  R * std::log(rss / R); // From http://www.answers.com/topic/schwarz-criterion
      const double sc = rss + (std::log((double)R) * d * this_k); // From Andrew Moore's slides

      std::cout << "k = " << this_k << "\trss = " << rss << "\tsic = " << sic << "\tsc = " << sc << std::endl;

//      if (sc < best_sc)
      if (sic < best_sic)
      {
        best_k = this_k;
        best_rss = rss;
        best_sic = sic;
        best_sc = sc;

        // Probably want to keep the best cluster centers and assignments
        best_cluster_centers = (*cluster_centers);
        best_cluster_assignment = (*cluster_assignment);
      }
      else
        // we have started increasing, so must have found the best
        found_best = true;
    }

    // Copy the best values back..
    (*cluster_centers) = best_cluster_centers;
    (*cluster_assignment) = best_cluster_assignment;

    return best_k;
}


}; //  End of the utility namespace


#endif
