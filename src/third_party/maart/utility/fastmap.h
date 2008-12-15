/**
 * @file fastmap.h
 *
 * Implements the FastMap data dimension reduction algorithm.
 *
 * Copyright 2003-2006 Steven Blackburn, http://www.beeka.org
 *
 * This file is part of the Music and Audio Retreival Tools (MaART).
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef maart_fastmap_h
#define maart_fastmap_h

#include <cassert>
#include <vector>

#include "utility/distance.h"


namespace utility {


/**
 * The class which implements the FastMap algorithm.
 *
 * @param T The object type. For example, vector<float>, if each object is a
 *          feature vector of floating-point numbers.
 *
 * Example:
 *   fastmap< vector<float> > fm;
 *
 * @note The FastMap algorithm is published in a paper by Christos Faloutsos 
 *       and King-Ip (David) Lin. "FastMap: A Fast Algorithm for Indexing, 
 *       Data-Mining and Visualization of Traditional and Multimedia Datasets", 
 *       SIGMOD Record (ACM Special Interest Group on Management of Data), 
 *       24, 2, issn:0163-5808, pp. 163--174, June, 1995.
 */
template <class T>
class fastmap {
  public:
    typedef T object_t;
    typedef std::vector< object_t > objects_t; // Define what a collection of objects is made of

    typedef float scoreT; // The result of a distance comparison

    typedef void (*progress_fn_t)(size_t total, size_t current);

    typedef std::vector< std::pair<unsigned int, unsigned int> > pivot_array_t;

  private:
    abstract_distance_functor<T, scoreT> const * distance_fn;
    progress_fn_t progress_fn;

    // Stores the ids of the pivot objects - one pair per recursive call [2 x k pivot array PA[]]
    pivot_array_t pivot_array;

    // At the end of the algorithm, the i-th row is the image of the i-th object [N x k array X[]]
    objects_t X;

  protected:
    const float fm_dist(const objects_t &objects, unsigned int a, unsigned int b, unsigned int p);
    void choose_distant_objects(const objects_t &objects, unsigned int *a, unsigned int *b, unsigned int p);
    void do_map(unsigned int k, const objects_t &objects, unsigned int column = 0);

  public:
    fastmap() : distance_fn(NULL), progress_fn(NULL) { };

    void set_distance_function(abstract_distance_functor<T, scoreT> const * fn) { distance_fn = fn; };
    const progress_fn_t set_progress_function(progress_fn_t fn) { progress_fn_t old = progress_fn; progress_fn = fn; return old; };

    void make_map(unsigned int k, const objects_t &objects);
    bool save_pivot_ids(const std::string &filename);
    bool load_pivot_ids(const std::string &filename);

    const std::vector<float> map_object(const object_t &object);

    double evaluate_stress(const objects_t &objects, unsigned int k);

    const objects_t &get_map() const { return X; };
};


/** 
 * Returns the disimilarity between the two objects, including the projection (p) 
 */
template <class T>
const float fastmap<T>::fm_dist(const objects_t &objects,  unsigned int a, unsigned int b, unsigned int p)
{
    if (p == 0)
      return (*distance_fn)(objects[a], objects[b]);
    else
    {
      // p must be > 0, so xa and xb must be in the p-1 column
      const float d=fm_dist(objects, a, b, p-1);
      const float xa=X[a][p-1];
      const float xb=X[b][p-1];
      return sqrt( (d*d) - ((xa-xb)*(xa-xb)) );
    }
}


/**
 * Choose distant objects.
 *  1) Choose arbitrarily an object, and let it be the second pivot object Ob
 *  2) let Oa = (the object that is farthest apart from Ob) (according to the distance function dist ())
 *  3) let Ob = (the object that is farthest apart from 0a)
 *  4) report the objects 0a and ob as the desired pair of objects.
 *
 * @todo Add a exit the iterations if the solution doesn't change
 */
template <class T>
void fastmap<T>::choose_distant_objects(const objects_t &objects, unsigned int *a, unsigned int *b, unsigned int p)
{

#if DEBUG_LEVEL > 0
   	DEBUG_STREAM << "Choosing distant objects for projection " << p << std::endl;
#endif

    /** The number of iterations to find the most distant objects */
    const unsigned int num_iterations = 5;

    //
    // Choose arbitrarily an object, and let it be the second pivot object Ob
    //
    *b = 0; // Start with the first object, to avoid randomness, apart from anything else

    float last_distance = 0.0;  // A note of the distance this iteration has to beat
    for (unsigned int iteration = 0; iteration < num_iterations; iteration++)
    {
      //
      // let Oa = (the object that is farthest apart from Ob) (according to the distance function dist())
      //
      *a = *b;
      float max_distance = 0.0;
      for (unsigned int n=0; n<objects.size(); n++)
      {
        const float distance = fm_dist(objects, *b, n, p);
      	if (distance > max_distance)
      	{
      		*a = n;
      		max_distance = distance;
      	}
      }


      //
      // let Ob = (the object that is farthest apart from 0a)
      //
      *b = *a;
      max_distance = 0.0;
      for (unsigned int n=0; n<objects.size(); n++)
      {
        const float distance = fm_dist(objects, *a, n, p);
      	if (distance > max_distance)
      	{
      		*b = n;
      		max_distance = distance;
      	}
      }


      //
      // Ensure each iteration is increasing the distance, stop if it isn't.
      // If this happens, it is probably due to the same two objects being
      // selected each time round the loop.
      //
      if (max_distance > last_distance)
        last_distance = max_distance;
      else
      {
#if DEBUG_LEVEL > 0
      	DEBUG_STREAM << "Stopping distance object search after iteration " << iteration << std::endl;
#endif
        break;
      }
    }
}


/**
 * FastMap. Initialise the map for k and then call do_map().
 *
 * @param k The desired number of resulting dimensions.
 * @param objects The array of objects.
 */
template <class T>
void fastmap<T>::make_map(unsigned int k, const objects_t &objects)
{
    assert(distance_fn != NULL);

#if DEBUG_LEVEL > 0
    DEBUG_STREAM << "Initialisting the output map (X) for " << objects.size() << " rows" << std::endl;
#endif
    // Initialise the array of mapped features
    X.resize(objects.size());
    for (unsigned int n=0; n<objects.size(); n++)
      X[n].resize(k);

    pivot_array.clear();
    pivot_array.reserve(k);

    do_map(k, objects);

#if DEBUG_LEVEL > 0
    DEBUG_STREAM << "Pivot objects indexes:" << std::endl;
    for (unsigned int n=0; n < pivot_array.size(); n++)
      DEBUG_STREAM << "   " << n << " : " << pivot_array[n].first << " --- " << pivot_array[n].second << std::endl;
#endif
}


/**
 * FastMap recursive function. This function is called by fastmap::make_map.
 *
 * @param k The desired number of resulting dimensions.
 *
 * @param objects The array of objects.
 *
 * @param column Indicates the column of the X[] array currently being updated
 */
template <class T>
void fastmap<T>::do_map(unsigned int k, const objects_t &objects, unsigned int column)
{
#if DEBUG_LEVEL > 0
    DEBUG_STREAM << "fastmap: k = " << k << ", column = " << column << std::endl;
#endif

    // 1) Check for the end of the recursion. (Could probably use k to imply column)
    if (k == 0)
      return;

    // 2) Choose pivot objects
    unsigned int a, b;
    choose_distant_objects(objects, &a, &b, column);

#if DEBUG_LEVEL > 1
    DEBUG_STREAM << "Pivot objects choosen: a = " << a << ", b = " << b << std::endl;
    DEBUG_STREAM << "Distance = " << fm_dist(objects, a, b, column) << std::endl;
#endif

    // 3) Record the ids of the pivot objects
    pivot_array.push_back(std::pair<unsigned int, unsigned int>(a, b));


    // 4)
    if ( fm_dist(objects, a, b, column) == 0.0)
    {
      // set X[ i, col#] =0 for every i and return
      // since all inter-object distances are 0

      for (unsigned int i=0; i<objects.size(); i++)
        for (unsigned int n=0; n<k; n++)
          X[i][column+n] = 0.0;

#if DEBUG_LEVEL > 0
      DEBUG_STREAM << "No further collapsing possible" << std::endl;
#endif

      return;
    }

    // Project objects on line (Oa, Ob)
    const float dab = fm_dist(objects, a, b, column);
    for (unsigned int i=0; i<objects.size(); i++)
    {
      const float dai = fm_dist(objects, a, i, column);
      const float dbi = fm_dist(objects, b, i, column);

      X[i][column] = ((dai*dai) + (dab*dab) - (dbi*dbi)) / (dab * 2.0);
    }

    // 6) Consider the projections of the objects on a hyper-plane
    //    perpendicular to the line (Oa, Ob); the distance function D’()
    //    between two projections is given by Eq. 4
    do_map(k-1, objects, column+1);
}


/**
 * Saves the pivot object ids
 */
template <class T>
bool fastmap<T>::save_pivot_ids(const std::string &filename)
{
    std::ofstream output(filename.c_str());
    if (!output)
      return false;

    for (pivot_array_t::const_iterator p=pivot_array.begin(); p!=pivot_array.end(); p++)
      output << p->first << ' ' << p->second << std::endl;

    return true;
}


/**
 * Load the pivot object ids
 */
template <class T>
bool fastmap<T>::load_pivot_ids(const std::string &filename)
{
    std::ifstream input(filename.c_str());
    if (!input)
      return false;

    pivot_array.clear();

    unsigned int a, b;
    do
    {
    	input >> a >> b;

    	if (input)
        pivot_array.push_back(std::pair<unsigned int, unsigned int>(a, b));
    } while (input);

    return true;
}


/**
 * Calculates the stress introduced by the mapping for a value of k.
 * The 'stress' function gives the relative error that the distances in k-d space
 * suffer from, on average.
 *
 * \verbatim
 *              [ sum_for_each_pair [ (new_distance - old_distance)^2 ] ]
 * stress = sqrt[ ----------------------------------------------------- ]
 *              [         sum_for_each_pair [ (old_distance)^2 ]        ]
 * \endverbatim
 *
 * The new_distance is the Eucildean distance between the FastMapped 'images'.
 * The old_distance is the one given by the distance function.
 *
 * @param objects The collection of vectors to compare with. This *must* be 
 *                the same as the set used to fastmap, as each is compared 
 *                with the mapping that is (currently) held internally.
 *
 * @param k The number of dimensions of the mapped images to consider in the 
 *          stress calculations.
 *                
 * @todo Could be clever about symetric distances and halve the calculations needed.
 */
template <class T>
double fastmap<T>::evaluate_stress(const objects_t &objects, unsigned int k)
{
#if DEBUG_LEVEL > 0
    DEBUG_STREAM << "Calculating stress, which takes " << objects.size() * objects.size() << " comparisons" << std::endl;
#endif

    assert(k <= X[0].size());

    double over = 0.0;
    double under = 0.0;

    for (unsigned int i=0; i<objects.size(); i++)
      for (unsigned int j=0; j<objects.size(); j++)
      {
        // Calculate the euclidean distance between the two images to this point
        float diff = 0.0;
        for (unsigned int n=0; n < k; n++)
          diff += (X[i][n] - X[j][n]) * (X[i][n] - X[j][n]);
        const double d1 = (diff == 0.0) ? 0.0 : sqrt(diff);

        const double d2 = (*distance_fn)(objects[i], objects[j]);

        over += (d1-d2) * (d1-d2);
        under += d1*d1;
      }

#if DEBUG_LEVEL > 1
    DEBUG_STREAM << "Calculating stress: over = " << over << ", under = " << under << std::endl;
#endif

    return sqrt(over/under);
}


/**
 * Map an object into k dimensions using the pivot objects from the current map
 */
template <class T>
const std::vector<float> fastmap<T>::map_object(const object_t &object)
{
    std::vector<float> mapped;
    return mapped;
}


/**
 * Map a query into the n-dimensions using a previous FastMap analysis.
 *
 * @todo Might want to make this note if this query would have changed
 *       the mapping
 *
void project_query(const string &queryfile, const objects_t &objects, const pivot_array_t &pivot_array, const string &querymap)
{
    const unsigned int k = pivot_array.size();

    float low, high;

    objects_t query;
    if (!load(queryfile, &query, &low, &high, false))
    {
      cerr << "Unable to open text file " << queryfile << endl;
//      return 1;
    }

    objects_t queryX;
    // Initialise the array of mapped features
    queryX.resize(query.size());
    for (unsigned int n=0; n<query.size(); n++)
      queryX[n].resize(k);

    for (unsigned int n=0; n<k; n++)
    {
    	const unsigned int a = pivot_array[n].first;
    	const unsigned int b = pivot_array[n].second;

      const float dab = fm_dist(objects, dist, a, b, n);
      for (unsigned int i=0; i<query.size(); i++)
      {
        // Project objects on line (Oa, Ob)
//        const float dai = fm_dist(objects, dist, a, i, n);  // TODO - this needs to use queryX for i
//        const float dbi = fm_dist(objects, dist, b, i, n);  // TODO - this needs to use queryX for i
        const float dai = dist(objects[a], query[i]);
        const float dbi = dist(objects[b], query[i]);

        queryX[i][n] = ((dai*dai) + (dab*dab) - (dbi*dbi)) / (dab * 2.0);
      }
    }

    if (!save(querymap, queryX))
    {
      cerr << "Unable to open text file for writing: " << querymap << endl;
//      return 1;
    }
}
*/


}; //  End of the utility namespace


#endif
