/* * Copyright (C) 2003  Columbia Genome Center * All Rights Reserved. *
* Mutual_Op.cc -- Implementation of the different operations. *
* $Id: Matrix_Op.cpp,v 1.3 2005/10/04 18:24:12 kw2110 Exp $ * */

#include <iostream>
#include <sstream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "Matrix_Op.h"

#ifdef __BCPLUSPLUS__   // For Borland Compiler
  #include <stlport/hash_map>
using namespace std;
#else
  #ifdef __GNUC__        // For GNU gcc compiler
    #if __GNUC__ < 3
      #include <hash_map.h>
      #define HASH_MAP hash_map
using namespace std;

    #elif __GNUC__ >= 4 && __GNUC_MINOR__ >= 3
      #include <unordered_map>
      #define HASH_MAP unordered_map

    #else
      #include <ext/hash_map>
      #define HASH_MAP hash_map

      #if __GNUC_MINOR__ == 0
using namespace std; // GCC 3.0
      #else
using namespace __gnu_cxx; // GCC 3.1 and later
      #endif
    #endif
  #endif
#endif


namespace aracne {

  using namespace std;

  typedef HASH_MAP< int, Node >::const_iterator const_line_iterator;

  double const MatrixOp::NO_CONNECTION = 0.0;

  /* Create the edge matrix
   * data       - microarray dataset
   * mutual     - mutual information calculator
   * matrix     - the edge matrix
   * threshld   - mutual informaiton threshold
   * controlId  - conditional gene id; when controlId = -1, no contraint (i.e. use
   *              all arrays to compute the mutual information)
   * ids        - vector of gene id for which edge matrix will be computed; if ids
   *              is empty, all genes will be computed
   * v          - a pointer to the vector which contains the array ids used for
   *              mutual information computation
   */
  void MatrixOp::createEdgeMatrix( Microarray_Set & data, Mutual_Info & mutual, Matrix & matrix,
      double threshold, int controlId, double noise2, vector < int > & ids, vector < int > * v )
  {
    time_t t1, t2;
    // Get a timestamp to assess how much time this will take
    time( & t1 );

    int markerNo = data.Get_Num_Markers();
    if ( ids.empty() )
    {
      size_t step = size_t(ceil(0.1*markerNo));
      createEntries( matrix, markerNo );
      // Loop over all the probes
      for ( int i = 0; i < markerNo; i++ ) {
        if ( ( i != controlId ) && ( isMarkerEnabled( data, i ) ) ) {
          // Process only markers that are active (i.e. ones that are not filtered out)
          computeOneRow( data, mutual, matrix, threshold, i, markerNo, controlId, v, true, true, noise2 );
        }
        if ( ( (i+1)%step ) == 0 ) {
          time( & t2 );
          cout << 10*(i+1)/step << "%, time: " << difftime( t2, t1 ) << endl << flush;
        }
      }
      time( & t2 );
      std::cout << "Gene: " << markerNo << " Time: " << difftime( t2, t1 ) << std::endl;
    }
    else
    {
      size_t step = size_t(ceil(0.1*ids.size()));
      for ( size_t i = 0; i < ids.size(); i++ ) {
        if ( ids[i] == controlId ) {
          continue;
        }

        while ( matrix.size() <= ( unsigned int ) ids[i] ) {
          matrix.push_back( HASH_MAP < int, Node > () );
        }
        computeOneRow( data, mutual, matrix, threshold, ids[i], markerNo, controlId, v, false, false, noise2 );
        if ( ( (i+1)%step ) == 0 ) {
          time( & t2 );
          cout << 10*(i+1)/step << "%, time: " << difftime( t2, t1 ) << endl << flush;
        }
      }
      time( & t2 );
      std::cout << "Gene: " << ids.size() << " Time: " << difftime( t2, t1 ) << std::endl;
    }
  }


  /* this function computes one row of the adjacency matrix. it is called by
   * the function createEdgeMatrix. Note that since the adjacency matrix is symmetric,
   * only upper right triangle is computed; so here this function computes only the * "upper right triangle" of the row.
   */
  void MatrixOp::computeOneRow( Microarray_Set & data, Mutual_Info & mutual, Matrix & matrix, double threshold, int row_idx,
      int markerNo, int controlId, vector < int > * v, bool half_matrix, bool symmetric, double noise2 )
  {
    int j = 0;
    if (half_matrix) {
      j = row_idx + 1;
    }

    for ( ; j < markerNo; j++ ) {
      if ( ( j != controlId ) && ( isMarkerEnabled( data, j ) ) ) {
        double edge = calculateMI( data, mutual, row_idx, j, threshold, noise2, v );
        if ( edge != MatrixOp::NO_CONNECTION ) {
          addNode( matrix, row_idx, j, edge, symmetric );
        }
      }
    }
  }


  /* removes non contiguous edges by using the triangle inequality
   * data       - microarray dataset
   * mutual     - mutual information calculator
   * matrix     - the edge matrix
   * threshld   - mutual informaiton threshold
   * epsilon    - DPI tolerance
   * ids        - vector of gene id for which edge matrix will be computed; if ids
   *              is empty, all genes will be computed
   * v          - a pointer to the vector which contains the array ids used for
   *              mutual information computation
   * constrained- a boolean to determine whether new MI need to be computed. For '-l'
   *              option, we don't want to compute any new MI other than those already
   *              computed for probes provided. For '-h' option, we first compute
   *              only the row in the adjacency matrix corresponding to the probe of
   *              interest; if DPI need to be performed, we will further compute more
   *              MIs as triplets can be formed.
   */
  void MatrixOp::reduce( Microarray_Set & data, Matrix & matrix, Mutual_Info & mutual, double epsilon, double threshold,
      double noise2, vector < int > & ids, vector < int > * v , map < size_t, size_t > & transfac )
  {
    time_t t1; // start time
    time_t t2; // end time
    time( & t1 ); // get starting time

    if ( ids.empty() )
    {
      // apply DPI on all the nodes
      size_t markerNo = matrix.size();
      for ( size_t i = 0; i < markerNo; i++ ) {
        reduceOneNode( i, data, matrix, mutual, epsilon, threshold, noise2, v, transfac );
      }
    }
    else {
      for ( size_t i = 0; i < ids.size(); i++ ) {
        reduceOneNode( ids[i], data, matrix, mutual, epsilon, threshold, noise2, v, transfac );
      }
    }
    time( & t2 );
    std::cout << "DPI running time is: " << difftime( t2, t1 ) << std::endl;
  }


  void MatrixOp::reduceOneNode( unsigned int row_idx, Microarray_Set & data, Matrix & matrix,
      Mutual_Info & mutual, double epsilon, double threshold, double noise2, vector < int > * v, map <size_t, size_t> & transfac )
  {
    SortDecreasing_ArrayValuePair sorter;
    HASH_MAP < int, Node > & m = matrix[row_idx];
    vector < ArrayValuePair > miVector;
    // put all the values of the i-th row in the matrix on a vector for sorting
    // the ArrayValuePair structure is borrowed from the Microarray_Set class,
    // it should actually be called EdgeMIPair here. They share the same
    // structure consisting of two fields: a int and a double. The assoiciated
    // sorting class will sort by the double fields.
    for ( HASH_MAP < int, Node >::const_iterator it = m.begin(); it != m.end(); it++ ) {
      miVector.push_back( ArrayValuePair( it->first, it->second.Get_Mutualinfo() ) );
    }
    // sort the vector
    sort( miVector.begin(), miVector.end(), sorter );

    // for each value in the vector (these are genes directly connected to row_idx.
    // i.e. geneId1 is connected to i)
    for ( size_t j = 0; j < miVector.size(); j++ ) {
      // (geneId1, valueAB) are the (id, value) of the selected gene, starting with the largest mi
      size_t geneId1 = miVector[j].getId();
      double valueAB = miVector[j].getValue();

      // Set the limits
      double minMI = valueAB / ( 1.0 - epsilon );
      // Loop on all the genes that have an equal or higher mutual information
      // (these genes are also connected to i) bool connected = true;
      for ( size_t k = 0; k < j; k++ ) {
        size_t geneId2 = miVector[k].getId();

        double valueAC = miVector[k].getValue();
        // We know that valueAB is smaller than valueAC by default; then tolerance is always positive
        if ( valueAC <= minMI ) {
          // we can exit from this loop because all future values are smaller
          break;
        }
        // We know that valueAB is smaller than valueAC by default.
        // Therefore, if it is also smaller than valueBC, then
        // geneId1 is connected to gene i via geneId2
        double valueBC = getNodeMI( matrix, geneId1, geneId2 );
        if ( valueBC == -1 ) {
          // MI need to be computed as necessary
          valueBC = calculateMI( data, mutual, geneId1, geneId2, threshold, noise2, v );
        }

        if ( valueBC > minMI ) {
          // if TF annotation information is provided, the triangle will be
          // broken only if certein logic is met
          if (not(transfac.empty())) {
            if ( protectedByTFLogic(transfac, row_idx, geneId1, geneId2) ) {
              continue;
            }
          }

          Node & node = m[geneId1];
          node.Set_Intermediate( geneId2 );
          break;
        }
      }
    }
  }


  /* compute mutual information between two gene expression vectors */
  double MatrixOp::calculateMI( Microarray_Set & data, Mutual_Info & mutual, int probeId1, int probeId2,
      double threshold, double noise2, vector < int > * v )
  {
    double edge = MatrixOp::NO_CONNECTION;
    if ( isSameGene( data, probeId1, probeId2 ) ) {
      return edge;
    }

    Pair_Vector pairs;
    fillPair( mutual.Get_Microarray_Num(), data, probeId1, probeId2, pairs, v );

    double mi = mutual.Compute_Pairwise_MI( pairs );
    if (mi >= threshold) {
      edge = mi;
      // correct for bias due to array measurement noise
      if (noise2 > 0) {
        double v1 = data.getMarkerVariance(probeId1);
        double v2 = data.getMarkerVariance(probeId2);
        double lambda = (v1/(v1-noise2))*(v2/(v2-noise2));
        edge = mi + 0.5*log( 1 + (exp(2*mi)-1)*(1-1/lambda) );
      }
    }
    return edge;
  }


  /* fill the pair vector with gene expression values for MI computation. For the
   * unconditional MI calculation (array = NULL), expresssion values from all
   * arrays are used; for the conditional case (array not NULL), only microarrays * in the array vector are used.
   *
   */
  void MatrixOp::fillPair( unsigned maNum, const Microarray_Set & matrix, int ID1, int ID2, Pair_Vector & pv,
      const std::vector < int > * arrays )
  {
    if ( arrays )
    {
      for ( unsigned int i = 0; i < maNum; i++ )
      {
        Gene_Pair pair;
        Set_Pair( pair, matrix.Get_Value( ( * arrays ) [i], ID1 ), matrix.Get_Value( ( * arrays ) [i], ID2 ), i );
        pv.push_back( pair );
      }
    }
    else
    {
      for ( unsigned int i = 0; i < maNum; i++ )
      {
        Gene_Pair pair;
        Set_Pair( pair, matrix.Get_Value( i, ID1 ), matrix.Get_Value( i, ID2 ), i );
        pv.push_back( pair );
      }
    }
  }


  inline void MatrixOp::Set_Pair( Gene_Pair & p, double xId, double yId, int maId )
  {
    p.Set_X( xId );
    p.Set_Y( yId );
    p.Set_MaID( maId );
  }


  inline bool MatrixOp::isMarkerEnabled( const Microarray_Set & data, int i )
  {
    return ( data.Get_Marker( i ).Enabled() );
  }


  inline bool MatrixOp::isSameGene( const Microarray_Set & data, int i, int j)
  {
    return ( (i == j) ||
             ( ( data.Get_Marker(i).Get_Label() == data.Get_Marker(j).Get_Label() ) &&
               ( data.Get_Marker(i).Get_Label() != "---" )
             )
           );
  }


  inline void MatrixOp::createEntries( Matrix & matrix, unsigned entry )
  {
    for ( unsigned i = 0; i < entry; i++ ) {
      matrix.push_back( HASH_MAP < int, Node > () );
    }
  }


  inline void MatrixOp::addNode( Matrix & matrix, int i, int j, double edgeValue, bool symmetric )
  {
    Node node( 0.0 );
    node.Set_Mutualinfo( edgeValue );
    matrix[i] [j] = node;
    if (symmetric) {
      matrix[j] [i] = node;
    }
  }


  /* the return of this function can be 3 cases:
   * case I: a postive real number - the node has been computed and survived the
             thresholding, so the MI value is returned
     case II: 0 - the node has been computed but did not survive the thresholding
     case III: -1 - the node has not yet been computed
   */
  double MatrixOp::getNodeMI(Matrix & matrix, unsigned int geneId1, unsigned int geneId2)
  {
    double mi = -1;
    if ( matrix.size() > geneId1 ) {
      HASH_MAP < int, Node > & m12 = matrix[geneId1];
      if ( not(m12.empty()) ) {
        // m12 can be empty for two reasons:
        // 1) it was simply pushed onto the hash as it is on the way for other rows
        //    to be computed
        // 2) it was actually computed, but not a single MI estimation pass the threshold.
        //    this case is possible but less likely to happen
        if ( m12.count( geneId2 ) > 0 ) {
          Node & n12 = m12[geneId2];
          mi = n12.Get_Mutualinfo();
        }
        else {
          // now that there are other entries on this row that have been computed,
          // this entry being absent is due to the fact that it is less than the threshold.
          // in this case, it obviously can not break the triangle
          mi = 0;
        }
      }
    }

    // for full network reconstruction, since the computations were symmetric, knowning
    // that n12 does not exist or has been thresholded will be sufficient.
    // however, for cases where we only compute a portion of the full adjacency matrix,
    // we have to check n21 to make sure.
    if ( (mi == -1) && (matrix.size() > geneId2) )
    {
      HASH_MAP < int, Node > & m21 = matrix[geneId2];
      if ( not(m21.empty()) )
      {
        if ( m21.count( geneId1 ) > 0 ) {
          Node & n21 = m21[geneId1];
          mi = n21.Get_Mutualinfo();
        }
        else {
          mi = 0;
        }
      }
    }
    return mi;
  }


  bool MatrixOp::protectedByTFLogic( map <size_t, size_t> & transfac, size_t geneId1, size_t geneId2, size_t geneId3 )
  {
    bool isA = ( transfac.find(geneId1) != transfac.end() );
    bool isB = ( transfac.find(geneId2) != transfac.end() );
    bool isC = ( transfac.find(geneId3) != transfac.end() );

    if ( (isA && isB) || not(isA || isB) ) {
      return false;
    }

    if ( isC ) {
      return false;
    }
    return true;
  }

} // aracne
