/* * Copyright (C) 2003  Columbia Genome Center * All Rights Reserved. *
* Matrix_Op.h -- Operation that can be done on a matrix. * Currently support only one operation at a time! *
* $Id: Matrix_Op.h,v 1.3 2005/10/04 18:24:12 kw2110 Exp $ */

#ifndef MATRIX_OP_H__
  #define MATRIX_OP_H__
  #include <string>
  #include <vector>
  #include "Microarray_Set.h"
  #include "Matrix.h"
  #include "MutualInfo.h"

namespace aracne {

  /**
   * Class MatrixOp
   *
   * This class implements all operations on an adjacency matrix, e.g. create an
   * adjacency matrix from data, thresholding, DPI, etc.
   */
  class MatrixOp
  {
  private:
    MatrixOp( const MatrixOp & other );
    MatrixOp & operator = ( const MatrixOp & other );

    /* Create entries for all the rows in the matrix */
    inline void createEntries( Matrix & matrix, unsigned entry );
    /* Add a node to the matrix */
    inline void addNode( Matrix & matrix, int i, int j, double edgeValue, bool symmetric );
    /* Checks that marker i is enabled. */
    inline bool isMarkerEnabled( const Microarray_Set & set, int i );
    /* Checks that marker i is enabled. */
    inline bool isSameGene( const Microarray_Set & data, int i, int j);
    /** Operation stop request */
    bool stopped;
    /* * The percent completed */
    double completion;

  public:
    static double const NO_CONNECTION;

    MatrixOp() : stopped( false ), completion( 0.0 )
    {
    };

    void computeOneRow( Microarray_Set & data, Mutual_Info & mutual, Matrix & matrix, double threshold, int row_idx,
         int markerNo, int controlId, vector < int > * v, bool half_matrix, bool symmetric, double noise2 );
    /* create egdes matrix */
    void createEdgeMatrix( Microarray_Set & data, Mutual_Info & mutual, Matrix & matrix, double threshold, int controlId,
         double noise2, vector < int > & ids, std::vector < int > * v );
    /* compute mutual information between two gene expression vectors */
    double calculateMI( Microarray_Set & data, Mutual_Info & mutual, int probeId1, int probeId2,
         double threshold, double noise2, std::vector < int > * v );
    /* removes non contiguous edges by using the triangle inequality * all edges in the adjacency matrix is processed */
    void reduce( Microarray_Set & data, Matrix & matrix, Mutual_Info & mutual, double epsilon, double threshold, double noise2,
         vector < int > & ids, vector < int > * v, map <size_t, size_t> & transfac );
    void reduceOneNode( unsigned int row_idx, Microarray_Set & data, Matrix & matrix, Mutual_Info & mutual,
         double epsilon, double threshold, double noise2, vector < int > * v,  map <size_t, size_t> & transfac);
    /* fill the pair vector with gene expression values for MI computation. For the
     * unconditional MI calculation (array = NULL), expresssion values from all
     * arrays are used; for the conditional case (array not NULL), only microarrays * in the array vector are used. */
    void fillPair( unsigned, const Microarray_Set &, int, int, Pair_Vector &, const std::vector < int > * );

    inline void Set_Pair( Gene_Pair &, double, double, int );

    double getNodeMI(Matrix & matrix, unsigned int geneId1, unsigned int geneId2);

    bool protectedByTFLogic( map <size_t, size_t> & transfac, size_t geneId1, size_t geneId2, size_t geneId3 );

  };

} // aracne

#endif
