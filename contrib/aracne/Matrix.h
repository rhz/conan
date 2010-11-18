/*
 * Copyright (C) 2003  Columbia Genome Center
 * All Rights Reserved.
 *
 * Matrix.h -- Adjacency Matrix Data Structure
 *
 * $Id: Matrix.h,v 1.2 2005/10/04 18:24:12 kw2110 Exp $
 */

/******************************************************************************
* Class Matrix
*
* This class represents the adjacency matrix data structure that consists a list
* of genes and their connections to other genes.
*
* One can think of an adjacency matrix as a two dimensional array: each entry
* in the array is called a node; a node[i,j] will present if there is a
* connection between gene i and gene j, absent otherwise. The matrix is
* symmetric.
*
* Inner Class Node - represent an edge between two genes
*
*******************************************************************************
*/


#ifndef MATRIX_H__
  #define MATRIX_H__

  #include <istream>
  #include <string>
  #include <vector>
  #include <map>
  #include "Microarray_Set.h"
  #include "param.h"

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
using namespace std;

      #else
        #include <ext/hash_map>
        #define HASH_MAP hash_map

        #if __GNUC__ == 3 && __GNUC_MINOR__ == 0
using namespace std; // GCC 3.0
        #else
using namespace __gnu_cxx; // GCC 3.1 and later
        #endif
      #endif
    #endif
  #endif

namespace aracne {

  /****************************************************************************
  * Node (Inner) Class
  *
  * A node is one entry in the adjacency matrix that represents an edge between
  * two genes.
  * */
  class Node
  {
    double mutinfo;
    int intermediate;

  public:
    Node( double mi = 0.0 )
    {
      mutinfo = mi;
      intermediate = -1;
    }

    bool Is_Active() const
    {
      return ( intermediate < 0 );
    }

    void Set_Intermediate( int id )
    {
      intermediate = id;
    }

    int Get_Intermediate() const
    {
      return ( intermediate );
    }

    void Set_Mutualinfo( double m )
    {
      mutinfo = m;
    }

    double Get_Mutualinfo() const
    {
      return ( mutinfo );
    }

    friend std::ostream & operator << ( std::ostream & out_file, const Node & n )
    {
      out_file << "(" << n.mutinfo << ", " << n.intermediate << ")";
      return ( out_file );
    }
  };


  /****************************************************************************
  * Class Matrix
  *
  * Conceptually a square matrix storing the connections between a set of genes
  *
  * */
  class Matrix : public std::vector < HASH_MAP < int, Node > >
  {
    bool writeTriangular;
    bool writeReduced;
    bool writeEmptyGenes;

  public:
    Matrix()
    {
      writeTriangular = false;  // if true, only triangular half of the adjacency matrix will be written
      writeReduced = false;     // if true, intermediate nodes will be written
      writeEmptyGenes = false;  // if true, empty lines will be written
    }

    ~Matrix()
    {
    }

    void read( Microarray_Set & data, Parameter & p ) throw(std::runtime_error);
    void read( std::istream& in, Microarray_Set & data, Parameter & p ) throw(std::runtime_error);
    void writeGeneLine( std::ostream & out, Microarray_Set & data, int geneId );
    void writeGeneList( Microarray_Set & data, std::string name, int probeId );
    void write( std::ostream & out, Microarray_Set & data, vector < int > ids );
    void write( Microarray_Set & data, vector < int > ids, Parameter & p, bool writeFull = true );
    bool hasNode( int i, int j );

    inline void setWriteTriangular( bool status )
    {
      writeTriangular = status;
    }

    inline void setWriteReduced( bool status )
    {
      writeReduced = status;
    }

    inline void setWriteEmptyGenes( bool status )
    {
      writeEmptyGenes = status;
    }
  };

} // aracne

#endif
