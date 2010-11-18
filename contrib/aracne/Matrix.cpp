//
//
// Copyright (C) 2003 Columbia Genome Center
// All Rights Reserved.
//
// Matrix.cc -- Create a matrix.
// vi: set ts=4
//
// $Id: Matrix.cpp,v 1.2 2005/10/04 18:24:12 kw2110 Exp $
//

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include "Matrix.h"

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


  //
  // reads an adjacency matrix from an input stream
  //
  void
  Matrix::read(Microarray_Set & data, Parameter & p) throw(std::runtime_error)
  {
    ifstream in(p.adjfile.c_str());
    if(!in.good())
    {
      cerr << "Bad File!" << endl;
      throw(std::runtime_error("Problem with reading the file: "+p.adjfile+"."));
    }
    read(in, data, p);
  }


  void
  Matrix::read(std::istream& in, Microarray_Set & data, Parameter & p) throw(std::runtime_error)
  {
    string line;
    string label;
    string value;
    getline( in, line );
    // by pass the lines starting with ">"
    while (strncmp(line.c_str(), ">", 1) == 0) {
      getline(in, line);
    }

    while ( in.good() ) {
      istringstream sin( line );

      getline( sin, label, '\t' );
      int geneId1 = data.getProbeId( label );
      if (geneId1 == -1) {
        throw(runtime_error("Cannot find marker: "+label+" in the ADJ file!"));
      }

      data.Get_Marker( geneId1 ).Enable();

      // Add enough entries to the vector of hashes to store data for geneId1
      while ( int(size()) <= geneId1 )
        push_back( HASH_MAP < int, Node > () );

      getline( sin, label, '\t' );
      while ( sin.good() ) {
        getline( sin, value, '\t' );

        double mi = atof( value.c_str() );
        if ( mi >= p.threshold ) {
          int geneId2 = data.getProbeId( label );
          if (geneId2 == -1) {
            throw(runtime_error("Cannot find marker: "+label+" in the ADJ file!"));
          }

          // Add enough entries to the vector of hashes to store data for geneId2
          while ( int(size()) <= geneId2 ) {
            push_back( HASH_MAP < int, Node > () );
          }

          Node & n = ( * this ) [geneId1] [geneId2];
          n.Set_Mutualinfo( mi );
        }
        getline( sin, label, '\t' );
      }
      getline( in, line );
    }
  }


  bool Matrix::hasNode( int i, int j )
  {
    size_t ii = max( i, j );
    size_t jj = min( i, j );
    if ( ii == jj )
      return ( false );
    if ( size() <= ii )
      return ( false );
    HASH_MAP < int, Node > & m = ( * this ) [i];
    if ( m.find( j ) == m.end() )
      return ( false );
    return ( true );
  }


  void Matrix::write( Microarray_Set & data, vector < int > ids, Parameter & p, bool writeFull )
  {
    if ( writeFull )
    {
      std::fstream output( p.outfile.c_str(), std::ios::out );
      cout << "Writing matrix: " << p.outfile << endl;

      output << ">  Input file      " << p.infile << endl;
      output << ">  ADJ file        " << p.adjfile << endl;
      output << ">  Output file     " << p.outfile << endl;
      output << ">  Algorithm       " << p.algorithm << endl;
      output << ">  Kernel width    " << p.sigma << endl;
      output << ">  No. bins        " << p.miSteps << endl;
      output << ">  MI threshold    " << p.threshold << endl;
      output << ">  MI P-value      " << p.pvalue << endl;
      output << ">  DPI tolerance   " << p.eps << endl;
      output << ">  Correction      " << p.correction << endl;
      output << ">  Subnetwork file " << p.subnetfile << endl;
      output << ">  Hub probe       " << p.hub << endl;
      output << ">  Control probe   " << p.controlId << endl;
      output << ">  Condition       " << p.condition << endl;
      output << ">  Percentage      " << p.percent << endl;
      output << ">  TF annotation   " << p.annotfile << endl;
      output << ">  Filter mean     " << p.mean << endl;
      output << ">  Filter CV       " << p.cv << endl;

      write( output, data, ids );
      output.flush();
      output.close();
    }
  }

  //
  // writes an adjacency matrix to an output stream. The format is
  //
  void Matrix::write( ostream & out, Microarray_Set & data, vector < int > ids )
  {
    if (ids.empty())
    {
      for ( size_t i = 0; i < size(); i++ )
      {
        writeGeneLine( out, data, i );
      }
      out.flush();
    }
    else
    {
      for ( size_t i = 0; i < ids.size(); i++ )
      {
        writeGeneLine( out, data, ids[i] );
      }
    }
  }


  //
  // Write a line in the matrix to an output file.
  //
  void Matrix::writeGeneLine( ostream & out, Microarray_Set & data, int Id )
  {
    HASH_MAP < int, Node > & m = ( * this ) [Id];
    Marker & marker = data.Get_Marker( Id );
    const string & label = marker.Get_Accnum();

    if ( ( m.size() == 0 ) && ( !writeEmptyGenes || ( !marker.Enabled() && !marker.Is_Control() ) ) )
      return;

    out << label;

    vector < int > sorted_keys;
    for ( const_line_iterator it = m.begin(); it != m.end(); ++it )
    {
      sorted_keys.push_back( it->first );
    }

    sort( sorted_keys.begin(), sorted_keys.end() );

    vector < int >::const_iterator viter;
    for ( viter = sorted_keys.begin(); viter != sorted_keys.end(); ++viter )
    {
      int Gene_Id = * viter;
      const Node & node = m[* viter];
      if ( writeTriangular && Gene_Id <= Id )
        continue;

      Marker & m = data.Get_Marker(Gene_Id);
      const string & affyid = m.Get_Accnum();
      if ( writeReduced )
      {
        out << '\t' << affyid;
        if ( node.Get_Intermediate() >= 0 )
          cout << '.' << node.Get_Intermediate();
        out << '\t' << node.Get_Mutualinfo();
      }
      else
      {
        if ( node.Get_Intermediate() >= 0 )
          continue;
        double mi = node.Get_Mutualinfo();
        out << '\t' << affyid << "\t" << mi;
      }
    }
    out << endl;
  }


  void Matrix::writeGeneList( Microarray_Set & data, std::string name, int probeId )
  {
    std::string filename = name + ".adj";
    std::fstream output( filename.c_str(), std::ios::out );
    cout << "Writing gene list: "<< filename << endl;

    HASH_MAP < int, Node > & m = ( * this ) [probeId];
    Marker & marker = data.Get_Marker( probeId );

    if ( ( m.size() == 0 ) && ( !writeEmptyGenes || ( !marker.Enabled() && !marker.Is_Control() ) ) )
    {
      return;
    }

    HASH_MAP < int, Node >::const_iterator it = m.begin();
    while ( it != m.end() )
    {
      Node const & n = it->second;
      int j = it->first;
      if ( j != probeId )
      {
        output << j << "\t" << data.Get_Marker( j ).Get_Accnum() << "\t" << n.Get_Mutualinfo() << endl;
      }
      it++;
    }

    output.flush();
    output.close();
  }

} // aracne
