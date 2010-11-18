//
// Copyright (C) 2003  Columbia Genome Center
// All Rights Reserved.
//
// Microarray_Set.h -- Microarray Data structures.
//
// $Id: Microarray_Set.h,v 1.3 2005/10/06 02:49:20 kw2110 Exp $
//


/******************************************************************************
* Class Microarray_Set
*
* A Microarray dataset consists of multiple arrays, each of which measures the
* expression of thousands of "Markers". Here we represent each of these
* measurement using a "Probe".
*
* Inner Class Marker - represent markers on each array
*
* Inner Class Probe - represent a measurement on a marker in an array
*
* Inner Class ArrayValuePair
*             SortDecreasing_ArrayValuePair
*             SortIncreasing_ArrayValuePair
*   - datastructure for sorting the arrays be the
*   the expression values of a particular marker
*
*******************************************************************************
*/


#ifndef DAT_H__
  #define DAT_H__

  #include <istream>
  #include <string>
  #include <vector>
  #include <map>
  #include "TypeManip.h"

  #ifdef __BCPLUSPLUS__   // For Borland Compiler
    #include <stlport/hash_map>
using namespace std;
  #else
    #ifdef __GNUC__        // For GNU gcc compiler
      #if __GNUC__ < 3
        #include <hash_map.h>
using namespace std;

      #elif __GNUC__ >= 4 && __GNUC_MINOR__ >= 3
        #include <unordered_map>

      #else
        #include <ext/hash_map>
        #if __GNUC_MINOR__ == 0
using namespace std; // GCC 3.0
        #else
using namespace __gnu_cxx; // GCC 3.1 and later
        #endif
      #endif
    #endif
  #endif

  #include <cstring>
  #include <stdexcept>

namespace aracne {

  using namespace std;

  /***********************************************************
   * Class Marker
   *
   * a marker describes a gene in a detailed way, including
   * attributes not useful in other places.
   */
  class Marker
  {
    int idnum;
    std::string accnum;
    std::string label;
    double var;
    bool isActive;
    bool isControl;

  public:
    Marker( int id = 0, const string & a = "", const string & l = "", double sigma = 0.0,
         bool active = true, bool control = false ) : idnum( id ), accnum( a ), label( l ), var( sigma ), isActive( active ),
         isControl( control )
         {
    }

    int Get_ID() const
    {
      return idnum;
    }

    const std::string & Get_Accnum() const
    {
      return accnum;
    }

    const std::string & Get_Label() const
    {
      return label;
    }

    void Set_ID( int id )
    {
      idnum = id;
    }

    void Set_Accnum( std::string & a )
    {
      accnum = a;
    }

    void Set_Label( std::string & l )
    {
      label = l;
    }

    void Set_Var ( double v )
    {
      var = v;
    }

    double Get_Var()
    {
      return var;
    }

    bool hasAccession( std::string & a )
    {
      // return ( accnum.find( a ) != accnum.npos );
      return ( accnum.compare(a)==0 );
    }

    void Enable()
    {
      isActive = true;
    }

    void Disable()
    {
      isActive = false;
    }

    bool Enabled() const
    {
      return ( isActive );
    }

    void set_Control()
    {
      isControl = true;
    }

    void set_No_Control()
    {
      isControl = false;
    }

    bool Is_Control()
    {
      return ( isControl );
    }

    friend std::ostream & operator << ( std::ostream & out_file, const Marker & m )
    {
      out_file << "(";
      out_file << m.label << "; " << m.accnum << "; ";
      out_file << m.idnum << "; ";
      if ( m.isActive )
        out_file << "Active; ";
      else
        out_file << "Not Active; ";
      if ( m.isControl )
        out_file << "Control)";
      else
        out_file << "Not Control)";

      return ( out_file );
    }
  };


  /*
   * A Marker_Set is an array of Markers, indexed by the
   * ID numbers contained within the marker.
   */
  typedef std::vector < Marker > Marker_Set;


  /***********************************************************
   * Class Probe
   *
   * Microarrays are really silicon chips that contain probes
   * to detect genes.  Each probe is a tuple containing values
   * and pvalues, and there are n such tuples arranged in an
   * array, indexed by the same numbers which indexes Markers
   * in a set.  Since Markers describe genes, the correspondence
   * is natural.  Associated with each set of Microarrays is
   * also a Marker_Set object.  If the cardinality of the marker
   * set is m, then each Microarray will have m probes
   *  associated with it.
   */
  class Probe
  {
    double val, pval;

  public:
    /** deliminator btw probes when they are written to a stream */
    static const char * DELIMINATOR;

    Probe( double v = 0.0, double pv = 0.0 )
    {
      Set_Value( v );
      Set_PValue( pv );
    }

    double Get_Value() const
    {
      return val;
    }

    double Get_PValue() const
    {
      return pval;
    }

    void Set_Value( double v )
    {
      val = v;
    }

    void Set_PValue( double pv )
    {
      pval = pv;
    }


  };

  std::ostream & operator << ( std::ostream & out, const Probe & p );


  typedef std::vector < Probe > Microarray;
  typedef std::vector < Microarray > Microarray_Vector;


  /*************************************************************
   * Class Microarray_Set
   *
   * Note that a Microarray_Set contains one or more microarrays
   * *and* has associated with it a set of markers.
   */
  class Microarray_Set
  {
    Marker_Set markerset;
    Microarray_Vector uarrays;
    std::map < std::string, int > markerMap;
    std::vector < std::string > header;
    unsigned int readMarkerWithPvalue( std::istream & in, const unsigned arrno ) throw( std::runtime_error );
    unsigned int readMarkerNoPvalue( std::istream & in, const unsigned arrno ) throw( std::runtime_error );
    unsigned int readHeader( std::istream & in ) throw( std::runtime_error );
  public:
    static const char * HEADER_DELIMINATOR;

    Microarray_Set()
    {
    }

    int Get_Num_Microarrays() const
    {
      return uarrays.size();
    }

    int Get_Num_Markers() const
    {
      return markerset.size();
    }

    int Get_Num_Active_Markers();

    std::string Get_Array_Header(int i) {
        return header[i];
    }

    Marker & Get_Marker( unsigned int i ) throw( std::runtime_error );
    const Marker & Get_Marker( unsigned int i ) const throw( std::runtime_error );
    const std::string Get_Marker_AffyId(unsigned int i) throw( std::runtime_error );

    int Get_MarkerId( const std::string & label )
    {
      if ( markerMap.find( label ) != markerMap.end() )
      {
        return markerMap[label];
      }
      return -1;
    }

    const std::vector < std::string > & Get_Header() const
    {
      return header;
    }

    void Set_ColHeader( unsigned int i, const string & hdr );
    const Probe & Get_Probe( unsigned int i, unsigned int m ) const throw( std::runtime_error );
    const Microarray & Get_Microarray( unsigned int i ) const throw( std::runtime_error );
    const double Get_Value( unsigned int i, unsigned int m ) const throw( std::runtime_error );
    const double Get_PValue( unsigned int i, unsigned int m ) const throw( std::runtime_error );
    void Set_Marker( unsigned int i, const Marker & m );
    void Set_Microarray( unsigned int i, const Microarray & a );
    void Set_Probe( unsigned int i, unsigned int j, Probe p );
    void Set_Value( unsigned int i, unsigned int m, double v ) throw( std::runtime_error );
    void Set_PValue( unsigned int i, unsigned int m, double v ) throw( std::runtime_error );
    int get_Id( std::string & a );
    /** Read a matrix from a file. */
    void read( const std::string & fileName ) throw( std::runtime_error );
    /** Read a matrix from an input stream */
    void read( std::istream & in ) throw( std::runtime_error );
    /* Used to filter all genes that have a mean lower than
     * minMean and a variance lower than mean * ratio. This limits
     * further analysis only to genes that have a reasonable dynamical range
     * */
    int filter( std::vector < int > & ids, double minMean = 50.0, double minSigma = 20.0, int ctlid = -1 );
    int filter( double minMean = 50.0, double sigma = 20.0, int ctlid = -1 );
    void computeMarkerVariance( vector < int > * arrays );
    double variance( unsigned int m, vector < int > * arrays );
    double getMarkerVariance( int m );
    int getProbeId( std::string label );
    void getHighLowPercent( double x, int mId, vector < int > & lower, vector < int > & upper );
    void getRandomSubsamples( double x, int mId, vector < int > & lower, vector < int > & upper );
    void bootStrap( vector < int > & boot, vector < int > * arrays );
    void addNoise();
    void shuffleGene(unsigned int idx);
  };

  std::ostream & operator << ( std::ostream & out, const Microarray_Set & p );


  class ArrayValuePair {
          int	arrayId;
          double	value;
  public:
          ArrayValuePair() { ; }
          ArrayValuePair(int Id, double v){ set(Id, v); }
          void	set(int iId, double v)	{ arrayId = iId; value = v; }
          int	getId()	const		{ return arrayId; }
          double	getValue() const	{ return value; }
  };

  class SortIncreasing_ArrayValuePair : std::binary_function<ArrayValuePair, ArrayValuePair, bool> {
  public:
          bool operator() (const ArrayValuePair& a, const ArrayValuePair& b) const {
                  return(a.getValue() < b.getValue());
          }
  };

  class SortDecreasing_ArrayValuePair : std::binary_function<ArrayValuePair, ArrayValuePair, bool> {
  public:
      bool operator() (const ArrayValuePair& a, const ArrayValuePair& b) const {
          return(a.getValue() > b.getValue());
      }
  };

} // aracne

#endif
