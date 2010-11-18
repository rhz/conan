//
// Copyright (C) 2003  Columbia Genome Center
// All Rights Reserved.
//
// MutualInfo.h -- Class and constant definitions for mutual
// information calculation progoram.
//
// $Id: MutualInfo.h,v 1.3 2005/09/29 20:45:49 kw2110 Exp $
//

#ifndef MUTUALINFO_H__
  #define MUTUALINFO_H__
  #include "TypeManip.h"
  #include <stdexcept>

  #ifndef M_SQRT2
    #define M_SQRT2 1.4142135623730950488016887
  #endif

namespace aracne {

  using namespace std;


  class Gene_Pair
  {

    double x;
    double y;
    int xi; // index of x
    int yi; // index y
    int maId;

  public:

    inline double Get_X() const { return x; }

    inline double Get_Y() const { return y; }

    inline int Get_XI() const { return xi; }

    inline int Get_YI() const { return yi; }

    inline int Get_MaID() const { return maId; }

    inline void Set_X( double X ) { x = X; }

    inline void Set_Y( double Y ) { y = Y; }

    inline void Set_XI( int XI ) { xi = XI; }

    inline void Set_YI( int YI ) { yi = YI; }

    inline void Set_MaID( int MaID ) { maId = MaID; }
  };


  typedef std::vector < Gene_Pair > Pair_Vector;


  class Sort_X : std::binary_function < Gene_Pair, Gene_Pair, bool >
  {

  public:

    bool operator() ( const Gene_Pair & a, const Gene_Pair & b ) const
    {
      if ( a.Get_X() != b.Get_X() )
      {
        return ( a.Get_X() < b.Get_X() );
      }
      else
      {
        return ( a.Get_MaID() < b.Get_MaID() );
      }
    }
  };


  class Sort_Y : std::binary_function < Gene_Pair, Gene_Pair, bool >
  {

  public:

    bool operator() ( const Gene_Pair & a, const Gene_Pair & b ) const
    {
      if ( a.Get_Y() != b.Get_Y() )
      {
        return ( a.Get_Y() < b.Get_Y() );
      }
      else
      {
        return ( a.Get_MaID() < b.Get_MaID() );
      }
    }
  };



  /****************************************************************************
  * Class MutualInfo
  *
  * The main MutualInfo class
  *
  **/
  class Mutual_Info
  {
    static const int MIBLOCKS = 2;
    static const int BINS = 1000;

    int     Microarray_Num;
    int     miSteps;
    int     ** MI_Space;
    double  ** MI_Prob;
    double  MA_Per_MI_Step;
    double  * histogram;
    double  * background;
    int     Max_Histo_Bin;
    int     Max_Background_Bin;

    //type of caluclation
    int type;

    // variance2 = 2*sigma*sigma, doing this just same some computation
    double variance2;
    double ** prob_table;
    double *  norm_1D_table;
    double ** norm_2D_table;
    bool Copula_Transform;

  public:

    static int const RANK_NONE = 0;
    static int const STD_REGRESSION = 1;
    static int const RANK_REGRESSION = 2;
    static int const MI_GAUSSIAN = 3;
    static int const MI_AVERAGE = 4;

    Mutual_Info( int, int, double, int );
    ~Mutual_Info();
    int Get_Type() const { return type; }
    void Set_Copula_Transform( bool t = true ) { Copula_Transform = t; }
    bool Is_Copula_Transform() const { return Copula_Transform; }
    double Compute_Pairwise_MI( Pair_Vector & );
    int Get_Microarray_Num() { return Microarray_Num; }

  private:

    void Initialize_Norm_Table( double );
    double Get_Kernel( Pair_Vector & pairs, unsigned int i );
    double Get_Mutual_Info_XY( Pair_Vector &, const Loki::Int2Type < Mutual_Info::MI_AVERAGE > & );
    double Get_Mutual_Info_XY( Pair_Vector &, const Loki::Int2Type < Mutual_Info::MI_GAUSSIAN > & );

  };

} // aracne

#endif
