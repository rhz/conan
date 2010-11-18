//
// Copyright (C) 2003 Columbia Genome Center
// All Rights Reserved.
//
// $Id: MutualInfo.cpp,v 1.4 2005/12/06 21:06:13 kw2110 Exp $
//

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include "MutualInfo.h"

namespace aracne {

  using namespace std;

  Mutual_Info::Mutual_Info( int maNo, int _miSteps, double sigma, int _type )
    : Microarray_Num( maNo ), type( _type ), variance2( 2*sigma*sigma )
  {
    switch ( type )
    {
      case MI_AVERAGE:
        miSteps = _miSteps;
        MA_Per_MI_Step = ( double )Microarray_Num / ( double )miSteps;
        MI_Space = new int * [miSteps];
        MI_Prob = new double * [miSteps];
        for ( int i = 0; i < miSteps; i++ )
        {
          MI_Space[i] = new int[miSteps];
          MI_Prob[i] = new double[miSteps];
          for ( int j = 0; j < miSteps; j++ )
          {
            MI_Space[i] [j] = 0;
            MI_Prob[i] [j] = 0.0;
          }
        }
      break;

      case MI_GAUSSIAN:

        // set coupla transform
        Set_Copula_Transform( true );

        // initialize the lookup table for computing Gaussian kernels
        prob_table = new double * [maNo];
        for ( int i = 0; i < maNo; i++ )
        {
          prob_table[i] = new double[maNo];
          for ( int j = 0; j < maNo; j++ )
          {
            prob_table[i][j] = -1.0;
          }
        }

        // initialize the lookup table for computing normalization factors
        norm_1D_table = new double[maNo];

        norm_2D_table = new double * [maNo];
        for ( int p = 0; p < maNo; p++ )
        {
          norm_2D_table[p] = new double[maNo];
        }

        Initialize_Norm_Table(sigma);

      break;

      default:
        throw std::runtime_error( "MI computation not supported" );
    }
  }

  void Mutual_Info::Initialize_Norm_Table( double sigma)
  {
    // initialize the 1D table
    for ( int i=0; i<Microarray_Num; i++ )
    {
      // array id is from 0 to M-1, after copula transform, it becomes from 0 to 1
      // with 1/(2M) at both sides
      double x = Copula_Transform ? ( ( double(i+1)-0.5 )/Microarray_Num ) : double( i );
      norm_1D_table[i] = 0.5*( erf( (1-x)/(sigma*M_SQRT2) ) - erf( (0-x)/(sigma*M_SQRT2) ) );
    }

    // initialize the 2D table
    for (int j=0; j<Microarray_Num; j++)
    {
      for ( int k=j; k<Microarray_Num; k++ )
      {
        norm_2D_table[j][k] = norm_1D_table[j]*norm_1D_table[k];
        norm_2D_table[k][j] = norm_2D_table[j][k];
      }
    }
  }

  /** Get Score for MI_GAUSSIAN */
  double Mutual_Info::Compute_Pairwise_MI( Pair_Vector & pairs )
  {
    double mi;
    switch ( type )
    {
      case Mutual_Info::MI_AVERAGE:
        mi = Get_Mutual_Info_XY( pairs, Loki::Int2Type < Mutual_Info::MI_AVERAGE > () );
      break;

      case Mutual_Info::MI_GAUSSIAN:
        mi = Get_Mutual_Info_XY( pairs, Loki::Int2Type < Mutual_Info::MI_GAUSSIAN > () );
      break;

      default:
        throw std::runtime_error( "MI computation not supported" );
    }
    return ( mi );
  }


  Mutual_Info::~Mutual_Info()
  {
    switch ( type )
    {
      case MI_AVERAGE:
        for ( int i = 0; i < miSteps; i++ )
        {
          delete[] MI_Space[i];
          delete[] MI_Prob[i];
        }
        delete[] MI_Space;
        delete[] MI_Prob;
      break;

      case MI_GAUSSIAN:
        for ( int i = 0; i < Microarray_Num; i++ )
        {
          delete[] prob_table[i];
          delete[] norm_2D_table[i];
        }

        delete[] prob_table;
        delete[] norm_2D_table;
        delete[] norm_1D_table;

      break;

      default:
        throw std::runtime_error( "MI computation not supported" );
    }
  }


  double Mutual_Info::Get_Mutual_Info_XY( Pair_Vector & pairs, const Loki::Int2Type < Mutual_Info::MI_GAUSSIAN > & rankType )
  {
    const size_t size = pairs.size();

    // rank and copula transformation
    Sort_X X_Sorter;
    sort( pairs.begin(), pairs.end(), X_Sorter );
    for ( size_t i = 0; i < size; i++ )
    {
      pairs[i].Set_XI(i);
      if (Copula_Transform)
      {
        pairs[i].Set_X( ( double(i+1)-0.5 )/size );
      }
    }

    Sort_Y Y_Sorter;
    sort( pairs.begin(), pairs.end(), Y_Sorter );
    for ( size_t i = 0; i < size; i++ )
    {
      pairs[i].Set_YI( i );
      if (Copula_Transform)
      {
        pairs[i].Set_Y( ( double(i+1)-0.5 )/size );
      }
    }

    double sum = 0.0;
    for ( size_t i=0; i<size; i++ )
    {
      double v = Get_Kernel( pairs, i );
      sum += log( v );
    }
    // double mi = sum / (static_cast <double>size);
    double mi = sum/(double)size;
    return ( std::max( mi, 0.0 ) );

  }

  double Mutual_Info::Get_Kernel( Pair_Vector & pairs, unsigned int i )
  {
    // mu_x, mu_y - the center of each gaussian kernel
    // ix, iy     - distance to the gaussian center in unit of index.
    // dx, dy     - the actual distance used to compute the kernel denstiy.
    //              if compula transformed, they are equal to ix and iy rescaled
    //              between 0 and 1;

    double fxy = 0.0;
    double fx = 0.0;
    double fy = 0.0;

    const size_t size = pairs.size();

    for ( size_t j = 0; j < size; j++ )
    {
      int mu_x = pairs[j].Get_XI();
      int mu_y = pairs[j].Get_YI();

      int ix = abs( pairs[i].Get_XI() - mu_x );
      int iy = abs( pairs[i].Get_YI() - mu_y );

      double dx = pairs[i].Get_X() - pairs[j].Get_X();
      double dy = pairs[i].Get_Y() - pairs[j].Get_Y();

      // compute the kernel density as necessary
      if ( prob_table[ix][iy] == -1.0 )
      {
        // if (ix, iy) is not computed, at least one of ix and iy is not computed
        // otherwise, all three should have been computed
        // when either ix=0 or iy=0, the trctdGssnBi.getProbability(dx, dy) is
        // equal to trctdGssnUni.getProbability(dx) or trctdGssnUni.getProbability(dy)
        if ( prob_table[ix][0] == -1.0 )
        {
          prob_table[ix][0] = std::exp(-(dx*dx)/variance2);
          prob_table[0][ix] = prob_table[ix][0];
        }

        if ( prob_table[0][iy] == -1.0 )
        {
          prob_table[0][iy] = std::exp(-(dy*dy)/variance2);
          prob_table[iy][0] = prob_table[0][iy];
        }

        prob_table[ix][iy] = prob_table[ix][0]*prob_table[0][iy];
        prob_table[iy][ix] = prob_table[ix] [iy];
      }

      fx += prob_table[ix][0]/norm_1D_table[mu_x];
      fy += prob_table[0][iy]/norm_1D_table[mu_y];
      fxy += prob_table[ix][iy]/norm_2D_table[mu_x][mu_y];
    }

    return ( ( fxy * size ) / ( fx * fy ) );
  }


  double Mutual_Info::Get_Mutual_Info_XY( Pair_Vector & pairs, const Loki::Int2Type < Mutual_Info::MI_AVERAGE > & rankType )
  {
    double H_xy = 0;
    double H_x = -log( 1 / ( double )miSteps );
    double H_y = H_x;
    Sort_X X_Sorter;
    Sort_Y Y_Sorter;

    sort( pairs.begin(), pairs.end(), X_Sorter );
    for ( size_t j = 0; j < pairs.size(); j++ )
    {
      pairs[j].Set_XI( j );
    }
    sort( pairs.begin(), pairs.end(), Y_Sorter );
    for ( size_t j = 0; j < pairs.size(); j++ )
    {
      pairs[j].Set_YI( j );
    }
    for ( int i = 0; i < miSteps; i++ )
    {
      for ( int j = 0; j < miSteps; j++ )
      {
        MI_Space[i] [j] = 0;;
      }
    }
    for ( size_t i = 0; i < pairs.size(); i++ )
    {
      int x = ( int )( ( double )pairs[i].Get_XI() / MA_Per_MI_Step );
      int y = ( int )( ( double )pairs[i].Get_YI() / MA_Per_MI_Step );
      MI_Space[x] [y] ++;
    }

    //
    // Compute the value for the initial block
    //
    int count = 0;
    for ( int di = 0; di < MIBLOCKS; di++ )
    {
      for ( int dj = 0; dj < MIBLOCKS; dj++ )
      {
        count += MI_Space[di] [dj];
      }
    }
    double sumP = 0;
    for ( int i = 0; i < miSteps; i++ )
    {
      int count1 = count;
      for ( int j = 0; j < miSteps; j++ )
      {
        double p = ( double )count1 / ( double )( pairs.size() * MIBLOCKS * MIBLOCKS );
        sumP += p;
        if ( p > 0 )
          H_xy -= p * log( p );
        MI_Prob[i] [j] = ( double )count1;
        for ( int off = 0; off < MIBLOCKS; off++ )
        {
          //
          // Wrap around the torus
          //
          int x = ( i + off ) % miSteps;
          int y = ( j + MIBLOCKS ) % miSteps;
          count1 -= MI_Space[x] [j];
          count1 += MI_Space[x] [y];
        }
      }
      for ( int off = 0; off < MIBLOCKS; off++ )
      {
        int x = ( i + MIBLOCKS ) % miSteps;
        count -= MI_Space[i] [off];
        count += MI_Space[x] [off];
      }
    }
    double mi = ( H_x + H_y - H_xy ) / ( H_x + H_y );
    return mi;
  }

} // aracne
