/* * Copyright (C) 2003  Columbia Genome Center * All Rights Reserved. *
* program_standalone.cc -- Class definitions for the standalone program. *
* $Id: standalone_program.cpp,v 1.4 2005/12/13 23:08:02 wl2131 Exp $ */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include "Microarray_Set.h"
#include "MutualInfo.h"
#include "Matrix_Op.h"
#include "parseargs.h"
#include "param.h"

using namespace std;
using namespace aracne;

/******************************************************************************/

void usage( const char * prog )
{
  ifstream infile;
  infile.open ("usage.txt", ifstream::in);

  if(!infile.good()) {
    throw runtime_error(string("Cannot find file \"usage.txt\"! Please make sure it is in the current working directory."));
  }

  while (infile.good() && infile.peek() != EOF) {
    cout << (char) infile.get();
  }
  infile.close();
  cout << endl << endl;
  exit( 1 );
}


void findKernelWidth(int n, Parameter & p)
{
  string filename = p.home_dir+"config_kernel.txt";
  ifstream infile(filename.c_str(), ifstream::in);

  if(!infile.good()) {
    throw runtime_error(string("Cannot find file: "+filename+"! Please use '-H' option to specify the ARACNE home directory."));
  }

  string line;
  double alpha;
  double beta;
  getline(infile, line);
  // by pass the lines starting with ">"
  while (strncmp(line.c_str(), ">", 1) == 0) {
    getline(infile, line);
  }

  istringstream sin(line);
  if ( sin.good() && (sin.peek() != EOF) ) {
    sin >> alpha >> beta;
  }
  else {
    throw runtime_error(string("Configuration file format error: "+filename));
  }
  infile.close();
  p.sigma = alpha*pow(n,beta);
}


void findThreshold(int n, Parameter & p)
{
  string filename = p.home_dir+"config_threshold.txt";
  ifstream infile(filename.c_str(), ifstream::in);

  if(!infile.good()) {
    throw runtime_error(string("Cannot find file: "+filename+"! Please use '-H' option to specify the ARACNE home directory."));
  }

  string line;
  double alpha;
  double beta;
  double gamma;
  getline(infile, line);
  // by pass the lines starting with ">"
  while (strncmp(line.c_str(), ">", 1) == 0) {
    getline(infile, line);
  }

  istringstream sin(line);
  if ( sin.good() && (sin.peek() != EOF) ) {
    sin >> alpha >> beta >> gamma;
  }
  else {
    throw runtime_error(string("Configuration file format error: "+filename));
  }
  infile.close();
  p.threshold = (alpha - log(p.pvalue))/((-beta)+(-gamma)*n);
}



/******************************************************************************/

Parameter parseParameter( int argc, char * argv[], int n )
{
  //initializes to default
  Parameter p;
  std::string progname = argv[0];
  std::string temp;

  // by pass the first n arguments
  int i = 0;
  while (i < n)
  {
    ++argv; --argc; i++;
  }

  // parse arguments
  ARGBEGIN
  {
    case 'i' :  // input file
      p.infile = ARGF();
    break;
    case 'o' :  // output file
      p.outfile = ARGF();
    break;
    case 'j' :  // adjacency matrix file
      p.adjfile = ARGF();
    break;
    case 'h' :  // hub gene
      p.hub = ARGF();
    break;
    case 's' :  // subset of probes
      p.subnetfile = ARGF();
    break;
    case 'l' :  // TF annotation file
      p.annotfile = ARGF();
    break;
    case 'c' :  // condition
      temp = ARGF();
      p.condition = temp.substr(0,1);
      p.controlId = temp.substr(1,temp.length()-1);
      p.percent = atof( ARGF() );
    break;
    case 'f' :  // filtering
      p.mean = atof( ARGF() );  // mean
      p.cv = atof( ARGF() );    // coefficient of variance
    break;
    case 't' : // mi threshold
         p.threshold = atof( ARGF() );
    break;
    case 'e' : // DPI tolerance
         p.eps = atof( ARGF() );
    break;
    case 'k' : // gaussian kernel width
         p.sigma = atof( ARGF() );
    break;
    case 'b' : // mi step size
         p.miSteps = atoi( ARGF() );
    break;
    case 'r' : // bootstrap sample number
         p.sample = atoi( ARGF() );
    break;
    case 'p' : // p-value
         p.pvalue = atof( ARGF() );
    break;
    case 'v' :  // verbose
         p.verbose = ARGF();
    break;
    case 'a' :  // algorithm
         p.algorithm = ARGF();
    break;
    case 'n' :  // correction for noise
         p.correction = atof( ARGF() );
    break;
    case 'H' : // ARACNE_HOME
         p.home_dir = ARGF();
    break;
    default :
         cout << "Unkown parameter: " << ARGC() << endl;
         cout << "Type \"--help\" for more information." << endl;
		 exit(1);
  }
  ARGEND;
  try {
    checkParameter(p);
  }
  catch ( string & s ) {
    std::cout << s << std::endl;
    exit( 1 );
  }
  displayParameter(p);
  return ( p );
}

/******************************************************************************/

int run( int argc, char * argv[] )
{
  if ( argc < 2 )
  {
    usage( argv[0] );
  }
  else if ( (strcmp(argv[1],"--help")==0) || (strcmp(argv[1],"--h")==0) )
  {
    usage( argv[0] );
  }
  /* ***************************************************************************
   * *          Calculate Mutual Information between two genes                 *
   * ***************************************************************************
   */
  else if ( strcmp( argv[1], "-cal" ) == 0 ) {
    if ( argc < 5 )
      usage( argv[0] );

    Parameter p;
    int nboot = 0;
    if ( strcmp(argv[4], "-boot")==0 )
    {
      nboot = atoi(argv[5]);
      std::cout << "No. of bootstrap: " << nboot << std::endl;
      p = parseParameter( argc, argv, 5 );
    }
    else
    {
      p = parseParameter( argc, argv, 3 );
    }

    Microarray_Set data;
    MatrixOp aracne;

    data.read( p.infile );

    std::string msg;
    int probeId1 = data.getProbeId( argv[2] );
    if (probeId1 == -1) {
      msg = "Cannot find marker: ";
      msg.append(argv[2]);
      throw runtime_error(msg+"!");
    }
    int probeId2 = data.getProbeId( argv[3] );
    if (probeId2 == -1) {
      msg = "Cannot find marker: ";
      msg.append(argv[3]);
      throw runtime_error(msg+"!");
    }

    int TYPE;
    if (cmpIgnoreCase( p.algorithm.c_str(), "accurate" )==0) {
      TYPE = Mutual_Info::MI_GAUSSIAN;
    }
    else  {
      TYPE = Mutual_Info::MI_AVERAGE;
    }

    vector < int > low;
    vector < int > high;
    vector < int > * v = NULL;
    int controlId = -1;
    int nsample = data.Get_Num_Microarrays();

    if (p.controlId != "") {
      controlId = data.getProbeId( p.controlId );
      bool isHigh = (p.condition == "+") ? true : false;
      data.getHighLowPercent( p.percent, controlId, low, high );
      v = ( isHigh ) ? & high : & low;
      nsample = v->size();
    }

    cout << "Marker No: " << data.Get_Num_Markers() << " (" << data.Get_Num_Active_Markers() << " active)"
         << ", Array No: " << nsample << endl;

    if (p.sigma == 99) {
      findKernelWidth( nsample, p );
      cout << "Kernel width determined for this dataset: " << p.sigma << endl;
    }

    if ( p.threshold == 0 && p.pvalue != 1 ) {
      findThreshold( nsample, p );
      cout << "MI threshold determined for p=" << p.pvalue << ": " << p.threshold << endl;
    }

    Mutual_Info mi(nsample , p.miSteps, p.sigma, TYPE );

    std::cout << "[NETWORK] Calculating MI ... " << std::endl;
    double mutualinfo = aracne.calculateMI( data, mi, probeId1, probeId2, p.threshold, p.correction, v );
    std::cout << "Probe 1: " << probeId1 << " (" << data.Get_Marker_AffyId(probeId1) << ")" << std::endl;
    std::cout << "Probe 2: " << probeId2 << " (" << data.Get_Marker_AffyId(probeId2) << ")" << std::endl;
    std::cout << "MI = " << mutualinfo << std::endl;

    // perform boostrapping
    if (nboot > 0)
    {
      std::cout << "BOOTSTRAPPING ... " << std::endl;
      if (p.outfile == "") {
          throw std::runtime_error("An output file must be specified!");
      }
      std::fstream output( p.outfile.c_str(), std::ios::out );
      output << mutualinfo << endl;

      vector < int > bs;
      vector < int > * u;
      int iter = 0;
      while (iter++ < nboot)
      {
        // sampling with replacement
        data.bootStrap(bs, v);
        u = & bs;
        output << aracne.calculateMI( data, mi, probeId1, probeId2, p.threshold, p.correction, u ) << endl;
        std::cout << ".";
        if (iter%80 == 0) {
          std::cout << std::endl;
        }
      }
      std::cout << std::endl;
      output.flush();
      output.close();
    }
  }
  /* ***************************************************************************
   * *                   Standard ARACNE operations                            *
   * ***************************************************************************
   */
  else
  {
    Parameter p = parseParameter( argc, argv, 0 );

    // read the input dataset
    Microarray_Set data;
    Matrix matrix;
    MatrixOp aracne;
    data.read( p.infile );

    if ( p.mean > 0 || p.cv > 0 ) {
      int ndisabled = data.filter( p.mean, p.cv );
      std::cout << ndisabled << " markers disabled due to lack of dynamic range." << std::endl << std::endl;
    }

    int TYPE;
    if (cmpIgnoreCase( p.algorithm.c_str(), "accurate" )==0) {
      TYPE = Mutual_Info::MI_GAUSSIAN;
    }
    else  {
      TYPE = Mutual_Info::MI_AVERAGE;
    }

    vector < int > ids;
    map < size_t, size_t > transfac;
    vector < int > low;
    vector < int > high;
    vector < int > bs;
    vector < int > * v = NULL;
    int controlId = -1;
    int nsample = data.Get_Num_Microarrays();

    if (p.controlId != "") {
      controlId = data.getProbeId( p.controlId );
      if (controlId == -1) {
        throw std::runtime_error(std::string("Cannot find marker: "+p.controlId+"!"));
      }
      bool isHigh = (p.condition == "+") ? true : false;
      data.getHighLowPercent( p.percent, controlId, low, high );
      v = ( isHigh ) ? & high : & low;
      nsample = v->size();
    }

    cout << "Marker No: " << data.Get_Num_Markers() << " (" << data.Get_Num_Active_Markers() << " active)"
         << ", Array No: " << nsample << endl;

    if (p.sigma == 99) {
      findKernelWidth( nsample, p );
      cout << "Kernel width determined for this dataset: " << p.sigma << endl;
    }

    if ( p.threshold == 0 && p.pvalue != 1 ) {
      findThreshold( nsample, p );
      cout << "MI threshold determined for p=" << p.pvalue << ": " << p.threshold << endl;
    }

    Mutual_Info mi(nsample , p.miSteps, p.sigma, TYPE );

    // to determine whether MI correction will be performed; if yes, prepare
    // the data to pre-compute the variance of each gene's expression in log scale
    // when MI is computed only within a subset of samples, variance is also computed
    // in the same subset.
    // Note, variance computation is the same no matter bootstrapping is being performed,
    // cause I don't think it make sense.
    if (p.correction != 0) {
      data.computeMarkerVariance( v );
    }

    // the parseParameter function should have made sure that p.hub and p.subnet
    // can not be supplied at the same time
    if (p.hub != "") {
      int hubId = data.getProbeId( p.hub );
      if (hubId == -1) {
        throw runtime_error(string("Cannot find the hub probe: " + p.hub + ", nothing to be computed!"));
      }
      ids.push_back(hubId);
    }

    if (not(p.subnet.empty())) {
      for (size_t i = 0; i < p.subnet.size(); i++) {
        int gid = data.getProbeId( p.subnet[i] );
        if (gid == -1) {
          cout << "Cannot find probe: " << p.subnet[i] << " in \"" << p.subnetfile
               << "\" ... ignored." << endl;
        }
        else {
          ids.push_back(gid);
        }
      }
    }

    if (not(p.tf_list.empty()))
    {
      for (size_t i = 0; i < p.tf_list.size(); i++)
      {
        int gid = data.getProbeId( p.tf_list[i] );
        if (gid == -1) {
          cout << "Cannot find probe: " << p.tf_list[i] << " in \"" << p.annotfile
               << "\" ... ignored." << endl;
        }
        else {
          transfac[gid] = 1;
        }
      }
    }

    // if existing adjacency matrix is provided, no MI computation needed
    if (p.adjfile != "")
    {
      matrix.read( data, p );
      if ( p.eps != 1 ) {
        cout << "[NETWORK] Applying DPI ..." << endl;
        aracne.reduce( data, matrix, mi, p.eps, p.threshold, p.correction, ids, v, transfac );
      }
    }
    else
    {
      // use resampling data
      if (p.sample > 0) {
        srand ( p.sample );
        data.bootStrap(bs, v);
        v = & bs;
      }
      // add noise to randomize order among identical expression values
      data.addNoise();

      // this is THE function that reconstruct the network!
      // 1) if conditional network is being reconstructed, 'controlId' will be the
      // gene id of the conditioning gene, and 'v' points to the set of array ids
      // in which the conditioning gene is certain percentage high/low.
      // 2) if 'ids' not empty, a subset of the network will be reconstructed
      // 3) if 'v' is non-empty while 'controlId' is -1, then bootstrapping is being performed.
      // 4) if MI correction is taken place, 'p.correction' will be a non-zero real number
      aracne.createEdgeMatrix(data, mi, matrix, p.threshold, controlId, p.correction, ids, v );

      if (p.eps != 1) {
        cout << "[NETWORK] Applying DPI ..." << endl;
        aracne.reduce( data, matrix, mi, p.eps, p.threshold, p.correction, ids, v, transfac );
      }
    }

    if ( p.outfile == "" ) {
      createOutfileName( p );
    }
    matrix.write(data, ids, p);
  }
  return 0;
}

/******************************************************************************/

int main( int argc, char * argv[] )
{
  try {
    return run( argc, argv );
  }
  catch ( string & s )
  {
    cout << s << endl;
  }
}
