/*
 * Conan - COmplex Network ANalisys
 * Copyright (C) 2008-2009  Ricardo Honorato Zimmer [rikardo.horo@gmail.com]
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef MUTUAL_INFORMATION_HPP
#define MUTUAL_INFORMATION_HPP
#include <conan/config.hpp>
#include <conan/inference/aracne.hpp>

#include <fstream>

#if 0
#define HAVE_INLINE
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics_double.h>
#endif

namespace conan { namespace inference { namespace mi {

  using namespace aracne;

  void findKernelWidth(int n, Parameter & p)
  {
    std::string filename = p.home_dir + "config_kernel.txt";
    std::ifstream infile(filename.c_str(), std::ifstream::in);

    if (!infile.good()) {
      throw std::runtime_error("Cannot find file: " + filename + "! Please use '-H' option to specify the ARACNE home directory.");
    }

    std::string line;
    double alpha;
    double beta;
    std::getline(infile, line);
    // by pass the lines starting with ">"
    while (strncmp(line.c_str(), ">", 1) == 0) {
      std::getline(infile, line);
    }

    std::istringstream sin(line);
    if ( sin.good() && (sin.peek() != EOF) ) {
      sin >> alpha >> beta;
    }
    else {
      throw std::runtime_error("Configuration file format error: " + filename);
    }

    infile.close();

    p.sigma = alpha * pow(n, beta);
  }


  void findThreshold(int n, Parameter & p)
  {
    std::string filename = p.home_dir + "config_threshold.txt";
    std::ifstream infile(filename.c_str(), ifstream::in);

    if (!infile.good()) {
      throw std::runtime_error("Cannot find file: " + filename + "! Please use '-H' option to specify the ARACNE home directory.");
    }

    std::string line;
    double alpha;
    double beta;
    double gamma;
    std::getline(infile, line);
    // by pass the lines starting with ">"
    while (strncmp(line.c_str(), ">", 1) == 0) {
      std::getline(infile, line);
    }

    std::istringstream sin(line);
    if ( sin.good() && (sin.peek() != EOF) ) {
      sin >> alpha >> beta >> gamma;
    }
    else {
      throw std::runtime_error("Configuration file format error: " + filename);
    }

    infile.close();

    p.threshold = (alpha - log(p.pvalue)) / ((-beta) + (-gamma) * n);
  }


  /**
   * @return Vector with computed MIs for the two genes. The first value (vector[0]) is the real value
   * for the MI and the rest are the bootstrapped ones.
   */
  std::vector<decimal> mutual_information_two_genes(
      //gsl_matrix * data_matrix,
      aracne::Parameter & p,
      int probeId1,
      int probeId2,
      size_t nboot
      )
  {
    using namespace aracne;
    using namespace std;

    try
    {
      checkParameter(p);
    }
    catch ( std::string & s )
    {
      std::cout << s << std::endl;
      exit( 1 );
    }
    displayParameter(p);

    std::vector<decimal> output_vec;
    Microarray_Set data;
    MatrixOp aracne;

    data.read( p.infile );
    //data.read( data_matrix );

#if 0
    std::string msg;
    int probeId1 = data.getProbeId( argv[2] );
    if (probeId1 == -1) {
      msg = "Cannot find marker: ";
      msg.append(argv[2]);
      throw (msg + "!");
    }
    int probeId2 = data.getProbeId( argv[3] );
    if (probeId2 == -1) {
      msg = "Cannot find marker: ";
      msg.append(argv[3]);
      throw (msg + "!");
    }
#endif

    int TYPE;
    if (cmpIgnoreCase( p.algorithm.c_str(), "accurate" ) == 0) {
      TYPE = Mutual_Info::MI_GAUSSIAN;
    }
    else {
      TYPE = Mutual_Info::MI_AVERAGE;
    }

    std::vector<int> low;
    std::vector<int> high;
    std::vector<int> * v = NULL;
    int controlId = -1;
    int nsample = data.Get_Num_Microarrays();

    if (p.controlId != "")
    {
      controlId = data.getProbeId( p.controlId );
      bool isHigh = (p.condition == "+") ? true : false;
      data.getHighLowPercent( p.percent, controlId, low, high );
      v = ( isHigh ) ? & high : & low;
      nsample = v->size();
    }

    std::cout << "Marker No: " << data.Get_Num_Markers() << " (" << data.Get_Num_Active_Markers() << " active)"
              << ", Array No: " << nsample << std::endl;

    if (p.sigma == 99) {
      findKernelWidth( nsample, p );
      std::cout << "Kernel width determined for this dataset: " << p.sigma << std::endl;
    }

    if ( p.threshold == 0 && p.pvalue != 1 ) {
      findThreshold( nsample, p );
      std::cout << "MI threshold determined for p=" << p.pvalue << ": " << p.threshold << std::endl;
    }

    Mutual_Info mi(nsample , p.miSteps, p.sigma, TYPE );

    std::cout << "[NETWORK] Calculating MI ... " << std::endl;
    double mutualinfo = aracne.calculateMI( data, mi, probeId1, probeId2, p.threshold, p.correction, v );
    std::cout << "Probe 1: " << probeId1 << " (" << data.Get_Marker_AffyId(probeId1) << ")" << std::endl;
    std::cout << "Probe 2: " << probeId2 << " (" << data.Get_Marker_AffyId(probeId2) << ")" << std::endl;
    std::cout << "MI = " << mutualinfo << std::endl;

    output_vec.push_back(mutualinfo);

    // perform boostrapping
    if (nboot > 0)
    {
      std::cout << "BOOTSTRAPPING ... " << std::endl;
      if (p.outfile == "") {
          throw std::runtime_error("An output file must be specified!");
      }
      //std::fstream output( p.outfile.c_str(), std::ios::out );
      //output << mutualinfo << endl;

      std::vector<int> bs;
      std::vector<int> * u;
      size_t iter = 0;
      while (iter++ < nboot)
      {
        // sampling with replacement
        data.bootStrap(bs, v);
        u = & bs;
        //output << aracne.calculateMI( data, mi, probeId1, probeId2, p.threshold, p.correction, u ) << endl;
        output_vec.push_back( aracne.calculateMI( data, mi, probeId1, probeId2, p.threshold, p.correction, u ) );
        std::cout << ".";
        if (iter % 80 == 0) {
          std::cout << std::endl;
        }
      }
      std::cout << std::endl;
      //output.flush();
      //output.close();
    }

    return output_vec;
  }


  template <class Graph>
  Graph read_aracne_output_file(
      std::string filename
      )
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_descriptor vertex;

    std::ifstream file(filename.c_str());

    file.exceptions ( std::ifstream::eofbit | std::ifstream::failbit | std::ifstream::badbit );

    std::string line;

    Graph g;

    while (true)
    {
      try
      {
        std::getline(file, line);

        if (line[0] == '>')
          continue;

        std::istringstream iss(line);

        std::string buf;

        iss >> buf; // current_vertex name
        vertex current_vertex = find_vertex_by_name(buf, g);

        if (current_vertex == GraphTraits::null_vertex())
        { // add the vertex if it doesn't exist yet
          current_vertex = boost::add_vertex(g);
          g[current_vertex].name = buf;
        }

        while ( ! (iss >> buf).eof() ) // target_vertex name
        {
          vertex target_vertex = find_vertex_by_name(buf, g);

          if (target_vertex == GraphTraits::null_vertex())
          {
            target_vertex = boost::add_vertex(g);
            g[target_vertex].name = buf;
          }

          iss >> buf; // weight

          if ( ! boost::edge(current_vertex, target_vertex, g).second ) // does not add two equal edges
            boost::add_edge(current_vertex, target_vertex, conan::from_string<decimal>(buf), g);
        }
      }
      catch (std::ios_base::failure & e)
      {
        if (file.good())
          std::cout << "Exception raised: " << e.what() << std::endl;

        break;
      }
    }

    return g;
  }


  /**
   *
   */
  template <class Graph>
  Graph mutual_information(
      //gsl_matrix * data_matrix,
      aracne::Parameter & p
      )
  {
    using namespace aracne;
    using namespace std;

    checkParameter(p);
    displayParameter(p);

    // read the input dataset
    Microarray_Set data;
    Matrix matrix;
    MatrixOp aracne;
    data.read( p.infile );
    //data.read( data_matrix );

    if ( p.mean > 0 || p.cv > 0 ) {
      int ndisabled = data.filter( p.mean, p.cv );
      std::cout << ndisabled << " markers disabled due to lack of dynamic range." << std::endl << std::endl;
    }

    int TYPE;
    if (cmpIgnoreCase( p.algorithm.c_str(), "accurate" ) == 0) {
      TYPE = Mutual_Info::MI_GAUSSIAN;
    }
    else  {
      TYPE = Mutual_Info::MI_AVERAGE;
    }

    std::vector<int> ids;
    std::map<size_t, size_t> transfac;
    std::vector<int> low;
    std::vector<int> high;
    std::vector<int> bs;
    std::vector<int> * v = NULL;
    int controlId = -1;
    int nsample = data.Get_Num_Microarrays();

    if (p.controlId != "")
    {
      controlId = data.getProbeId( p.controlId );

      if (controlId == -1) {
        throw std::runtime_error("Cannot find marker: " + p.controlId + "!");
      }

      bool isHigh = (p.condition == "+") ? true : false;
      data.getHighLowPercent( p.percent, controlId, low, high );
      v = ( isHigh ) ? & high : & low;
      nsample = v->size();
    }

    std::cout << "Marker No: " << data.Get_Num_Markers() << " (" << data.Get_Num_Active_Markers() << " active)"
              << ", Array No: " << nsample << std::endl;

    if ( p.sigma == 99 )
    {
      findKernelWidth( nsample, p );
      std::cout << "Kernel width determined for this dataset: " << p.sigma << std::endl;
    }

    if ( p.threshold == 0 && p.pvalue != 1 )
    {
      findThreshold( nsample, p );
      std::cout << "MI threshold determined for p=" << p.pvalue << ": " << p.threshold << std::endl;
    }

    Mutual_Info mi(nsample , p.miSteps, p.sigma, TYPE);

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
    if (p.hub != "")
    {
      int hubId = data.getProbeId( p.hub );

      if (hubId == -1) {
        throw std::runtime_error("Cannot find the hub probe: " + p.hub + ", nothing to be computed!");
      }

      ids.push_back(hubId);
    }

    if (not( p.subnet.empty() ))
    {
      for (size_t i = 0; i < p.subnet.size(); i++)
      {
        int gid = data.getProbeId( p.subnet[i] );

        if (gid == -1)
        {
          std::cout << "Cannot find probe: " << p.subnet[i] << " in \"" << p.subnetfile
                    << "\" ... ignored." << std::endl;
        }
        else
        {
          ids.push_back(gid);
        }
      }
    }

    if (not( p.tf_list.empty() ))
    {
      for (size_t i = 0; i < p.tf_list.size(); i++)
      {
        int gid = data.getProbeId( p.tf_list[i] );
        if (gid == -1) {
          std::cout << "Cannot find probe: " << p.tf_list[i] << " in \"" << p.annotfile
                    << "\" ... ignored." << std::endl;
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

      if ( p.eps != 1 )
      {
        std::cout << "[NETWORK] Applying DPI ..." << std::endl;
        aracne.reduce( data, matrix, mi, p.eps, p.threshold, p.correction, ids, v, transfac );
      }
    }
    else
    {
      // use resampling data
      if (p.sample > 0)
      {
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
        std::cout << "[NETWORK] Applying DPI ..." << std::endl;
        aracne.reduce( data, matrix, mi, p.eps, p.threshold, p.correction, ids, v, transfac );
      }
    }

    if ( p.outfile == "" ) {
      createOutfileName( p );
    }

    matrix.write(data, ids, p);
    return read_aracne_output_file<Graph>(p.outfile);
  }

}}} // conan::inference::mi

#endif // MUTUAL_INFORMATION_HPP
