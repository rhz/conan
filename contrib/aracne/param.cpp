/*
 * Copyright (C) 2003, 2004  Columbia Genome Center
 * All Rights Reserved.
 *
 * param.cc --
 *
 * $Id: param.cpp,v 1.3 2005/12/13 23:08:02 wl2131 Exp $
 */

#include <fstream>
#include <iostream>
#include <sstream>
#include "param.h"

namespace aracne {

  using namespace std;

  string getFileName(std::string matrixName)
  {
    unsigned int len = matrixName.length();
    string::size_type b   = matrixName.find_last_of("/");
    if(b == std::string::npos) {
      b = matrixName.find_last_of("\\");
    }
    string dirname("");
    string basename(matrixName);
    if ((b != matrixName.npos) && (b < len)) {
      // Extract the directory and the filename if path is included
      basename = matrixName.substr(b + 1, len);
      dirname = matrixName.substr(0, b) + std::string("/");
    }
    string::size_type c = basename.find_last_of(".");
    if(c != basename.npos) {
      basename = basename.substr(0,c);
    }
    string filename = dirname + basename;
    return filename;
  }

  void checkParameter(Parameter &p)
  {
    if (p.infile == "") {
      throw runtime_error(string("No input file specified!"));
    }

    if ( (p.subnetfile != "") && (p.hub != "") )
    {
      throw runtime_error(string("Either supply one hub gene by '-h' or multiple genes in a file by '-l', but not both!"));
    }

    if (not((p.condition == "+") || (p.condition == "-") || (p.condition == ""))) {
      throw runtime_error(string("Condition can be either '+' or '-', but no else!"));
    }

    if ((p.condition == "+") || (p.condition == "-")) {
      if (p.controlId == "") {
        throw runtime_error(string("Control gene ID must be specified using '-c'!"));
      }
    }

    if ( (p.sigma != 99) && (p.sigma <= 0 || p.sigma >= 1) ) {
      throw runtime_error(string("Kernel width '-k' must be within (0,1)!"));
    }

    if (p.threshold < 0) {
      throw runtime_error(string("MI threshold '-t' must be non-negative!"));
    }

    if ( p.threshold > 0 && p.pvalue != 1 ) {
      cout << "P-value will not be used, since a threshold has been specified." << endl;
    }

    if ( p.pvalue <= 0 || p.pvalue >1 ) {
      throw runtime_error(string("P-value '-p' must be in the range (0,1]!"));
    }

    if (p.eps < 0 || p.eps>1) {
      throw runtime_error(string("DPI tolerance '-e' must be within [0,1]!"));
    }

    if (p.percent <= 0 || p.percent >= 1) {
      throw runtime_error(string("Percentage microarray must be within (0,1)!"));
    }

    if (p.miSteps <= 0) {
      throw runtime_error(string("Number of bins for fast method '-b' must be postive!"));
    }

    if (p.mean < 0) {
      throw runtime_error(string("Gene filter mean must be non-negative!"));
    }

    if (p.cv < 0) {
      throw runtime_error(string("Gene filter cv (coefficient of variance) must be non-negative!"));
    }

    if (p.correction < 0) {
      throw runtime_error(string("Array measurement noise level '-n' must be non-negative!"));
    }

    if (not(cmpIgnoreCase(p.verbose.c_str(),"on")==0 || cmpIgnoreCase(p.verbose.c_str(), "off")==0)) {
      throw runtime_error(string("Verbose '-v' can be either 'on' or 'off', but no else!"));
    }

    if (not(cmpIgnoreCase(p.algorithm.c_str(),"accurate")==0 || cmpIgnoreCase(p.algorithm.c_str(),"fast")==0)) {
      throw runtime_error(string("Supported algorithm '-a' is either 'accurate' or 'fast', but no else!"));
    }

    if (p.home_dir != "./") {
      unsigned int  len = p.home_dir.length();
      string::size_type b = p.home_dir.find_last_of("/");
      if ( b == p.home_dir.npos || b != (len-1) ) {
        p.home_dir = p.home_dir + string("/");
      }
    }
  }

  void createOutfileName(Parameter &p)
  {
    string filename = getFileName(p.infile);
    string temp;
    if (p.hub !="") {
      filename = filename + "_h" + p.hub;
    }
    if (p.controlId != "") {
      temp = (p.condition == "+") ? "H" : "L";
      filename = filename + "_c" + p.controlId + temp;
    }
    char kernel[10];
    sprintf(kernel, "%0.3g", p.sigma);
    filename = filename + "_k" + kernel;

    if (p.threshold > 0) {
      char threshold[10];
      sprintf(threshold, "%0.2g", p.threshold);
      filename = filename + "_t" + threshold;
    }

    if (p.eps < 1) {
      char eps[10];
      sprintf(eps, "%0.2g", p.eps);
      filename = filename + "_e" + eps;
    }

    if (p.sample > 0) {
      char sample[10];
      sprintf(sample, "%03i", p.sample);
      filename = filename + "_r" + sample;
    }
    filename = filename + ".adj";
    p.outfile = filename;
  }


  void displayParameter(Parameter &p)
  {
    cout << endl;
    if (p.adjfile != "") {
      cout << "[PARA] Pairwise MI will be read from: " << p.adjfile << endl;
      cout << "[PARA] Input file:    " << p.infile << endl;
      cout << "[PARA] Output file:   " << p.outfile << endl;
      if (p.threshold > 0) {
        cout << "[PARA] MI threshold:  " << p.threshold << endl;
      }
      else {
        cout << "[PARA] MI P-value:    " << p.pvalue << endl;
      }
      cout << "[PARA] DPI tolerance: " << p.eps << endl;
      return;
    }

    std::cout << "[PARA] Input file:    " << p.infile << endl;
    std::cout << "[PARA] Output file:   " << p.outfile << endl;
    if (cmpIgnoreCase(p.algorithm.c_str(),"accurate")==0)
    {
      cout << "[PARA] Algorithm:     " << p.algorithm << endl;
      if (p.sigma == 99) {
        cout << "[PARA] Kernel width:  determined by program" << endl;
      }
      else {
        cout << "[PARA] Kernel width:  " << p.sigma << endl;
      }
    }
    else
    {
      cout << "[PARA] Algorithm:     " << p.algorithm << endl;
      cout << "[PARA] No. bins:      " << p.miSteps << endl;
    }
    if (p.threshold > 0) {
      cout << "[PARA] MI threshold:  " << p.threshold << endl;
    }
    else {
      cout << "[PARA] MI P-value:    " << p.pvalue << endl;
    }

    cout << "[PARA] DPI tolerance: " << p.eps << endl;

    if (p.correction > 0) {
      cout << "[PARA] Correction for MI estimation (array noise level: " << p.correction << ")" << endl;
    }

    if (p.subnetfile != "")
    {
      int cnt = readProbeList(p.subnetfile, p.subnet);
      cout << "[PARA] Subset of probes to reconstruct: " << p.subnetfile << " (" << cnt << ")" << endl;
    }
    if (p.hub != "")
    {
      cout << "[PARA] Hub probe to reconstruct: " << p.hub << endl;
    }
    if (p.controlId != "")
    {
      cout << "[PARA] Control gene:  " << p.controlId << endl;
      cout << "[PARA] Condition:     " << p.condition << endl;
      cout << "[PARA] Percentage:    " << p.percent << endl;
    }
    if (p.annotfile != "")
    {
      int cnt = readProbeList(p.annotfile, p.tf_list);
      cout << "[PARA] TF annotation list: " << p.annotfile << " (" << cnt << ")" << endl;
    }
    if (p.mean != 0 || p.cv != 0)
    {
      cout << "[PARA] Filter mean:   " << p.mean << endl;
      cout << "[PARA] Filter CV:     " << p.cv << endl;
    }
  }


  /* read the list of nodes to be included in constructing a subnetwork from a
   * file.
   */
  int readProbeList(string infile, vector<string> & probe_list)
  {
    ifstream in(infile.c_str());
    if(!in.good()) {
      throw runtime_error(string("Problem with reading the file: " + infile));
    }

    string line;
    size_t lnum = 0;
    while ( in.good() && (in.peek() != EOF) && (in.peek() != '\012') ) {
      getline(in, line);
      istringstream sin(line);
      string gid;
      sin >> gid;
      probe_list.push_back(gid);
      lnum++;
    }
    return lnum;
  }


  int cmpIgnoreCase(const char* a, const char* b)
  {
    int i = 0;
    for(; toupper(a[i])==toupper(b[i]); ++i){
      if(a[i] == '\0') return 0;
    }
    return(a[i]-b[i]);
  }


  // Default parameters
  // mi threshold
  const double Parameter::default_threshold  = 0.0;
  // p-value for mi threshold
  const double Parameter::default_pvalue  = 1.0;
  // DPI tolerance
  const double Parameter::default_eps = 1.0;
  // kernel width
  const double Parameter::default_sigma   = 99;
  // mi step size
  const int Parameter::default_miSteps    = 6;
  // sample number
  const int Parameter::default_sample   = 0;
  // high/low percentage
  const double Parameter::default_percent = 0.35;
  // filter mean
  const double Parameter::default_mean    = 0.0;
  // filter standard deviation
  const double Parameter::default_cv   = 0.0;
  // array measurement noise level
  const double Parameter::default_correction = 0.0;

} // aracne
