/*
 * Copyright (C) 2003, 2004  Columbia Genome Center
 * All Rights Reserved.
 *
 * param.h --
 *
 * $Id: param.h,v 1.3 2005/12/13 23:08:02 wl2131 Exp $
 */
#ifndef PARAM_H__
#define PARAM_H__

#include <string>
#include <vector>
#include <stdexcept>

namespace aracne {

  using namespace std;

  /**
    * Parameters for the algorithm
    */
  struct Parameter{
    // Default parameters
    static const double default_threshold;
    static const double default_pvalue;
    static const double default_eps;
    static const double default_sigma;
    static const int default_miSteps;
    static const int default_sample;
    static const double default_percent;
    static const double default_mean;
    static const double default_cv;
    static const double default_correction;

    double threshold;   // mi threshold
    double pvalue;
    double eps;         // DPI tolerance
    double sigma;       // gaussian kernel width
    int miSteps;        // mi step size
    int sample;         // sample number
    double percent;     // high/low percentage for the conditional analysis
    double mean;        // filter mean
    double cv;          // filter coefficient of variance
    double correction;  // coorection for noise
    string verbose;
    string algorithm;
    string infile;
    string outfile;
    string adjfile;
    string hub;
    string subnetfile;
    string annotfile;
    string controlId;
    string condition;
    string home_dir;
    vector < string > subnet;
    vector < string > tf_list;

    Parameter():threshold(default_threshold), pvalue(default_pvalue), eps(default_eps),
                  sigma(default_sigma), miSteps(default_miSteps),
                  sample(default_sample), percent(default_percent),
                  mean(default_mean), cv(default_cv), correction(default_correction),
                  verbose("off"), algorithm("accurate"), infile(""), outfile(""),
                  adjfile(""), hub(""), subnetfile(""), annotfile(""), controlId(""),
                  condition(""), home_dir("./")
      {
      }
  };

  // service functions
  string getFileName(std::string matrixName);
  void checkParameter(Parameter &p);
  void displayParameter(Parameter &p);
  int readProbeList(std::string subnetfile, std::vector<std::string> & subnet);
  int cmpIgnoreCase(const char* a, const char* b);
  void createOutfileName(Parameter &p);

} // aracne

#endif
