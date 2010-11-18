// Copyright (C) 2003 Columbia Genome Center
// All Rights Reserved.
//
// $Id: Microarray_Set.cpp,v 1.3 2005/10/06 02:49:20 kw2110 Exp $
//

#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "Microarray_Set.h"

namespace aracne {

  using namespace std;

  const char* Probe::DELIMINATOR = "\t";
  const char* Microarray_Set::HEADER_DELIMINATOR = "\t";

  /**
    * Dump a single probe into a stream
    */
  std::ostream& operator << (std::ostream& out, const Probe& p)
  {
    out << "(" << p.Get_Value() << ", " << p.Get_PValue() << ")";
    return (out);
  }

  /**
    * Dumps a microarray set into a stream
    */
  std::ostream& operator << (std::ostream& out, const Microarray_Set& ms)
  {
    // first dump headers;
    const char* delim = Microarray_Set::HEADER_DELIMINATOR;
    unsigned int mrkNo = ms.Get_Num_Markers();
    const vector<string>& hdr = ms.Get_Header();
    if (0 != hdr.size())
    {
      std::copy(hdr.begin(), hdr.end(),
          std::ostream_iterator<string>(out,delim));
    }
    else
    {
      // make up a col header
      out << "Id" << delim << "Desc" << delim;
      unsigned int col = 0;
      for (; col< mrkNo-1 ; ++col) {
        out << "exp" << col << delim;
      }
      out << "exp" << col;
    }
    out << "\n";

    // dump markers
    for (unsigned int i = 0; i < mrkNo; ++i)
    {
      const Marker& m = ms.Get_Marker(i);
      out << m.Get_Accnum() << "\t" << m.Get_Label();
      for (int j = 0; j < ms.Get_Num_Microarrays(); j++)
      {
        const Probe& p = ms.Get_Probe(j, i);
        out << Probe::DELIMINATOR << p.Get_Value();
      }

      if (i != mrkNo - 1)
        out << "\n";
    }

    return(out);
  }

  /**
   * Read a matrix data file into a Microarray_Set.
   * The file is consists of sets of lines, that start with a description of a
   * markers (including accession number and descriptive label). The remainder of
   * the lines are probe values (the values only, or (value, pvalue) pairs) for
   * each microarray for the marker that started the line. We assign each
   * marker an ID number according to the line number it appeared on.
   * Conceptually, the file can be thought of as a matrix of probe data in column
   * major order.
   *
   * Kai Wang: currently, the read function will handle the following situation
   * correctly:
   * 1) DOS format vs. UNIX format
   * 2) A blank line at the end of the file
   * 3) A TAB at the end of each line
   *
   * Kai Wang: the program can also automatically tell whether the input file
   * consisits of (value, pvalue) pairs, or single expression values. It will
   * read the data and set up the structure correctly according the format of the
   * file. (If the data contains only expression measurements, a p-value of 0.0
   * will be assigned to each expression value in the data.
   */
  void
  Microarray_Set::read(const string &file) throw(std::runtime_error)
  {
    ifstream in(file.c_str());
    if (!in.good())
    {
      cerr << "Bad File" << endl;
      throw(std::runtime_error("Problem with reading the file."));
    }
    read(in);
  }

  /**
    * Read a matrix from a stream.
    */
  void
  Microarray_Set::read(std::istream& in) throw(std::runtime_error)
  {
    if (!in.good()) {
      throw(std::runtime_error("Bad stream."));
    }

    ios_base::iostate oldState = in.exceptions();
    in.exceptions(ios_base::badbit | ios_base::failbit);

    try
    {
      //read header
      string line;
      getline(in, line);
      istringstream sin(line);
      // std::cout << "Reading the header line ..." << std::endl;

      // newly added begins here
      unsigned int arrayNo = readHeader(sin);
      // std::cout << "no. arrays: " << arrayNo << std::endl;

      /* we need to decide whether the input file contain only expression
       * values or (value, pvalue) pairs; to do this, we read the next
       * line, count how many entries it has.
       * after we decide the file format, we need to read the contents
       * of this line we just extracted again according to the format we
       * determined.
       * So, we did getline(in, line) here once, but we will use the
       * string line twice
       */
      getline(in, line);  // read only once
      // by pass the "Description" lines
      unsigned int bypass_line_cnt = 0;
      while (strncmp(line.c_str(), "Description", 11) == 0) {
        bypass_line_cnt++;
        getline(in, line);
      }
      std::cout << "\n[READ] " << bypass_line_cnt << " Description lines bypassed." << std::endl;

      istringstream pin(line);  // first time used to decide format
      std::vector < std::string > firstprobe;
      std::string token;

      do
      {
        // pin >> token;
        getline(pin, token, '\t');
        firstprobe.push_back(token);
      }
      while ( pin.good() && (pin.peek() != '\015') && (pin.peek() != EOF) );
      // while (pin.good() && (pin.get() != '\015') && ( (pin.peek() != EOF) && (pin.peek() != '\015') ) );

      unsigned int valueNo = firstprobe.size() - 2;
      // std::cout << "no. entries in the data: " << valueNo << std::endl;

      if (arrayNo == valueNo)
      {
        // single expression values
        std::cout << "[READ] P-value columns not found." << std::endl;
        unsigned int proben = 0;  // probe number

        istringstream fin(line);  // second time used to read the content
        string accnum, label;

        getline(fin, accnum, '\t');
        getline(fin, label, '\t');
        // std::cout << accnum << '\t' << label << std::endl;
        markerMap[accnum] = proben;
        Set_Marker(proben, Marker(proben, accnum, label));
        readMarkerNoPvalue(fin, proben);
        // std::cout << std::endl;
        proben++;

        // enter a loop to read the rest of probes (i.e. lines)
        // make sure the next line is not a blank line
        while ( in.good() && (in.peek() != EOF) )
        {
          getline(in, line);
          istringstream sin(line);
          getline(sin, accnum, '\t');
          getline(sin, label, '\t');
          // std::cout << accnum << '\t' << label << '\t';
          markerMap[accnum] = proben;
          Set_Marker(proben, Marker(proben, accnum, label));
          if ( readMarkerNoPvalue(sin, proben) != arrayNo )
          {
            std::ostringstream s;
            s << "Incorrect data format at line no: " << (proben+2);
            throw std::runtime_error(s.str());
          }
          // std::cout << std::endl;
          proben++;
        }
      }
      else if (arrayNo * 2 == valueNo)
      {
        // (value, pvalue) pairs
        std::cout << "[READ] (value, p-value) pairs found." << std::endl;
        unsigned int proben = 0;  // probe number

        istringstream fin(line);  // second time used to read the content
        string accnum, label;

        getline(fin, accnum, '\t');
        getline(fin, label, '\t');
        // std::cout << accnum << '\t' << label << std::endl;
        markerMap[accnum] = proben;
        Set_Marker(proben, Marker(proben, accnum, label));
        readMarkerWithPvalue(fin, proben);
        // std::cout << std::endl;
        proben++;

        // enter a loop to read the rest of probes (i.e. lines)
        // make sure the next line is not a blank line
        while ( in.good() && (in.peek() != EOF) )
        {
          getline(in, line);
          istringstream sin(line);
          getline(sin, accnum, '\t');
          getline(sin, label, '\t');
          // std::cout << accnum << '\t' << label << std::endl;
          markerMap[accnum] = proben;
          Set_Marker(proben, Marker(proben, accnum, label));
          if ( readMarkerWithPvalue(sin, proben) != arrayNo )
          {
            std::ostringstream s;
            s << "Incorrect data format at line no: " << (proben+2);
            throw std::runtime_error(s.str());
          }
          // std::cout << std::endl;
          proben++;
        }
      }
      else
      {
        // wrong format
        throw std::runtime_error("Incorrect file format: header line doesn't match the rest data.");
      }
    }
    catch (ios_base::failure& f)
    {
      throw std::runtime_error("Could not read data (last line empty?).");
    }

    in.exceptions(oldState);	//reset excep state
  }

  /**
    * Read a marker line with (value, pvalue) format
    */
  unsigned int
  Microarray_Set::readMarkerWithPvalue(std::istream& sin,
      const unsigned proben) throw(std::runtime_error)
  {
    if (!sin) {
      throw std::runtime_error("Bad stream (end of stream?).");
    }

    std::string sPVal;
    double val, pval;
    unsigned int markern = 0; // number of expression values of each probe
    sin.exceptions(ios_base::badbit | ios_base::failbit);
    try
    {
      do
      {
        sin >> val >> sPVal;
        switch(sPVal[0])
        {
          case 'A':
            pval = 0.7;
            break;
          case 'M':
            pval = 0.5;
            break;
          case 'P':
            pval = 0.1;
            break;
          default:
            pval = atof(sPVal.c_str());
            break;
        }
        Set_Probe(markern, proben, Probe(val, pval));
        markern++;
        // std::cout << "*" << val << '\t' << "*" << sPVal << '\t';
      }
      while ( sin.good() && (sin.get() != '\015') && ( (sin.peek() != EOF) && (sin.peek() != '\015') ) );
      /* the get() function mainly get rid of the TAB, since they will
       * be bypassed by >> operator anyway.
       * the peek() function is to see whether the next charaacter is
       * carriage return (for DOS format) or EOF
       *
       * end of line format under DOS: '\015' '\012'
       * end of line format under UNIX: '\012'
       * so, when reading a DOS file under UNIX, an extra '\015' will
       * occur at the end of each line by getline().
       * '\015' - Carriage return
       * '\012' - line feed
       */

    }
    catch (ios_base::failure& f)
    {
      std::ostringstream s;
      s << "Could not read data at line no: " << (proben+2);
      std::cout << s.str() << std::endl;
      throw std::runtime_error(s.str());
    }

    return markern;
  }

  /**
    * Read a marker line with no pvalues
    */
  unsigned int
  Microarray_Set::readMarkerNoPvalue(std::istream& sin,
                const unsigned proben) throw(std::runtime_error)
  {
    if (!sin) {
      throw std::runtime_error("Bad stream. (end of stream?)");
    }

    double val;
    unsigned int markern = 0; // number of expression values of each probe
    sin.exceptions(ios_base::badbit | ios_base::failbit);

    try
    {
      do
      {
        sin >> val;
        Set_Probe(markern, proben, Probe(val, 0.0));
        markern++;
        // std::cout << "*" << val << '\t';
      }
      while ( sin.good() && (sin.get() != '\015') && ( (sin.peek() != EOF) && (sin.peek() != '\015') ) );
      // the get() function mainly get rid of the TAB, since they will be bypassed by >> operator anyway.
      // the peek() function is to see whether the next charaacter is carriage return (for DOS format)
      // or EOF
    }
    catch (ios_base::failure& f)
    {
      std::ostringstream s;
      s << "Could not read data at line no: " << proben;
      std::cout << s.str() << std::endl;
      throw std::runtime_error(s.str());
    }

    return markern;
  }


  unsigned int
  Microarray_Set::readHeader(std::istream& in)throw(std::runtime_error)
  {
    if (!in) {
      throw std::runtime_error("Bad stream. (end of stream?)");
    }

    in.exceptions(ios_base::badbit | ios_base::failbit);
    try
    {
      do
      {
        std::string hdr;
        // here has to be implemented this way, instead of the way
        // as in the readMarker function, because some headers
        // contains white space characters in them.
        getline(in, hdr, '\t');
        Set_ColHeader(header.size(), hdr);
      }
      while ( in.good() && (in.peek() != '\015') && (in.peek() != EOF) );
      // std::cout << std::endl;
    }
    catch (ios_base::failure& f)
    {
      throw std::runtime_error("Error while reading file headers(win/*nix/mac end of line?).");
    }

    return header.size() - 2;
  }


  void
  Microarray_Set::Set_ColHeader(unsigned int i, const string& hdr)
  {
    while (i >= header.size()) {
      header.push_back("");
    }
    header[i] = hdr;
  }

  /**
    * filters the Microarray_set by disabling all the markers whose mean
    * is smaller than minMean and standard deviation is smaller than the mean * ratio
    */
  int
  Microarray_Set::filter(vector<int> &ids, double min_mean, double min_sigma, int ctlid)
  {
    double iNo = double(uarrays.size());
    int	ndisabled = markerset.size();

    for (int i = 0; i < ndisabled; i++)
      markerset[i].Disable();

    for (size_t m = 0; m < markerset.size(); m++)
    {
      if ((int)m == ctlid)
        continue;

      double nx = 0;
      double nxx = 0;

      for (vector<int>::iterator i = ids.begin(); i != ids.end(); i++)
      {
        double v = Get_Value(*i, m);
        nx += v;
        nxx += v * v;
      }

      double mean = nx / iNo;
      // double stdev = sqrt((nxx - mean * mean) / (iNo - 1));
      double stdev = sqrt( (iNo*nxx - nx*nx)/(iNo*iNo) );

      if ((mean >= min_mean) && (stdev >= mean * min_sigma))
      {
        markerset[m].Enable();
        ndisabled--;
      }
    }

    return(ndisabled);
  }

  //
  // This version generates the control set.
  //
  int Microarray_Set::filter(double min_mean, double min_sigma, int ctlid)
  {
    vector<int> v;
    unsigned int i;
    int r;

    for (i = 0; i < uarrays.size(); i++)
      v.push_back(i);

    r = filter(v, min_mean, min_sigma, ctlid);

    return(r);
  }

  void Microarray_Set::computeMarkerVariance(vector < int > * arrays)
  {
    for (size_t i=0; i < markerset.size(); i++)
    {
      double var = variance( i, arrays );
      // std::cout << var << std::endl;
      markerset[i].Set_Var( var );
    }
  }

  // compute the variance of a gene expression vector in log-space
  double Microarray_Set::variance(unsigned int m, vector < int > * arrays)
  {
    double var;

    if ( arrays ) // compute variance of a gene's expression only within selectd samples
    {
      size_t n = arrays->size();
      double s = 0;
      double ss = 0;

      for (size_t i = 0; i < n; i++)
      {
        double v = log(max(Get_Value((* arrays)[i], m), 0.1));
        s += v;
        ss += v*v;
      }
      var = (ss - s * s / n) / (n - 1);
    }
    else // compute variance across all arrays
    {
      size_t n = uarrays.size();
      double s = 0;
      double ss = 0;

      for (size_t i = 0; i < n; i++)
      {
        double v = log(max(Get_Value(i, m), 0.1));
        s += v;
        ss += v*v;
      }
      var = (ss - s * s / n) / (n - 1);
    }

    return var;
  }

  double Microarray_Set::getMarkerVariance(int i)
  {
    return markerset[i].Get_Var();
  }


  void Microarray_Set::Set_Marker(unsigned int i, const Marker& m)
  {
    while (i >= markerset.size()) {
      Marker v;
      markerset.push_back(v);
    }
    markerset[i] = m;
  }

  void Microarray_Set::Set_Microarray(unsigned i, const Microarray& a)
  {
    while (i >= uarrays.size()) {
      Microarray v;
      uarrays.push_back(v);
    }
    uarrays[i] = a;
  }


  void Microarray_Set::Set_Probe(unsigned int i, unsigned int j, Probe p)
  {
    while (i >= uarrays.size()) {
      Microarray v;
      uarrays.push_back(v);
    }
    while (j >= uarrays[i].size()) {
      Probe v;
      uarrays[i].push_back(v);
    }
    uarrays[i] [j] = p;
  }

  void Microarray_Set::Set_Value(unsigned int i, unsigned int m, double v) throw(std::runtime_error)
  {
    while (i >= uarrays.size()) {
      Microarray v;
      uarrays.push_back(v);
    }
    while (m >= uarrays[i].size()) {
      Probe v;
      uarrays[i].push_back(v);
    }
    uarrays[i] [m].Set_Value(v);
  }

  void Microarray_Set::Set_PValue(unsigned int i, unsigned int m, double v) throw(std::runtime_error)
  {
    while (i >= uarrays.size()) {
      Microarray v;
      uarrays.push_back(v);
    }
    while (m >= uarrays[i].size()) {
      Probe v;
      uarrays[i].push_back(v);
    }
    uarrays[i] [m].Set_PValue(v);
  }

  int Microarray_Set::get_Id(std::string & a)
  {
    for (unsigned i = 0; i < markerset.size(); i++) {
      if (markerset[i].hasAccession(a))
        return i;
    }
    return -1;
  }

  const Marker& Microarray_Set::Get_Marker(unsigned int i)const throw(std::runtime_error)
  {
    if (i >= markerset.size())
      throw(std::runtime_error("Out of bounds!"));
    return markerset[i];
  }

  Marker& Microarray_Set::Get_Marker(unsigned int i) throw(std::runtime_error)
  {
    if (i >= markerset.size())
      throw(std::runtime_error("Out of bounds!"));
    return markerset[i];
  }

  const std::string Microarray_Set::Get_Marker_AffyId(unsigned int i) throw(std::runtime_error)
  {
    if (i >= markerset.size())
      throw(std::runtime_error("Out of bounds!"));
    return markerset[i].Get_Accnum();
  }


  int Microarray_Set::Get_Num_Active_Markers()
  {
    int num_active_markers = 0;
    for (unsigned int i = 0; i < markerset.size(); i++)
    {
      if ( markerset[i].Enabled() )
      {
        num_active_markers++;
      }
    }
    return num_active_markers;
  }


  const Probe & Microarray_Set::Get_Probe(unsigned int i, unsigned int m) const throw(std::runtime_error)
  {
    if (i >= uarrays.size())
      throw(std::runtime_error("Out of bounds!"));
    if (m >= uarrays[i].size())
      throw(std::runtime_error("Out of bounds!"));

    return uarrays[i][m];
  }

  const Microarray & Microarray_Set::Get_Microarray(unsigned int i) const throw(std::runtime_error)
  {
    if (i >= uarrays.size())
        throw(std::runtime_error("Out of bounds!"));
    return uarrays[i];
  }

  const double Microarray_Set::Get_Value(unsigned int i, unsigned int m) const throw(std::runtime_error)
  {
    if (i >= uarrays.size())
      throw(std::runtime_error("Out of bounds!"));
    if (m >= uarrays[i].size())
      throw(std::runtime_error("Out of bounds!"));

    return uarrays[i][m].Get_Value();
  }

  const double Microarray_Set::Get_PValue(unsigned int i, unsigned int m) const throw(std::runtime_error)
  {
    if (i >= uarrays.size())
      throw(std::runtime_error("Out of bounds!"));
    if (m >= uarrays[i].size())
      throw(std::runtime_error("Out of bounds!"));
    return uarrays[i] [m].Get_PValue();
  }

  int Microarray_Set::getProbeId(std::string label)
  {
    int controlId;// = -1;
    bool isNumber = true;
    for (unsigned int i = 0; i < label.length(); i++)
    {
      if (!isdigit(label.at(i))) {
        isNumber = false;
        break;
      }
    }

    if(isNumber)
      controlId = atoi(label.c_str());
    else
     controlId = get_Id(label);

    return controlId;
  }

  void Microarray_Set::getHighLowPercent( double x, int mId, vector < int > & lower, vector < int > & upper )
  {
    vector < ArrayValuePair > sortArray;
    SortIncreasing_ArrayValuePair sorter;
    int idNo = uarrays.size();
    for ( int iId = 0; iId < idNo; iId++ )
    {
      ArrayValuePair pair;
      pair.set( iId, Get_Value( iId, mId ) );
      sortArray.push_back( pair );
    }
    sort( sortArray.begin(), sortArray.end(), sorter );
    int idPercNo = ( int )( ( double )idNo * x );
    for ( int iId = 0; iId < idPercNo; iId++ )
    {
      lower.push_back( sortArray[iId].getId() );
      upper.push_back( sortArray[idNo - idPercNo - 1 + iId].getId() );
    }
  }

  void Microarray_Set::getRandomSubsamples( double x, int mId, vector < int > & lower, vector < int > & upper )
  {
    // mId is not always used. in some case it is needed
    // to make sure that the high and low subsets are the same for the
    // high and low runs
    // srand ( mId );
    // srand ( time(NULL) );

    vector < ArrayValuePair > sortArray;
    SortIncreasing_ArrayValuePair sorter;
    int idNo = uarrays.size();

    for ( int iId = 0; iId < idNo; iId++ )
    {
      ArrayValuePair pair;
      // rand() generates a random number between 0 and 32767
      pair.set( iId, rand() );
      sortArray.push_back( pair );
    }
    sort( sortArray.begin(), sortArray.end(), sorter );
    int idPercNo = ( int )( ( double )idNo * x );
    for ( int iId = 0; iId < idPercNo; iId++ )
    {
      lower.push_back( sortArray[iId].getId() );
      upper.push_back( sortArray[idNo - idPercNo - 1 + iId].getId() );
    }
  }


  void Microarray_Set::bootStrap( vector < int > & boot, vector < int > * arrays )
  {
    if (not(boot.empty()))
    {
      boot.clear();
    }
    // if a subset of arrays is selected, bootstrap only from this subet of arrays
    if ( arrays )
    {
      int idNo = arrays->size();
      for ( int iId = 0; iId < idNo; iId++ )
      {
        // a random number between [0, idNo-1]
        int r = rand()%idNo;
        boot.push_back( (* arrays)[r] );
      }
    }
    else
    {
      int idNo = uarrays.size();
      for ( int iId = 0; iId < idNo; iId++ )
      {
        // a random number between [0, idNo-1]
        int r = rand()%idNo;
        boot.push_back( r );
      }
    }
  }

  void Microarray_Set::addNoise()
  {
     int idNo = uarrays.size();
     int mNo = markerset.size();
     for ( int id = 0; id < idNo; id++)
     {
       for ( int mid = 0; mid < mNo; mid++)
       {
          double noise = (((double)rand()) / (double)RAND_MAX) * 1e-7;
          double value = Get_Value(id, mid);
          double changedValue = value + noise;
          Set_Value(id, mid, changedValue);
       }
     }
  }

  void Microarray_Set::shuffleGene(unsigned int idx)
  {
    for (size_t i=uarrays.size(); i>1; i--)
    {
      // pick a random sample between [0, i-1]
      size_t rind = rand()%i;
      // swap probe rind and (i-1)
      double temp = Get_Value(i-1, idx);
      Set_Value(i-1, idx, Get_Value(rind, idx));
      Set_Value(rind, idx, temp);
    }
  }
  /* proof: suppose a vector with n elements, with indices 1,2,3, ..., n,
   * the the element at i is shuffled to poistion j can happen only if it is not
   * chosen on the first n - j steps of the algorithm and is chosen on the step
   * n - j + 1. Therefore, the probability of this event is:
   * prob = [(n-1)/n]*[(n-2)/(n-1)]*...*[j/(j+1)]*(1/j) = 1/n
   */

} // aracne
