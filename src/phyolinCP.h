/*!
 * @author Leah L. Weber
 * @version 1.0
 * @date June 14, 2021
 */

#ifndef PHYOLINCP_H
#define PHYOLINCP_H


#include <ilcp/cp.h>
#include "singlecellmatrix.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

class PhyolinCP
{
public:
  /*!
   * @brief constructor for an instance of Phyolin
   * @param B a single cell matrix object of dimension cells x mutations containing
   * the binary input data
   * @param fp the false positive rate for the input data
   * @param time the time limit for the CP solver
   * @param threshold 
   */
   
  PhyolinCP(const SingleCellMatrix B, double fp, int time, double threshold);
   
  /*! 
   * @return get the number of flips made by Phyolin
   */
  int getObjective() {
      return _objValue;
  }

  /*! 
   * @return a boolean indicating if a solution has been found
   */
  int getSolveStatus() {
      return _solve;
  }

  /*!
   * @brief write the flipped output data to a file
   * @param filename the string path to the file where the output file should be written
   * @param delim a string delimiter for the output file
   */
  void write_csv(std::string filename, std::string delim);
  
  /*!
   * @brief 
   * @param filename the string path to the file where the output file should be written
   * @param delim a string delimiter for the output file
   */
  void write_counts(std::string filename);


  double getLikelihood();
   
  
private:

  void init();
  
  SingleCellMatrix _B; /*!< a single cell matrix object (shape: cells x mutations) with input data */
 
  IloEnv _env; /*!< environment variable */

  IloModel _model; /*!< model variable */

  IloCP _cp; /*!< solver  variable */

  IloArray <IloIntVarArray> _x; /*!< decision variable indicating value after flipping shape( cells x mutations) */

  IloIntVarArray _c; /*!< decision variable for permutation of columns */

  IloIntVar _y; /*!< decision variable for the number of false positive flips */

  bool _solve; /*!< value of the solve status of the model */

  double _objValue; /*!< the number of flips made by Phyolin */

  double _fp; /*!< false positive rate */

  int _time; /*!< time limit for the solver */

  double _threshold;

  int _fpCounts;
  
  int _cells;

  int _sites;

  double _estFP;

  double _estFN;

  int _inputZeros;

  int _inputOnes;

  int _inputMissing;

  int _outputZeros;
  int _outputOnes;
  int _outputMissing;
  int N_01;
  int N_11;
  int N_00;
  int N_10;
  // double _doublet;


  std::vector<std::vector<int>> _Bout;

};

#endif // SETCOVERILP_H
