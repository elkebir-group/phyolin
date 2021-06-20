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
   * @param threshold the estimated false negative rate
   */

  PhyolinCP(const SingleCellMatrix B, double fp, int time, double threshold);

  /*! 
   * @return get the number of flips made by Phyolin
   */
  int getObjective()
  {
    return _objValue;
  }

  /*! 
   * @return get Phyolin estimate false negative rate
   */
  double getEstFN()
  {
    return _estFN;
  }

  /*! 
   * @return a boolean indicating if a solution has been found
   */
  int getSolveStatus()
  {
    return _solve;
  }

  /*!
   * @brief write the flipped output data to a file
   * @param filename the string path to the file where the output file should be written
   * @param delim a string delimiter for the output file
   */
  void write_csv(std::string filename, std::string delim);

  /*!
   * @brief write a report on the input parameters and output values to a file
   * @param filename the string path to the file where the output file should be written
   */
  void write_counts(std::string filename);

  /*!
   * @brief Computes the likelihood of the observed data under a linear phylogeny
   */
  double getLikelihood();

private:
  /*!
 *
 * @brief initalize and solve a constraint programming model to find the minimum
 * number of flips such that the input single-cell matrix represents a linear phylogeny
 * 
 */
  void init();

  SingleCellMatrix _B; /*!< a single cell matrix object (shape: cells x mutations) with input data */

  IloEnv _env; /*!< environment variable */

  IloModel _model; /*!< model variable */

  IloCP _cp; /*!< solver  variable */

  IloArray<IloIntVarArray> _x; /*!< decision variable indicating value after flipping shape( cells x mutations) */

  IloIntVarArray _c; /*!< decision variable for permutation of columns */

  IloIntVar _y; /*!< decision variable for the number of false positive flips */

  bool _solve; /*!< value of the solve status of the model */

  double _objValue; /*!< the number of flips made by Phyolin */

  double _fp; /*!< false positive rate */

  int _time; /*!< time limit for the solver */

  double _threshold;

  int _fpCounts; /*!< the number of false positives afters flipping */

  int _cells; /*!< the number of cells in the input data */

  int _sites; /*!< the number of cells in the input data */

  double _estFP; /*!< the estimated false positive rate after flipping */

  double _estFN; /*!< the estimated false negative rate after flipping */

  int _inputZeros; /*!< the total number of 0's in the input data */

  int _inputOnes; /*!< the total number of 1's in the input data */

  int _inputMissing; /*!< the total number of -1's in the input data representing missing data */

  int _outputZeros; /*!< the total number of 0's after flipping */

  int _outputOnes; /*!< the total number of 1's after flipping */

  int _outputMissing; /*!< the total number of -1's after flipping (should equal number of _inputMissing) */

  int N_01; /*!< the total number of flips from 0 to 1 (false negatives) */

  int N_11; /*!< the total number of non-flips from 1 to 1 (true positives) */

  int N_00; /*!< the total number of non-flips from 0 to 0  (true negatives) */

  int N_10; /*!< the total number of non-flips from 1 to 0  (false positives) */

  std::vector<std::vector<int> > _Bout; /*!< the flipped Phyolin output matrix (shape: cells x sites) */
};

#endif // PHYOLINCP_H
