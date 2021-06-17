/*!
 * @author Leah L. Weber
 * @version 1.0
 * @date June 14, 2021
 */

#include "phyolinCP.h"
#include <vector>
#include <array>
#include <cstdlib>
#include <lemon/arg_parser.h>
#include <math.h>

PhyolinCP::PhyolinCP(const SingleCellMatrix B, double fp, int time, double threshold)

    : _env(), _model(_env), _cp(_model), _fp(fp), _B(B), _cells(B._cells), _sites(B._sites), N_11(0), N_00(0), N_01(0), N_10(0), _Bout(_cells, std::vector<int>(_sites)), _inputOnes(0), _inputZeros(0), _outputZeros(0), _outputOnes(0), _inputMissing(0), _outputMissing(0), _estFN(0), _estFP(0), _fpCounts(0), _time(time), _threshold(threshold), _solve(false)
{
  init();
}

void PhyolinCP::init()
{

  std::cout << "Initializing Phyolin" << std::endl;

  IloInt rows = _cells;
  IloInt cols = _sites;
  std::cout << "cells: " << rows << std::endl;
  std::cout << "mutations: " << cols << std::endl;

  /*! Initialize the decision variables */
  typedef IloArray<IloIntVarArray> IloIntVarArray2;
  IloIntVarArray2 _x(_env, rows);

  for (IloInt i = 0; i < rows; i++)
  {
    _x[i] = IloIntVarArray(_env, cols, -1, 1);
  }
 


  _c = IloIntVarArray(_env, cols, 0, cols - 1);

  /*! Add constraints to the model */
  for (int i = 0; i < rows; i++)
  {

    for (int j = 0; j < cols; j++)

      if (_B.is(i, j) == 1)
      {
        _inputOnes++;
        _model.add(_x[i][j] >= 0); /*!< constraint ensures input 1s must be either 0 or 1*/
      }
      else if (_B.is(i, j) == -1)
      {
        _inputMissing++;
        _model.add(_x[i][j] == -1); /*!< constraint prevents missing data from being modified*/
      }
      else
      {
        _model.add(_x[i][j] >= 0); /*!< constraint ensures input 0s must be either 0 or 1*/
        _inputZeros++;
      }
  }

   /*! Calculate the z value to determine the budget for false positive flips
   */

  int total_fps = ceil(_fp * _inputOnes);
  _y = IloIntVar(_env, 0, total_fps);

  // int total_fps = ceil(_fp * _inputZeros / (1 - _fp));

  std::cout << "input zeros:" << _inputZeros << std::endl;
  std::cout << "input ones:" << _inputOnes << std::endl;
  std::cout << "input missing:" << _inputMissing << std::endl;
  std::cout << "total fp flips allowed:" << total_fps << std::endl;

  //constraint to ensure the matrix _x represent a linear phylogeny

  for (IloInt k = 0; k < rows; k++)
  {

    for (IloInt i = 0; i < cols; i++)
    {

      for (IloInt j = 0; j < cols; j++)
      {

        if (_B.is(k, j) != -1 && _B.is(k, i) != -1)
        {
          _model.add(IloIfThen(_env, _c[i] < _c[j], _x[k][j] <= _x[k][i]));
        }
      }
    }
  }

  //the number of false positives must be less than integer variable y,
  // which has an upper bound of allowed false positives
  IloIntExpr OnetoZeroFlips(_env, 0);

  for (int i = 0; i < rows; i++)
  {
    for (int j = 0; j < cols; j++)
    {
      if (_B.is(i, j) == 1)
      {
        OnetoZeroFlips += _x[i][j] == 0;
      }
    }
  }
  _model.add(OnetoZeroFlips <= _y); /*!< limits the number of flips from 1 to 0 (bounded by domain of _y)  */

  _model.add(IloAllDiff(_env, _c)); /*!< constraint to ensure every column in matrix has a different value in ordering */

  /*! Add objective to the model that minimizes the number of flips from 0 to 1 (false negatives) */
  std::cout << "adding objective" << std::endl;

  IloIntExpr obj(_env, 0);

  for (IloInt i = 0; i < rows; i++)
  {
    for (IloInt j = 0; j < cols; j++)
    {
      if (_B.is(i, j) == 0)
      {
        obj += _x[i][j];
      }
    }
  }

  _model.add(IloMinimize(_env, obj));

  /*! Set the parameters and search phases  for the constraint programming solver*/
  _cp.setParameter(IloCP::TimeLimit, _time);
  _cp.setParameter(IloCP::LogPeriod, 10000);
  _cp.setParameter(IloCP::Workers, 1);
  _cp.setSearchPhases(IloSearchPhase(_env, _c));
  _cp.startNewSearch();

  /*! Call the solver and compute the number of each flip type */
  if (_cp.solve())
  {
    _solve = true;

    _objValue = _cp.getObjValue();
    _fpCounts = _cp.getValue(_y);

    for (IloInt i = 0; i < rows; i++)
    {
      for (IloInt j = 0; j < cols; j++)
      {
        _Bout[i][j] = _cp.getValue(_x[i][j]);
      }
    }

    _outputOnes = 0;
    _outputZeros = 0;
    _outputMissing = 0;

    for (int i = 0; i < _Bout.size(); i++)
    {
      for (int j = 0; j < _Bout[i].size(); j++)
      {

        if (_Bout[i][j] == 1)
        {
          _outputOnes++;
        }
        else if (_Bout[i][j] == -1)
        {
          if (_B.is(i, j) != -1)
          {
            std::cout << "warning: missing data flipped at " << i << "," << j << std::endl;
          }
          _outputMissing++;
        }
        else
        {
          _outputZeros++;
        }
      }
    }

    _estFP = _fpCounts * 1.0 / _outputZeros;
    _estFN = _objValue * 1.0 / _outputOnes;
  }
  _cp.end();

  std::cout << "solve:" << _solve << std::endl;
}

void PhyolinCP::write_counts(std::string filename)
{

  std::ofstream myfile(filename);

  std::string pred = _estFN <= _threshold ? "linear" : "branched";
  myfile << "variable"
         << ","
         << "value" << std::endl;
  myfile << "prediction"
         << "," << pred << std::endl;
  myfile << "cells"
         << "," << _cells << std::endl;
  myfile << "sites"
         << "," << _sites << std::endl;
  myfile << "inputFP"
         << "," << _fp << std::endl;
  myfile << "estFP"
         << "," << _estFP << std::endl;
  myfile << "inputFN"
         << "," << _threshold << std::endl;
  myfile << "estFN"
         << "," << _estFN << std::endl;
  myfile << "input1"
         << "," << _inputOnes << std::endl;
  myfile << "input0"
         << "," << _inputZeros << std::endl;
  // myfile << "loglikelihood" << "," << getLikelihood() << std::endl;
  myfile << "N_00"
         << "," << N_00 << std::endl;
  myfile << "N_10"
         << "," << N_10 << std::endl;
  myfile << "N_11"
         << "," << N_11 << std::endl;
  myfile << "N_01"
         << "," << N_01 << std::endl;
  myfile << "N_missing"
         << "," << _inputMissing << std::endl;

  myfile.close();
}

double PhyolinCP::getLikelihood()
{

  int missing = 0;

  for (int i = 0; i < _cells; i++)
  {

    for (int j = 0; j < _sites; j++)
    {
      if (_B.is(i, j) == 0 && _Bout[i][j] == 0)
      {
        N_00++;
      }
      else if (_B.is(i, j) == 0 && _Bout[i][j] == 1)
      {
        N_01++;
      }
      else if (_B.is(i, j) == 1 && _Bout[i][j] == 1)
      {
        N_11++;
      }
      else if (_B.is(i, j) == 1 && _Bout[i][j] == 0)
      {

        N_10++;
      }
      else
      {
        missing++;
      }
    }
  }

  double log_likelihood = N_00 * log(1 - _fp) + N_11 * log(1 - _threshold) + N_10 * log(_fp) + N_01 * log(_threshold);

  return (log_likelihood);
}

void PhyolinCP::write_csv(std::string filename, std::string delim)
{

  std::vector<std::string> colnames = _B._mutIDs;

  std::ofstream myFile(filename);

  /*! write the column names to the output file 
  */
  for (int j = 0; j < colnames.size(); ++j)
  {
    myFile << colnames[j];
    if (j != colnames.size() - 1)
    {
      myFile << delim;
    }
  }
  if (colnames.size() > 0)
  {
    myFile << "\n";
  }

  /*! write the flipped data to the output file 
     */
  for (int i = 0; i < _Bout.size(); ++i)
  {
    for (int j = 0; j < _Bout[i].size(); ++j)
    {
      myFile << _Bout[i][j];
      if (j != _Bout[i].size() - 1)
      {

        myFile << delim;
      }
    }
    myFile << "\n";
  }

  myFile.close();
}

int main(int argc, char **argv)
{

  /*! Initialize input variables and default parameters 
  */
  std::string inputFile;
  std::string sep(",");
  std::string outFile;
  std::string countFile;
  double fp = 0;
  bool headers = false;
  int time = 300;
  double threshold = 0.25;

  lemon::ArgParser ap(argc, argv);

  ap.refOption("input", "Binary Single Cell Matrix File", inputFile, true);
  ap.refOption("output", "Output filename for output linear phylogeny", outFile);
  ap.refOption("counts", "Filename recording the number of false positive and negative flips", countFile);
  ap.refOption("fp", "Number of False Positives to Allow", fp);
  ap.refOption("fn", "False Negative Threshold", threshold);
  ap.refOption("headers", "Input File has Headers", headers);
  ap.refOption("sep", "Character separator in input and output files, default is ','", sep);
  ap.refOption("time", "Max runtime in seconds of the solver", time);
  ap.run();

  /*! create a single cell matrix object from the input data */
  SingleCellMatrix Bin(inputFile, ",", headers);

  /*! create a Phyolin object and solve the constrain programming model */
  PhyolinCP phy(Bin, fp, time, threshold);

  if (phy.getSolveStatus())
  {
    std::cout << "Flips: " << phy.getObjective() << std::endl;
    double like = phy.getLikelihood();
    std::cout << "Log Likelihood: " << like << std::endl;
    phy.write_csv(outFile, sep);
    phy.write_counts(countFile);
  }
  else
  {
    std::cout << "No Solution" << std::endl;
  }

  return 0;
}
