

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
  PhyolinCP(const SingleCellMatrix B, double fp, int time, double threshold);
  
   int getObjective() {
      return _objValue;
    }

    int getSolveStatus() {
      return _solve;
    }
  

   void write_csv(std::string filename, 
                   std::string delim);
  
   void write_counts(std::string filename);


  double getLikelihood();
   
  
private:
  void init();
  
private:
  SingleCellMatrix _B;
  /// Number of features in distinguishing feature set
  //const int _k;
  /// Environment
 

  IloEnv _env;
  /// CPlex model
  IloModel _model;
  /// Solver
  IloCP _cp;
  /// Cover variables
  IloArray <IloIntVarArray> _x;

  IloIntVarArray _c;

  IloIntVar _y;

  bool _solve;

  // IloIntVarArray _d;
  /// Minimum weight variable
  //IloNumVar _z;
  /// Objective value
  double _objValue;

  double _fp;

  int _time;

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
