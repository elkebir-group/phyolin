/*
 * setcoverilp.h
 *
 *  Created on: 5-dec-2019
 *      Author: M. El-Kebir
 */

#ifndef PHYOLINCP_H
#define PHYOLINCP_H


#include <ilcp/cp.h>
#include <vector>
//#include <fstream>

class PhyolinCP
{
public:
  PhyolinCP(const std::vector<std::vector<int>> B);
  
  double solve();
  //void printSolutions(std::string, std::string, int ,const std::vector<std::vector<int>> );


  
private:
  void init(const std::vector<std::vector<int>> B);
  
private:
  /// Number of features in distinguishing feature set
  //const int _k;
  /// Environment
 
  IloEnv _env;
  /// CPlex model
  IloModel _model;
  /// Solver
  IloCP _cp;
  /// Cover variables
  IloArray<IloIntVarArray> _x;

  IloIntVarArray _c;
  /// Minimum weight variable
  //IloNumVar _z;
  /// Objective value
  double _objValue;
};

#endif // SETCOVERILP_H
