/*
 * phyolinCP.h
 *
 *  created on 30 Apr 2020
 *      Author: Leah Weber
 */

#ifndef PHYOLINCP_H
#define PHYOLINCP_H


#include <ilcp/cp.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

class PhyolinCP
{
public:
  PhyolinCP(const std::vector<std::vector<int>> B, double fn);
  
  double hypoTest();

  void solve();
  

   void write_csv(std::string filename, std::vector<std::string> colnames, 
                   std::string delim);

   bool _solve;

   double _flips;



   
  

  
private:
  //input matrix B
 std::vector<std::vector<int>> _B;

  /// Environment
  IloEnv _env;

  /// CP model
  IloModel _model;

  /// Solver
  IloCP _cp;
  
  /// decision variablse
  IloArray <IloIntVarArray> _x;

  IloIntVarArray _c;
  
  //false negative rate
  double _fnr;

  //B' output matrix after flipping
  std::vector<std::vector<int>> _Bout;

  //time limit for the solver
  int _time;

};

#endif 
