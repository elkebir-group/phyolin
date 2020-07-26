/*
 * phyolinCP.cpp
 *
 *  created on 30 Apr 2020
 *  refactored on 26 July 2020
 *   Leah Weber
 */

#include "phyolinCP.h"
#include <vector>
#include <array>
#include <cstdlib>
// #include "utils.h"

PhyolinCP::PhyolinCP( const std::vector<std::vector<int>> B, double fn)
                         
  : _env()
  , _model(_env)
  , _cp(_model)
  ,_fnr(fn)
  ,_solve(false)
  ,_flips()
  ,_time(500)
  ,_B(B)
{
  std::cout << "Phyolin instance created!" << std::endl; 
}

/*
*  intializes as instance of Phyolin with a give input matrix B
*  and create the CP model to be solved, then calls the solver
*/
void PhyolinCP::solve()
{

  //get the size of the input matrix B
  IloInt rows = _B.size();
  IloInt cols = _B[0].size();

  //determine the number of zeros in B
  int zero_count = 0;
  for(int i=0; i < rows; i++){
    for(int j=0; j < cols; j++)
      if(_B[i][j]==0){
        zero_count = zero_count + 1;
      }
  }

  // int total_flips_allowed = zero_count * _fnr;

  //initialize variables 
  typedef IloArray<IloIntVarArray> IloIntVarArray2;
  IloIntVarArray2 _x(_env, rows);

  //add x variables to the model evironment with domain 0-1
  for (IloInt i = 0; i < rows; i++) {
      _x[i]  = IloIntVarArray(_env, cols, 0, 1);
  
  }

  //add permutation column variables to the model
  _c = IloIntVarArray(_env, cols, 0, cols-1);
 

  
 //counts the total number of constraints added to the model
 int total_const = 0;

//adds the constraint that forces the flipped matrix to represent a linear phylogeny
 for(IloInt k=0; k < rows; k++){
  for(IloInt i=0; i < cols; i++){
    for(IloInt j=0; j < cols; j++){
      total_const += 1;
      _model.add(IloIfThen(_env, _c[i] < _c[j], _x[k][j] <= _x[k][i]));
    }
  }
}

  
 
  
  //add the set of constraints to the model that prevent 1's from being flipped to 0
  for(int i=0; i < rows; i++){
    for(int j=0; j < cols; j++){

      if(_B[i][j]==1){
        total_const = total_const + 1;
       
        _model.add(_x[i][j]==1);
      }
    }
  }
 
  /*
  *Generate the objective value (the total number of flips) which consists of two terms
  * The first term is the sum of _x variables or B' after flipping
  * and the second term subtracts the total number of 1's in B
  *  flips = sum(_x) - sum(B)
  */
  IloIntExpr obj(_env,0);
  IloIntExpr Bsum(_env);
  for(IloInt i=0; i < rows; i++){
     for(IloInt j =0; j < cols; j++){
        Bsum += -1*_B[i][j];
    }
  }
 

  for(IloInt i=0; i < rows; i++){
      obj += IloSum(_x[i]); 
  }
  obj += Bsum;

  //we want to minimize the total number of flips 
  _model.add(IloMinimize(_env, obj));

 //for the constraint satisfaction version of the problem
 // _model.add(obj <= total_flips_allowed);  
  //_model.add(obj == total);
     
  //global constraint to ensure each column takes a different value in the permutation
  _model.add(IloAllDiff(_env, _c));
  total_const = total_const + 1;


   
  
  

  //set solver parameters

    _cp.setParameter(IloCP::TimeLimit, _time);
    _cp.setParameter(IloCP::LogVerbosity, IloCP::Quiet);
  
  
  /*
  * Solve the cp model
  */
  if(_cp.solve()){
     _solve = true;

     //set the minimum number of flips found in order for B to represent a linear phylogeny
    _flips = _cp.getObjValue();

  
  //Find the values of B' (the ouput matrix)    
    for(IloInt i = 0; i < rows; i++){
      std::vector<int> v;
      for(IloInt j =0; j < cols; j++){
        v.push_back(_cp.getValue(_x[i][j]));
      }
      _Bout.push_back(v);
  
    }
  }
      
      
  
}  //end solve

//https://www.gormanalysis.com/blog/reading-and-writing-csv-files-with-cpp/
void PhyolinCP::write_csv(std::string filename, std::vector<std::string> colnames, 
                   std::string delim){
 
    
    // Create an output filestream object
    std::ofstream myFile(filename);
    
    // Send column names to the stream
    for(int j = 0; j < colnames.size(); ++j)
    {
        myFile << colnames[j];
        if(j != colnames.size() - 1) myFile << delim; // No comma at end of line
    }
    if(colnames.size() > 0){
        myFile << "\n";
    }
    
    
    // Send data to the stream
    for(int i = 0; i < _Bout.size(); ++i)
    {
        for(int j = 0; j < _Bout[i].size(); ++j)
        {
            myFile << _Bout[i][j];
            if(j != _Bout[i].size() - 1) myFile << delim; // No comma at end of line
        }
        myFile << "\n";
    }
    
    // Close the file
    myFile.close();
}

/*
* Perform the hypothesis test 
*
*/
double PhyolinCP::hypoTest()
{
  int r = _Bout.size();
  int c = _Bout[0].size();
  int all_ones = 0;
  for(int i =0; i < _Bout.size(); i++){
    for(int j=0; j < _Bout[i].size(); j++){
      all_ones += _Bout[i][j];
    }
  }
  double estFn = (double) _flips / all_ones;
  return estFn;
}

int main(int argc, char** argv){
  
  //read in path to csv file
  //convert to B input
  std::string inputFile  = argv[1];
  //std::cout << inputFile << std::endl;
  std::string outputFile = argv[2];
  std::vector<std::vector<int>> Bin;
  std::vector<std::string> colnames;
  std::ifstream myfile(inputFile);
  
  //myfile.ignore(1000,"\n");
  if(myfile.is_open()){
    //bool firstline = true;
    int n =0;
    while(!myfile.eof()){
      
      std::string line;

      std::getline(myfile, line);
      if(line == "") break;
      n = n + 1;
      
    
        std::istringstream input(line);
        
        if(n > 1){
        std::vector<int> bits;
        
          std::string s;
          for (std::array<char, 4> a; input.getline(&a[0], 4, ','); ) {
                  bits.push_back(std::stoi(&a[0]));
          }

          Bin.push_back(bits);
      }else{
        
          std::string s;
          for (std::array<char, 1000> a; input.getline(&a[0], 1000, ','); ) {
                  colnames.push_back(&a[0]);
          }
      }   
      
    }
      
  }



  double fn = strtod(argv[3], NULL);
  std::cout << std::endl;
  std::cout << "Starting Phyolin...." <<std::endl;
  PhyolinCP phy(Bin, fn);
  float estFN = 0;
 
   phy.solve();
  
  if(phy._solve){
    std::cout << "Solve complete, writing output..." << std::endl;
    std::cout << "Total Flips: " << phy._flips << std::endl;
    phy.write_csv(outputFile, colnames, ",");
    estFN = phy.hypoTest();
    std::cout << "Phyolin estimated false negative rate: " << estFN << std::endl;
    if(estFN <= fn){
      std::cout << "Topology is Linear" << std::endl;
    }else{
        std::cout << "Topology is Branched" << std::endl;
    }
  }else{
    std::cout << "Error! Model could not be solved" << std::endl;
  }
  
    std::cout << std::endl;
  
  return 0;
}

