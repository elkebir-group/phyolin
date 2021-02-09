/*
 * phyolinCP.cpp
 *
 *  created on 30 Apr 2020
 *   Leah Weber
 */

#include "phyolinCP.h"
#include <vector>
#include <array>
#include <cstdlib>
#include <lemon/arg_parser.h>
#include <math.h>



PhyolinCP::PhyolinCP( const SingleCellMatrix B, double fp, int time, double threshold)
                         
  : _env()
  , _model(_env)
  , _cp(_model)
  ,_fp(fp)
  , _B(B)
  , _cells(B._cells)
  , _sites(B._sites)
  , N_11(0)
  , N_00(0)
  , N_01(0)
  , N_10(0)
  , _Bout(_cells, std::vector<int>(_sites))
  , _inputOnes(0)
  , _inputZeros(0)
  , _outputZeros(0)
  , _outputOnes(0)
  , _inputMissing(0)
  , _outputMissing(0)
  , _estFN(0)
  , _estFP(0)
  , _fpCounts(0)
  , _time(time)
  , _threshold(threshold)

  ,_solve(false)
{
  init();
}

void PhyolinCP::init()
{

  
  IloInt rows = _B._cells;
  IloInt cols = _B._sites;
  std::cout << "rows: " << rows << std::endl;
  std::cout << "cols: " << cols << std::endl;

  //initialize variables 
  typedef IloArray<IloIntVarArray> IloIntVarArray2;
  IloIntVarArray2 _x(_env, rows);

  //add x variables to the model with domain -1, 0,1
  for (IloInt i = 0; i < rows; i++) {
      _x[i]  = IloIntVarArray(_env, cols, -1, 1);

      
  }
 
  
  

  for(int i=0; i < rows; i++){
    for(int j=0; j < cols; j++)
      if(_B.is(i,j)==1){
        _inputOnes++;
        _model.add(_x[i][j] >= 0);
        //_model.add(_x[i][j]==1
      }else if(_B.is(i,j)==-1){
        _inputMissing++;
             //add a constraint so missing data doesn't change
        _model.add(_x[i][j]== -1);
      }else{
         _model.add(_x[i][j] >= 0);
         _inputZeros++;
      }
  }



   int total_fps = ceil(_fp * _inputZeros/(1- _fp));

  // int total_fps = one_count * _fp;
  std::cout << "input zeros:" << _inputZeros << std::endl;
  std::cout << "input ones:" << _inputOnes << std::endl;
  std::cout << "input missing:" << _inputMissing << std::endl;
  std::cout << "total fp flips allowed:" << total_fps << std::endl;

  
  //_y = IloIntVarArray(_env, one_count, 0, 1);

  _y = IloIntVar(_env, 0, total_fps);

  //add permutation column variables to the model
  _c = IloIntVarArray(_env, cols, 0, cols-1);
 
 
 
 
//adding constraints

//constraint to ensure the matrix _x represent a linear phylogeny
 int total_const = 0;
 for(IloInt k=0; k < rows; k++){
  for(IloInt i=0; i < cols; i++){
    for(IloInt j=0; j < cols; j++){
      total_const += 1;
      if(_B.is(k,j) != -1 && _B.is(k,i) != -1){
        _model.add(IloIfThen(_env, _c[i] < _c[j], _x[k][j] <= _x[k][i]));
      }
      
    }
  }
}


 
 //the number of false positives must be less than integer variable y,
 // which has an upper bound of allowed false positives
  IloIntExpr OnetoZeroFlips(_env, 0);

 for(int i=0; i< rows; i++){
   for(int j = 0; j < cols; j++){
     if(_B.is(i,j)==1){
       OnetoZeroFlips += _x[i][j] ==0;
       
     }

   }
 }
 _model.add(OnetoZeroFlips <= _y);
 
 //every column has a different value in the ordering
 _model.add(IloAllDiff(_env, _c));






  // std::cout << "Total FPS " << total_fps << std::endl;

  std::cout << "adding objective" <<std::endl;  
  IloIntExpr obj(_env,0);
  IloIntExpr denom(_env,0);
  //IloIntExpr Bsum(_env);
  for(IloInt i=0; i < rows; i++){
     for(IloInt j =0; j < cols; j++){
        if(_B.is(i,j)==0){
          obj +=  _x[i][j];
        }
        
    }
  }
 

     
  
 

  
  _model.add(IloMinimize(_env, obj));
  //_model.add(obj/denom <= _threshold);
   
 

  
 


  _cp.setParameter(IloCP::TimeLimit,_time);
  _cp.setParameter(IloCP::LogPeriod, 10000);
  _cp.setParameter(IloCP::Workers, 1);
  _cp.setSearchPhases(IloSearchPhase(_env, _c));
    //_cp.setParameter(IloCP::SearchType, "Multipoint")
  _cp.startNewSearch();

  if(_cp.solve()){
     _solve = true;
    
  


  
      _objValue = _cp.getObjValue();
      _fpCounts = _cp.getValue(_y);

      
      for(IloInt i = 0; i < rows; i++){
        for(IloInt j =0; j < cols; j++){
         _Bout[i][j]= _cp.getValue(_x[i][j]);
        }
      
      }
      
      _outputOnes = 0;
      _outputZeros = 0;
      _outputMissing=0;
      
      for(int i=0; i < _Bout.size(); i++){
        for(int j = 0; j < _Bout[i].size(); j++){
    
          if(_Bout[i][j] == 1){
            _outputOnes++;
          }else if(_Bout[i][j]==-1){
            if(_B.is(i,j) != -1){
              std::cout << "warning: missing data flipped at " << i << "," << j << std::endl;
            }
            _outputMissing++;
          }else{
            _outputZeros++;
          }
        }
      }

      _estFP = _fpCounts*1.0/ _outputZeros;
      _estFN = _objValue* 1.0/ _outputOnes;
 
      // if(_estFN <= _threshold){
           
      //   break;
      // }
  }
  _cp.end();



  std::cout << "solve:" << _solve << std::endl;  

       //std::cout << _cp.getValue(_x[0][0]) << std::endl;
      // for(IloInt i=0; i < 4; i++){
      //   for(IloInt j=0; j < 5; j++){
      //       std::cout << "x_" << i << "," << j << "=" << _cp.getValue(_x[i][j]) << std::endl;
      //   }
      // }
 
  
}

void PhyolinCP::write_counts(std::string filename){
  
  std::ofstream myfile(filename);

  // myfile << "prediction\tcells\tloci\tinputFP\testFP\tflipsFP\tinptFN\testFN\tflipsFN\tinput0\tinput1\tinputMissing\toutput0\toutput1\toutputMissing" << std::endl;
  
  std::string pred = _estFN <= _threshold ? "linear": "branched";
  myfile << "variable" << "," << "value" << std::endl;
  myfile << "prediction" << "," << pred << std::endl;
  myfile << "cells" << "," << _cells << std::endl;
   myfile << "sites" << "," << _sites << std::endl;
  myfile << "inputFP" << "," << _fp << std::endl;
  myfile << "estFP" << "," << _estFP << std::endl;
  myfile << "inputFN" << "," << _threshold << std::endl;
  myfile << "estFN" << "," << _estFN << std::endl;
  myfile << "input1" << "," << _inputOnes << std::endl;
  myfile << "input0" << "," << _inputZeros << std::endl;
  // myfile << "loglikelihood" << "," << getLikelihood() << std::endl;
  myfile << "N_00" << "," << N_00 << std::endl;
  myfile << "N_10" << "," << N_10 << std::endl;
  myfile << "N_11" << "," << N_11 << std::endl;
  myfile << "N_01" << "," << N_01 << std::endl;
  myfile << "N_missing" << "," << _inputMissing << std::endl;

  // myfile  << _fp << "\t" ;
  // myfile << _estFP <<  "\t";
  // myfile   <<_fpCounts << "\t" ;
  // myfile  << _threshold << "\t" ;
  // myfile  << _estFN << "\t" ;  
  // myfile  << _objValue << "\t";
  // myfile << _inputZeros << "\t";
  // myfile << _inputOnes << "\t";
  // myfile << _inputMissing << "\t";
  // myfile << _outputZeros << "\t";
  // myfile << _outputOnes << "\t";
  // myfile << _outputMissing << "\t";
  

  myfile.close();
}

double PhyolinCP::getLikelihood(){
 int missing = 0;
  for(int i=0; i < _cells; i++){
    for(int j=0; j < _sites; j++){
      if(_B.is(i,j)==0 && _Bout[i][j]==0){
        N_00++;
      }else if(_B.is(i,j)==0 && _Bout[i][j]==1){
        N_01++;
      }else if(_B.is(i,j)==1 && _Bout[i][j]==1){
        N_11++;
      }else if(_B.is(i,j)==1 && _Bout[i][j]==0){
        N_10++;
      }else{
        missing++;
      }
    }
  }
  double log_likelihood = N_00*log(1-_fp) +  N_11*log(1-_threshold) + N_10*log(_fp) + N_01*log(_threshold);

  return(log_likelihood);
  
}
//https://www.gormanalysis.com/blog/reading-and-writing-csv-files-with-cpp/
void PhyolinCP::write_csv(std::string filename, 
                   std::string delim){
    std::vector<std::string> colnames = _B._mutIDs;
    
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

// double PhyolinCP::solve()
// {
//   int r = _Bout.size();
//   int c = _Bout[0].size();
//   std::cout << "objective" << _objValue << std::endl;
//   double estFn = (double) _objValue / (r*c);
//   return estFn;
// }



int main(int argc, char** argv){
  
  

 
  std::string inputFile;
  std::string sep (",");
  std::string outFile;
  std::string countFile;
  double fp =0;
  bool headers = false;
  int time = 300;
  double threshold = 0.25;
  
  
  lemon::ArgParser ap(argc, argv);

  ap.refOption("input", "Binary Single Cell Matrix File", inputFile, true);
  ap.refOption("output", "Output filename for output linear phylogeny", outFile);
  ap.refOption("counts", "Filename recording the number of false positive and negative flips", countFile);
  ap.refOption("fp" ,"Number of False Positives to Allow", fp );
  ap.refOption("fn" ,"False Negative Threshold", threshold );
  ap.refOption("headers", "Input File has Headers", headers);
  ap.refOption("sep", "Character separator in input and output files, default is ','", sep);
  ap.refOption("time" ,"Max runtime in seconds of the solver", time );
  ap.run();
  

  SingleCellMatrix Bin(inputFile, ",", headers);

  PhyolinCP phy(Bin, fp, time, threshold);
  
 
  if(phy.getSolveStatus()){
    std::cout << "Flips: " << phy.getObjective() << std::endl;
    double like = phy.getLikelihood();
    std::cout << "Log Likelihood: " << like << std::endl;
    phy.write_csv(outFile,  sep);
    phy.write_counts(countFile);
  }else{
    std::cout << "No Solution" << std::endl;
  }


  
  

  
  return 0;
}

