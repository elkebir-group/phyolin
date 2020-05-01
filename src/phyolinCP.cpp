/*
 * phyolinCP.cpp
 *
 *  created on 30 Apr 2020
 *   Leah Weber
 */

#include "phyolinCP.h"
#include <vector>
// #include "utils.h"

PhyolinCP::PhyolinCP( const std::vector<std::vector<int>> B)
                         
  : _env()
  , _model(_env)
  , _cp(_model)
{
  init(B);
}

void PhyolinCP::init(const std::vector<std::vector<int>> B)
{
  char buf[1024];
  
  IloInt rows = B.size();
  IloInt cols = B[0].size();

  std::cout << cols << std::endl;
  

typedef IloArray<IloIntVarArray> IloIntVarArray2;

  //  for k in range(ROWS):
  //       for i in range(COLS):
  //           for j in range(COLS):
  //               mdl.add(mdl.if_then(c[i] < c[j], x[k][j] <= x[k][i]))
  
  // Initialize variables
  
 //public IloIfThen(const IloEnv env, const IloConstraint left, const IloConstraint right, const char * name=0)

    IloIntVarArray2 _x(_env, rows);

    for (IloInt i = 0; i < rows; i++) {
      _x[i]  = IloIntVarArray(_env, cols, 0, 1);
      
    }
  _c = IloIntVarArray(_env, cols, 0, cols-1);
  
  //_z = IloNumVar(_env, -IloInfinity, 0, "z");
  
  // // Initialize constraints
  // IloExpr sum(_env);
  // sum = 0;

  
  
  IloIntVar obj(_env, 0, rows*cols);

    
  IloIntExpr total(_env); 
  for(IloInt i=0; i < rows; i++){
      total += IloSum(_x[i]);
    
  }

 for(IloInt k=0; k < rows; k++){
  for(IloInt i=0; i < cols; i++){
    for(IloInt j=0; j < cols; j++){
      _model.add(IloIfThen(_env, _c[i] < _c[j], _x[k][j] <= _x[k][i]));
    }
  }
}
  _model.add(obj == total);
  std::cout << "here" << std::endl;
  
  for(int i=0; i < rows; i++){
    for(int j=0; j < cols; j++){
    
      if(B[i][j]==1){
        _model.add(_x[i][j]==1);
      }
    }
  }
  
     
     
     
  
  _model.add(IloAllDiff(_env, _c));


  _model.add(IloMinimize(_env, obj));

   
  // // there are _k features in the cover
  // for (int i = 0; i < nrFeatures; ++i)
  // {
  //   sum += _x[i];
  // }
  // _model.add(sum == _k);
  // sum.clear();
  
  

  // for( std::vector<int> sol: solution_set){
  //   int card = sol.size();
  //   for( int s: sol){
  //     sum += _x[s];
  //   }
  //   _model.add(sum < card);
  //   sum.clear();
  // }
 
  // //create the objective function
  // for (int i = 0; i < nrFeatures; ++i)
  // {
  //   //assert(_freqMap.count(_mutationVector[i]) == 1);
  //   //double freq = log(_freqMap.find(_mutationVector[i])->second);
  //   obj += _x[i];
  // }
  
  
}


double PhyolinCP::solve()
{
  //std::vector<int> new_sol;
  //_cp.setParam(IloCP::Threads, 1);
    _cp.setParameter(IloCP::TimeLimit, 20);
    _cp.setParameter(IloCP::LogPeriod, 10000);
  _cp.solve();
  // if (!_cp.solve())
  // {
  //   new_sol.push_back(-1);
  //   return new_sol;
  // }



  
  _objValue = _cp.getObjValue();
  std::cout << "Obj value: " << _objValue << std::endl;
 
  // const int nrFeatures = _featureVector.size();
  // for (int i = 0; i < nrFeatures; ++i)
  // {
  //   if (_cplex.getValue(_x[i]) > 0.4)
  //   {
      
  //     new_sol.push_back(i);
  //     std::cout << "Featurette:";
  //     for(std::string s: _featureVector[i]){
  //        std::cout  << s << "->";
  //     }

        

        
        
        
  //   }

  //    std::cout << "" << std::endl;
      
      
      //std::cout << _mutationVector[i] << " " << _freqMap.find(_mutationVector[i])->second << std::endl;
    
  
  
  return _objValue;
}

int main(int argc, char** argv){
  
  //read in path to csv file
  //convert to B input
  
  std::vector<int> v(4,0);
  std::vector<std::vector<int>> Bin(1000,v);
  Bin[0][0] =1;
  Bin[3][1] = 1;
  Bin[4][1] = 1;
  std::cout << "launched" << std::endl;
  PhyolinCP phy(Bin);
  std::cout << "init" << std::endl;
  phy.solve();

  
  return 0;
}

// void SetCoverIlp::printSolutions(std::string origFile, std::string subdir, 
//                                   int treeNum,
//                                   const std::vector<std::vector<int>> solutions)
// {
//   std::string fname =  subdir + "/" + origFile + "_T" + 
//                       std::to_string(treeNum) + "_DF.csv";
//   std::ofstream myfile;
//   myfile.open (fname);
//    if(!myfile) 
//    { 
//        std::cout<<"Error in creating file!!!"; 
       
//    } 

//   std::string dfFam = "";
//   //int nr = _featureVector.size();
//   std::vector<StringSet> features = _featureVector;
//   for(const std::vector<int> df: solutions){
//     for(auto it = df.begin(); it != df.end(); ++it){
//       dfFam += "'";
//       StringSet featurettes = features[(*it)];
//       for( auto feat = featurettes.begin(); feat != featurettes.end(); ++feat ){
//         dfFam += (*feat);
//         if(std::next(feat) == featurettes.end()){
//           dfFam += "'";
//         }else{
//           dfFam += " ";
//         }
//       }
//       if(std::next(it) == df.end()){
//         dfFam += "\n";

//       }else{
//         dfFam += ",";
//       }
//     }
//   }
//     myfile << dfFam;
//     myfile.close();
// }


