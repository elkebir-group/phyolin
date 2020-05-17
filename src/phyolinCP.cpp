/*
 * phyolinCP.cpp
 *
 *  created on 30 Apr 2020
 *   Leah Weber
 */

#include "phyolinCP.h"
#include <vector>
#include <array>
// #include "utils.h"

PhyolinCP::PhyolinCP( const std::vector<std::vector<int>> B, double fn)
                         
  : _env()
  , _model(_env)
  , _cp(_model)
  ,_fnr(fn)
  ,_solve(false)
{
  init(B);
}

void PhyolinCP::init(const std::vector<std::vector<int>> B)
{

  
  IloInt rows = B.size();
  IloInt cols = B[0].size();
  //std::cout << "rows: " << rows << std::endl;
  //std::cout << "cols: " << cols << std::endl;
  
  int zero_count = 0;
  for(int i=0; i < rows; i++){
    for(int j=0; j < cols; j++)
      if(B[i][j]==0){
        zero_count = zero_count + 1;
      }
  }

  int total_flips_allowed = zero_count * _fnr;

  //initialize variables 
  typedef IloArray<IloIntVarArray> IloIntVarArray2;
  IloIntVarArray2 _x(_env, rows);

  //add x variables to the model with domain 0-1
  for (IloInt i = 0; i < rows; i++) {
      _x[i]  = IloIntVarArray(_env, cols, 0, 1);

      
  }



  //add permutation column variables to the model
  _c = IloIntVarArray(_env, cols, 0, cols-1);
 

  
  
 int total_const = 0;
 for(IloInt k=0; k < rows; k++){
  for(IloInt i=0; i < cols; i++){
    for(IloInt j=0; j < cols; j++){
      total_const += 1;
      _model.add(IloIfThen(_env, _c[i] < _c[j], _x[k][j] <= _x[k][i]));
    }
  }
}
 std::cout << "Constraints added: " << total_const << std::endl;
  
 

  int tc = 0;
  for(int i=0; i < rows; i++){
    for(int j=0; j < cols; j++){

      if(B[i][j]==1){
        tc += 1;
       
        _model.add(_x[i][j]==1);
      }
    }
  }
   std::cout << "adding domain constraint: " << tc << std::endl;
  // IloIntExpr InitTotal(_env, Bsum);
  // IloIntExpr total(_env);
  std::cout << "adding objective" <<std::endl;  
  IloIntExpr obj(_env,0);
  IloIntExpr Bsum(_env);
  for(IloInt i=0; i < rows; i++){
     for(IloInt j =0; j < cols; j++){
        Bsum += -1*B[i][j];
    }
  }
 
   
  std::cout << "init obj" <<std::endl; 
  for(IloInt i=0; i < rows; i++){
      //std::cout << i << std::endl;
      obj += IloSum(_x[i]);
    
  }
  obj += Bsum;
  
 // _model.add(obj <= total_flips_allowed);

  std::cout << "adding all different" <<std::endl;
     
  //_model.add(obj == total);
     
  
  _model.add(IloAllDiff(_env, _c));

   std::cout << "obj" <<std::endl; 
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
  
  std::cout << "model initialized" << std::endl;

    _cp.setParameter(IloCP::TimeLimit, 500);
    _cp.setParameter(IloCP::LogPeriod, 10000);
  if(_cp.solve()){
     _solve = true;
  // if (!_cp.solve())
  // {
  //   new_sol.push_back(-1);
  //   return new_sol;
  // }



  
      _objValue = _cp.getObjValue();
      std::cout << "Obj value: " << _objValue  << std::endl;
      
      for(IloInt i = 0; i < rows; i++){
        std::vector<int> v;
        for(IloInt j =0; j < cols; j++){
          v.push_back(_cp.getValue(_x[i][j]));
        }
        _Bout.push_back(v);
      
      }
  }
      
    
       

       //std::cout << _cp.getValue(_x[0][0]) << std::endl;
      // for(IloInt i=0; i < 4; i++){
      //   for(IloInt j=0; j < 5; j++){
      //       std::cout << "x_" << i << "," << j << "=" << _cp.getValue(_x[i][j]) << std::endl;
      //   }
      // }
 
  
}

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

double PhyolinCP::solve()
{
  int r = _Bout.size();
  int c = _Bout[0].size();
  double estFn = (double) _objValue / (r*c);
  return estFn;
}

int main(int argc, char** argv){
  
  //read in path to csv file
  //convert to B input

  std::string inputFile  = argv[1];
  std::cout << inputFile << std::endl;
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
      //std::cout << n << std::endl;
      std::getline(myfile, line);
      if(line == "") break;
      n = n + 1;
      
        //std::cout << line << std::endl;
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


  // std::vector<int> v(4,0);
  // std::vector<std::vector<int>> Bin(1000,v);
  // Bin[0][0] =1;
  // Bin[3][1] = 1;
  // Bin[4][1] = 1;
  // std::cout << "launched" << std::endl;
  //std::cout << Bin.size() << std::endl;
  double fn = std::stoi(argv[3])/ (double) 100;
  PhyolinCP phy(Bin, fn);
  std::cout << "Flips: " << phy.solve() << std::endl;
  //std::cout << "Phylogeny is Linear:" << phy._solve << std::endl;
  if(phy._solve){
    phy.write_csv(outputFile, colnames, ",");
  }
  
  //std::cout << "init" << std::endl;
  //phy.solve();

  
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


