

#include "singlecellmatrix.h"
#include <sstream>
#include<iostream>


SingleCellMatrix::SingleCellMatrix(int n, int m)
:_cellIDs()
,_mutIDs()
,_sscm(n*m)
{
    _cells = n;
    _sites = m;
    for(int i=0; i < _cells; ++i){
      _cellIDs.push_back("c" + std::to_string(i));
    }

    for(int i=0; i < _sites; ++i){
      _mutIDs.push_back("m" + std::to_string(i));
    }
}


SingleCellMatrix::SingleCellMatrix(std::vector<std::vector<int>> sscm)
:_cells(sscm.size())
,_sites(sscm[0].size())
,_sscm()
, _cellIDs()
, _mutIDs()
{
    int matsize = _cells* _sites;
    _sscm.resize(matsize);
    // std::cout << _cells << ":" << _sites << ":" <<  _sscm.size()<< std::endl;


      
  for(int i=0; i < sscm.size(); i++){
    for(int j=0; j < sscm[i].size();j++){
      _sscm[i*_sites + j] = sscm[i][j];
   }
  }
 
    //   std::vector<std::vector<double>> sscm;
    //   sscm.reserve(ilil.size());
    //     for (auto&& il : ilil)
    //     {
    //         sscm.emplace_back(il);
    //     }
    //  _sites = sscm[0].size();

  

     for(int i=0; i < _cells; ++i){
      _cellIDs.push_back("c" + std::to_string(i));
    }

    for(int i=0; i < _sites; ++i){
      _mutIDs.push_back("m" + std::to_string(i));
    }
}


SingleCellMatrix::SingleCellMatrix(std::string filename, std::string sep, bool header):
_sscm()
{
     std::cout << "reading  from file " << filename <<std::endl;
     readfromfile(filename, sep, header);
   
   if(!header){

  
    for(int i=0; i < _cells; ++i){
      _cellIDs.push_back("c" + std::to_string(i));
    }

    for(int i=0; i < _sites; ++i){
      _mutIDs.push_back("m" + std::to_string(i));
    }
   }

}
void SingleCellMatrix::readfromfile(std::string inputFile, 
                                                         std::string sep, bool header)
{

//   std::vector<std::string> colnames;

  std::vector<std::vector<double>> temp;
  std::ifstream myfile(inputFile);
  

  if(myfile.is_open()){
          
      std::string line;
      if(header){
        std::getline(myfile, line);
        std::stringstream input(line);
        int mutname =0;
        while(input.good()){
            std::string call;
            if(std::getline(input, call, ',')){
                  if(mutname > 0){
                    _mutIDs.push_back(call);
                  }
                  
            }
            mutname++;

        }
      }
    while(!myfile.eof()){
      std::getline(myfile, line);
      
      

      std::vector<double> cellmuts;
      if(line == "") break;
              

        std::stringstream input(line);
        int cellidcol =1;
        while(input.good()){
            std::string call;
            if(std::getline(input, call, ',') && (!header || cellidcol > 0)){
                
                   cellmuts.push_back(std::stoi(call));
            }else{
               
                  _cellIDs.push_back(call);
            }
            cellidcol++;
       
        } 
         
          temp.push_back(cellmuts);
       
    }
  }
  _sites = temp[0].size();
  _cells = temp.size();
  _sscm.resize(_sites*_cells);
  for(int i=0; i < temp.size(); i++){
    for(int j=0; j < temp[i].size();j++){
      _sscm[i*_sites + j] = temp[i][j];
   }
  }
    myfile.close();


}

double SingleCellMatrix::is(const int i, const int j){
  return _sscm[i*_sites + j];
}
void SingleCellMatrix::writetofile(std::string filename, 
                                                        std::string delim)
{
  
 
    
    // Create an output filestream object
    std::ofstream myFile(filename);
    
    // Send column names to the stream
    for(int j = 0; j < _mutIDs.size(); ++j)
    {
        myFile << _mutIDs[j];
        if(j != _mutIDs.size() - 1) myFile << delim; // No comma at end of line
    }
    if(_mutIDs.size() > 0){
        myFile << "\n";
    }
    
    
    // Send data to the stream
    for(int i = 0; i < _cells; ++i)
    {
        for(int j = 0; j < _sites; ++j)
        {
            myFile << is(i,j);
            if(j != _sites - 1) myFile << delim; // No comma at end of line
        }
        myFile << "\n";
    }
    
    // Close the file
    myFile.close();
}

