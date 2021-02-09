

#ifndef SINGLECELLMATRIX_H
#define SINGLECELLMATRIX_H

#include <fstream>
#include <vector>


class SingleCellMatrix
{
   public:

   SingleCellMatrix(int n, int m);

   SingleCellMatrix(std::string filename, std::string sep, bool header);
   
   SingleCellMatrix(std::vector<std::vector<int>>);

   std::vector<int> _sscm;

   int _cells;

   int _sites;

   std::vector<std::string> _cellIDs;

   std::vector<std::string> _mutIDs;

   double is(const int i, const int j);

   void readfromfile(std::string filename, std::string delim, bool header);

   void writetofile(std::string filename, std::string delim);

};

#endif // SINGLECELLMATRIX_H