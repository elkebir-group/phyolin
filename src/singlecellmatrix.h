

#ifndef SINGLECELLMATRIX_H
#define SINGLECELLMATRIX_H

#include <fstream>
#include <vector>
/*!
 * @brief A binary matrix containing single-cell sequencing data (shape: cells x SNVs) which
 * an entry of 1 indicates that an SNV is present in a cell and 0 otherwise.
 */
class SingleCellMatrix
{
public:
    /*!
    * @brief Constructor to intialize an empty single cell matrix of dimension n x m
    * @param n the number of desired rows
    * @param m the number of desired columns
    * */
   SingleCellMatrix(int n, int m);

    /*!
    * @brief Constructor to intialize a single cell matrix of dimension from a specified input file
    * @param filename the string path to the location of the output file
    * @param sep the delimiter of the input file (i.e "," for csv or "\t" for tsv files)
    * @param header a boolean indicating if the input file has column names
    * */
   SingleCellMatrix(std::string filename, std::string sep, bool header);

   /*!
    * @brief Constructor to intialize a single cell matrix of dimension from 2-d vector
    * @param sscm a 2-d vector of integers containing the single-cell matrix input data
    * */
   SingleCellMatrix(std::vector<std::vector<int> > sscm);

   std::vector<int> _sscm;  /*!< a one-dimension vector containing the single-cell data */

   int _cells; /*!< the number of cells (rows) contained in the matrix */

   int _sites; /*!< the number of mutations (columns) contained in the matrix */

   std::vector<std::string> _cellIDs; /*!< the labels of the cells */

   std::vector<std::string> _mutIDs; /*!< the labels of the mutations */
   /*!
    * @brief gets the value of the matrix at row i and column j
    * @param i the row of the desired entry
    * @param j the column of the desired entry 
    * */
   double is(const int i, const int j);

   /*!
    * @brief read an input matrix from a file into a SingleCellMatrix
    * @param filename the string path to file
    * @param delim the delimiter denoting entries in the input file (i.e "," for csv or "\t" for tsv files)
    * @param header a boolean indicating if the input file has column names
    */
   void readfromfile(std::string filename, std::string delim, bool header);

    /*!
    * @brief writes the Single Cell Matrix to a file
    * @param filename the string path to the location of the output file
    * @param delim the delimiter denoting entries in the input file (i.e "," for csv or "\t" for tsv files)
    */
   void writetofile(std::string filename, std::string delim);
};

#endif // SINGLECELLMATRIX_H