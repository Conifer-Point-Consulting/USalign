#ifndef _USalignLib_h_
#define _USalignLib_h_ 1

#include "USalign.h"
#include "basic_fun.h"
#include <string>
#include <vector>

int doAlignment(std::string pdb1, std::string pdb2, TMatrix &tmatrix);
TMatrix getMatrix(std::string pdb_str1, std::string pdb_str2);

#endif
