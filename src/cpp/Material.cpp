/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Material.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;



//	Write material data to Stream
bool CQ4Material::Read(ifstream& Input)
{
	Input >> nset;
	Input >> E >> Nu >> thickness;
	return true;
}

void CQ4Material::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << Nu << setw(16) << thickness << endl;
}
