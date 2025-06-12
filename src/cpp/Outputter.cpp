/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include <ctime>

#include "Domain.h"
#include "Outputter.h"
#include "SkylineMatrix.h"

using namespace std;

//	Output current time and date
void COutputter::PrintTime(const struct tm* ptm, COutputter &output)
{
	const char* weekday[] = {"Sunday", "Monday", "Tuesday", "Wednesday",
							 "Thursday", "Friday", "Saturday"};
	const char* month[] = {"January", "February", "March", "April", "May", "June",
						   "July", "August", "September", "October", "November", "December"};

	output << "        (";
	output << ptm->tm_hour << ":" << ptm->tm_min << ":" << ptm->tm_sec << " on ";
	output << month[ptm->tm_mon] << " " << ptm->tm_mday << ", " << ptm->tm_year + 1900 << ", "
		   << weekday[ptm->tm_wday] << ")" << endl
		   << endl;
}

COutputter* COutputter::_instance = nullptr;

//	Constructor
COutputter::COutputter(string FileName)
{
	OutputFile.open(FileName);

	if (!OutputFile)
	{
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}
}

//	Return the single instance of the class
COutputter* COutputter::GetInstance(string FileName)
{
	if (!_instance)
		_instance = new COutputter(FileName);
    
	return _instance;
}

//	Print program logo
void COutputter::OutputHeading()
{
	CDomain* FEMData = CDomain::GetInstance();

	*this << "TITLE : " << FEMData->GetTitle() << endl;

	time_t rawtime;
	struct tm* timeinfo;

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	PrintTime(timeinfo, *this);
}

//	Print nodal data
void COutputter::OutputNodeInfo()
{
	CDomain* FEMData = CDomain::GetInstance();

	CNode* NodeList = FEMData->GetNodeList();

	*this << "C O N T R O L   I N F O R M A T I O N" << endl
		  << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	unsigned int NUMNP = FEMData->GetNUMNP();
	unsigned int NUMEG = FEMData->GetNUMEG();
	unsigned int NLCASE = FEMData->GetNLCASE();
	unsigned int MODEX = FEMData->GetMODEX();

	*this << "      NUMBER OF NODAL POINTS . . . . . . . . . . (NUMNP)  =" << setw(6) << NUMNP << endl;
	*this << "      NUMBER OF ELEMENT GROUPS . . . . . . . . . (NUMEG)  =" << setw(6) << NUMEG << endl;
	*this << "      NUMBER OF LOAD CASES . . . . . . . . . . . (NLCASE) =" << setw(6) << NLCASE << endl;
	*this << "      SOLUTION MODE  . . . . . . . . . . . . . . (MODEX)  =" << setw(6) << MODEX << endl;
	*this << "         EQ.0, DATA CHECK" << endl
		  << "         EQ.1, EXECUTION" << endl
		  << endl;

	*this << " N O D A L   P O I N T   D A T A" << endl << endl;
	*this << "    NODE       BOUNDARY                         NODAL POINT" << endl
		  << "   NUMBER  CONDITION  CODES                     COORDINATES" << endl;

	for (unsigned int np = 0; np < NUMNP; np++)
		NodeList[np].Write(*this);

	*this << endl;
}

//	Output equation numbers
void COutputter::OutputEquationNumber()
{
	CDomain* FEMData = CDomain::GetInstance();
	unsigned int NUMNP = FEMData->GetNUMNP();

	CNode* NodeList = FEMData->GetNodeList();

	*this << " EQUATION NUMBERS" << endl
		  << endl;
	*this << "   NODE NUMBER   DEGREES OF FREEDOM" << endl;
	*this << "        N           X    Y    Z" << endl;

	for (unsigned int np = 0; np < NUMNP; np++) // Loop over for all node
		NodeList[np].WriteEquationNo(*this);

	*this << endl;
}

//	Output element data
void COutputter::OutputElementInfo()
{
	//	Print element group control line

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NUMEG = FEMData->GetNUMEG();

	*this << " E L E M E N T   G R O U P   D A T A" << endl
		  << endl
		  << endl;

	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		*this << " E L E M E N T   D E F I N I T I O N" << endl
			  << endl;

		ElementTypes ElementType = FEMData->GetEleGrpList()[EleGrp].GetElementType();
		unsigned int NUME = FEMData->GetEleGrpList()[EleGrp].GetNUME();

		*this << " ELEMENT TYPE  . . . . . . . . . . . . .( NPAR(1) ) . . =" << setw(5)
			  << ElementType << endl;
		*this << "     EQ.1, TRUSS ELEMENTS" << endl
			  << "     EQ.2, ELEMENTS CURRENTLY" << endl
			  << "     EQ.3, NOT AVAILABLE" << endl
			  << endl;

		*this << " NUMBER OF ELEMENTS. . . . . . . . . . .( NPAR(2) ) . . =" << setw(5) << NUME
			  << endl
			  << endl;

		switch (ElementType)
		{
		case ElementTypes::Q4:
			OutputQ4Elements(EleGrp);
			break;
		default:
			*this << ElementType << " has not been implemented yet." << endl;
			break;
		}
	}
}

//      Output Q4 element data
void COutputter::OutputQ4Elements(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::GetInstance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		<< endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		<< endl
		<< endl;

	*this << "  SET       YOUNG'S        POISSON" << endl
		<< " NUMBER     MODULUS          RATIO      THICK" << endl
		<< "               E              NU            T" << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	for (unsigned int mset = 0; mset < NUMMAT; mset++)
	{
		*this << setw(5) << mset + 1;
		ElementGroup.GetMaterial(mset).Write(*this);
	}

	*this << endl << endl
		<< " E L E M E N T   I N F O R M A T I O N" << endl;

	*this << " ELEMENT     NODE     NODE     NODE     NODE     MATERIAL" << endl
		<< " NUMBER-N      1        2        3        4   SET NUMBER" << endl;

	unsigned int NUME = ElementGroup.GetNUME();

	for (unsigned int Ele = 0; Ele < NUME; Ele++)
	{
		*this << setw(5) << Ele + 1;
		ElementGroup[Ele].Write(*this);
	}

	*this << endl;
}

//	Print load data
void COutputter::OutputLoadInfo()
{
	CDomain* FEMData = CDomain::GetInstance();

	for (unsigned int lcase = 1; lcase <= FEMData->GetNLCASE(); lcase++)
	{
		CLoadCaseData* LoadData = &FEMData->GetLoadCases()[lcase - 1];

		*this << setiosflags(ios::scientific);
		*this << " L O A D   C A S E   D A T A" << endl
			  << endl;

		*this << "     LOAD CASE NUMBER . . . . . . . =" << setw(6) << lcase << endl;
		*this << "     NUMBER OF CONCENTRATED LOADS . =" << setw(6) << LoadData->nloads << endl
			  << endl;
		*this << "    NODE       DIRECTION      LOAD" << endl
			  << "   NUMBER                   MAGNITUDE" << endl;

		LoadData->Write(*this);

		*this << endl;
	}
}

//	Print nodal displacement
void COutputter::OutputNodalDisplacement()
{
	CDomain* FEMData = CDomain::GetInstance();
	CNode* NodeList = FEMData->GetNodeList();
	double* Displacement = FEMData->GetDisplacement();

	*this << setiosflags(ios::scientific);

	*this << " D I S P L A C E M E N T S" << endl
		  << endl;
	*this << "  NODE           X-DISPLACEMENT    Y-DISPLACEMENT    " << endl;

	for (unsigned int np = 0; np < FEMData->GetNUMNP(); np++)
		NodeList[np].WriteNodalDisplacement(*this, Displacement);

	*this << endl;
}


//	Calculate stresses
void COutputter::OutputElementStress()
{
	CDomain* FEMData = CDomain::GetInstance();

	double* Displacement = FEMData->GetDisplacement();

	unsigned int NUMEG = FEMData->GetNUMEG();

	for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++)
	{
		*this << " S T R E S S  C A L C U L A T I O N S  F O R  E L E M E N T  G R O U P" << setw(5)
			  << EleGrpIndex + 1 << endl
			  << endl;

		CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];
		unsigned int NUME = EleGrp.GetNUME();
		ElementTypes ElementType = EleGrp.GetElementType();

		switch (ElementType)
		{
		

		case ElementTypes::Q4:
			*this << "  ELEMENT       SIGMA-X       SIGMA-Y         TAU-XY" << endl;

			for (unsigned int Ele = 0; Ele < NUME; Ele++)
			{
				CElement& Element = EleGrp[Ele];
				double stress[3];
				Element.ElementStress(stress, Displacement);
				*this << setw(5) << Ele + 1 << setw(15) << stress[0]
					<< setw(15) << stress[1] << setw(15) << stress[2] << endl;
			}

			*this << endl;

			break;

			default: // Invalid element type
				cerr << "*** Error *** Elment type " << ElementType
					<< " has not been implemented.\n\n";
		}
	}
}

//	Print total system data
void COutputter::OutputTotalSystemData()
{
	CDomain* FEMData = CDomain::GetInstance();

	*this << "	TOTAL SYSTEM DATA" << endl
		  << endl;

	*this << "     NUMBER OF EQUATIONS . . . . . . . . . . . . . .(NEQ) = " << FEMData->GetNEQ()
		  << endl
		  << "     NUMBER OF MATRIX ELEMENTS . . . . . . . . . . .(NWK) = " << FEMData->GetStiffnessMatrix()->size()
		  << endl
		  << "     MAXIMUM HALF BANDWIDTH  . . . . . . . . . . . .(MK ) = " << FEMData->GetStiffnessMatrix()->GetMaximumHalfBandwidth()
		  << endl
		  << "     MEAN HALF BANDWIDTH . . . . . . . . . . . . . .(MM ) = " << FEMData->GetStiffnessMatrix()->size() / FEMData->GetNEQ() << endl
		  << endl
		  << endl;
}

#ifdef _DEBUG_

//	Print column heights for debuging
void COutputter::PrintColumnHeights()
{
	*this << "*** _Debug_ *** Column Heights" << endl;

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NEQ = FEMData->GetNEQ();
	CSkylineMatrix<double> *StiffnessMatrix = FEMData->GetStiffnessMatrix();
	unsigned int* ColumnHeights = StiffnessMatrix->GetColumnHeights();

	for (unsigned int col = 0; col < NEQ; col++)
	{
		if (col + 1 % 10 == 0)
		{
			*this << endl;
		}

		*this << setw(8) << ColumnHeights[col];
	}

	*this << endl
		  << endl;
}

//	Print address of diagonal elements for debuging
void COutputter::PrintDiagonalAddress()
{
	*this << "*** _Debug_ *** Address of Diagonal Element" << endl;

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NEQ = FEMData->GetNEQ();
	CSkylineMatrix<double> *StiffnessMatrix = FEMData->GetStiffnessMatrix();
	unsigned int* DiagonalAddress = StiffnessMatrix->GetDiagonalAddress();

	for (unsigned int col = 0; col <= NEQ; col++)
	{
		if (col + 1 % 10 == 0)
		{
			*this << endl;
		}

		*this << setw(8) << DiagonalAddress[col];
	}

	*this << endl
		  << endl;
}

//	Print banded and full stiffness matrix for debuging
void COutputter::PrintStiffnessMatrix()
{
	*this << "*** _Debug_ *** Banded stiffness matrix" << endl;

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NEQ = FEMData->GetNEQ();
	CSkylineMatrix<double> *StiffnessMatrix = FEMData->GetStiffnessMatrix();
	unsigned int* DiagonalAddress = StiffnessMatrix->GetDiagonalAddress();

	*this << setiosflags(ios::scientific) << setprecision(5);

	for (unsigned int i = 0; i < DiagonalAddress[NEQ] - DiagonalAddress[0]; i++)
	{
		*this << setw(14) << (*StiffnessMatrix)(i);

		if ((i + 1) % 6 == 0)
		{
			*this << endl;
		}
	}

	*this << endl
		  << endl;

	*this << "*** _Debug_ *** Full stiffness matrix" << endl;

	for (int I = 1; I <= NEQ; I++)
	{
		for (int J = 1; J <= NEQ; J++)
		{
			int J_new = (J > I) ? J : I;
			int I_new = (J > I) ? I : J;
			int H = DiagonalAddress[J_new] - DiagonalAddress[J_new - 1];
			if (J_new - I_new - H >= 0)
			{
				*this << setw(14) << 0.0;
			}
			else
			{
				*this << setw(14) << (*StiffnessMatrix)(I_new, J_new);
			}
		}

		*this << endl;
	}

	*this << endl;
}

void COutputter::OutputVTKFile(const string& filename)
{
	CDomain* FEMData = CDomain::GetInstance();

	ofstream vtkFile(filename);

	if (!vtkFile)
	{
		cerr << "*** Error *** Cannot open VTK file: " << filename << endl;
		return;
	}

	vtkFile << "# vtk DataFile Version 3.0\n";
	vtkFile << "STAP++ Q4 patch output\n";
	vtkFile << "ASCII\n";
	vtkFile << "DATASET UNSTRUCTURED_GRID\n";

	// Write points
	CNode* NodeList = FEMData->GetNodeList();
	unsigned int NUMNP = FEMData->GetNUMNP();

	vtkFile << "POINTS " << NUMNP << " float\n";

	for (unsigned int np = 0; np < NUMNP; np++)
	{
		vtkFile << NodeList[np].XYZ[0] << " "
			<< NodeList[np].XYZ[1] << " "
			<< 0.0 << "\n";
	}

	// Write cells
	unsigned int NUMEG = FEMData->GetNUMEG();
	unsigned int total_cells = 0;
	for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++)
	{
		total_cells += FEMData->GetEleGrpList()[EleGrpIndex].GetNUME();
	}

	vtkFile << "CELLS " << total_cells << " " << total_cells * 5 << "\n";

	for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++)
	{
		CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];
		unsigned int NUME = EleGrp.GetNUME();

		for (unsigned int Ele = 0; Ele < NUME; Ele++)
		{
			CElement& Element = EleGrp[Ele];
			CNode** nodes = Element.GetNodes();

			vtkFile << "4 "
				<< nodes[0]->NodeNumber - 1 << " "
				<< nodes[1]->NodeNumber - 1 << " "
				<< nodes[2]->NodeNumber - 1 << " "
				<< nodes[3]->NodeNumber - 1 << "\n";
		}
	}

	// Write cell types (VTK_QUAD = 9)
	vtkFile << "CELL_TYPES " << total_cells << "\n";
	for (unsigned int i = 0; i < total_cells; i++)
	{
		vtkFile << "9\n";
	}

	// Write point data (displacement)
	double* Displacement = FEMData->GetDisplacement();

	vtkFile << "POINT_DATA " << NUMNP << "\n";
	vtkFile << "VECTORS displacement float\n";

	for (unsigned int np = 0; np < NUMNP; np++)
	{
		double u = (NodeList[np].bcode[0] == 0) ? 0.0 : Displacement[NodeList[np].bcode[0] - 1];
		double v = (NodeList[np].bcode[1] == 0) ? 0.0 : Displacement[NodeList[np].bcode[1] - 1];

		vtkFile << u << " " << v << " " << 0.0 << "\n";
	}

	// Write cell data (stress)
	vtkFile << "CELL_DATA " << total_cells << "\n";
	vtkFile << "SCALARS sigma_x float 1\n";
	vtkFile << "LOOKUP_TABLE default\n";

	for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++)
	{
		CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];
		unsigned int NUME = EleGrp.GetNUME();
		ElementTypes ElementType = EleGrp.GetElementType();

		for (unsigned int Ele = 0; Ele < NUME; Ele++)
		{
			CElement& Element = EleGrp[Ele];
			double stress[3];
			Element.ElementStress(stress, FEMData->GetDisplacement());
			vtkFile << stress[0] << "\n"; // Sigma-x
		}
	}

	vtkFile << "SCALARS sigma_y float 1\n";
	vtkFile << "LOOKUP_TABLE default\n";

	for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++)
	{
		CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];
		unsigned int NUME = EleGrp.GetNUME();

		for (unsigned int Ele = 0; Ele < NUME; Ele++)
		{
			CElement& Element = EleGrp[Ele];
			double stress[3];
			Element.ElementStress(stress, FEMData->GetDisplacement());
			vtkFile << stress[1] << "\n"; // Sigma-y
		}
	}

	vtkFile << "SCALARS tau_xy float 1\n";
	vtkFile << "LOOKUP_TABLE default\n";

	for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++)
	{
		CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];
		unsigned int NUME = EleGrp.GetNUME();

		for (unsigned int Ele = 0; Ele < NUME; Ele++)
		{
			CElement& Element = EleGrp[Ele];
			double stress[3];
			Element.ElementStress(stress, FEMData->GetDisplacement());
			vtkFile << stress[2] << "\n"; // Tau-xy
		}
	}

	vtkFile.close();

	cout << "*** VTK file written to " << filename << endl;
}


//	Print displacement vector for debuging
void COutputter::PrintDisplacement()
{
	*this << "*** _Debug_ *** Displacement vector" << endl;

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NEQ = FEMData->GetNEQ();
	double* Force = FEMData->GetForce();

	*this << setiosflags(ios::scientific) << setprecision(5);

	for (unsigned int i = 0; i < NEQ; i++)
	{
		if ((i + 1) % 6 == 0)
		{
			*this << endl;
		}

		*this << setw(14) << Force[i];
	}

	*this << endl
		  << endl;
}

#endif
