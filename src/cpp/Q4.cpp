#include "Q4.h"
#include "Material.h"
#include <cmath>
#include <iomanip>

using namespace std;

CQ4::CQ4()
{
    NEN_ = 4;
    nodes_ = new CNode * [NEN_];
    ND_ = 8; // 4 nodes * 2 dofs per node
    LocationMatrix_ = new unsigned int[ND_];
    ElementMaterial_ = nullptr;
}

CQ4::~CQ4()
{
}

bool CQ4::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
    unsigned int MSet, N1, N2, N3, N4;
    Input >> N1 >> N2 >> N3 >> N4 >> MSet;
    ElementMaterial_ = dynamic_cast<CQ4Material*>(MaterialSets) + MSet - 1;
    nodes_[0] = &NodeList[N1 - 1];
    nodes_[1] = &NodeList[N2 - 1];
    nodes_[2] = &NodeList[N3 - 1];
    nodes_[3] = &NodeList[N4 - 1];
    return true;
}

void CQ4::Write(COutputter& output)
{
    output << setw(11) << nodes_[0]->NodeNumber
        << setw(9) << nodes_[1]->NodeNumber
        << setw(9) << nodes_[2]->NodeNumber
        << setw(9) << nodes_[3]->NodeNumber
        << setw(12) << ElementMaterial_->nset << endl;
}

static void ShapeFunc(double xi, double eta, double N[4])
{
    N[0] = 0.25 * (1 - xi) * (1 - eta);
    N[1] = 0.25 * (1 + xi) * (1 - eta);
    N[2] = 0.25 * (1 + xi) * (1 + eta);
    N[3] = 0.25 * (1 - xi) * (1 + eta);
}

static void ShapeDerivs(double xi, double eta, double dNdxi[4], double dNdeta[4])
{
    dNdxi[0] = -0.25 * (1 - eta);
    dNdxi[1] = 0.25 * (1 - eta);
    dNdxi[2] = 0.25 * (1 + eta);
    dNdxi[3] = -0.25 * (1 + eta);

    dNdeta[0] = -0.25 * (1 - xi);
    dNdeta[1] = -0.25 * (1 + xi);
    dNdeta[2] = 0.25 * (1 + xi);
    dNdeta[3] = 0.25 * (1 - xi);
}

void CQ4::ElementStiffness(double* Matrix)
{
    // zero matrix (upper triangular storage later)
    double full[8][8];
    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 8; j++)
            full[i][j] = 0.0;

    CQ4Material* mat = dynamic_cast<CQ4Material*>(ElementMaterial_);
    double E = mat->E;
    double nu = mat->Nu;
    double t = mat->thickness;

    // Constitutive matrix for plane stress
    double coeff = E / (1.0 - nu * nu);
    double D[3][3] = { {coeff, coeff * nu, 0},
                       {coeff * nu, coeff, 0},
                       {0,0, coeff * (1 - nu) / 2.0} };

    // 2x2 Gauss points
    double gp[2] = { -sqrt(1.0 / 3.0), sqrt(1.0 / 3.0) };

    for (int i1 = 0; i1 < 2; i1++)
        for (int i2 = 0; i2 < 2; i2++)
        {
            double xi = gp[i1];
            double eta = gp[i2];
            double dNdxi[4], dNdeta[4];
            ShapeDerivs(xi, eta, dNdxi, dNdeta);

            // Jacobian
            double J[2][2] = { {0,0},{0,0} };
            for (int k = 0; k < 4; k++) {
                J[0][0] += dNdxi[k] * nodes_[k]->XYZ[0];
                J[0][1] += dNdxi[k] * nodes_[k]->XYZ[1];
                J[1][0] += dNdeta[k] * nodes_[k]->XYZ[0];
                J[1][1] += dNdeta[k] * nodes_[k]->XYZ[1];
            }
            double detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
            double invJ[2][2];
            invJ[0][0] = J[1][1] / detJ;
            invJ[0][1] = -J[0][1] / detJ;
            invJ[1][0] = -J[1][0] / detJ;
            invJ[1][1] = J[0][0] / detJ;

            double B[3][8];
            for (int col = 0; col < 8; col++)
                for (int row = 0; row < 3; row++)
                    B[row][col] = 0.0;

            for (int a = 0; a < 4; a++)
            {
                double dNdx = invJ[0][0] * dNdxi[a] + invJ[0][1] * dNdeta[a];
                double dNdy = invJ[1][0] * dNdxi[a] + invJ[1][1] * dNdeta[a];
                B[0][2 * a] = dNdx;
                B[1][2 * a + 1] = dNdy;
                B[2][2 * a] = dNdy;
                B[2][2 * a + 1] = dNdx;
            }

            // compute stiffness k += B^T * D * B * detJ * t * w
            double BtD[8][3];
            for (int i = 0; i < 8; i++)
                for (int j = 0; j < 3; j++) {
                    BtD[i][j] = 0.0;
                    for (int k = 0; k < 3; k++)
                        BtD[i][j] += B[k][i] * D[k][j];
                }
            double k[8][8];
            for (int i = 0; i < 8; i++)
                for (int j = 0; j < 8; j++) {
                    k[i][j] = 0.0;
                    for (int m = 0; m < 3; m++)
                        k[i][j] += BtD[i][m] * B[m][j];
                    k[i][j] *= detJ * t * 1.0; // weight =1 for each gp
                }
            // accumulate into element stiffness matrix
            int map[8] = { 0,1,3,4,6,7,9,10 };
            for (int i = 0; i < 8; i++)
                for (int j = 0; j < 8; j++)
                    full[i][j] += k[i][j];
        }

    // copy to upper triangular packed storage
    int index = 0;
    for (int j = 0; j < 8; j++) {
        for (int i = 0; i <= j; i++) {
            Matrix[index++] = full[i][j];
        }
    }
}

void CQ4::ElementStress(double* stress, double* Displacement)
{
    // stress vector [sx, sy, txy] at element center
    CQ4Material* mat = dynamic_cast<CQ4Material*>(ElementMaterial_);
    double E = mat->E;
    double nu = mat->Nu;
    double coeff = E / (1 - nu * nu);
    double D[3][3] = { {coeff, coeff * nu, 0},
                       {coeff * nu, coeff, 0},
                       {0,0, coeff * (1 - nu) / 2.0} };
    double dNdxi[4], dNdeta[4];
    ShapeDerivs(0, 0, dNdxi, dNdeta);
    double J[2][2] = { {0,0},{0,0} };
    for (int k = 0; k < 4; k++) {
        J[0][0] += dNdxi[k] * nodes_[k]->XYZ[0];
        J[0][1] += dNdxi[k] * nodes_[k]->XYZ[1];
        J[1][0] += dNdeta[k] * nodes_[k]->XYZ[0];
        J[1][1] += dNdeta[k] * nodes_[k]->XYZ[1];
    }
    double detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    double invJ[2][2];
    invJ[0][0] = J[1][1] / detJ;
    invJ[0][1] = -J[0][1] / detJ;
    invJ[1][0] = -J[1][0] / detJ;
    invJ[1][1] = J[0][0] / detJ;

    double B[3][8];
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 8; c++)
            B[r][c] = 0.0;
    for (int a = 0; a < 4; a++) {
        double dNdx = invJ[0][0] * dNdxi[a] + invJ[0][1] * dNdeta[a];
        double dNdy = invJ[1][0] * dNdxi[a] + invJ[1][1] * dNdeta[a];
        B[0][2 * a] = dNdx;
        B[1][2 * a + 1] = dNdy;
        B[2][2 * a] = dNdy;
        B[2][2 * a + 1] = dNdx;
    }

    double u[8];
    for (int i = 0; i < 8; i++) {
        int lm = LocationMatrix_[i];
        u[i] = (lm == 0) ? 0.0 : Displacement[lm - 1];
    }
    double strain[3] = { 0,0,0 };
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 8; j++)
            strain[i] += B[i][j] * u[j];
    for (int i = 0; i < 3; i++) {
        stress[i] = 0.0;
        for (int j = 0; j < 3; j++)
            stress[i] += D[i][j] * strain[j];
    }
}

