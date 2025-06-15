#include "Q4.h"
#include "Material.h"
#include <cmath>
#include <iomanip>

CQ4::CQ4()
{
    NEN_ = 4;
    nodes_ = new CNode*[NEN_];
    ND_ = 8; // 2 DOF per node
    LocationMatrix_ = new unsigned int[ND_];
    ElementMaterial_ = nullptr;
}

void CQ4::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN_; ++N)
        for (unsigned int d = 0; d < 2; ++d)
            LocationMatrix_[i++] = nodes_[N]->bcode[d];
}

bool CQ4::Read(ifstream & Input, CMaterial * MaterialSets, CNode * NodeList)
{
    unsigned int MSet;
    unsigned int N1, N2, N3, N4;
    Input >> N1 >> N2 >> N3 >> N4 >> MSet;
    ElementMaterial_ = dynamic_cast<CQ4Material*>(MaterialSets) + MSet - 1;
    nodes_[0] = &NodeList[N1 - 1];
    nodes_[1] = &NodeList[N2 - 1];
    nodes_[2] = &NodeList[N3 - 1];
    nodes_[3] = &NodeList[N4 - 1];
    return true;
}

void CQ4::Write(COutputter & output)
{
    output << std::setw(11) << nodes_[0]->NodeNumber
           << std::setw(9) << nodes_[1]->NodeNumber
           << std::setw(9) << nodes_[2]->NodeNumber
           << std::setw(9) << nodes_[3]->NodeNumber
           << std::setw(12) << ElementMaterial_->nset << std::endl;
}

void CQ4::ElementStiffness(double * Matrix)
{
    clear(Matrix, SizeOfStiffnessMatrix());

    CQ4Material* material = dynamic_cast<CQ4Material*>(ElementMaterial_);
    double E = material->E;
    double nu = material->nu;
    double t = material->thickness;

    double x[4], y[4];
    for (unsigned int i = 0; i < 4; i++)
    {
        x[i] = nodes_[i]->XYZ[0];
        y[i] = nodes_[i]->XYZ[1];
    }

    double D[3][3];
    double coeff = E / (1.0 - nu * nu);
    D[0][0] = coeff; D[0][1] = coeff * nu; D[0][2] = 0.0;
    D[1][0] = coeff * nu; D[1][1] = coeff; D[1][2] = 0.0;
    D[2][0] = 0.0; D[2][1] = 0.0; D[2][2] = coeff * (1.0 - nu) / 2.0;

    double gp[2] = { -std::sqrt(1.0 / 3.0), std::sqrt(1.0 / 3.0) };
    double full[8][8] = { 0.0}
    ;

    for (int ii = 0; ii < 2; ++ii)
    {
        for (int jj = 0; jj < 2; ++jj)
        {
            double xi = gp[ii];
            double eta = gp[jj];

            double dNdxi[4] = { -(1 - eta) / 4.0, (1 - eta) / 4.0, (1 + eta) / 4.0, -(1 + eta) / 4.0 };
            double dNdeta[4] = { -(1 - xi) / 4.0, -(1 + xi) / 4.0, (1 + xi) / 4.0, (1 - xi) / 4.0 };

            double J[2][2] = { { 0.0,0.0},{ 0.0,0.0} }
            ;
            for (int k = 0; k < 4; k++)
            {
                J[0][0] += dNdxi[k] * x[k];
                J[0][1] += dNdxi[k] * y[k];
                J[1][0] += dNdeta[k] * x[k];
                J[1][1] += dNdeta[k] * y[k];
            }
            double detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
            double invJ[2][2];
            invJ[0][0] = J[1][1] / detJ;
            invJ[0][1] = -J[0][1] / detJ;
            invJ[1][0] = -J[1][0] / detJ;
            invJ[1][1] = J[0][0] / detJ;

            double B[3][8];
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 8; j++)
                    B[i][j] = 0.0;

            for (int k = 0; k < 4; k++)
            {
                double dNdx = invJ[0][0] * dNdxi[k] + invJ[0][1] * dNdeta[k];
                double dNdy = invJ[1][0] * dNdxi[k] + invJ[1][1] * dNdeta[k];
                B[0][2 * k] = dNdx;
                B[1][2 * k + 1] = dNdy;
                B[2][2 * k] = dNdy;
                B[2][2 * k + 1] = dNdx;
            }

            double BtD[8][3];
            for (int i = 0; i < 8; i++)
                for (int j = 0; j < 3; j++)
                {
                    BtD[i][j] = 0.0;
                    for (int k = 0; k < 3; k++)
                        BtD[i][j] += B[k][i] * D[k][j];
                }

            for (int i = 0; i < 8; i++)
                for (int j = 0; j < 8; j++)
                {
                    double val = 0.0;
                    for (int k = 0; k < 3; k++)
                        val += BtD[i][k] * B[k][j];
                    full[i][j] += val * detJ * t;
                }
        }
    }

    int index = 0;
    for (int i = 0; i < 8; i++)
        for (int j = 0; j <= i; j++)
            Matrix[index++] = full[i][j];
}

void CQ4::ElementStress(double * stress, double * Displacement)
{
    CQ4Material* material = dynamic_cast<CQ4Material*>(ElementMaterial_);
    double E = material->E;
    double nu = material->nu;

    double x[4], y[4];
    for (unsigned int i = 0; i < 4; i++)
    {
        x[i] = nodes_[i]->XYZ[0];
        y[i] = nodes_[i]->XYZ[1];
    }

    // compute strain at center (xi=0, eta=0)
    double dNdxi[4] = { -0.25, 0.25, 0.25, -0.25 };
    double dNdeta[4] = { -0.25, -0.25, 0.25, 0.25 };

    double J[2][2] = { { 0.0,0.0},{ 0.0,0.0} }
    ;
    for (int k = 0; k < 4; k++)
    {
        J[0][0] += dNdxi[k] * x[k];
        J[0][1] += dNdxi[k] * y[k];
        J[1][0] += dNdeta[k] * x[k];
        J[1][1] += dNdeta[k] * y[k];
    }
    double detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    double invJ[2][2];
    invJ[0][0] = J[1][1] / detJ;
    invJ[0][1] = -J[0][1] / detJ;
    invJ[1][0] = -J[1][0] / detJ;
    invJ[1][1] = J[0][0] / detJ;

    double B[3][8] ={ 0}
    ;
    for (int k = 0; k < 4; k++)
    {
        double dNdx = invJ[0][0] * dNdxi[k] + invJ[0][1] * dNdeta[k];
        double dNdy = invJ[1][0] * dNdxi[k] + invJ[1][1] * dNdeta[k];
        B[0][2 * k] = dNdx;
        B[1][2 * k + 1] = dNdy;
        B[2][2 * k] = dNdy;
        B[2][2 * k + 1] = dNdx;
    }

    double strain[3] = { 0 };
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 8; j++)
            if (LocationMatrix_[j])
                strain[i] += B[i][j] * Displacement[LocationMatrix_[j] - 1];

    double D[3][3];
    double coeff = E / (1.0 - nu * nu);
    D[0][0] = coeff; D[0][1] = coeff * nu; D[0][2] = 0.0;
    D[1][0] = coeff * nu; D[1][1] = coeff; D[1][2] = 0.0;
    D[2][0] = 0.0; D[2][1] = 0.0; D[2][2] = coeff * (1.0 - nu) / 2.0;

    for (int i = 0; i < 3; i++)
    {
        stress[i] = 0.0;
        for (int j = 0; j < 3; j++)
            stress[i] += D[i][j] * strain[j];
    }
}