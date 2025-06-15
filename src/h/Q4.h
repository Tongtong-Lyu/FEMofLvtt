#pragma once
#include "Element.h"

//! 4-node quadrilateral plane element
class CQ4 : public CElement
{
public:
    CQ4();
    ~CQ4() = default;

    bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList) override;
    void Write(COutputter& output) override;
    void ElementStiffness(double* Matrix) override;
    void ElementStress(double* stress, double* Displacement) override;

    void GenerateLocationMatrix() override;
};