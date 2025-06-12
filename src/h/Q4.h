#pragma once

#include "Element.h"

//! 4 node quadrilateral element (plane stress/strain)
class CQ4 : public CElement
{
public:
    CQ4();
    ~CQ4();

    //! Read element data from stream Input
    virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

    //! Write element data to stream
    virtual void Write(COutputter& output);

    //! Calculate element stiffness matrix
    virtual void ElementStiffness(double* Matrix);

    //! Calculate element stress at the element center
    virtual void ElementStress(double* stress, double* Displacement);
};
