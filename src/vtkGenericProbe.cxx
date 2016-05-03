#include "vtkGenericProbe.h"

#include "vtkObjectFactory.h"
#include "vtkDataSet.h"
#include "vtkMath.h"
#include "vtkGenericCell.h"
#include "vtkPlane.h"
#include "vtkTriangle.h"
#include "vtkDoubleArray.h"
#include "vtkStringArray.h"
#include "vtkCellLocator.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include <iostream>
#include <map>
#include <set>
#include <math.h> 

#if VTK_MAJOR_VERSION <= 5
vtkCxxRevisionMacro(vtkGenericProbe, "$Revision: 0.0$");
#endif
vtkStandardNewMacro(vtkGenericProbe);

void GetPcoordsTriangle(double x[3],double pcoords[3], vtkGenericCell* cell)
{

  int i, j;
  double pt1[3], pt2[3], pt3[3], n[3], fabsn;
  double rhs[2], c1[2], c2[2];
  double det;
  double maxComponent;
  int idx=0, indices[2];
  double cp[3];

  cell->Points->GetPoint(1, pt1);
  cell->Points->GetPoint(2, pt2);
  cell->Points->GetPoint(0, pt3);

  vtkTriangle::ComputeNormalDirection(pt1, pt2, pt3, n);

  // Project point to plane
  //
  vtkPlane::GeneralizedProjectPoint(x,pt1,n,cp);

  // Construct matrices.  Since we have over determined system, need to find
  // which 2 out of 3 equations to use to develop equations. (Any 2 should
  // work since we've projected point to plane.)
  //
  for (maxComponent=0.0, i=0; i<3; i++)
    {
    // trying to avoid an expensive call to fabs()
    if (n[i] < 0)
      {
      fabsn = -n[i];
      }
    else
      {
      fabsn = n[i];
      }
    if (fabsn > maxComponent)
      {
      maxComponent = fabsn;
      idx = i;
      }
    }
  for (j=0, i=0; i<3; i++)
    {
    if ( i != idx )
      {
      indices[j++] = i;
      }
    }

  for (i=0; i<2; i++)
    {
    rhs[i] = cp[indices[i]] - pt3[indices[i]];
    c1[i] = pt1[indices[i]] - pt3[indices[i]];
    c2[i] = pt2[indices[i]] - pt3[indices[i]];
    }

  pcoords[0] = vtkMath::Determinant2x2(rhs,c2) / det;
  pcoords[1] = vtkMath::Determinant2x2(c1,rhs) / det;
  pcoords[2] = 0.0;
};

void GetPcoordsTetra(double x[3],double pcoords[3],vtkGenericCell* cell)
{
  double pt1[3], pt2[3], pt3[3], pt4[3];
  int i;
  double rhs[3], c1[3], c2[3], c3[3];
  double det, p4;

  pcoords[0] = pcoords[1] = pcoords[2] = 0.0;

  cell->Points->GetPoint(1, pt1);
  cell->Points->GetPoint(2, pt2);
  cell->Points->GetPoint(3, pt3);
  cell->Points->GetPoint(0, pt4);

  for (i=0; i<3; i++)
    {
    rhs[i] = x[i] - pt4[i];
    c1[i] = pt1[i] - pt4[i];
    c2[i] = pt2[i] - pt4[i];
    c3[i] = pt3[i] - pt4[i];
    }

  pcoords[0] = vtkMath::Determinant3x3 (rhs,c2,c3) / det;
  pcoords[1] = vtkMath::Determinant3x3 (c1,rhs,c3) / det;
  pcoords[2] = vtkMath::Determinant3x3 (c1,c2,rhs) / det;
};

vtkGenericProbe::vtkGenericProbe(){
  this->locator = vtkCellLocator::New();
  this->dataset = 0;
};
vtkGenericProbe::~vtkGenericProbe(){
  this->locator->Delete();
};

void vtkGenericProbe::SetDataSet(vtkDataSet* data){
  this->dataset = data;
  this->locator->SetDataSet(this->dataset);
  this->locator->BuildLocator();
};


vtkDoubleArray* vtkGenericProbe::GetValues(double x[3],vtkStringArray* value_namelist,
			  vtkStringArray* gradient_namelist){

  this->locator->BuildLocatorIfNeeded();
  
  double closestPoint[3];
  vtkGenericCell* cell=vtkGenericCell::New();
  vtkIdType cellId;
  int subId;
  double dist2;

  this->locator->FindClosestPoint(x,closestPoint,cell, cellId, subId, dist2);

  vtkDoubleArray* output=vtkDoubleArray::New();

  double* shapes = new double [cell->GetNumberOfPoints()];

  double pcoords[3];

  switch (cell->GetCellType()) 
    {
    case VTK_TRIANGLE:
      GetPcoordsTriangle(x,pcoords,cell);
      break;
    case VTK_QUADRATIC_TRIANGLE:
      GetPcoordsTriangle(x,pcoords,cell);
      break;
    case VTK_TETRA:
      GetPcoordsTetra(x,pcoords,cell);
      break;
    case VTK_QUADRATIC_TETRA:
      GetPcoordsTetra(x,pcoords,cell);
      break;
    }
  
  cell->InterpolateFunctions(pcoords, shapes);

  int index;
  vtkDataArray* pdata;

  for (int i=0; i<value_namelist->GetNumberOfValues(); i++) {
    pdata = this->dataset->GetPointData()->GetArray(value_namelist->GetValue(i), index);
    for  (int j=0; j<pdata->GetNumberOfComponents(); j++) {
      double ddata=0.0;
      for  (int k=0; k<cell->GetNumberOfPoints(); k++) {
	ddata=ddata+shapes[k]*pdata->GetComponent(cell->GetPointIds()->GetId(k),j);
      };
      output->InsertNextValue(ddata);
    }
  }

  delete [] shapes;
  cell->Delete();

  return output;
};



