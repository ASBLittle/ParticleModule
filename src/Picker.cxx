#include "Picker.h"
#include "vtkCellLocator.h"
#include "vtkUnstructuredGrid.h"
#include "vtkGenericCell.h"
#include "vtkPointData.h"
#include "vtkCellData.h"

double weights[10];

void find_cell(vtkCellLocator *locator, double* x, vtkIdType &cellId,
		 double* pcoords, double tol=1.0e-6)
{ 
  vtkGenericCell * cell = vtkGenericCell::New();
  // cellId = locator->FindCell(x,1.0e-16, cell, pcoords, weights);
  int subId;
  double dist2;
  double cp[3];
  locator->FindClosestPoint(x,cp, cell, cellId, subId, dist2);
  cell->EvaluatePosition(cp,NULL,subId, pcoords, dist2, weights);
  if (dist2>tol) cellId =-1;
  cell->Delete();

  return;
}

bool evaluate_field(vtkUnstructuredGrid *ugrid, vtkCellLocator *locator, double* x, 
		    char* name, double* output, double tol=1.0e-8)
{ 
  vtkGenericCell * cell = vtkGenericCell::New();
  // cellId = locator->FindCell(x,1.0e-16, cell, pcoords, weights);
  int subId;
  vtkIdType cellId;
  double dist2;
  double cp[3], pcoords[3];
  locator->FindClosestPoint(x,cp, cell, cellId, subId, dist2);
  cell->EvaluatePosition(cp,NULL,subId, pcoords, dist2, weights);

  vtkDataArray* data = ugrid->GetPointData()->GetArray(name);
  if(!data) data = ugrid->GetCellData()->GetArray(name);

  if (data) {
    for (int i=0; i<data->GetNumberOfTuples(); ++i)
      {
	output[i] = 0.0;
	
	for (int j=0; j <cell->GetNumberOfPoints(); ++j) 
	  {
	    output[i] += weights[j]*data->GetComponent(cell->GetPointIds()->GetId(j), i);
	  }
      }
  }
  return data != NULL && dist2 >= tol;
}
