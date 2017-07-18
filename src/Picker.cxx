#include "Picker.h"
#include "vtkCellLocator.h"
#include "vtkGenericCell.h"

double weights[10];

void find_cell(vtkCellLocator *locator, double* x, vtkIdType &cellId, double* pcoords)
{ 
  vtkGenericCell * cell = vtkGenericCell::New();
  // cellId = locator->FindCell(x,1.0e-16, cell, pcoords, weights);
  int subId;
  double dist2;
  double cp[3];
  locator->FindClosestPoint(x,cp, cell, cellId, subId, dist2);
  cell->EvaluatePosition(cp,NULL,subId, pcoords, dist2, weights);
  if (dist2>1.0e-8) cellId =-1;
  cell->Delete();

  return;
}
