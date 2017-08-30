#include "vtkCellLocator.h"
#include "vtkUnstructuredGrid.h"

void find_cell(vtkCellLocator *, double*, vtkIdType&, double*, double);
bool evaluate_field(vtkUnstructuredGrid*, vtkCellLocator*, double*, char*, double*, double);

