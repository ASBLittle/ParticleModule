#include "vtkCellLocator.h"
#include "vtkDataArray.h"
#include "vtkUnstructuredGrid.h"
#include "vtkGenericCell.h"

void find_cell(vtkAbstractCellLocator *, double*, vtkIdType&, double*, double, vtkGenericCell*);
bool evaluate_field(vtkUnstructuredGrid*, vtkAbstractCellLocator*, double*, char*, double*, double, vtkGenericCell*);
bool evaluate_field(vtkDataArray*, vtkAbstractCellLocator*, double*, double*, double, vtkGenericCell*);

