/* @(#)PickerObject.h
 */

#ifndef _PICKEROBJECT_H
#define _PICKEROBJECT_H 1

#include <Python.h>

#include "vtkCellLocator.h"
#include "vtkDataArray.h"
#include "vtkUnstructuredGrid.h"
#include "vtkGenericCell.h"

extern "C" {

  typedef struct {
    PyObject_HEAD
    vtkCellLocator* locator;
    vtkDataSet* ugrid;
    char data_name[255];
    vtkDoubleArray* data;
    vtkGenericCell* cell;
    vtkIdType cell_index;
    double tol2, pos[3], pcoords[3], weights[10];
  } vtk_extrasPicker;

  PyObject* vtk_extrasPicker_NEW();
  int vtk_extrasPicker_Check(PyObject*);
  PyTypeObject* vtk_extrasPicker_Type();

};

#endif /* _PICKEROBJECT_H */

