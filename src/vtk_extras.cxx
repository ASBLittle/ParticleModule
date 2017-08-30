#include "vtkPythonArgs.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkCellLocator.h"
#include "stdio.h"
#include "numpy/arrayobject.h"

#include "vtkExtrasErrors.h"
#include "BoundingSurface.h"
#include "Picker.h"

extern "C" {

  static PyObject *extras_find_cell(PyObject *self, PyObject *args) {

    vtkPythonArgs argument_parser(args, "extras_bounding_surface");
    vtkCellLocator *locator;
    double x[3];

    if (!argument_parser.GetVTKObject(locator, "vtkCellLocator")) {
      PyErr_SetString(PyExc_TypeError, "Need VTK unstructured grid as first argument");
      return NULL;
    }
    argument_parser.GetArray(x,3);
    
    // apply our function
    npy_intp dims[1]={3};
    PyObject* pcoords = PyArray_SimpleNew(1,dims,NPY_DOUBLE);
    vtkIdType cellId;
    find_cell(locator, x, cellId, (double*) PyArray_GETPTR1(pcoords,0), 1.0e-6);


    // Now back to Python
    PyObject* out= Py_BuildValue("OO",PyInt_FromLong(cellId),pcoords);
    return out;
  }
  

  static PyObject *extras_bounding_surface(PyObject *self, PyObject *args) {

    vtkPythonArgs argument_parser(args, "extras_bounding_surface");
    vtkUnstructuredGrid *input;    

    if (!argument_parser.GetVTKObject(input, "vtkUnstructuredGrid")) {
      PyErr_SetString(PyExc_TypeError, "Need VTK unstructured grid as first argument");
      return NULL;
    }
    
    // apply our function

    vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::New();
    BoundingSurface* boundary = BoundingSurface::New();
    boundary->GetSurface(input, ugrid);
    boundary->Delete();

    // The object below is what we'll return (this seems to add a reference)
    PyObject* pyugrid = vtkPythonUtil::GetObjectFromPointer(ugrid);

    // Clean up our spare reference now (or you could use smart pointers)
    ugrid->Delete();

    // Now back to Python
    return pyugrid;
  }

  static PyObject *extras_evaluate_field(PyObject *self, PyObject *args) {

    vtkPythonArgs argument_parser(args, "extras_evaluate_field");
    vtkCellLocator *locator;
    vtkUnstructuredGrid *ugrid;
    double x[3];
    char* name;

    if (!argument_parser.GetVTKObject(ugrid, "vtkUnstructuredGrid")) {
      PyErr_SetString(PyExc_TypeError, "Need VTK unstructured grid as first argument");
      return NULL;
    }

    if (!argument_parser.GetVTKObject(locator, "vtkCellLocator")) {
      PyErr_SetString(PyExc_TypeError, "Need VTK cell locator as second argument");
      return NULL;
    }

    argument_parser.GetArray(x,3);
    argument_parser.GetValue(name);
    
    npy_intp dims[1] = {0};
    
    if (ugrid->GetPointData()->HasArray(name)) {
      dims[1] = ugrid->GetPointData()->GetArray(name)->GetNumberOfComponents();
    } else if (ugrid->GetCellData()->HasArray(name)) {
      dims[1] = ugrid->GetCellData()->GetArray(name)->GetNumberOfComponents();
    }

    PyObject* output = PyArray_SimpleNew(1,dims,NPY_DOUBLE);

    // apply our function
    evaluate_field(ugrid, locator, x, name, (double*) PyArray_GETPTR1(output,0), 1.0e-6);


    // Now back to Python
    return output;
  }

  char bounding_surface_docstring[] = "ReadGmsh(vtkUnstructuredGrid) -> vtkUnstructuredGrid\n\n Extract the boundary from a VTK unstructured grid object.";  

  static PyMethodDef extrasMethods[] = {
    { (char *)"BoundingSurface", (PyCFunction) extras_bounding_surface, METH_VARARGS, bounding_surface_docstring},
    { (char *)"FindCell", (PyCFunction) extras_find_cell, METH_VARARGS, bounding_surface_docstring},
    { (char *)"EvaluateField", (PyCFunction) extras_evaluate_field, METH_VARARGS, bounding_surface_docstring},    
    { NULL, NULL, 0, NULL }
  };

  static char extrasDocString[] = "Module collecting python wrappers to VTK stuff.";

  PyMODINIT_FUNC initvtk_extras() {
    
    PyObject* m = Py_InitModule3("vtk_extras", extrasMethods, extrasDocString);
    if (m == NULL) return;
    import_array();

    PyObject* vtk = PyImport_ImportModule("vtk");

    Py_INCREF(vtk);
    PyModule_AddObject(m,"vtk",vtk);
    
    
  }

}
