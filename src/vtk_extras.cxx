#include "Python.h"

#define PY_ARRAY_UNIQUE_SYMBOL PickerObject_ARRAY_API
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"

#include "vtkPythonArgs.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkGenericCell.h"
#include "vtkCellData.h"
#include "vtkCellLocator.h"
#include "vtkIdList.h"
#include "stdio.h"

#include "vtkExtrasErrors.h"
#include "BoundingSurface.h"
#include "Picker.h"
#include "PickerObject.h"

extern "C" {

  static PyObject* pyVersion = PyString_FromString(GIT_VERSION_STRING);
  
  static vtkGenericCell* cell;

  static PyObject *extras_find_cell(PyObject *self, PyObject *args) {

    vtkPythonArgs argument_parser(args, "extras_bounding_surface");
    vtkCellLocator *locator;
    double x[3];

    if (!argument_parser.GetVTKObject(locator, "vtkAbstractCellLocator")) {
      PyErr_SetString(PyExc_TypeError, "Need VTK unstructured grid as first argument");
      return NULL;
    }
    argument_parser.GetArray(x,3);
    
    // apply our function
    npy_intp dims[1]={3};
    PyObject* pcoords = PyArray_SimpleNew(1,dims,NPY_DOUBLE);
    vtkIdType cellId;
    find_cell(locator, x, cellId, (double*) PyArray_GETPTR1((PyArrayObject*)pcoords,0), 1.0e-6, cell);


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

  static PyObject *extras_vInterpolate(PyObject *self, PyObject *args) {

    vtkPythonArgs argument_parser(args, "extras_vInterpolate");
    vtkDoubleArray *data;
    vtkIdList *ids;
    PyObject* weights;

    if (!argument_parser.GetVTKObject(data, "vtkDoubleArray")) {
      PyErr_SetString(PyExc_TypeError, "Need VTK double array as first argument");
      return NULL;
    }

    if (!argument_parser.GetVTKObject(ids, "vtkIdList")) {
      PyErr_SetString(PyExc_TypeError, "Need VTK id list as second argument");
      return NULL;
    }

    weights = PyTuple_GetItem(args,2);

    npy_intp dims[1] = {0};

    dims[0] = data->GetNumberOfComponents();

    PyObject* output = PyArray_ZEROS(1,dims,NPY_DOUBLE,0);

    double *doutput = (double*) PyArray_GETPTR1((PyArrayObject*)output,0);
    double *dweights = (double*) PyArray_GETPTR1((PyArrayObject*)weights,0);

    for (int c=0; c<ids->GetNumberOfIds(); ++c) {
      for (int i=0; i<dims[0]; ++i) {
	doutput[i] = doutput[i] 
	  + dweights[c] * data->GetComponent(ids->GetId(c), i);
      }	
    }



    return output;

  }

  char vInterpolate_docstring[] = "vInterpolate(vtkDoubleArray, vtkIdList, ndarray weights) -> ndarray result\n\n Interpolate values from a VTK double array.";  


  static PyObject *extras_evaluate_field(PyObject *self, PyObject *args) {

    vtkPythonArgs argument_parser(args, "extras_evaluate_field");
    vtkAbstractCellLocator *locator;
    vtkObject *tmp;
    vtkDataArray* data = NULL;
    double x[3];
    char* name;

    if (!argument_parser.GetVTKObject(tmp, "vtkObject")) {
      PyErr_SetString(PyExc_TypeError, "Need VTK unstructured grid or Data Array as first argument");
      return NULL;
    }

    if (!argument_parser.GetVTKObject(locator, "vtkAbstractCellLocator")) {
      PyErr_SetString(PyExc_TypeError, "Need VTK cell locator as second argument");
      return NULL;
    }

    argument_parser.GetArray(x,3);
    argument_parser.GetValue(name);

    if (tmp->IsA("vtkUnstructuredGrid")) {
      data = ((vtkUnstructuredGrid*)tmp)->GetPointData()->GetArray(name);
      if (!data) data = ((vtkUnstructuredGrid*)tmp)->GetCellData()->GetArray(name);
    } else if (tmp->IsA("vtkDataArray")) {
      data = (vtkDataArray*) tmp;
    } else {
      PyErr_SetString(PyExc_TypeError, "Need VTK unstructured grid or Data Array as first argument");
      return NULL;
    }
    if (!data) {
      PyErr_SetString(PyExc_TypeError, "Can't find data!");
      return NULL;
    }
    
    npy_intp dims[1] = { data->GetNumberOfComponents() };
    PyObject* output = PyArray_SimpleNew(1,dims,NPY_DOUBLE);

    // apply our function
    evaluate_field(data, locator, x, (double*) PyArray_GETPTR1((PyArrayObject*)output,0), 1.0e-6, cell);

    // Now back to Python
    return output;
  }

  char bounding_surface_docstring[] = "ReadGmsh(vtkUnstructuredGrid) -> vtkUnstructuredGrid\n\n Extract the boundary from a VTK unstructured grid object.";    

  static PyMethodDef extrasMethods[] = {
    { (char *)"BoundingSurface", (PyCFunction) extras_bounding_surface, METH_VARARGS, bounding_surface_docstring},
    { (char *)"FindCell", (PyCFunction) extras_find_cell, METH_VARARGS, bounding_surface_docstring},
    { (char *)"EvaluateField", (PyCFunction) extras_evaluate_field, METH_VARARGS, bounding_surface_docstring},  
    { (char *)"vInterpolate", (PyCFunction) extras_vInterpolate, METH_VARARGS, vInterpolate_docstring},    
    { NULL, NULL, 0, NULL }
  };

  static char extrasDocString[] = "Module collecting python wrappers to VTK stuff.";

#if PY_MAJOR_VERSION >= 3
    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "vtk_extras",     /* m_name */
        extrasDocString,  /* m_doc */
        -1,                  /* m_size */
        extrasMethods,    /* m_methods */
        NULL,                /* m_reload */
        NULL,                /* m_traverse */
        NULL,                /* m_clear */
        NULL,                /* m_free */
    };
#endif

#if PY_MAJOR_VERSION >= 3
  PyMODINIT_FUNC PyInit_vtk_extras() {

    PyObject* m = PyModule_Create(&moduledef);
    if (m == NULL) return NULL;
    import_array();

    PyModule_AddObject(m, "__version__", pyVersion);

    if (PyType_Ready(vtk_extrasPicker_Type()) < 0)  return NULL;

    cell = vtkGenericCell::New();

    PyObject* vtk = PyImport_ImportModule("vtk");

    Py_INCREF(vtk);
    PyModule_AddObject(m,"vtk",vtk);
    Py_INCREF(vtk_extrasPicker_Type());
    PyModule_AddObject(m,"Picker", (PyObject*)vtk_extrasPicker_Type());

    return m;
  }

#else
  PyMODINIT_FUNC initvtk_extras() {


    PyObject* m = Py_InitModule3("vtk_extras", extrasMethods, extrasDocString);
    if (m == NULL) return;
    import_array();

    PyModule_AddObject(m, "__version__", pyVersion);

    if (PyType_Ready(vtk_extrasPicker_Type()) < 0)  return;

    cell = vtkGenericCell::New();

    PyObject* vtk = PyImport_ImportModule("vtk");

    Py_INCREF(vtk);
    PyModule_AddObject(m,"vtk",vtk);
    Py_INCREF(vtk_extrasPicker_Type());
    PyModule_AddObject(m,"Picker", (PyObject*)vtk_extrasPicker_Type());
  }
#endif
}
