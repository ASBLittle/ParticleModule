#include <Python.h>

#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL PickerObject_ARRAY_API
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"

#include <iostream>
#include <vector>
#include <numeric>

#include "vtkPythonUtil.h"
#include "vtkPythonArgs.h"
#include "vtkDoubleArray.h"
#include "PickerObject.h"


extern "C" {

  static PyObject* vtk_extrasPicker_new(PyTypeObject* type, PyObject *args, PyObject *kw)
  {
    vtk_extrasPicker *self = (vtk_extrasPicker *) type->tp_alloc(type, 0);
    return (PyObject *) self;

  }

  static PyObject* get_locator(PyObject* pyself, void *closure) {
    vtk_extrasPicker *self = (vtk_extrasPicker *)pyself;
    return vtkPythonUtil::GetObjectFromPointer(self->locator);
  }

  static int set_locator(PyObject* pyself, PyObject* o, void *closure) {
    vtk_extrasPicker *self = (vtk_extrasPicker *)pyself;
    if (self->locator) {
      self->locator->Delete();
      self->locator=NULL;
    }
    self->locator = (vtkCellLocator*) vtkPythonUtil::GetPointerFromObject(o,"vtkCellLocator");
    self->locator->BuildLocatorIfNeeded();
    self->locator->Register(NULL);
    return 0;
  }

  static PyObject* get_grid(PyObject* pyself, void *closure) {
    vtk_extrasPicker *self = (vtk_extrasPicker *)pyself;
    return vtkPythonUtil::GetObjectFromPointer(self->ugrid);
  }

  static int set_grid(PyObject* pyself, PyObject* o, void *closure) {
    vtk_extrasPicker *self = (vtk_extrasPicker *)pyself;
    if (self->ugrid) {
      self->ugrid->Delete();
      self->ugrid=NULL;
    }
    self->ugrid = (vtkDataSet*) vtkPythonUtil::GetPointerFromObject(o,"vtkDataSet");
    return 0;
  }

  static PyObject* get_cell_index(PyObject* pyself, void *closure) {
    vtk_extrasPicker *self = (vtk_extrasPicker *)pyself;
    if (self->cell_index >= 0) {
      return PyInt_FromLong(self->cell_index);
    } else {
      Py_RETURN_NONE;
    }				 
  }

  static PyObject* get_cell(PyObject* pyself, void *closure) {
    vtk_extrasPicker *self = (vtk_extrasPicker *)pyself;
    return vtkPythonUtil::GetObjectFromPointer(self->cell);
  }

  static PyObject* get_pos(PyObject* pyself, void *closure) {
    vtk_extrasPicker *self = (vtk_extrasPicker *)pyself;
    npy_intp dims[1] = {3};
    return PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, self->pos);
  }

  static int set_pos(PyObject* pyself, PyObject* o, void *closure) {
    vtk_extrasPicker *self = (vtk_extrasPicker *)pyself;
    PyObject* pypos = PyArray_FromObject(o,NPY_DOUBLE,1,1);
    double* pos=(double*) PyArray_GETPTR1((PyArrayObject*)pypos,0);
    for (int i=0; i<3; ++i) {
      self->pos[i] = pos[i];
    }
    self->cell_index = self->locator->FindCell(self->pos, self->tol2, self->cell, self->pcoords, self->weights);
    return 0;
  }

  static int vtk_extrasPicker_init(PyObject *pyself, PyObject *args, PyObject *kw)
  { 
    vtkPythonArgs argument_parser(args, "vtk_extrasPicker_init");
    vtkCellLocator* locator=NULL;
    
    if (!argument_parser.GetVTKObject(locator, "vtkCellLocator")) {
      PyErr_SetString(PyExc_TypeError, "Need VTK cell locator as first argument");
      return -1;
    }
 
    if (!pyself) return -1;
    vtk_extrasPicker* self = (vtk_extrasPicker*)pyself;
						     
    self->locator = locator;
    self->locator->Register(NULL);
    self->cell = vtkGenericCell::New();
    self->tol2=1.0e-8;
    for (int i=0; i<3; ++i) {
      self->pos[i] =0.0;
    }
    if (self->locator->GetDataSet()) {
      self->locator->BuildLocator();
      self->cell_index = self->locator->FindCell(self->pos, self->tol2, self->cell, self->pcoords, self->weights);
      self->ugrid = self->locator->GetDataSet();
      self->ugrid->Register(NULL);
    } else {
      self->cell_index = -1;
      self->ugrid = NULL;
    }
    return 0;
  }

  static void vtk_extrasPicker_dealloc(PyObject *self)
  {
    vtk_extrasPicker *p = (vtk_extrasPicker *)self;
    p->locator->Delete();
    p->cell->Delete();
    Py_TYPE(self)->tp_free(self);
  }

  static bool do_update(vtk_extrasPicker* picker, double x[3]) {
    for (int i=0; i<3; ++i) {
	if (x[i] != picker->pos[i]) return true;
      }
    return false;
  } 

  PyObject* vtk_extrasPicker_call(PyObject *self, PyObject *args, PyObject *kw)
  {

    vtkPythonArgs argument_parser(args, "vtk_extrasPicker_call");
    vtk_extrasPicker *p= (vtk_extrasPicker*) self;
    PyObject *pyobj, *pypos;

    if (!argument_parser.GetPythonObject(pyobj)) {
      PyErr_SetString(PyExc_TypeError, "Issue with first argument.");
      return NULL;
    }

    vtkDoubleArray* data;
    if (!argument_parser.GetVTKObject(data, "vtkDoubleArray")|| !data) {
      PyErr_SetString(PyExc_TypeError, "Need VTK double data array as second argument");
      return NULL;
    }
    
    pypos = PyArray_FromObject(pyobj,NPY_DOUBLE,1,1);
    double* pos=(double*) PyArray_GETPTR1((PyArrayObject*)pypos,0);
    if (!pos) {
      PyErr_SetString(PyExc_TypeError, "Cannot convert first argument to numpy array.");
      return NULL;
    }
    
    if (!p->locator->GetDataSet()) {
      p->cell_index = -1;
    } else if (p->cell_index == -1 || do_update(p, pos)) {
      p->locator->BuildLocatorIfNeeded();
      p->cell_index = p->locator->FindCell(pos, p->tol2, p->cell, p->pcoords, p->weights);
    }

    if (p->cell_index == -1) {
      Py_RETURN_NONE;
    }

    int subid=0;
    double x[3];
    if (p->locator->GetDataSet() != p->ugrid) {
      p->ugrid->GetCell(p->cell_index, p->cell);
      p->cell->EvaluateLocation(subid, p->pcoords, x, p->weights);
    }

    Py_XDECREF(pypos);

    npy_intp dims[1] = { data->GetNumberOfComponents() };
    PyObject* output = PyArray_SimpleNew(1,dims,NPY_DOUBLE);
    double* cout = (double*) PyArray_GETPTR1((PyArrayObject*)output,0);

    for (int j=0; j<dims[0]; ++j) {
      cout[j] = 0.0;
      for (int i=0; i<p->cell->GetNumberOfPoints(); ++i) {
	cout[j] = cout[j] + p->weights[i]*data->GetComponent(p->cell->GetPointId(i),j) ;
      }
    }

    return (PyObject*) output;
    
  }

  static PyObject* vtk_extrasPicker_repr(PyObject* self) {
    vtk_extrasPicker *p = (vtk_extrasPicker *)self;
    char buf[500];
    sprintf(buf, "%s ", Py_TYPE(self)->tp_name) ;
    return PyUnicode_FromString(buf);
  }

  static PyGetSetDef vtk_extrasPicker_getset[] = {
    { (char*)"locator", get_locator, set_locator, (char*)"Get/set picker locator.", NULL},
    { (char*)"grid", get_grid, set_grid, (char*)"Get/set picker grid.", NULL},
    { (char*)"cell_index", get_cell_index, NULL, (char*)"Get picker cell index.", NULL},
    { (char*)"cell", get_cell, NULL, (char*)"Get picker cell.", NULL},
    { (char*)"pos", get_pos, (setter)set_pos, (char*)"Get/set picker position.", NULL},
    { NULL}
  };

  static PyTypeObject vtk_extrasPickerType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "vtk_extras.Picker",            /*tp_name*/
    sizeof(vtk_extrasPicker),       /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    vtk_extrasPicker_dealloc,       /*tp_dealloc*/
    0,                         /*tp_print PSEUDODEPRECATED*/
    0,                         /*tp_getattr DEPRECATED*/
    0,                         /*tp_setattr DEPRECATED*/
    0,                         /*tp_compare*/
    vtk_extrasPicker_repr,          /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    vtk_extrasPicker_call,          /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,        /* tp_flags */
    "Picker object to interogate ugrid data field by field.", /* tp_doc */
    0,  /* tp_traverse */
    0,  /* tp_clear */
    0,  /* tp_richcompare */
    0,  /* tp_weaklistoffset */
    0,  /* tp_iter: __iter__() method */
    0,  /* tp_iternext: next() method */
    0,                    /* tp_methods */
    0,                    /* tp_members */
    vtk_extrasPicker_getset, /* tp_getset */
    0, /* tp_base */
    0, /* tp_dict */
    0, /* tp_descr_get */
    0, /* tp_descr_set */
    0, /* tp_dictoffset */
    vtk_extrasPicker_init, /* tp_init */
    0, /* tp_alloc */
    vtk_extrasPicker_new, /* tp_new */
    0, /* tp_free Low-level free-memory routine */
    0, /* tp_is_gc For PyObject_IS_GC */
  };

  PyObject* vtk_extrasPicker_NEW(){
    vtk_extrasPicker* object = PyObject_NEW(vtk_extrasPicker, &vtk_extrasPickerType);
    return (PyObject*) object;
  }

  int vtk_extrasPicker_Check(PyObject* ob)
  {
    return PyObject_IsInstance(ob, (PyObject*) &vtk_extrasPickerType);
  }
  
  PyTypeObject* vtk_extrasPicker_Type() {
    return &vtk_extrasPickerType;
  }
}
