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
#include "vtkPointData.h"
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
    vtkCellLocator* locator =  (vtkCellLocator*) vtkPythonUtil::GetPointerFromObject(o,"vtkCellLocator");
    if (self->locator!=locator) {
      if(self->locator) self->locator->Delete();
      self->locator = locator;
      if (locator) {
	self->locator->Register(NULL);
	self->locator->BuildLocator();
      }
    }
    return 0;
  }

  static PyObject* get_grid(PyObject* pyself, void *closure) {
    vtk_extrasPicker *self = (vtk_extrasPicker *)pyself;
    return vtkPythonUtil::GetObjectFromPointer(self->ugrid);
  }

  static int set_grid(PyObject* pyself, PyObject* o, void *closure) {
    vtk_extrasPicker *self = (vtk_extrasPicker *)pyself;
    vtkDataSet* ugrid = (vtkDataSet*) vtkPythonUtil::GetPointerFromObject(o,"vtkDataSet");
    if (self->ugrid!=ugrid) {
      if(self->ugrid) self->ugrid->Delete();
      self->ugrid = ugrid;
      if (ugrid) {
	self->ugrid->Register(NULL);
        if(self->data_name) {
	  self->data = (vtkDoubleArray*) self->ugrid->GetPointData()->GetArray(self->data_name);
	}
      } else return -1;
    }
    return 0;
  } 

  static PyObject* get_name(PyObject* pyself, void *closure) {
    vtk_extrasPicker *self = (vtk_extrasPicker *)pyself;
    if (self->data_name) {
      return PyString_FromString(self->data_name);
    } else {
      Py_RETURN_NONE;
    }
  }

  static int set_name(PyObject* pyself, PyObject* o, void *closure) {
    vtk_extrasPicker *self = (vtk_extrasPicker *)pyself;
    strncpy(self->data_name, PyString_AsString(o), 255);
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
    if (!pyself) return -1;
    vtk_extrasPicker* self = (vtk_extrasPicker*)pyself;

    vtkPythonArgs argument_parser(args, "vtk_extrasPicker_init");
    vtkCellLocator* locator=NULL;

    if(PyObject_Length(args)) {
      if (!argument_parser.GetVTKObject(locator, "vtkCellLocator")) {
	PyErr_SetString(PyExc_TypeError, "Need VTK cell locator as first argument");
	return -1;
      }						     
      self->locator = locator;
      self->locator->Register(NULL);
    } else {
      self->locator = NULL;
    }
    self->cell = vtkGenericCell::New();
    self->tol2=1.0e-32;
    for (int i=0; i<3; ++i) {
      self->pos[i] = 0.0;
    }
    if (locator && self->locator->GetDataSet()) {
      self->locator->BuildLocator();
      self->cell_index = self->locator->FindCell(self->pos, self->tol2, self->cell, self->pcoords, self->weights);
      self->ugrid = self->locator->GetDataSet();
      self->ugrid->Register(NULL);
      self->data=NULL;
      strcpy(self->data_name, "Velocity");
    } else {
      self->cell_index = -1;
      self->ugrid = NULL;
      self->data = NULL;
      strcpy(self->data_name, "");
    }
    return 0;
  }

  static void vtk_extrasPicker_dealloc(PyObject *self)
  {
    vtk_extrasPicker *p = (vtk_extrasPicker *)self;
    if (p->locator) p->locator->Delete();
    if (p->cell) p->cell->Delete();
    if (p->ugrid) p->ugrid->Delete();
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

#if VTK_MAJOR_VERSION<6
    void * tmp;
    argument_parser.GetValue(tmp);
    pyobj = PyTuple_GetItem(args, 0) ;
#else
    if (!argument_parser.GetPythonObject(pyobj)) {
      PyErr_SetString(PyExc_TypeError, "Issue with first argument.");
      return NULL;
    }
#endif

    vtkDoubleArray* data;
    if (PyObject_Length(args)>1) {
      if (!argument_parser.GetVTKObject(data, "vtkDoubleArray")|| !data) {
	PyErr_SetString(PyExc_TypeError, "Need VTK double data array as second argument");
	return NULL;
      }
    } else {
      data = p->data;
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

    npy_intp dims[1] = { 3 };
    PyObject* output = PyArray_SimpleNew(1,dims,NPY_DOUBLE);
    double* cout = (double*) PyArray_GETPTR1((PyArrayObject*)output,0);

    for (int j=0; j<3; ++j) {
      cout[j] = 0.0;
    }
    for (int j=0; j<data->GetNumberOfComponents(); ++j) {
      for (int i=0; i<p->cell->GetNumberOfPoints(); ++i) {
	cout[j] = cout[j] + p->weights[i]*data->GetComponent(p->cell->GetPointId(i),j) ;
      }
    }

    return (PyObject*) output;
    
  }

  PyObject* vtk_extrasPickerNearest(PyObject *self, PyObject *args, PyObject *kw)
  {

    vtkPythonArgs argument_parser(args, "vtk_extrasPicker_call");
    vtk_extrasPicker *p= (vtk_extrasPicker*) self;
    PyObject *pyobj, *pypos;

#if VTK_MAJOR_VERSION<6
    void * tmp;
    argument_parser.GetValue(tmp);
    pyobj = PyTuple_GetItem(args, 0) ;
#else
    if (!argument_parser.GetPythonObject(pyobj)) {
      PyErr_SetString(PyExc_TypeError, "Issue with first argument.");
      return NULL;
    }
#endif

    vtkDoubleArray* data;
    if (PyObject_Length(args)>1) {
      if (!argument_parser.GetVTKObject(data, "vtkDoubleArray")|| !data) {
	PyErr_SetString(PyExc_TypeError, "Need VTK double data array as second argument");
	return NULL;
      }
    } else {
      data = p->data;
    }
    
    pypos = PyArray_FromObject(pyobj,NPY_DOUBLE,1,1);
    double* pos=(double*) PyArray_GETPTR1((PyArrayObject*)pypos,0);
    if (!pos) {
      PyErr_SetString(PyExc_TypeError, "Cannot convert first argument to numpy array.");
      return NULL;
    }
    
    double closestPoint[3];
    int subId;
    double dist2;
    if (!p->locator->GetDataSet()) {
      p->cell_index = -1;
    } else if (p->cell_index == -1 || do_update(p, pos)) {
      p->locator->BuildLocatorIfNeeded();
      p->locator->FindClosestPoint(pos, closestPoint, p->cell, p->cell_index, subId, dist2);
    }

    if (p->cell_index == -1 || dist2>p->tol2) {
      Py_RETURN_NONE;
    }

    int subid=0;
    double x[3];
    if (p->locator->GetDataSet() != p->ugrid) {
      p->cell->EvaluatePosition(pos, NULL, subId, p->pcoords, dist2, p->weights);
    }

    Py_XDECREF(pypos);

    npy_intp dims[1] = { 3 };
    PyObject* output = PyArray_SimpleNew(1,dims,NPY_DOUBLE);
    double* cout = (double*) PyArray_GETPTR1((PyArrayObject*)output,0);

    for (int j=0; j<3; ++j) {
      cout[j] = 0.0;
    }
    for (int j=0; j<data->GetNumberOfComponents(); ++j) {
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

  static PyMethodDef vtk_extrasPicker_methods[] = {
    { (char*)"nearest", (PyCFunction) vtk_extrasPickerNearest, METH_VARARGS|METH_KEYWORDS, (char*) "Get value at nearest point in domain."},
    {NULL}
  };

  static PyGetSetDef vtk_extrasPicker_getset[] = {
    { (char*)"locator", get_locator, set_locator, (char*)"Get/set picker locator.", NULL},
    { (char*)"grid", get_grid, set_grid, (char*)"Get/set picker grid.", NULL},
    { (char*)"name", get_name, set_name, (char*)"Get/set data name.", NULL},
    { (char*)"cell_index", get_cell_index, NULL, (char*)"Get picker cell index.", NULL},
    { (char*)"cell", get_cell, NULL, (char*)"Get picker cell.", NULL},
    { (char*)"pos", get_pos, (setter)set_pos, (char*)"Get/set picker position.", NULL},
    { NULL}
  };

  static PyTypeObject vtk_extrasPickerType = {
    PyVarObject_HEAD_INIT(NULL, 0)
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
    vtk_extrasPicker_methods,                    /* tp_methods */
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
    0 /* tp_is_gc For PyObject_IS_GC */
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
