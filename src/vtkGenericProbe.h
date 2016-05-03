#include "vtkObject.h"

#include "vtkDoubleArray.h"
#include "vtkStringArray.h"
#include "vtkDataSet.h"
#include "vtkCellLocator.h"

class vtkGenericProbe : public vtkObject
{
  public: 
  static vtkGenericProbe *New();
#if VTK_MAJOR_VERSION <= 5
  vtkTypeRevisionMacro(vtkGenericProbe,vtkObject);
#else
  vtkTypeMacro(vtkGenericProbe,vtkObject);
#endif
  void SetDataSet(vtkDataSet* data);
  vtkDoubleArray* GetValues(double x[3],vtkStringArray* value_namelist,
			    vtkStringArray* gradient_namelist);

  protected:

  vtkGenericProbe();
  ~vtkGenericProbe();

  vtkCellLocator* locator;
  vtkDataSet* dataset;

};
