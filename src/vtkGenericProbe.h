#include "vtkObject.h"

#include "vtkDoubleArray.h"
#include "vtkStringArray.h"
#include "vtkDataSet.h"
#include "vtkCellLocator.h"

class vtkGenericProbe : public vtkObject
{
  public: 
  static vtkGenericProbe *New(); 
  vtkTypeRevisionMacro(vtkGenericProbe,vtkObject);

  void SetDataSet(vtkDataSet* data);
  vtkDoubleArray* GetValues(double x[3],vtkStringArray* value_namelist,
			    vtkStringArray* gradient_namelist);

  protected:

  vtkGenericProbe();
  ~vtkGenericProbe();

  vtkCellLocator* locator;
  vtkDataSet* dataset;

};
