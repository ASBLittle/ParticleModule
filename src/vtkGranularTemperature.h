#include "vtkPolyDataAlgorithm.h"

class vtkGranularTemperature : public vtkPolyDataAlgorithm
{
public:
  static vtkGranularTemperature *New();
#if VTK_MAJOR_VERSION <= 5
  vtkTypeRevisionMacro(vtkGranularTemperature,vtkPolyDataAlgorithm);
#else
  vtkTypeMacro(vtkGranularTemperature,vtkPolyDataAlgorithm);
#endif

 protected:

  vtkGranularTemperature();
  ~vtkGranularTemperature();

  virtual int RequestData(
			  vtkInformation* request,
			  vtkInformationVector** InputVector,
			  vtkInformationVector* outputVector);
  virtual int RequestUpdateExtent(
			  vtkInformation* request,
			  vtkInformationVector** InputVector,
			  vtkInformationVector* outputVector);
  virtual int FillInputPortInformation(int port,vtkInformation *info);




};
