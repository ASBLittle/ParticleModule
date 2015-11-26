#include "vtkPolyDataAlgorithm.h"

class vtkGranularTemperature : public vtkPolyDataAlgorithm
{
public:
  static vtkGranularTemperature *New();
  vtkTypeRevisionMacro(vtkGranularTemperature,vtkPolyDataAlgorithm);

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
