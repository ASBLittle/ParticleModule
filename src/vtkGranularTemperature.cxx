#include "vtkGranularTemperature.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkCell.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkDataObject.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkPointLocator.h"
#include "vtkMath.h"
#include "vtkDoubleArray.h"
#include <iostream>
#include <map>
#include <set>
#include <math.h> 


vtkCxxRevisionMacro(vtkGranularTemperature, "$Revision: 0.0$");
vtkStandardNewMacro(vtkGranularTemperature);

double radial_distribution(double volfrac){
  
  double ALPHA_MAX = 0.59999;
  double ALPHA0    = 0.6;

  if ( volfrac> ALPHA_MAX ) {volfrac = ALPHA_MAX; }
  
  return 1.0/(1.0-pow(volfrac/ALPHA0,1.0/3.0));
    }

vtkGranularTemperature::vtkGranularTemperature(){
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
};
vtkGranularTemperature::~vtkGranularTemperature(){};

int vtkGranularTemperature::RequestData(
		      vtkInformation* vtkNotUsed(request),
		      vtkInformationVector **inputVector,
		      vtkInformationVector* outputVector )
{

#ifndef NDEBUG
  this->DebugOn();
#endif
  // output0 will be a copy of the high dimensional cells and 
  // output1 will be the rest
  vtkInformation* outInfo=outputVector->GetInformationObject(0);
  vtkPolyData* output= vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT() ) );

  vtkInformation *inInfo=inputVector[0]->GetInformationObject(0);
  vtkPolyData* input= vtkPolyData::GetData(inputVector[0]);

  vtkDebugMacro(<<"In Granular temperature code");

  // first copy the input

  output->CopyStructure(input);
  output->SetPoints(input->GetPoints());
  output->DeepCopy(input);

  vtkDebugMacro(<<"Data copied");

  vtkDebugMacro(<< output->GetNumberOfPoints());

  vtkPointLocator* locator = vtkPointLocator::New();
  locator->SetDataSet(input);
  locator->BuildLocator();

  vtkDebugMacro(<<"Locator built");

  double LENGTH = 0.03;
  double MODIFIER = 2.e3;

  vtkDoubleArray* volume = vtkDoubleArray::New();
  vtkDoubleArray* velocity = vtkDoubleArray::New();
  vtkDoubleArray* temperature = vtkDoubleArray::New();
  vtkDoubleArray* pressure = vtkDoubleArray::New();
  vtkDoubleArray* grad_pres = vtkDoubleArray::New();

  volume->SetName("SolidVolumeFraction");
  volume->Allocate(output->GetNumberOfPoints());
  volume->SetNumberOfComponents(1);
  volume->SetNumberOfTuples(output->GetNumberOfPoints());

  velocity->SetName("SolidVolumeVelocity");
  velocity->Allocate(3*output->GetNumberOfPoints());
  velocity->SetNumberOfComponents(3);
  velocity->SetNumberOfTuples(output->GetNumberOfPoints());

  temperature->SetName("GranularTemperature");
  temperature->Allocate(output->GetNumberOfPoints());
  temperature->SetNumberOfComponents(1);
  temperature->SetNumberOfTuples(output->GetNumberOfPoints());

  pressure->SetName("SolidPressure");
  pressure->Allocate(output->GetNumberOfPoints());
  pressure->SetNumberOfComponents(1);
  pressure->SetNumberOfTuples(output->GetNumberOfPoints());

  grad_pres->SetName("SolidPressureGradient");
  grad_pres->Allocate(3*output->GetNumberOfPoints());
  grad_pres->SetNumberOfComponents(3);
  grad_pres->SetNumberOfTuples(output->GetNumberOfPoints());

  vtkDebugMacro(<< velocity->GetNumberOfTuples());

  vtkDebugMacro(<<"Entering main loop");

  for (vtkIdType i=0;i<output->GetNumberOfCells();i++) {
    vtkCell* cell=output->GetCell(i);
    vtkIdList* point_list = vtkIdList::New();

    double diameter=1e-3;
    double density=2.5e3;

    locator->FindPointsWithinRadius(LENGTH,
				    cell->GetPoints()->GetPoint(0),
				    point_list);

    volume->SetValue(i,0.0);
    velocity->SetTuple3(i,0.0,0.0,0.0);
    grad_pres->SetTuple3(i,0.0,0.0,0.0);
    temperature->SetValue(i,0.0);
    pressure->SetValue(i,0.0);
    


    for (vtkIdType j=0;j<point_list->GetNumberOfIds();j++) {

      double rad2=vtkMath::Distance2BetweenPoints(output->GetPoints()->GetPoint(point_list->GetId(j)),
						  output->GetPoints()->GetPoint(i));

      rad2 = rad2 / LENGTH;

      double gamma = 1.0/6.0*vtkMath::Pi()*pow(diameter,3)*exp(-rad2)*MODIFIER / (0.5 * LENGTH * LENGTH*(1.0-exp(-1.0)));   

      // calculate the solid volume fraction:

      volume->SetValue(i, volume->GetValue(i) + gamma );  

      double lvel[3]={0,0,0};

      vtkMath::Add(output->GetPointData()->GetVectors("Particle Velocity")->GetTuple(point_list->GetId(j)),lvel,lvel);
      vtkMath::MultiplyScalar(lvel,gamma);
      vtkMath::Add(velocity->GetTuple3(i),lvel,lvel);


    // calculate solid velocity:
	velocity->SetTuple3(i, lvel[0], lvel[1],lvel[2] ); 
	
  }
    double vel[3]={0,0,0};
    vtkMath::Add(velocity->GetTuple3(i),vel,vel);
    vtkMath::MultiplyScalar(vel,1.0/volume->GetValue(i));
    velocity->SetTuple3(i, vel[0], vel[1],vel[2] );


      
    for (vtkIdType j=0;j<point_list->GetNumberOfIds();j++) {
      	
      double rad2=vtkMath::Distance2BetweenPoints(output->GetPoints()->GetPoint(point_list->GetId(j)),
						  output->GetPoints()->GetPoint(i));

      rad2 = rad2 / LENGTH;

      double gamma = 1.0/6.0*vtkMath::Pi()*pow(diameter,3)*exp(-rad2)*MODIFIER / (0.5 * LENGTH * LENGTH*(1.0-exp(-1.0)));   
      // Now form the granular temperature and its gradient:
      temperature->SetValue(i,temperature->GetValue(i)+
			      vtkMath::Distance2BetweenPoints(velocity->GetTuple3(i), output->GetPointData()->GetVectors("Particle Velocity")->GetTuple(point_list->GetId(j)))*gamma);							     	
    }

    for (vtkIdType j=0;j<point_list->GetNumberOfIds();j++) {
      	
      double rad2=vtkMath::Distance2BetweenPoints(output->GetPoints()->GetPoint(point_list->GetId(j)),
						  output->GetPoints()->GetPoint(i));

      rad2 = rad2 / LENGTH;

      double gamma = 1.0/6.0*vtkMath::Pi()*pow(diameter,3)*exp(-rad2)*MODIFIER / (0.5 * LENGTH * LENGTH*(1.0-exp(-1.0))); 

      double lvel[3]={0,0,0};

      vtkMath::Add(output->GetPoints()->GetPoint(point_list->GetId(j)),
		        lvel, lvel);
      vtkMath::Subtract(output->GetPoints()->GetPoint(i),
		       lvel, lvel);
      vtkMath::MultiplyScalar(lvel,density*radial_distribution(volume->GetValue(i))*gamma);
      vtkMath::Add(grad_pres->GetTuple3(i),lvel,lvel);
      // Now form the granular temperature and its gradient:
      grad_pres->SetTuple3(i,lvel[0],lvel[1],lvel[2]);
    }

    temperature->SetValue(i, temperature->GetValue(i)/volume->GetValue(i));

    pressure->SetValue(i,density*volume->GetValue(i)*radial_distribution(volume->GetValue(i))*temperature->GetValue(i));

   
    point_list->Delete();
  }

  vtkDebugMacro(<<"Main Loop done");   
  
  output->GetPointData()->AddArray(volume);
  output->GetPointData()->AddArray(velocity);
  output->GetPointData()->AddArray(temperature);
  output->GetPointData()->AddArray(pressure);
  output->GetPointData()->AddArray(grad_pres);
  locator->Delete();

  return 1;
  }

int vtkGranularTemperature::RequestUpdateExtent(
			  vtkInformation* request,
			  vtkInformationVector** inputVector,
			  vtkInformationVector* outputVector)
 {

  vtkInformation* outInfo=outputVector->GetInformationObject(0);
  vtkInformation *inInfo=inputVector[0]->GetInformationObject(0);

  int piece, numPieces, ghostLevels;
  
  piece=outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  numPieces=outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
  ghostLevels=outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS());

  if (numPieces > 1)
  {
    ++ghostLevels;
  }

  vtkDebugMacro(<<"Running Update Extent"<<piece<<numPieces<<ghostLevels);

  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(),piece);
  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(),numPieces);
  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(),ghostLevels);
  inInfo->Set(vtkStreamingDemandDrivenPipeline::EXACT_EXTENT(),1);


  return 1;

};


int vtkGranularTemperature::FillInputPortInformation(int,vtkInformation *info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(),"vtkPolyData");
  return 1;
}
