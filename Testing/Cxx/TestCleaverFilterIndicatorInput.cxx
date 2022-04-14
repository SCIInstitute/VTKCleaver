
// First include the required header files for the VTK classes we are using.
#include "vtkActor.h"
#include "vtkAppendFilter.h"
#include "vtkAxes.h"
#include "vtkCamera.h"
#include "vtkCellData.h"
#include "vtkConeSource.h"
#include "vtkCubeSource.h"
#include "vtkCullerCollection.h"
#include "vtkDataSetMapper.h"
#include "vtkElevationFilter.h"
#include "vtkFloatArray.h"
#include "vtkGeometryFilter.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkLight.h"
#include "vtkLogger.h"
#include "vtkLookupTable.h"
#include "vtkPNGWriter.h"
#include "vtkPointDataToCellData.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRegressionTestImage.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkTransform.h"
#include "vtkTransformFilter.h"
#include "vtkTubeFilter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkWindowToImageFilter.h"

#include "vtkCleaverImageToUnstructuredGridFilter.h"

const int TEST_IMAGE_SIZE = 24;

// Test the use of indicator functions as input. Indicator functions are
// increasingly positive where the label/region is present, and more negative
// further from the label/region.
int TestCleaverFilterIndicatorInput(int argc, char* argv[])
{
  vtkNew<vtkImageData> indicatorImage[3];
  indicatorImage[0]->SetDimensions(TEST_IMAGE_SIZE, TEST_IMAGE_SIZE, TEST_IMAGE_SIZE);
  indicatorImage[0]->AllocateScalars(VTK_DOUBLE, 1);
  indicatorImage[1]->SetDimensions(TEST_IMAGE_SIZE, TEST_IMAGE_SIZE, TEST_IMAGE_SIZE);
  indicatorImage[1]->AllocateScalars(VTK_DOUBLE, 1);
  indicatorImage[2]->SetDimensions(TEST_IMAGE_SIZE, TEST_IMAGE_SIZE, TEST_IMAGE_SIZE);
  indicatorImage[2]->AllocateScalars(VTK_DOUBLE, 1);

  // indicator functions based on axes
  for (int i = 0; i < TEST_IMAGE_SIZE; ++i)
  {
    for (int j = 0; j < TEST_IMAGE_SIZE; ++j)
    {
      for (int k = 0; k < TEST_IMAGE_SIZE; ++k)
      {
        // cylinder, axis at zero
        // double val = - (i * i + j * j) + TEST_IMAGE_SIZE ;
        // linear ramp crossing zero on one axis.
        double val = i - TEST_IMAGE_SIZE * 0.66;
        // curve, no zeros
        double val2 = (j - k * i + 0.3);
        // constant background
        indicatorImage[0]->SetScalarComponentFromDouble(i, j, k, 0, 0.0);
        indicatorImage[1]->SetScalarComponentFromDouble(i, j, k, 0, val2);
        indicatorImage[2]->SetScalarComponentFromDouble(i, j, k, 0, val);
      }
    }
  }

  vtkNew<vtkCleaverImageToUnstructuredGridFilter> cleaverFilter;

  cleaverFilter->SetInputData(0, indicatorImage[0]);
  cleaverFilter->AddInputData(0, indicatorImage[1]);
  cleaverFilter->AddInputData(0, indicatorImage[2]);
  cleaverFilter->SetSamplingRate(1.0);
  cleaverFilter->SetRateOfChange(0.1);
  // cleaverFilter->SetFeatureScaling(0.5);
  cleaverFilter->Update();

  // Output checking
  vtkLogIf(ERROR, cleaverFilter->GetOutput()->GetNumberOfCells() != 4152,
    "Unexpected number of cells " << cleaverFilter->GetOutput()->GetNumberOfCells());
  vtkLogIf(ERROR, cleaverFilter->GetOutput()->GetNumberOfPoints() != 864,
    "Unexpected number of points " << cleaverFilter->GetOutput()->GetNumberOfPoints());

  vtkNew<vtkDataSetMapper> allMapper;
  allMapper->SetInputConnection(cleaverFilter->GetOutputPort());
  allMapper->ScalarVisibilityOn();
  allMapper->SetScalarModeToUseCellData();
  allMapper->SetColorModeToMapScalars();
  allMapper->SetScalarRange(cleaverFilter->GetOutput()->GetCellData()->GetScalars()->GetRange());

  vtkNew<vtkLookupTable> rainbowBlueRedLut;
  rainbowBlueRedLut->SetNumberOfColors(256);
  rainbowBlueRedLut->SetHueRange(0.667, 0.0);
  rainbowBlueRedLut->Build();
  allMapper->SetLookupTable(rainbowBlueRedLut);

  vtkNew<vtkActor> allActor;
  allActor->SetMapper(allMapper);
  allActor->GetProperty()->SetEdgeVisibility(true);

  vtkNew<vtkRenderer> ren1;
  // turn off all cullers
  ren1->GetCullers()->RemoveAllItems();
  ren1->AddActor(allActor);
  ren1->SetBackground(0.1, 0.2, 0.4);
  ren1->GetActiveCamera()->SetPosition(-30.0, -20.0, 100.0);
  ren1->GetActiveCamera()->SetFocalPoint(10.0, 10.0, 10.0);
  ren1->GetActiveCamera()->SetViewAngle(45.0);

  vtkNew<vtkRenderWindow> renWin;
  renWin->AddRenderer(ren1);
  renWin->SetSize(600, 600);

  renWin->Render();

  vtkNew<vtkRenderWindowInteractor> iren;
  iren->SetRenderWindow(renWin);

  vtkNew<vtkInteractorStyleTrackballCamera> style;
  //  iren->SetInteractorStyle(style);

  iren->Initialize();

  int retVal = vtkRegressionTestImage(renWin);

  if (retVal == vtkRegressionTester::DO_INTERACTOR)
  {
    iren->Start();
  }

  return 0;
}
