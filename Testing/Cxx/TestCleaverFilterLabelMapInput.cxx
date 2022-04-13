
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
// Test the use of a labelmap as input.
int TestCleaverFilterLabelMapInput(int argc, char* argv[])
{
  vtkNew<vtkImageData> labelImage;
  labelImage->SetDimensions(TEST_IMAGE_SIZE, TEST_IMAGE_SIZE, TEST_IMAGE_SIZE);
  labelImage->AllocateScalars(VTK_DOUBLE, 1);

  // label map with 2 spheres
  double cX = 20.0;
  double cY = 20.0;
  double cZ = 20.0;
  double rad = 5.0;
  for (int i = 0; i < TEST_IMAGE_SIZE; ++i)
  {
    for (int j = 0; j < TEST_IMAGE_SIZE; ++j)
    {
      for (int k = 0; k < TEST_IMAGE_SIZE; ++k)
      {
        double val =
          (cX - i) * (cX - i) + (cY - j) * (cY - j) + (cZ - k) * (cZ - k) < rad * rad ? 1.0 : 0.0;
        // second sphere at zero
        if (i * i + j * j + k * k < TEST_IMAGE_SIZE * TEST_IMAGE_SIZE * 0.25)
        {
          val = 2.0;
        }
        labelImage->SetScalarComponentFromDouble(i, j, k, 0, val);
      }
    }
  }

  vtkNew<vtkCleaverImageToUnstructuredGridFilter> cleaverFilter;

  cleaverFilter->SetInputData(labelImage);
  cleaverFilter->SetRateOfChange(0.19);
  cleaverFilter->SetSamplingRate(1.0);
  cleaverFilter->Update();

  // Output checking
  vtkLogIf(ERROR, cleaverFilter->GetOutput()->GetNumberOfCells() != 1170,
    "Unexpected number of cells " << cleaverFilter->GetOutput()->GetNumberOfCells());
  vtkLogIf(ERROR, cleaverFilter->GetOutput()->GetNumberOfPoints() != 281,
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
