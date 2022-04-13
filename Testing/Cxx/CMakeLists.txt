# Test has a baseline, and writes a file.
vtk_add_test_cxx(vtkCleaverImageToUnstructuredGridFilterCxxTests tests
  TestCleaverFilterLabelMapInput.cxx,NO_DATA
  TestCleaverFilterIndicatorInput.cxx,NO_DATA
  )

vtk_test_cxx_executable(vtkCleaverImageToUnstructuredGridFilterCxxTests tests
  RENDERING_FACTORY
  DISABLE_FLOATING_POINT_EXCEPTIONS
  )
