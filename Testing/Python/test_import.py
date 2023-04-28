from vtk import vtkCleaver


def test_basic_import():
    # Make sure it has some of the attributes we are expecting
    assert hasattr(vtkCleaver, 'vtkCleaverImageToUnstructuredGridFilter')
    assert hasattr(vtkCleaver.vtkCleaverImageToUnstructuredGridFilter, 'GetInputIsIndicatorFunction')
