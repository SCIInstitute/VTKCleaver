/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkCleaverImageToUnstructuredGridFilter.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#ifndef vtkCleaverImageToUnstructuredGridFilter_h
#define vtkCleaverImageToUnstructuredGridFilter_h

#include "vtkCleaverModule.h" // For export macro

#include "vtkUnstructuredGridAlgorithm.h"

class vtkInformation;
class vtkInformationVector;
class vtkDataSet;

/// Convert an image containing a labelmap to a mesh, using Cleaver
class VTKCLEAVER_EXPORT vtkCleaverImageToUnstructuredGridFilter
  : public vtkUnstructuredGridAlgorithm
{
public:
  static vtkCleaverImageToUnstructuredGridFilter* New();
  vtkTypeMacro(vtkCleaverImageToUnstructuredGridFilter, vtkUnstructuredGridAlgorithm);

  void PrintSelf(ostream& os, vtkIndent indent) override;

  ///@{
  /** Is the input image a label image or an indicator function? This is only used if
   * there is only one input. Otherwise, indicator functions are assumed.
   */
  vtkSetMacro(InputIsIndicatorFunction, bool);
  vtkGetMacro(InputIsIndicatorFunction, bool);
  vtkBooleanMacro(InputIsIndicatorFunction, bool);
  ///@}

  ///@{
  /** Blending function sigma for input(s) to remove alias artifacts. */
  vtkSetMacro(Sigma, double);
  vtkGetMacro(Sigma, double);
  ///@}

  ///@{
  /** Sizing field sampling rate. The sampling rate of the input indicator functions or calculated
   * indicator functions from segmentation files.
   * The default sample rate will be the dimensions of the volume. Smaller sampling creates coarser
   * meshes.
   * Adjusting this parameter will also affect Cleaver’s runtime, with smaller values running
   * faster. */
  vtkSetMacro(SamplingRate, double);
  vtkGetMacro(SamplingRate, double);
  ///@}

  ///@{
  /** Sizing field rate of change. the maximum rate of change of element size throughout a mesh.
   * Helpful for meshes with high and low curvature.
   * Will have no effect on meshes with constant element sizing methods. */
  vtkSetMacro(RateOfChange, double);
  vtkGetMacro(RateOfChange, double);
  ///@}

  ///@{
  /** Sizing field feature scaling. Scales features of the mesh effecting element size. Higher
   * feature scaling creates coarser meshes. */
  vtkSetMacro(FeatureScaling, double);
  vtkGetMacro(FeatureScaling, double);
  ///@}

  ///@{
  /** Sizing field padding. Adds a volume buffer around the data. Useful when volumes intersect near
   * the boundary. */
  vtkSetMacro(Padding, int);
  vtkGetMacro(Padding, int);
  ///@}

  ///@{
  vtkSetMacro(Alpha, double);
  vtkGetMacro(Alpha, double);
  ///@}

protected:
  vtkCleaverImageToUnstructuredGridFilter();
  ~vtkCleaverImageToUnstructuredGridFilter() override;

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillInputPortInformation(int port, vtkInformation* info) override;

  bool InputIsIndicatorFunction{ false };
  double Alpha{ 0.4 };
  double SamplingRate{ 1.0 };
  double RateOfChange{ 0.2 };
  double FeatureScaling{ 1.0 };
  int Padding{ 0 };
  double Sigma{ 1.0 };

private:
  vtkCleaverImageToUnstructuredGridFilter(const vtkCleaverImageToUnstructuredGridFilter&) = delete;
  void operator=(const vtkCleaverImageToUnstructuredGridFilter&) = delete;
};

#endif
