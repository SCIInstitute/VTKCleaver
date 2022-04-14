/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkCleaverImageToUnstructuredGridFilter.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkCleaverImageToUnstructuredGridFilter.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkImageCast.h"
#include "vtkImageData.h"
#include "vtkImageEuclideanDistance.h"
#include "vtkImageGaussianSmooth.h"
#include "vtkImageIterator.h"
#include "vtkImageMathematics.h"
#include "vtkImageThreshold.h"
#include "vtkInformation.h"
#include "vtkInformationIterator.h"
#include "vtkInformationStringKey.h"
#include "vtkInformationVector.h"
#include "vtkIntArray.h"
#include "vtkLogger.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkProperty.h"
#include "vtkSmartPointer.h"
#include "vtkStdString.h"
#include "vtkTetra.h"
#include "vtkUnstructuredGrid.h"

#include <cleaver/Cleaver.h>
#include <cleaver/CleaverMesher.h>
#include <cleaver/InverseField.h>
#include <cleaver/SizingFieldCreator.h>
#include <cleaver/TetMesh.h>

#include <cassert>
#include <cmath>
#include <list>
#include <sstream>

namespace
{
// NRRDTools is at https://github.com/SCIInstitute/Cleaver2/blob/master/src/lib/nrrd2cleaver/itk/NRRDTools.cpp
// from NRRDTools::checkImageSize
bool checkImageSize(vtkImageData* inputImg, double sigma)
{
  int dims[3];
  inputImg->GetDimensions(dims);
  auto* spacing = inputImg->GetSpacing();
  std::vector<double> imageSize{ dims[0] * spacing[0], dims[1] * spacing[1], dims[2] * spacing[2] };
  double imageSizeMin = *(std::min_element(std::begin(imageSize), std::end(imageSize)));

  return (sigma / imageSizeMin) >= 0.1;
}

// Follow the pattern of NRRDTools::segmentationToIndicatorFunctions but
// substitute VTK filters.
std::vector<cleaver::AbstractScalarField*> segmentationToIndicatorFunctions(
  vtkImageData* image, double sigma)
{
  std::vector<cleaver::AbstractScalarField*> fields;
  double range[2];
  image->GetScalarRange(range);
  auto minLabel = static_cast<size_t>(range[0]);
  auto maxLabel = static_cast<size_t>(range[1]);
  vtkLog(INFO, "Range " << minLabel << " " << maxLabel);
  // setup two pipelines which are combined with subtraction below.
  vtkNew<vtkImageThreshold> threshold[2];
  vtkNew<vtkImageGaussianSmooth> blur[2];
  vtkNew<vtkImageEuclideanDistance> distance[2];
  vtkNew<vtkImageThreshold> thresholdDistance[2];
  vtkNew<vtkImageMathematics> imgdiff;
  for (auto idx = 0; idx < 2; ++idx)
  {
    // setup to create a binary mask from each label
    threshold[idx]->ReplaceInOn();
    threshold[idx]->ReplaceOutOn();
    threshold[idx]->SetInputData(image);
    // make sure blur input is floating point.
    threshold[idx]->SetOutputScalarTypeToFloat();
    // bluring
    blur[idx]->SetInputConnection(threshold[idx]->GetOutputPort());
    blur[idx]->SetStandardDeviation(sigma);
    // second threshold
    thresholdDistance[idx]->ReplaceInOn();
    thresholdDistance[idx]->SetInValue(0.0);
    thresholdDistance[idx]->SetOutputScalarTypeToFloat();
    thresholdDistance[idx]->SetInputConnection(blur[idx]->GetOutputPort());
    // distance map, output always doubles
    distance[idx]->SetInputConnection(thresholdDistance[idx]->GetOutputPort());
  }
  // Create an indicator function which is zero at the edges of a region,
  // has increasing positive values based on distance inside the region,
  // and negative values for distance outside.
  // Mimic itkApproximateSignedDistanceMapImageFilter, which produces the
  // inverse of an indicator function, by combining two
  // threshold + distance that are the inverse of each other
  threshold[0]->SetInValue(0.0);
  threshold[0]->SetOutValue(1.0);
  threshold[1]->SetInValue(1.0);
  threshold[1]->SetOutValue(0.0);
  // extract images from each label for an indicator function
  for (size_t ii = minLabel, num = 0; ii <= maxLabel; ii++, num++)
  {
    for (auto idx = 0; idx < 2; ++idx)
    {
      threshold[idx]->ThresholdBetween(
        static_cast<double>(ii) - 0.001, static_cast<double>(ii) + 0.001);
      blur[idx]->Update();
      auto* img = blur[idx]->GetOutput();

      img->GetScalarRange(range);
      // get the midpoint of the range
      auto md = (range[0] + range[1]) / 2.f;
      vtkLog(INFO, "Md " << md);

      // ITK itkApproximateSignedDistanceMapImageFilter has inside/outside values.
      // For VTK, anything non-zero is considered NOT an object. Threshold again to
      // set 0.0 where the object is present. threshold[1] treats the background
      // as the object, effectively computing distances inside the object.
      thresholdDistance[idx]->ThresholdByLower(md);

      distance[idx]->Update();
    }
    // combine the two distances, so distance inside is negative, distance outside
    // is positive, consistent with the ITK filter.
    imgdiff->SetOperationToSubtract();
    imgdiff->SetInputConnection(0, distance[0]->GetOutputPort());
    imgdiff->SetInputConnection(1, distance[1]->GetOutputPort());
    imgdiff->Update();
    // this is now the inverse of an indicator function.
    auto* img = imgdiff->GetOutput();

    int dims[3];
    img->GetDimensions(dims);
    size_t numPixel = dims[0] * dims[1] * dims[2];
    float* data = new float[numPixel];
    fields.push_back(new cleaver::FloatField(data, dims[0], dims[1], dims[2]));
    std::string name("SegmentationLabel");
    std::stringstream ss;
    ss << name << ii;
    fields[num]->setName(ss.str());
    bool warning = checkImageSize(img, sigma);
    fields[num]->setWarning(warning);

    // One field for each label in the input labelmap image.
    vtkImageIterator<double> it(img, img->GetExtent());
    size_t pixel = 0;
    float min = static_cast<float>(*(it.BeginSpan()));
    float max = min;
    auto* spacing = img->GetSpacing();
    std::string error = "none";
    while (!it.IsAtEnd())
    {
      double* valIt = it.BeginSpan();
      double* valEnd = it.EndSpan();
      while (valIt != valEnd)
      {
        float val = static_cast<float>(*valIt++);
        // ((cleaver::FloatField*)fields[num])->data()[pixel++] = -val;
        data[pixel++] = -val;
        // Error checking
        if (std::isnan(val) && error.compare("none") == 0)
        {
          error = "nan";
        }
        else if (val < min)
        {
          min = val;
        }
        else if (val > max)
        {
          max = val;
        }
      }
      it.NextSpan();
    }
    if ((min >= 0 || max <= 0) && (error.compare("none") == 0))
    {
      error = "maxmin";
    }

    std::cout << "Field " << fields[num]->name() << " " << error << std::endl;
    fields[num]->setError(error);
    ((cleaver::FloatField*)fields[num])
      ->setScale(cleaver::vec3(spacing[0], spacing[1], spacing[2]));
  }
  return fields;
}

// Follow the pattern of NRRDTools::loadNRRDFiles
std::vector<cleaver::AbstractScalarField*> loadIndicatorFunctions(
  const std::vector<vtkImageData*> &images, double sigma)
{
  std::vector<cleaver::AbstractScalarField*> fields;
  vtkNew<vtkImageCast> makeFloat;
  makeFloat->SetOutputScalarTypeToFloat();
  vtkNew<vtkImageGaussianSmooth> blur;
  blur->SetStandardDeviation(sigma);
  blur->SetDimensionality(3);
  for(auto image: images)
  {
    makeFloat->SetInputData(image);
    // blur, input should be floats.
    blur->SetInputConnection(makeFloat->GetOutputPort());
    blur->Update();
    auto* img = blur->GetOutput();

    int dims[3];
    img->GetDimensions(dims);
    size_t numPixel = dims[0] * dims[1] * dims[2];
    float* data = new float[numPixel];
    size_t num = fields.size();
    fields.push_back(new cleaver::FloatField(data, dims[0], dims[1], dims[2]));
    std::string name = image->GetObjectName();
    if (name.empty())
    {
      std::stringstream ss;
      ss << "IndicatorFunction" << num;
      fields[num]->setName(ss.str());
    }
    else
    {
      fields[num]->setName(name);
    }
    bool warning = checkImageSize(img, sigma);
    fields[num]->setWarning(warning);

    // One field per input image.
    vtkImageIterator<float> it(img, img->GetExtent());
    size_t pixel = 0;
    float min = static_cast<float>(*(it.BeginSpan()));
    float max = min;
    auto* spacing = img->GetSpacing();
    std::string error = "none";
    while (!it.IsAtEnd())
    {
      float* valIt = it.BeginSpan();
      float* valEnd = it.EndSpan();
      while (valIt != valEnd)
      {
        float val = static_cast<float>(*valIt++);
        // ((cleaver::FloatField*)fields[num])->data()[pixel++] = val;
        data[pixel++] = val;
        // Error checking
        if (std::isnan(val) && error.compare("none") == 0)
        {
          error = "nan";
        }
        else if (val < min)
        {
          min = val;
        }
        else if (val > max)
        {
          max = val;
        }
      }
      it.NextSpan();
    }
    if ((min >= 0 || max <= 0) && (error.compare("none") == 0))
    {
      error = "maxmin";
    }

    std::cout << "Field " << fields[num]->name() << " " << error << std::endl;
    fields[num]->setError(error);
    ((cleaver::FloatField*)fields[num])
      ->setScale(cleaver::vec3(spacing[0], spacing[1], spacing[2]));
  }
  return fields;
}

// Support for fillUnstructuredGrid(), copied from Cleaver TetMesh.cpp
using vec3 = cleaver::vec3;
using TetMesh = cleaver::TetMesh;
using Vertex = cleaver::Vertex;
using Tet = cleaver::Tet;
class vec3_compare
{
public:
  bool operator()(const vec3& a, const vec3& b) const
  {
    if ((a.x < b.x) && (b.x - a.x) > 1e-9)
      return true;
    else if ((a.x > b.x) && (a.x - b.x) > 1e-9)
      return false;
    if ((a.y < b.y) && (b.y - a.y) > 1e-9)
      return true;
    else if ((a.y > b.y) && (a.y - b.y) > 1e-9)
      return false;
    if ((a.z < b.z) && (b.z - a.z) > 1e-9)
      return true;
    else if ((a.z > b.z) && (a.z - b.z) > 1e-9)
      return false;
    return false;
  }
};
typedef std::map<const vec3, size_t, vec3_compare> VertMap;

// given a cleaver TetMesh, fill in the UnstructuredGrid output. Based on
// TetMesh::writeVtkUnstructuredGrid()
// https://github.com/SCIInstitute/Cleaver2/blob/master/src/lib/cleaver/TetMesh.cpp#L1358
void fillUnstructuredGrid(vtkUnstructuredGrid* ugrid, TetMesh* tetMesh)
{
  //-----------------------------------
  //         Create Pruned Vertex List
  //-----------------------------------
  VertMap vert_map;
  std::vector<vec3> pruned_verts;
  size_t pruned_pos = 0;
  auto numTets = tetMesh->tets.size();
  for (size_t t = 0; t < numTets; t++)
  {
    Tet* tet = tetMesh->tets[t];

    Vertex* v1 = tet->verts[0];
    Vertex* v2 = tet->verts[1];
    Vertex* v3 = tet->verts[2];
    Vertex* v4 = tet->verts[3];

    vec3 p1 = v1->pos();
    vec3 p2 = v2->pos();
    vec3 p3 = v3->pos();
    vec3 p4 = v4->pos();

    if (!vert_map.count(p1))
    {
      vert_map.insert(std::pair<vec3, size_t>(p1, pruned_pos));
      pruned_pos++;
      pruned_verts.push_back(p1);
    }
    if (!vert_map.count(p2))
    {
      vert_map.insert(std::pair<vec3, size_t>(p2, pruned_pos));
      pruned_pos++;
      pruned_verts.push_back(p2);
    }
    if (!vert_map.count(p3))
    {
      vert_map.insert(std::pair<vec3, size_t>(p3, pruned_pos));
      pruned_pos++;
      pruned_verts.push_back(p3);
    }
    if (!vert_map.count(p4))
    {
      vert_map.insert(std::pair<vec3, size_t>(p4, pruned_pos));
      pruned_pos++;
      pruned_verts.push_back(p4);
    }
  }

  //-----------------------------------
  //         Write Vertex List
  //-----------------------------------
  vtkPoints* points = ugrid->GetPoints();
  if (points)
  {
    points->Reset();
  }
  else
  {
    vtkNew<vtkPoints> pts;
    vtkNew<vtkDoubleArray> coords;
    coords->SetNumberOfComponents(3);
    pts->SetData(coords);
    ugrid->SetPoints(pts);
    points = pts;
  }
  for (size_t i = 0; i < pruned_verts.size(); i++)
  {
    points->InsertNextPoint(pruned_verts[i].x, pruned_verts[i].y, pruned_verts[i].z);
  }

  //-----------------------------------
  //         Write Cell/Face List
  //-----------------------------------
  // \todo make writing background optional
  ugrid->AllocateExact(numTets, numTets * 4);
  vtkNew<vtkTetra> tet;
  for (size_t f = 0; f < numTets; f++)
  {
    Tet* t = tetMesh->tets.at(f);

    for (unsigned int r = 0; r < 4; ++r)
    {
      Vertex* v1 = t->verts[r];
      size_t idx = vert_map.find(v1->pos())->second;
      tet->GetPointIds()->SetId(r, idx);
    }
    ugrid->InsertNextCell(VTK_TETRA, tet->GetPointIds());
  }

  vtkIntArray* labels = vtkIntArray::SafeDownCast(ugrid->GetCellData()->GetArray("labels"));
  ;
  if (!labels)
  {
    vtkNew<vtkIntArray> lbl;
    labels = lbl;
    lbl->SetName("labels");
    ugrid->GetCellData()->AddArray(lbl);
  }
  labels->SetNumberOfTuples(numTets);
  //-----------------------------------
  //         Write Labels
  //-----------------------------------
  for (size_t f = 0; f < numTets; ++f)
  {
    Tet* t = tetMesh->tets.at(f);
    labels->SetValue(f, static_cast<int>(t->mat_label));
  }
  ugrid->GetCellData()->SetActiveScalars("labels");
}
}

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkCleaverImageToUnstructuredGridFilter);

//----------------------------------------------------------------------------
vtkCleaverImageToUnstructuredGridFilter::vtkCleaverImageToUnstructuredGridFilter() = default;
//----------------------------------------------------------------------------
vtkCleaverImageToUnstructuredGridFilter::~vtkCleaverImageToUnstructuredGridFilter() = default;

//----------------------------------------------------------------------------
int vtkCleaverImageToUnstructuredGridFilter::FillInputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  // allow more than one input for indicator functions instead of a single
  // label map.
  info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
  return 1;
}

//----------------------------------------------------------------------------
int vtkCleaverImageToUnstructuredGridFilter::RequestData(vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  const auto element_sizing_method = cleaver::Adaptive;
  const bool strip_exterior = true;
  const bool fix_tets = true;
  const bool verbose = true;
  const bool segmentation = this->GetNumberOfInputConnections(0) == 1 && !this->GetInputIsIndicatorFunction();
  // get the output info object
  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  // get the output
  vtkUnstructuredGrid* output =
    vtkUnstructuredGrid::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);

  // get the first input
  vtkImageData* inputImg = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

  // calculate indicator functions from label map.
  std::vector<cleaver::AbstractScalarField*> fields;
  if (segmentation)
  {
    fields = segmentationToIndicatorFunctions(inputImg, this->Sigma);
  }
  else
  {
    // inputs are indicator functions directly.
    std::vector<vtkImageData *> images;
    for (auto i = 0; i < this->GetNumberOfInputConnections(0); ++i)
    {
      inInfo = inputVector[0]->GetInformationObject(i);
      inputImg = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
      if (inputImg)
      {
        images.push_back(inputImg);
      }
    }
    fields = loadIndicatorFunctions(images, this->Sigma);
  }
  bool simple = false;
  cleaver::TetMesh* bgMesh = nullptr;
  cleaver::Volume* volume = new cleaver::Volume(fields);
  cleaver::CleaverMesher mesher(simple);
  mesher.setVolume(volume);
  mesher.setAlphaInit(this->Alpha);

  // Create the sizing field
  std::vector<cleaver::AbstractScalarField*> sizingField;
  sizingField.push_back(cleaver::SizingFieldCreator::createSizingFieldFromVolume(volume,
    (float)(1.0 / this->RateOfChange), (float)this->SamplingRate, (float)this->FeatureScaling,
    (int)this->Padding, (element_sizing_method != cleaver::Constant), verbose));

  // Set Sizing Field on Volume
  volume->setSizingField(sizingField[0]);

  // Construct Background Mesh
  mesher.setConstant(false);
  bgMesh = mesher.createBackgroundMesh(verbose);

  // Apply Mesh Cleaving
  mesher.buildAdjacency(verbose);
  mesher.sampleVolume(verbose);
  mesher.computeAlphas(verbose);
  mesher.computeInterfaces(verbose);
  mesher.generalizeTets(verbose);
  mesher.snapsAndWarp(verbose);
  mesher.stencilTets(verbose);

  cleaver::TetMesh* mesh = mesher.getTetMesh();

  // Strip Exterior Tets
  if (strip_exterior)
  {
    cleaver::stripExteriorTets(mesh, volume, verbose);
  }

  // Compute Quality If Havn't Already
  mesh->computeAngles();

  // Fix jacobians if requested.
  if (fix_tets)
  {
    mesh->fixVertexWindup(verbose);
  }

  // cleaver::MeshFormat output_format = cleaver::VtkUSG;
  // std::string output_path = "./";
  // std::string output_name = "mesh";
  // mesh->writeMesh(output_path + output_name, output_format, verbose);
  // mesh->writeInfo(output_path + output_name, verbose);
  fillUnstructuredGrid(output, mesh);

  return 1;
}

void vtkCleaverImageToUnstructuredGridFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  Superclass::PrintSelf(os, indent);
  os << "InputIsIndicatorFunction: " << this->InputIsIndicatorFunction << std::endl;
  os << "Alpha: " << this->Alpha << std::endl;
  os << "SamplingRate: " << this->SamplingRate << std::endl;
  os << "RateOfChange: " << this->RateOfChange << std::endl;
  os << "FeatureScaling: " << this->FeatureScaling << std::endl;
  os << "Padding: " << this->Padding << std::endl;
  os << "Sigma: " << this->Sigma << std::endl;
}
