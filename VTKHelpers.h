/*=========================================================================
 *
 *  Copyright David Doria 2011 daviddoria@gmail.com
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#ifndef VTKHelpers_H
#define VTKHelpers_H

// VTK
class vtkImageData;
class vtkPoints;
class vtkPolyData;

// STL
#include <string>

namespace VTKHelpers
{

/** The values that VTK interprets in the 4th channel of an image as opaque and transparent. */
enum PixelValuesEnum {TRANSPARENT_PIXEL = 0, OPAQUE_PIXEL = 255};

/** Set the size of 'output' to that of 'input'. */
void SetImageSizeToMatch(vtkImageData* const input, vtkImageData* const output);

/** Compute how many unique points there are in a vtkPoints object. This is a special case where the points are all along a line in order. */
unsigned int NumberOfUniquePoints(vtkPoints* const points, const float tolerance);

/** Convert an ordered list of 0D 'points' into a path (i.e. add 1D topology). */
void PathFromPoints(vtkPoints* const points, vtkPolyData* const path);

/** Write a vtkPolyData to a .vtp file. */
void WritePolyData(vtkPolyData* const polyData, const std::string& fileName);

/** Write a vtkImageData to a .vti file. */
void WriteImageData(vtkImageData* const imageData, const std::string& fileName);

/** Get the center of a cell. */
void GetCellCenter(vtkImageData* const imageData, const unsigned int cellId, double center[3]);

/** Set the center pixel of an 'image' to the specified 'color'. The image is assumed to have odd dimensions. */
void SetImageCenterPixel(vtkImageData* const image, const unsigned char color[3]);

/** Set an image to black except for its border, which is set to 'color'. */
void BlankAndOutlineImage(vtkImageData* const image, const unsigned char color[3]);

/** Set an image to black except for its border, which is set to 'color'. */
void OutlineImage(vtkImageData* const image, const unsigned char color[3]);

/** Set an image to black. */
void ZeroImage(vtkImageData* const image, const unsigned int channels);

/** Extract the non-zero pixels of a "vector image" and convert them to vectors in a vtkPolyData.
 * This is useful because glyphing a vector image is too slow to use as a visualization,
 * because it "draws" the vectors, even if they are zero length. In this code we are often
 * interested in displaying vectors along a contour, so this is a very very small subset of a whole vector image.
 */
void KeepNonZeroVectors(vtkImageData* const image, vtkPolyData* const output);

/** Make pixels with value 'value' transparent */
void MakeValueTransparent(vtkImageData* const inputImage, const unsigned char value[3]);

/** Make an entire image transparent. */
void MakeImageTransparent(vtkImageData* const image);

/** Output the names of all arrays in the PointData of a PolyData */
void OutputAllArrayNames(vtkPolyData* const polyData);

/** Scale 'image' to the range [0,255] */
void ScaleImage(vtkImageData* const image);

void MaskImage(vtkImageData* const VTKImage, vtkImageData* const VTKSegmentMask, vtkImageData* const VTKMaskedImage);

} // end namespace

#endif
