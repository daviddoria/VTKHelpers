Usage
-----
Everything in this repository is in the VTKHelpers namespace.
After including VTKHelpers.h, the functions can be used as follows:

vtkPolyData* polyData = ...;
VTKHelpers::WritePolyData(polyData, "output.vtp");
