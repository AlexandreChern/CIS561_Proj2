/*=========================================================================

  Program:   Visualization Toolkit
  Module:    SpecularSpheres.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//
// This examples demonstrates the effect of specular lighting.
//
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkJPEGReader.h"
#include "vtkImageData.h"

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>

#include <algorithm>

// ****************************************************************************
//  Function: GetNumberOfPoints
//
//  Arguments:
//     dims: an array of size 3 with the number of points in X, Y, and Z.
//           2D data sets would have Z=1
//
//  Returns:  the number of points in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfPoints(const int *dims)
{
    // 3D
    //return dims[0]*dims[1]*dims[2];
    // 2D
    return dims[0]*dims[1];
}

// ****************************************************************************
//  Function: GetNumberOfCells
//
//  Arguments:
//
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the number of cells in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfCells(const int *dims)
{
    // 3D
    //return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    return (dims[0]-1)*(dims[1]-1);
}


// ****************************************************************************
//  Function: GetPointIndex
//
//  Arguments:
//      idx:  the logical index of a point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1]
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the point index
//
// ****************************************************************************

int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    return idx[1]*dims[0]+idx[0];
}


// ****************************************************************************
//  Function: GetCellIndex
//
//  Arguments:
//      idx:  the logical index of a cell.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the cell index
//
// ****************************************************************************

int GetCellIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)*idx[0];
    // 2D
    return idx[1]*(dims[0]-1)+idx[0];
}

// ****************************************************************************
//  Function: GetLogicalPointIndex
//
//  Arguments:
//      idx (output):  the logical index of the point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1] 
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      pointId:  a number between 0 and (GetNumberOfPoints(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
    // 3D
    // idx[0] = pointId%dim[0];
    // idx[1] = (pointId/dims[0])%dims[1];
    // idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    idx[0] = pointId%dims[0];
    idx[1] = pointId/dims[0];
}


// ****************************************************************************
//  Function: GetLogicalCellIndex
//
//  Arguments:
//      idx (output):  the logical index of the cell index.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      cellId:  a number between 0 and (GetNumberOfCells(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
    // 3D
    // idx[0] = cellId%(dims[0]-1);
    // idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    // idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    idx[0] = cellId%(dims[0]-1);
    idx[1] = cellId/(dims[0]-1);
}

// ****************************************************************************
//  Function: LocateCell (self-defined)
//
//  Arguments:
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     cellID: an array storing the id for the cell which contains the point where the field is to be evaluated
//     pt: a two-dimensional location
//     dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument cellID is output)
//
// ****************************************************************************

void LocateCell(const float *X, const float *Y, int *cellID, const float *pt, const int *dims ){
    // 2D
    // Arguments: pt: 2d array of point locations;
    cellID[0] = 0;
    cellID[1] = 0;
    for (int i=0; i<=dims[0]-1; i++){
        if (X[i] <= pt[0] && X[i+1] > pt[0]){
            cellID[0] = i;
        }
    }

    for (int j=0; j<=dims[1]-1; j++){
        if (Y[j] <= pt[1] && Y[j+1] > pt[1]){
            cellID[1] = j;
        }
    }
}


// ****************************************************************************
//  Function: EvaluateFieldAtLocation
//
//  Arguments:
//     pt: a two-dimensional location
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     F: a scalar field defined on the mesh.  Its size is dims[0]*dims[1].
//
//   Returns: the interpolated field value. 0 if the location is out of bounds.
//
// ****************************************************************************

float
EvaluateFieldAtLocation(const float *pt, const int *dims, 
                        const float *X, const float *Y, const float *F)
{
    // if (pt[0] < X[0] || pt[0] > X[dims[0]-1] || pt[1] < Y[0] || Y[dims[1]-1])
    if (pt[0] < X[0] || pt[1] < Y[0] || pt[0] > X[dims[0]-1] || pt[1] > Y[dims[1]-1])
        return 0; // IMPLEMENT ME!!
    else {
        int cellID[2];
        cellID[0] = 0;
        cellID[1] = 0;
        LocateCell(X,Y,cellID,pt,dims);
        std::cout << cellID[0] << '\t' << cellID[1] << std::endl;
        int pt1[2];
        int pt2[2];
        int pt3[2];
        int pt4[2];

        int pt1_real = GetPointIndex(cellID,dims);
        int pt2_real = pt1_real + 1;
        int pt3_real = pt1_real + dims[0];
        int pt4_real = pt1_real + dims[0] + 1;


        pt1[0] = cellID[0];
        pt1[1] = cellID[1];
        pt2[0] = cellID[0] +1;
        pt2[1] = cellID[1];
        pt3[0] = cellID[0];
        pt3[1] = cellID[1]+1;
        pt4[0] = cellID[0]+1;
        pt4[1] = cellID[1]+1;

        


        float prop_x = float(pt[0] - X[pt1[0]])/float(X[pt1[0]+1] - X[pt1[0]]);
        float prop_y = float(pt[1] - Y[pt1[1]])/float(Y[pt1[1]+1] - Y[pt1[1]]);
        
        float field_1 = F[pt1_real] + prop_x*(F[pt2_real] - F[pt1_real]);
        float field_2 = F[pt3_real] + prop_x*(F[pt4_real] - F[pt3_real]);

        float final_field_1 = field_1 + prop_y*(field_2 - field_1);

        // float field_3 = F[pt1_real] + prop_y*(F[pt3_real] - F[pt1_real]);
        // float field_4 = F[pt2_real] + prop_y*(F[pt4_real] - F[pt2_real]);

        // float final_field_2 = field_3 + prop_x*(field_4-field_3);
        // // return prop_y;
        // float difference = final_field_1 - final_field_2;
        return final_field_1;
    }
    // cellID[1] = 0;
    // LocateCell(X,Y, cellid, pt, dims);
}



// ****************************************************************************
//  Function: BoundingBoxForCell
//
//  Arguments:
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     cellId: a cellIndex (I.e., between 0 and GetNumberOfCells(dims))
//     bbox (output): the bounding box of cellId.  Format should be
//                     bbox[0]: the minimum X value in cellId.
//                     bbox[1]: the maximum X value in cellId.
//                     bbox[2]: the minimum Y value in cellId.
//                     bbox[3]: the maximum Y value in cellId.
//
//  Returns:  None (argument bbox is output)
//
// ****************************************************************************

void
BoundingBoxForCell(const float *X, const float *Y, const int *dims,
                   int cellId, float *bbox)
{
    //bbox[0] = -100;    
    //bbox[1] = +100;
    //bbox[2] = -100;
    //bbox[3] = +100;
    // IMPLEMENT ME!
    int logical_index[3];
    GetLogicalCellIndex(logical_index, cellId, dims);
    bbox[0] = X[logical_index[0]];
    bbox[1] = X[logical_index[0] + 1];
    bbox[2] = Y[logical_index[1]];
    bbox[3] = Y[logical_index[1] + 1];
}

// ****************************************************************************
//  Function: CountNumberOfStraddingCells
//
//  Arguments:
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     F: a scalar field defined on the mesh.  Its size is dims[0]*dims[1].
//
//  Returns:  the number of cells that straddle 0, i.e., the number of cells
//            that contains points who have F>0 and also have points with F<0.
//
// ****************************************************************************

int
CountNumberOfStraddlingCells(const float *X, const float *Y, const int *dims,
                             const float *F)
{   
    int NumberOfStraddlingCells=0;
    int NumberOfPositiveCells = 0;
    int NumberOfNegativeCells = 0;
    int NumberOfTotalCells = (dims[0]-1)*(dims[1]-1);
    int logical_index[3];
    for (int cellId = 0; cellId< NumberOfTotalCells; cellId++){
        GetLogicalCellIndex(logical_index, cellId, dims);
        int pt1 = GetPointIndex(logical_index,dims); // pt1 being the left-bottom point
        int pt2 = pt1+1;
        int pt3 = pt1+dims[0];
        int pt4 = pt1+dims[0]+1;
        if (F[pt1] >= 0 && F[pt2] >= 0 && F[pt3] >= 0 && F[pt4] >= 0){
            NumberOfPositiveCells += 1;
        }
        else if (F[pt1] <= 0 && F[pt2] <= 0 && F[pt3] <= 0 && F[pt4] <= 0){
            NumberOfNegativeCells += 1;
        }
        else{
            NumberOfStraddlingCells += 1;
        }
    }
    return NumberOfStraddlingCells;
    //return NumberOfTotalCells;
    // IMPLEMENT ME!
}

int main()
{
    int  i;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj2_data.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    
    int numCells = CountNumberOfStraddlingCells(X, Y, dims, F);
    cerr << "The number of cells straddling zero is " << numCells << endl;

    float bbox[4];
    const int ncells = 5;
    int cellIds[ncells] = { 0, 50, 678, 1000, 1200 };
    for (i = 0 ; i < ncells ; i++)
    {
        BoundingBoxForCell(X, Y, dims, cellIds[i], bbox);
        cerr << "The bounding box for cell " << cellIds[i] << " is " 
             << bbox[0] << "->" << bbox[1] << ", " << bbox[2] << "->" << bbox[3]
             << endl;
    }

    const int npts = 10;
    // const int npts = 11;
    float pt[npts][3] = 
         {
            // {0.01, 0.02, 0}, // test point remove later
            {1.01119, 0.122062, 0},
            {0.862376, 1.33839, 0},
            {0.155026, 0.126123, 0},
            {0.69736, 0.0653565, 0},
            {0.2, 0.274117, 0},
            {0.893699, 1.04111, 0},
            {0.608791, -0.0533753, 0},
            {1.00543, 0.138024, 0},
            {0.384128, -0.0768977, 0},
            {0.666757, 0.60259, 0},
         };

    

    for (i = 0 ; i < npts ; i++)
    {
        float f = EvaluateFieldAtLocation(pt[i], dims, X, Y, F);
        cerr << "Evaluated field at (" << pt[i][0] <<"," << pt[i][1] << ") as "
             << f << endl;
    }
    
   
    cerr << "Infinite loop here, else Windows people may have the terminal "
         << "disappear before they see the output."
         << " Remove these lines if they annoy you." << endl;
    cerr << "(press Ctrl-C to exit program)" << endl;
    while (1) ; 
}




