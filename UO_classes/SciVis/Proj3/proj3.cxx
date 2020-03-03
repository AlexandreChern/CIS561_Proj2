#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDataSetReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>


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
    //return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
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


void
WriteImage(vtkImageData *img, const char *filename)
{
    std::string full_filename = filename;
    full_filename += ".png";
    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(full_filename.c_str());
    writer->Write();
    writer->Delete();
}

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *image = vtkImageData::New();
    image->SetDimensions(width, height, 1);
    //image->SetWholeExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetUpdateExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetNumberOfScalarComponents(3);
    image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
    //image->AllocateScalars();

    return image;
}

// ****************************************************************************
//  Function: ApplyBlueHotColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using the blue 
//     hot color map.
//
//     The blue hot color map has:
//        F=0: (0,0,128) 
//        F=1: (255,255,255) 
//       and smooth interpolation in between
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
// ****************************************************************************

void
ApplyBlueHotColorMap(float F, unsigned char *RGB)
{
    for (int i=0; i <= 1; i++){
        RGB[i] = F*255;
    }
    RGB[2] = 128 + F*(255-128);
}


// ****************************************************************************
//  Function: ApplyDifferenceColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using a divergent colormap
//
//     The divergent color map has:
//        F=0: (0,0,128) 
//        F=0.5: (255,255,255) 
//        F=1: (128, 0, 0)
//       and smooth interpolation in between
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
// ****************************************************************************
void
ApplyDifferenceColorMap(float F, unsigned char *RGB)
{
        if (F>=0 && F<= 0.5){
            float ratio_1 = (F-0.0)/0.5;
            RGB[0] = ratio_1*(255-0);
            RGB[1] = ratio_1*(255-0);
            RGB[2] = 128 + ratio_1*(255-128);
        }
        else{
            float ratio_2 = (F-0.5)/0.5;
            RGB[0] = 255 + ratio_2*(128-255);
            RGB[1] = 255 + ratio_2*(0-255);
            RGB[2] = 255 + ratio_2*(0-255);
        }
    
}

// ****************************************************************************
//  Function: ApplyBHSVColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using an HSV rainbow colormap
//
//     The rainbow colormap uses a saturation =1.0, value = 1.0, 
//     and interpolates hue from 0 to 360 degrees 
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
// ****************************************************************************
void
ApplyHSVColorMap(float F, unsigned char *RGB)
{
    float hue = F*360.0f;
    float saturation = 1.0f;
    float value = 1.0f;

    hue /= 60.f;

    int i = floor(hue);
    float f = hue - i;

    float r,g,b;
    float v = 1.0f;

    float p = value * (1 - saturation);
    float q = value * (1 - saturation * f);
    float t = value * (1 - saturation * (1 - f));

    switch(i) 
    {
        case 0: r = v; g = t; b = p;
        break;
        case 1: r = q; g = v; b = p;
        break;
        case 2: r = p; g = v; b = t;
        break;
        case 3: r = p; g = q; b = v;
        break;
        case 4: r = t; g = p; b = v;
        break;
        case 5: r = v; g = p; b = q;
        break;
    }

    RGB[0] = r * 255;
    RGB[1] = g * 255;
    RGB[2] = b * 255;
}


int main()
{
    int  i, j;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj3_data.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    
    int nx = 500;
    int ny = 500;

    vtkImageData *images[3];
    unsigned char *buffer[3];
    for (i = 0 ; i < 3 ; i++)
    {
        images[i] = NewImage(nx, ny);
        buffer[i] = (unsigned char *) images[i]->GetScalarPointer(0,0,0);
    }

    for (i = 0 ; i < 3*nx*ny ; i++)
        for (j = 0 ; j < 3 ; j++)
            buffer[j][i] = 0;

    for (i = 0 ; i < nx ; i++)
        for (j = 0 ; j < ny ; j++)
        {
            // ITERATE OVER PIXELS
            float pt[2];
            pt[0] = -9 + float(i)/float(nx-1)*18.0;
            pt[1] = -9 + float(j)/float(ny-1)*18.0;
            float f = EvaluateFieldAtLocation(pt, dims, X, Y, F);
            float normalizedF = float(f - 1.2)/float(5.02-1.2); //...; see step 5 re 1.2->5.02
            
            // I TAKE OVER HERE
            int offset = 3*(j*nx+i);
            ApplyBlueHotColorMap(normalizedF, buffer[0]+offset);
            ApplyDifferenceColorMap(normalizedF, buffer[1]+offset);
            ApplyHSVColorMap(normalizedF, buffer[2]+offset);
        }

    WriteImage(images[0], "bluehot");
    WriteImage(images[1], "difference");
    WriteImage(images[2], "hsv");
}
