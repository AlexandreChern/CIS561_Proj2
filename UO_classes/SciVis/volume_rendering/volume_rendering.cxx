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
#include <vtkPNGWriter.h>

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkRectilinearGrid.h>
#include <vtkDataSetWriter.h>
#include <vtkRectilinearGridToTetrahedra.h>
#include <vtkUnstructuredGrid.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <sys/time.h>


// #include <iostream.h>


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
    return dims[0]*dims[1]*dims[2];
    // 2D
    //return dims[0]*dims[1];
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
    return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    //return (dims[0]-1)*(dims[1]-1);
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
    return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    //return idx[1]*dims[0]+idx[0];
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
    return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
    // 2D
    //return idx[1]*(dims[0]-1)+idx[0];
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
     idx[0] = pointId%dims[0];
     idx[1] = (pointId/dims[0])%dims[1];
     idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    // idx[0] = pointId%dims[0];
    // idx[1] = pointId/dims[0];
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
    idx[0] = cellId%(dims[0]-1);
    idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    //idx[0] = cellId%(dims[0]-1);
    //idx[1] = cellId/(dims[0]-1);
}



struct Camera
{
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];
};


struct TransferFunction
{
    double          min;
    double          max;
    int             numBins;
    unsigned char  *colors;  // size is 3*numBins
    double         *opacities; // size is numBins

    // Take in a value and applies the transfer function.
    // Step #1: figure out which bin "value" lies in.
    // If "min" is 2 and "max" is 4, and there are 10 bins, then
    //   bin 0 = 2->2.2
    //   bin 1 = 2.2->2.4
    //   bin 2 = 2.4->2.6
    //   bin 3 = 2.6->2.8
    //   bin 4 = 2.8->3.0
    //   bin 5 = 3.0->3.2
    //   bin 6 = 3.2->3.4
    //   bin 7 = 3.4->3.6
    //   bin 8 = 3.6->3.8
    //   bin 9 = 3.8->4.0
    // and, for example, a "value" of 3.15 would return the color in bin 5
    // and the opacity at "opacities[5]".

    int GetBin(double value){
        double range_of_bins = (max - min)/numBins;
        return (int)((value - min)/range_of_bins);
    }

    void ApplyTransferFunction(double value, unsigned char *RGB, double &opacity)
    {
        int bin = GetBin(value);
        RGB[0] = colors[3*bin+0];
        RGB[1] = colors[3*bin+1];
        RGB[2] = colors[3*bin+2];
        opacity = opacities[bin];
    }
};

TransferFunction
SetupTransferFunction(void)
{
    int  i;

    TransferFunction rv;
    rv.min = 10;
    rv.max = 15;
    rv.numBins = 256;
    rv.colors = new unsigned char[3*256];
    rv.opacities = new double[256];
    unsigned char charOpacity[256] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 13, 14, 14, 14, 14, 14, 14, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 5, 4, 3, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 17, 17, 17, 17, 17, 17, 16, 16, 15, 14, 13, 12, 11, 9, 8, 7, 6, 5, 5, 4, 3, 3, 3, 4, 5, 6, 7, 8, 9, 11, 12, 14, 16, 18, 20, 22, 24, 27, 29, 32, 35, 38, 41, 44, 47, 50, 52, 55, 58, 60, 62, 64, 66, 67, 68, 69, 70, 70, 70, 69, 68, 67, 66, 64, 62, 60, 58, 55, 52, 50, 47, 44, 41, 38, 35, 32, 29, 27, 24, 22, 20, 20, 23, 28, 33, 38, 45, 51, 59, 67, 76, 85, 95, 105, 116, 127, 138, 149, 160, 170, 180, 189, 198, 205, 212, 217, 221, 223, 224, 224, 222, 219, 214, 208, 201, 193, 184, 174, 164, 153, 142, 131, 120, 109, 99, 89, 79, 70, 62, 54, 47, 40, 35, 30, 25, 21, 17, 14, 12, 10, 8, 6, 5, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        };

    for (i = 0 ; i < 256 ; i++)
        rv.opacities[i] = charOpacity[i]/255.0;
    const int numControlPoints = 8;
    unsigned char controlPointColors[numControlPoints*3] = { 
           71, 71, 219, 0, 0, 91, 0, 255, 255, 0, 127, 0, 
           255, 255, 0, 255, 96, 0, 107, 0, 0, 224, 76, 76 
       };
    double controlPointPositions[numControlPoints] = { 0, 0.143, 0.285, 0.429, 0.571, 0.714, 0.857, 1.0 };
    for (i = 0 ; i < numControlPoints-1 ; i++)
    {
        int start = controlPointPositions[i]*rv.numBins;
        int end   = controlPointPositions[i+1]*rv.numBins+1;
cerr << "Working on " << i << "/" << i+1 << ", with range " << start << "/" << end << endl;
        if (end >= rv.numBins)
            end = rv.numBins-1;
        for (int j = start ; j <= end ; j++)
        {
            double proportion = (j/(rv.numBins-1.0)-controlPointPositions[i])/(controlPointPositions[i+1]-controlPointPositions[i]);
            if (proportion < 0 || proportion > 1.)
                continue;
            for (int k = 0 ; k < 3 ; k++)
                rv.colors[3*j+k] = proportion*(controlPointColors[3*(i+1)+k]-controlPointColors[3*i+k])
                                 + controlPointColors[3*i+k];
        }
    }    

    return rv;
}

Camera
SetupCamera(void)
{
    Camera rv;
    rv.focus[0] = 0;
    rv.focus[1] = 0;
    rv.focus[2] = 0;
    rv.up[0] = 0;
    rv.up[1] = -1;
    rv.up[2] = 0;
    rv.angle = 30;
    rv.near = 7.5e+7;
    rv.far = 1.4e+8;
    rv.position[0] = -8.25e+7;
    rv.position[1] = -3.45e+7;
    rv.position[2] = 3.35e+7;

    return rv;
}

void cross_product(const double *vec_A, const double *vec_B, double *cross_product){
    cross_product[0] = (vec_A[1]*vec_B[2] - vec_A[2]*vec_B[1]);
    cross_product[1] = (vec_A[2]*vec_B[0] - vec_A[0]*vec_B[2]);
    cross_product[2] = (vec_A[0]*vec_B[1] - vec_A[1]*vec_B[0]);
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


// ********************************
// Function: LocateCell
// ********************************
void LocateCell(const float *X, const float *Y, const float *Z, int *cellID, const double *pt, const int *dims ){
    // 3D
    // Arguments: pt: 3D array of point locations;
    cellID[0] = 0;
    cellID[1] = 0;
    cellID[2] = 0;
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

    for (int k=0;k<=dims[2]-1;k++){
        if (Z[k] <= pt[2] && Z[k+1] > pt[2]){
            cellID[2] = k;
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

double
EvaluateFieldAtLocation(const double *pt, const int *dims, 
                        const float *X, const float *Y, const float *Z, const float *F)
{
    // if (pt[0] < X[0] || pt[0] > X[dims[0]-1] || pt[1] < Y[0] || Y[dims[1]-1])
    if (pt[0] < X[0] || pt[1] < Y[0] || pt[0] > X[dims[0]-1] || pt[1] > Y[dims[1]-1])
        return 0; // IMPLEMENT ME!!
    else {
        int cellID[3];
        cellID[0] = 0;
        cellID[1] = 0;
        cellID[2] = 0;
        LocateCell(X,Y,Z,cellID,pt,dims);
        std::cout << cellID[0] << '\t' << cellID[1] << '\t' << cellID[2] << std::endl;
        int pt1[3];
        int pt2[3];
        int pt3[3];
        int pt4[3];

        int pt5[3];
        int pt6[3];
        int pt7[3];
        int pt8[3];

        int pt1_real = GetPointIndex(cellID,dims);
        int pt2_real = pt1_real + 1;
        int pt3_real = pt1_real + dims[0];
        int pt4_real = pt1_real + dims[0] + 1;

        int pt5_real = pt1_real + dims[0]*dims[1];
        int pt6_real = pt2_real + dims[0]*dims[1];
        int pt7_real = pt3_real + dims[0]*dims[1];
        int pt8_real = pt4_real + dims[0]*dims[1];



        pt1[0] = cellID[0];
        pt1[1] = cellID[1];
        pt1[2] = cellID[2];

        pt2[0] = cellID[0] +1;
        pt2[1] = cellID[1];
        pt2[2] = cellID[2];

        pt3[0] = cellID[0];
        pt3[1] = cellID[1]+1;
        pt3[2] = cellID[2];

        pt4[0] = cellID[0]+1;
        pt4[1] = cellID[1]+1;
        pt4[2] = cellID[2];

        pt5[0] = pt1[0];
        pt5[1] = pt1[1];
        pt5[2] = pt1[2] + 1;

        pt6[0] = pt2[0];
        pt6[1] = pt2[1];
        pt6[2] = pt2[2] + 1;

        pt7[0] = pt3[0];
        pt7[1] = pt3[1];
        pt7[2] = pt3[2] + 1;

        pt8[0] = pt4[0];
        pt8[1] = pt4[1];
        pt8[2] = pt4[1] + 1;


        


        float prop_x = float(pt[0] - X[pt1[0]])/float(X[pt1[0]+1] - X[pt1[0]]);
        float prop_y = float(pt[1] - Y[pt1[1]])/float(Y[pt1[1]+1] - Y[pt1[1]]);
        float prop_z = float(pt[2] - Z[pt1[2]])/float(Y[pt1[2]+1] - Y[pt1[2]]);
        
        float field_1 = F[pt1_real] + prop_x*(F[pt2_real] - F[pt1_real]);
        float field_2 = F[pt3_real] + prop_x*(F[pt4_real] - F[pt3_real]);


        float final_field_1 = field_1 + prop_y*(field_2 - field_1);

        float field_3 = F[pt5_real] + prop_x*(F[pt6_real] - F[pt5_real]);
        float field_4 = F[pt7_real] + prop_x*(F[pt8_real] - F[pt7_real]);

        float final_field_2 = field_3 + prop_y*(field_4 - field_3);

        float final_field = final_field_1 + prop_z*(final_field_2 - final_field_1);

        // float field_3 = F[pt1_real] + prop_y*(F[pt3_real] - F[pt1_real]);
        // float field_4 = F[pt2_real] + prop_y*(F[pt4_real] - F[pt2_real]);

        // float final_field_2 = field_3 + prop_x*(field_4-field_3);
        // // return prop_y;
        // float difference = final_field_1 - final_field_2;
        return final_field;
    }
    // cellID[1] = 0;
    // LocateCell(X,Y, cellid, pt, dims);
}

// *************************************
// Function Euclidean Distance
// 
// Arguments:
//
// V: Vector in 3D
// 
// *************************************

double euclidean_distance(const double * V){
    double distance;
    // int length = sizeof(V)/sizeof(V[0]);
    for (int i=0;i<3;i++){
        distance += V[i]*V[i];
    }
    return sqrt(distance);
}




int main(){
    std::cout << "Hello" << endl;
    // double test_a[3] = {0,-1,0};
    // double test_b[3] = {0.863954,0.361286,-0.350814};
    // double test_c[3];
    // cross_product(test_a,test_b,test_c);
    // std::cout << test_c[0] << '\t' << test_c[1] << '\t' << test_c[2] << std::endl;

    // double A[3];
    // double B[3];
    // A[0] = 1;
    // A[1] = 0;
    // A[2] = 0;

    // B[0] = 0;
    // B[1] = 1;
    // B[2] = 0;

    // double C[3];
    // cross_product(A,B,C);
    // std::cout << C[0] << '\t' << C[1] << '\t' << C[2] << std::endl;
    // return 0;
    // okay cross_product works

    // Test euclidean

    // float v[3] = {0,2,2};
    // float distance = euclidean_distance(v);
    // std::cout <<"Euclidean distance of v is \t" << distance << endl;

    // Euclidean Works!

    int i, j, k;
    int x, y;
    int  offset;
    int sample_size = 100;

    Camera camera = SetupCamera();

    TransferFunction transferFunction = SetupTransferFunction();

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("astro64.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid * ) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *Z = (float *) rgrid->GetZCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);

    int nx = 500;
    int ny = 500;

    vtkImageData *images = NewImage(nx,ny);
    unsigned char *buffer = (unsigned char *) images->GetScalarPointer(0,0,0);

    for (i=0; i<3*nx*ny; i++){
        buffer[i]=0;
    }

    // Starting Ray Casting
    double look[3];
    double u[3],v[3];
    double delta_x[3];
    double delta_y[3];
    double euclidean_dist;

    std::cout << camera.up[0] << '\t' << camera.up[1] << '\t' << camera.up[2] << std::endl; 

    for (i=0;i<3;i++){
        look[i] = camera.focus[i] - camera.position[i];
    }
    euclidean_dist = euclidean_distance(look);
    for (i=0;i<3;i++){
        look[i] = look[i] / euclidean_dist;
    }
    std::cout << look[0] << '\t' << look[1] << '\t' << look[2] << std::endl;

    cross_product(look, camera.up, u);
    
    euclidean_dist = euclidean_distance(u);
    for (i=0;i<3;i++){
        u[i] = u[i] / euclidean_dist;
    }
    std::cout << u[0] << '\t' << u[1] << '\t' << u[2] << std::endl;
    

    cross_product(look,u,v);
    euclidean_dist = euclidean_distance(v);
     for (i=0;i<3;i++){
        v[i] = v[i] / euclidean_dist;
    }

    std::cout << v[0] << '\t' << v[1] << '\t' << v[2] << std::endl;

    for (i=0;i<3;i++){
        delta_x[i] = (2 * tan(camera.angle/360)*M_PI * u[i])/nx;
        delta_y[i] = (2 * tan(camera.angle/360)*M_PI * v[i])/ny;
    }

    unsigned char RGB[3];
    double opacity_1, opacity_2, pixel_RGB[3];

    struct timeval startTime;
    gettimeofday(&startTime,0);

    for (x=0; x< nx; x++){
        for (y=0; y<ny; y++){
            double ray[3];
            double sample_position[3];
            double eval_result;


            // Get the ray     
            for(i=0; i<3; i++){
                ray[i] = look[i] + ((2*x+1-nx)/2.0)*delta_x[i] + ((2*y+1-ny)/2.0)*delta_y[i];
            }
            euclidean_dist = euclidean_distance(ray);
            for (i=0; i<3; i++){
                ray[i] /= euclidean_dist;
            }
            
            // Get Intersection
            for (i=0;i<3;i++){
                sample_position[i] = camera.position[i] - (camera.near * ray[i]);
            }

            double step_size = (camera.far - camera.near)/(sample_size);

            eval_result = EvaluateFieldAtLocation(sample_position,dims,X,Y,Z,F);
            transferFunction.ApplyTransferFunction(eval_result,RGB,opacity_1);

            opacity_1 = 1.0 - pow((1.0 - opacity_1), (500/sample_size));

            for (i=0; i<3; i++){
                pixel_RGB[i] = (RGB[i]/255) * opacity_1;
            }

            for (i=0; i < sample_size; i++){
                for (j = 0; j < 3; j++){
                    sample_position[j] += (step_size * ray[j]);
                }
                eval_result = EvaluateFieldAtLocation(sample_position,dims,X,Y,Z,F);

                if (transferFunction.GetBin(eval_result) >= 0) {
                    transferFunction.ApplyTransferFunction(eval_result, RGB, opacity_2);
                    opacity_2 = 1.0 - pow((1.0 - opacity_2), (500/sample_size));
                }

                for (j = 0; j < 3; j++){
                    pixel_RGB[j] = pixel_RGB[j] + (1-opacity_1)*opacity_2*(RGB[j]/255.0);
                }
                opacity_1 = opacity_1 + (1 - opacity_1)*opacity_2;
            }

            offset = 3*(y*nx+x);
            buffer[offset] = (unsigned char)(pixel_RGB[0] * 255.0);
            buffer[offset+1] = (unsigned char)(pixel_RGB[1] * 255.0);
            buffer[offset+2] = (unsigned char)(pixel_RGB[2] * 255.0);
        }
    }



     


    /*****************
     * Implementing Ray Casting
     * STEP 1: Find RAY
     * STEP 2: Intersect Volume with the Ray
     * STEP 3: Calculate Color from intersection
     * STEP 4: Assign Color to that pixel
     ******************/

    WriteImage(images,"test_volume_rendering");

}