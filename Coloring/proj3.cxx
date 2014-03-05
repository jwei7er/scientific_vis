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
// 
// Modified by Jordan Weiler
// November 1, 2013
//
// This program creates three different images based on three different
// colormaps:
//
//     Blue hot colormap
//     Divergent colormap
//     HSV rainbow colormap
//
// ****************************************************************************

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

// ****************************************************************************
//  Function: GetProportion
//
//  Arguments:
//     x: value in question
//	   a: low boundary value
//     b: high boundary value
//
//  Returns: the proportional value of x between a and b
// ****************************************************************************

float GetProportion(float x, float a, float b)
{
    return (x-a) / (b-a);
}

// ****************************************************************************
//  Function: Interpolate
//
//  Arguments:
//     x: value to interpolate
//     a: low boundary value
//     f_a: F(a) value
//     b: high boundary value
//     f_b: F(b) value
//
//  Returns: the interpolated value F(x)
// ****************************************************************************

float Interpolate(float x, float a, float f_a, float b, float f_b)
{
    return f_a + GetProportion(x, a, b) * (f_b - f_a);
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

float EvaluateFieldAtLocation(const float *pt, const int *dims, const float *X, 
                              const float *Y, const float *F)
{
    float value = 0.0;
    
    int minX, maxX, minY, maxY;
    
    bool foundX = false;
    for (int i = 0; i < dims[0] - 1; i++) {
//     	cerr << "X i: " << i << " X[i]: " << X[i] << " pt[0]: " << pt[0] << " X[i+1]: " << X[i+1] << endl;
    
    	if (X[i] <= pt[0] && pt[0] < X[i+1]) {
    		minX = i;
    		maxX = i+1;
    		foundX = true;
    		break;
    	}
    }
    
    if (!foundX) {
    	return value;
    }
    
    bool foundY = false;
    for (int i = 0; i < dims[1] - 1; i++) {
    	if (Y[i] <= pt[1] && pt[1] < Y[i+1]) {
    		minY = i;
    		maxY = i+1;
    		foundY = true;
    		break;
    	}
    }
    
    if (!foundY) {
    	return value;
    }
    
    int lMinIndex = minY * dims[0] + minX;
	int lMaxIndex = maxY * dims[0] + minX;
	int rMinIndex = minY * dims[0] + maxX;
	int rMaxIndex = maxY * dims[0] + maxX;
	
	// interpolate the top corners and bottom corners of cell to X pt value
	// then interpolate from top interpolation to bottom interpolation to Y pt value
	float bottomF = Interpolate(pt[0], X[minX], F[lMinIndex], X[maxX], F[rMinIndex]);
	float topF = Interpolate(pt[0], X[minX], F[lMaxIndex], X[maxX], F[rMaxIndex]);
	value = Interpolate(pt[1], Y[minY], bottomF, Y[maxY], topF);
    
    return value;
}

// ****************************************************************************
//  Function: WriteImage
//
//  Purpose: 
//     Writes out a PNG image.
//
//  Arguments:
//       img (input):      vtkImageData used to create PNG image
//       filename (input): filename for PNG image
//      
// ****************************************************************************

void WriteImage(vtkImageData *img, const char *filename)
{
    std::string full_filename = filename;
    full_filename += ".png";
    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInput(img);
    writer->SetFileName(full_filename.c_str());
    writer->Write();
    writer->Delete();
}

// ****************************************************************************
//  Function: NewImage
//
//  Purpose: 
//     Creates and returns a vtkImageData image from given width height 
//     dimensions.
//
//  Arguments:
//       width (input):  image width
//       height (input): image height
//      
// ****************************************************************************

vtkImageData * NewImage(int width, int height)
{
    vtkImageData *image = vtkImageData::New();
    image->SetDimensions(width, height, 1);
    image->SetWholeExtent(0, width-1, 0, height-1, 0, 0);
    image->SetUpdateExtent(0, width-1, 0, height-1, 0, 0);
    image->SetNumberOfScalarComponents(3);
    image->SetScalarType(VTK_UNSIGNED_CHAR);
    image->AllocateScalars();

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

void ApplyBlueHotColorMap(float F, unsigned char *RGB)
{
    int rMin, rMax, gMin, gMax, bMin, bMax;
    rMin = 0;
    rMax = 255;
    gMin = 0;
    gMax = 255;
    bMin = 128;
    bMax = 255;
    
    RGB[0] = rMin + F * (rMax - rMin);
    RGB[1] = gMin + F * (gMax - gMin);
    RGB[2] = bMin + F * (bMax - bMin);
}

// ****************************************************************************
//  Function: ApplyDifferenceColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using a divergent 
//     colormap.
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

void ApplyDifferenceColorMap(float F, unsigned char *RGB)
{
    int rMin, rMax, gMin, gMax, bMin, bMax;
    float t = 0.0;
    if (F < 0.5) {
        rMin = 0;
        rMax = 255;
        gMin = 0;
        gMax = 255;
        bMin = 128;
        bMax = 255;
        t = F/0.5;
    } else if (F > 0.5) {
        rMin = 255;
        rMax = 128;
        gMin = 255;
        gMax = 0;
        bMin = 255;
        bMax = 0;
        t = (F-0.5)/(1.0 - 0.5);
    } else {
        rMin = 255;
        rMax = 255;
        gMin = 255;
        gMax = 255;
        bMin = 255;
        bMax = 255;
        t = 1.0;
    }
    
    RGB[0] = rMin + t * (rMax - rMin);
    RGB[1] = gMin + t * (gMax - gMin);
    RGB[2] = bMin + t * (bMax - bMin);
}

// ****************************************************************************
//  Function: ApplyBHSVColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using an HSV 
//     rainbow colormap.
//
//     The rainbow colormap uses a saturation =1.0, value = 1.0, 
//     and interpolates hue from 0 to 360 degrees 
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
// ****************************************************************************

void ApplyHSVColorMap(float F, unsigned char *RGB)
{
    float hue = 360 * F / 60.0;
    float saturation = 1.0;
    float value = 1.0;
    
    int i = floor(hue);
    float f = hue - i;
    
    // cerr << "F: " << F << " hue: " << hue << " i: " << i << endl;
    
    float p = 255 * value * (1 - saturation);
    float q = 255 * value * (1 - saturation * f);
    float t = 255 * value * (1 - saturation * (1 - f));
    
    switch(i) {
        case 0:
            RGB[0] = 255 * value;
            RGB[1] = t;
            RGB[2] = p;
            break;
        case 1:
            RGB[0] = q;
            RGB[1] = 255 * value;
            RGB[2] = p;
            break;
        case 2:
            RGB[0] = p;
            RGB[1] = 255 * value;
            RGB[2] = t;
            break;
        case 3:
            RGB[0] = p;
            RGB[1] = q;
            RGB[2] = 255 * value;
            break;
        case 4:
            RGB[0] = t;
            RGB[1] = p;
            RGB[2] = 255 * value;
            break;
        case 5:
            RGB[0] = 255 * value;
            RGB[1] = p;
            RGB[2] = q;
            break;
    }
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
            float A = -9.0;
            float B = 9.0;
            pt[0] = A + (i/(nx - 1.0)) * (B - A);
            pt[1] = A + (j/(ny - 1.0)) * (B - A);
            
			float ePt[1][3] = {{pt[0], pt[1], 0}};
			float f = EvaluateFieldAtLocation(ePt[0], dims, X, Y, F);

            float fMin = 1.2;
            float fMax = 5.02;
            float normalizedF = (f-fMin)/(fMax-fMin);
            
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
