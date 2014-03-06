/*=============================================================

  Written by Jordan Weiler

  This produces 4 renderers inside a render window:
    - bottom left: isosurface
    - top left: slices from a rectilinear grid
    - top right: streamlines
    - bottom right: hedgehog glyphs

=============================================================*/

#include <vtkPolyData.h>
#include <vtkDataSetReader.h>
#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSurfaceReconstructionFilter.h>
#include <vtkContourFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkReverseSense.h>
#include <vtkMarchingCubes.h>
#include <vtkRectilinearGrid.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkHedgeHog.h>
#include <vtkProperty.h>
#include <vtkArrowSource.h>
#include <vtkGlyph3D.h>
#include <vtkExtractVOI.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkMaskPoints.h>
#include <vtkPlane.h>
#include <vtkClipPolyData.h>
#include <vtkCutter.h>
#include <vtkAppendPolyData.h>
#include <vtkLookupTable.h>
#include <vtkStreamLine.h>
#include <vtkPlaneSource.h>
#include <vtkInitialValueProblemSolver.h>
#include <vtkRungeKutta4.h>

int main(int, char *[])
{
    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj7.vtk");
    rdr->Update();

    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();

    // hardyglobal scalars are used for render 1 and render 2
    rgrid->GetPointData()->SetActiveScalars("hardyglobal");

    // Start of Render 1
    vtkSmartPointer<vtkContourFilter> cf = vtkSmartPointer<vtkContourFilter>::New();
    cf->SetInputConnection(rgrid->GetProducerPort());
    cf->SetValue(0, 2.5);
    cf->SetValue(1, 5.0);
    cf->Update();

    vtkSmartPointer<vtkPolyDataMapper> map1 = vtkSmartPointer<vtkPolyDataMapper>::New();
    map1->SetInputConnection(cf->GetOutputPort());
    map1->SetScalarRange(rgrid->GetScalarRange());
    map1->Update();

    // Start of Render 2
    vtkSmartPointer<vtkPlane> plane1 = vtkSmartPointer<vtkPlane>::New();
    plane1->SetOrigin(0, 0, 0);
    plane1->SetNormal(1,0,0);

    vtkSmartPointer<vtkPlane> plane2 = vtkSmartPointer<vtkPlane>::New();
    plane2->SetOrigin(0, 0, 0);
    plane2->SetNormal(0,1,0);
    
    vtkSmartPointer<vtkPlane> plane3 = vtkSmartPointer<vtkPlane>::New();
    plane3->SetOrigin(0, 0, 0);
    plane3->SetNormal(0,0,1);
    
    vtkSmartPointer<vtkCutter> cutter1 = vtkSmartPointer<vtkCutter>::New();
    cutter1->SetInputConnection(rgrid->GetProducerPort());
    cutter1->SetCutFunction(plane1);

    vtkSmartPointer<vtkCutter> cutter2 = vtkSmartPointer<vtkCutter>::New();
    cutter2->SetInputConnection(rgrid->GetProducerPort());
    cutter2->SetCutFunction(plane2);

    vtkSmartPointer<vtkCutter> cutter3 = vtkSmartPointer<vtkCutter>::New();
    cutter3->SetInputConnection(rgrid->GetProducerPort());
    cutter3->SetCutFunction(plane3);

    vtkSmartPointer<vtkAppendPolyData> appender = vtkSmartPointer<vtkAppendPolyData>::New();
    appender->AddInputConnection(cutter1->GetOutputPort());
    appender->AddInputConnection(cutter2->GetOutputPort());
    appender->AddInputConnection(cutter3->GetOutputPort());

    vtkSmartPointer<vtkPolyDataMapper> cutterMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    cutterMapper->SetInputConnection(appender->GetOutputPort());
    cutterMapper->SetScalarRange(rgrid->GetScalarRange());
    cutterMapper->Update();


    // grad vectors are used for render 3 and render 4
    rgrid->GetPointData()->SetActiveVectors("grad");


    // Start of Render 3
    vtkSmartPointer<vtkMaskPoints> maskPoints = vtkSmartPointer<vtkMaskPoints>::New();
    maskPoints->SetOnRatio(131);
    maskPoints->SetInputConnection(rgrid->GetProducerPort());
    maskPoints->Update();
    
  
    vtkSmartPointer<vtkArrowSource> arrowSource = vtkSmartPointer<vtkArrowSource>::New();

    vtkSmartPointer<vtkGlyph3D> g3 = vtkSmartPointer<vtkGlyph3D>::New();
    g3->SetSourceConnection(arrowSource->GetOutputPort());
    g3->SetInputConnection(maskPoints->GetOutputPort());
    g3->SetScaleFactor(0.4);

    vtkSmartPointer<vtkPolyDataMapper> glyphMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    glyphMapper->SetInputConnection(g3->GetOutputPort());
    glyphMapper->SetScalarRange(rgrid->GetScalarRange());
    glyphMapper->Update();


    // Start of Render 4
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (int i = 0; i < 19; i++)
        points->InsertNextPoint(-9 + i, 0, 0);
   
    vtkSmartPointer<vtkPolyData> polyD = vtkSmartPointer<vtkPolyData>::New();
    polyD->SetPoints(points);

    vtkSmartPointer<vtkRungeKutta4> RungeKutta4 = vtkSmartPointer<vtkRungeKutta4>::New();

    vtkSmartPointer<vtkStreamLine> streamLine = vtkSmartPointer<vtkStreamLine>::New();
    streamLine->SetInputConnection(rgrid->GetProducerPort());
    streamLine->SetSource(polyD);
    streamLine->SetMaximumPropagationTime(150);
    streamLine->SetStepLength(.01);
    streamLine->SetNumberOfThreads(1);
    streamLine->SetIntegrator(RungeKutta4);

    vtkSmartPointer<vtkPolyDataMapper> streamMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    streamMapper->SetInputConnection(streamLine->GetOutputPort());
    streamMapper->SetScalarRange(rgrid->GetScalarRange());
    streamMapper->Update();


    // Set up actors
    vtkSmartPointer<vtkActor> ren1Actor = vtkSmartPointer<vtkActor>::New();
    vtkSmartPointer<vtkActor> ren2Actor = vtkSmartPointer<vtkActor>::New();
    vtkSmartPointer<vtkActor> ren3Actor = vtkSmartPointer<vtkActor>::New();
    vtkSmartPointer<vtkActor> ren4Actor = vtkSmartPointer<vtkActor>::New();

    // Set up mappers
    ren1Actor->SetMapper(map1);
    ren2Actor->SetMapper(cutterMapper);
    ren3Actor->SetMapper(glyphMapper);
    ren4Actor->SetMapper(streamMapper);

    // Set up render window and interactor
    vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
    vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renWin);

    // Render 1: Isosurface of hardyglobal with isovalues 2.5 and 5.0. Any color
    vtkSmartPointer<vtkRenderer> ren1 = vtkSmartPointer<vtkRenderer>::New();

    // Render 2: 3 slices of hardyglobal at x=0, y=0, and z=0. Rainbow colormap (default one).
    vtkSmartPointer<vtkRenderer> ren2 = vtkSmartPointer<vtkRenderer>::New();

    // Render 3: Hedgehog glyphs of the variable grad. Choose density and colors.
    vtkSmartPointer<vtkRenderer> ren3 = vtkSmartPointer<vtkRenderer>::New();

    // Render 4: Streamlines of the variable grad. Use RK4 for integration. Seed locations should be in a line
    // from (-9, 0, 0) to (9, 0, 0). Should be 19 total seeds, over each integer in the range.
    vtkSmartPointer<vtkRenderer> ren4 = vtkSmartPointer<vtkRenderer>::New();
    
    renWin->AddRenderer(ren1);
    renWin->AddRenderer(ren2);
    renWin->AddRenderer(ren3);
    renWin->AddRenderer(ren4);

    bool doRender1 = true;
    bool doRender2 = true;
    bool doRender3 = true;
    bool doRender4 = true;

    if (doRender1)
        ren1->AddActor(ren1Actor);
    if (doRender2)
        ren2->AddActor(ren2Actor);
    if (doRender3)
        ren3->AddActor(ren3Actor);
    if (doRender4)
        ren4->AddActor(ren4Actor);

    // Set viewports
    ren1->SetViewport(0.0, 0.0, 0.5, 0.5);
    ren2->SetViewport(0.0, 0.5, 0.5, 1.0); 
    ren3->SetViewport(0.5, 0.0, 1.0, 0.5);
    ren4->SetViewport(0.5, 0.5, 1.0, 1.0);
    
    renWin->SetSize(1000,1000);

    iren->Initialize();
    iren->Start();

    return EXIT_SUCCESS;
}


