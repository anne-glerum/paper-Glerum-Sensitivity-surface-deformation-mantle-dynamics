#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Partitioned Unstructured Grid Reader'
solution00000pvtu = XMLPartitionedUnstructuredGridReader(FileName=['/Users/visu/Glerum/01092016/aspect/EGU/K08_notomo_2/solution/solution-00000.pvtu'])
solution00000pvtu.PointArrayStatus = ['density', 'viscosity']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
solution00000pvtuDisplay = Show(solution00000pvtu, renderView1)
# trace defaults for the display properties.
solution00000pvtuDisplay.AmbientColor = [0.0, 0.0, 0.0]
solution00000pvtuDisplay.ColorArrayName = [None, '']
solution00000pvtuDisplay.EdgeColor = [0.0, 0.0, 0.0]
solution00000pvtuDisplay.GlyphType = 'Arrow'
solution00000pvtuDisplay.CubeAxesColor = [0.0, 0.0, 0.0]
solution00000pvtuDisplay.ScalarOpacityUnitDistance = 51152.954947465136
solution00000pvtuDisplay.SetScaleArray = ['POINTS', 'C_1']
solution00000pvtuDisplay.ScaleTransferFunction = 'PiecewiseFunction'
solution00000pvtuDisplay.OpacityArray = ['POINTS', 'C_1']
solution00000pvtuDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# reset view to fit data
renderView1.ResetCamera()

# create a new 'Plot Over Line'
plotOverLine1 = PlotOverLine(Input=solution00000pvtu,
    Source='High Resolution Line Source')

# init the 'High Resolution Line Source' selected for 'Source'
plotOverLine1.Source.Point1 = [3064970.1310329866, -2179010.3331278353, -2179010.3331278353]
plotOverLine1.Source.Point2 = [6371000.000000001, 2179010.3331278353, 2179010.3331278353]

# Properties modified on plotOverLine1
plotOverLine1.Tolerance = 2.22044604925031e-16

# Properties modified on plotOverLine1.Source
plotOverLine1.Source.Point1 = [3471000.0, 0.0, 0.0]
plotOverLine1.Source.Point2 = [6371000.000000001, 0.0, 0.0]
plotOverLine1.Source.Resolution = 200

# get layout
viewLayout1 = GetLayout()

# Create a new 'Line Chart View'
lineChartView1 = CreateView('XYChartView')
lineChartView1.ViewSize = [752, 818]

# place view in the layout
viewLayout1.AssignView(2, lineChartView1)

# show data in view
plotOverLine1Display = Show(plotOverLine1, lineChartView1)
# trace defaults for the display properties.
plotOverLine1Display.CompositeDataSetIndex = [0]
plotOverLine1Display.UseIndexForXAxis = 0
plotOverLine1Display.XArrayName = 'arc_length'
plotOverLine1Display.SeriesVisibility = ['density', 'viscosity']
plotOverLine1Display.SeriesLabel = ['arc_length', 'arc_length', 'density', 'density', 'viscosity', 'viscosity', 'vtkValidPointMask', 'vtkValidPointMask', 'Points_X', 'Points_X', 'Points_Y', 'Points_Y', 'Points_Z', 'Points_Z', 'Points_Magnitude', 'Points_Magnitude']
plotOverLine1Display.SeriesColor = ['arc_length', '0', '0', '0', 'density', '0.3', '0.69', '0.29', 'viscosity', '0', '0', '0', 'vtkValidPointMask', '0.89', '0.1', '0.11', 'Points_X', '0.22', '0.49', '0.72', 'Points_Y', '0.3', '0.69', '0.29', 'Points_Z', '0.6', '0.31', '0.64', 'Points_Magnitude', '1', '0.5', '0']
plotOverLine1Display.SeriesPlotCorner = ['arc_length', '0', 'density', '0', 'viscosity', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']
plotOverLine1Display.SeriesLineStyle = ['arc_length', '1', 'density', '1', 'viscosity', '1', 'vtkValidPointMask', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'Points_Magnitude', '1']
plotOverLine1Display.SeriesLineThickness = ['arc_length', '2', 'density', '2', 'viscosity', '2', 'vtkValidPointMask', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'Points_Magnitude', '2']
plotOverLine1Display.SeriesMarkerStyle = ['arc_length', '0', 'density', '0', 'viscosity', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']

# set active view
SetActiveView(renderView1)

# set active source
SetActiveSource(solution00000pvtu)

# set active view
SetActiveView(lineChartView1)

# save data
SaveData('/Users/visu/Glerum/01092016/aspect/EGU/K08_notomo_2/slices/profiles_K08_notomo_2__0.csv', proxy=plotOverLine1)

# create a new 'Plot Over Line'
plotOverLine2 = PlotOverLine(Input=solution00000pvtu,
    Source='High Resolution Line Source')

# Properties modified on plotOverLine2
plotOverLine2.Tolerance = 2.22044604925031e-16

# Properties modified on plotOverLine2.Source
# Nb
plotOverLine2.Source.Point2 = [6.07942e6, -79945.1, -1.90367e6]
plotOverLine2.Source.Point1 = [plotOverLine2.Source.Point2[0]*3471.0/6371.0,plotOverLine2.Source.Point2[1]*3471.0/6371.0,plotOverLine2.Source.Point2[2]*3471.0/6371.0]
plotOverLine2.Source.Resolution = 200

# get layout
viewLayout2 = GetLayout()

# Create a new 'Line Chart View'
lineChartView2 = CreateView('XYChartView')
lineChartView2.ViewSize = [752, 818]

# place view in the layout
viewLayout2.AssignView(2, lineChartView2)

# show data in view
plotOverLine2Display = Show(plotOverLine2, lineChartView2)
# trace defaults for the display properties.
plotOverLine2Display.CompositeDataSetIndex = [0]
plotOverLine2Display.UseIndexForXAxis = 0
plotOverLine2Display.XArrayName = 'arc_length'
plotOverLine2Display.SeriesVisibility = ['density', 'viscosity']
plotOverLine2Display.SeriesLabel = ['arc_length', 'arc_length', 'density', 'density', 'viscosity', 'viscosity', 'vtkValidPointMask', 'vtkValidPointMask', 'Points_X', 'Points_X', 'Points_Y', 'Points_Y', 'Points_Z', 'Points_Z', 'Points_Magnitude', 'Points_Magnitude']
plotOverLine2Display.SeriesColor = ['arc_length', '0', '0', '0', 'density', '0.3', '0.69', '0.29', 'viscosity', '0', '0', '0', 'vtkValidPointMask', '0.89', '0.1', '0.11', 'Points_X', '0.22', '0.49', '0.72', 'Points_Y', '0.3', '0.69', '0.29', 'Points_Z', '0.6', '0.31', '0.64', 'Points_Magnitude', '1', '0.5', '0']
plotOverLine2Display.SeriesPlotCorner = ['arc_length', '0', 'density', '0', 'viscosity', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']
plotOverLine2Display.SeriesLineStyle = ['arc_length', '1', 'density', '1', 'viscosity', '1', 'vtkValidPointMask', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'Points_Magnitude', '1']
plotOverLine2Display.SeriesLineThickness = ['arc_length', '2', 'density', '2', 'viscosity', '2', 'vtkValidPointMask', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'Points_Magnitude', '2']
plotOverLine2Display.SeriesMarkerStyle = ['arc_length', '0', 'density', '0', 'viscosity', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']

# save data
SaveData('/Users/visu/Glerum/01092016/aspect/EGU/K08_notomo_2/slices/profiles_K08_notomo_2__Nb.csv', proxy=plotOverLine2)

# create a new 'Plot Over Line'
plotOverLine3 = PlotOverLine(Input=solution00000pvtu,
    Source='High Resolution Line Source')

# Properties modified on plotOverLine2
plotOverLine3.Tolerance = 2.22044604925031e-16

# Properties modified on plotOverLine2.Source
# Nb
plotOverLine3.Source.Point2 = [6.08711e6, -373318, 1.78707e6]
plotOverLine3.Source.Point1 = [plotOverLine3.Source.Point2[0]*3471.0/6371.0,plotOverLine3.Source.Point2[1]*3471.0/6371.0,plotOverLine3.Source.Point2[2]*3471.0/6371.0]
plotOverLine3.Source.Resolution = 200

# get layout
viewLayout3 = GetLayout()

# Create a new 'Line Chart View'
lineChartView3 = CreateView('XYChartView')
lineChartView3.ViewSize = [752, 818]

# place view in the layout
viewLayout3.AssignView(2, lineChartView3)

# show data in view
plotOverLine3Display = Show(plotOverLine3, lineChartView3)
# trace defaults for the display properties.
plotOverLine3Display.CompositeDataSetIndex = [0]
plotOverLine3Display.UseIndexForXAxis = 0
plotOverLine3Display.XArrayName = 'arc_length'
plotOverLine3Display.SeriesVisibility = ['density', 'viscosity']
plotOverLine3Display.SeriesLabel = ['arc_length', 'arc_length', 'density', 'density', 'viscosity', 'viscosity', 'vtkValidPointMask', 'vtkValidPointMask', 'Points_X', 'Points_X', 'Points_Y', 'Points_Y', 'Points_Z', 'Points_Z', 'Points_Magnitude', 'Points_Magnitude']
plotOverLine3Display.SeriesColor = ['arc_length', '0', '0', '0', 'density', '0.3', '0.69', '0.29', 'viscosity', '0', '0', '0', 'vtkValidPointMask', '0.89', '0.1', '0.11', 'Points_X', '0.22', '0.49', '0.72', 'Points_Y', '0.3', '0.69', '0.29', 'Points_Z', '0.6', '0.31', '0.64', 'Points_Magnitude', '1', '0.5', '0']
plotOverLine3Display.SeriesPlotCorner = ['arc_length', '0', 'density', '0', 'viscosity', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']
plotOverLine3Display.SeriesLineStyle = ['arc_length', '1', 'density', '1', 'viscosity', '1', 'vtkValidPointMask', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'Points_Magnitude', '1']
plotOverLine3Display.SeriesLineThickness = ['arc_length', '2', 'density', '2', 'viscosity', '2', 'vtkValidPointMask', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'Points_Magnitude', '2']
plotOverLine3Display.SeriesMarkerStyle = ['arc_length', '0', 'density', '0', 'viscosity', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']

# save data
SaveData('/Users/visu/Glerum/01092016/aspect/EGU/K08_notomo_2/slices/profiles_K08_notomo_2__NEu.csv', proxy=plotOverLine3)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [4717985.065516494, 0.0, 13511161.326738741]
renderView1.CameraFocalPoint = [4717985.065516494, 0.0, 0.0]
renderView1.CameraParallelScale = 3496945.8728126283

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
