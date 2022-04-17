#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Partitioned Unstructured Grid Reader'
solution00000pvtu = XMLPartitionedUnstructuredGridReader(FileName=['/Users/visu/Glerum/01092016/aspect/EGU/K08_notomo_2/solution/solution-00000.pvtu'])

# Properties modified on solution00000pvtu
solution00000pvtu.PointArrayStatus = ['p', 'T', 'density', 'viscosity', 'nonadiabatic_temperature', 'thermal_expansivity', 'specific_heat']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1219, 822]

# show data in view
solution00000pvtuDisplay = Show(solution00000pvtu, renderView1)
# trace defaults for the display properties.
solution00000pvtuDisplay.ColorArrayName = [None, '']
solution00000pvtuDisplay.GlyphType = 'Arrow'
solution00000pvtuDisplay.ScalarOpacityUnitDistance = 38396.72469188691

# reset view to fit data
renderView1.ResetCamera()

# create a new 'Slice'
slice1 = Slice(Input=solution00000pvtu)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [4721612.875, 0.0, 0.0]

# Properties modified on slice1.SliceType
slice1.SliceType.Normal = [0.0, 1.0, 0.0]

# show data in view
slice1Display = Show(slice1, renderView1)
# trace defaults for the display properties.
slice1Display.ColorArrayName = [None, '']
slice1Display.GlyphType = 'Arrow'

# hide data in view
Hide(solution00000pvtu, renderView1)

# create a new 'Plot Over Line'
plotOverLine1 = PlotOverLine(Input=slice1,
    Source='High Resolution Line Source')

# init the 'High Resolution Line Source' selected for 'Source'
plotOverLine1.Source.Point1 = [3265797.25, 0.0, -2157930.0]
plotOverLine1.Source.Point2 = [6370500.0, 0.0, 2157930.0]

# Properties modified on plotOverLine1
plotOverLine1.Tolerance = 2.22044604925031e-16

# Properties modified on plotOverLine1.Source
plotOverLine1.Source.Point1 = [3352730.0, 0.0, -898361.0]
plotOverLine1.Source.Point2 = [6153910.0, 0.0, -1648940.0]

# get layout
viewLayout1 = GetLayout()

# Create a new 'Line Chart View'
lineChartView1 = CreateView('XYChartView')
lineChartView1.ViewSize = [605, 822]

# place view in the layout
viewLayout1.AssignView(2, lineChartView1)

# show data in view
plotOverLine1Display = Show(plotOverLine1, lineChartView1)
# trace defaults for the display properties.
plotOverLine1Display.CompositeDataSetIndex = [0]
plotOverLine1Display.UseIndexForXAxis = 0
plotOverLine1Display.XArrayName = 'arc_length'
plotOverLine1Display.SeriesVisibility = ['density', 'nonadiabatic_temperature', 'p', 'T', 'viscosity', 'thermal_expansivity', 'specific_heat']
#plotOverLine1Display.SeriesLabel = ['arc_length', 'arc_length', 'density', 'density', 'nonadiabatic_temperature', 'nonadiabatic_temperature', 'p', 'p', 'T', 'T', 'viscosity', 'viscosity', 'vtkValidPointMask', 'vtkValidPointMask', 'Points_X', 'Points_X', 'Points_Y', 'Points_Y', 'Points_Z', 'Points_Z', 'Points_Magnitude', 'Points_Magnitude']
#plotOverLine1Display.SeriesColor = ['arc_length', '0', '0', '0', 'density', '0.89', '0.1', '0.11', 'nonadiabatic_temperature', '0.22', '0.49', '0.72', 'p', '0.3', '0.69', '0.29', 'T', '0.6', '0.31', '0.64', 'viscosity', '1', '0.5', '0', 'vtkValidPointMask', '0.65', '0.34', '0.16', 'Points_X', '0', '0', '0', 'Points_Y', '0.89', '0.1', '0.11', 'Points_Z', '0.22', '0.49', '0.72', 'Points_Magnitude', '0.3', '0.69', '0.29']
#plotOverLine1Display.SeriesPlotCorner = ['arc_length', '0', 'density', '0', 'nonadiabatic_temperature', '0', 'p', '0', 'T', '0', 'viscosity', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']
#plotOverLine1Display.SeriesLineStyle = ['arc_length', '1', 'density', '1', 'nonadiabatic_temperature', '1', 'p', '1', 'T', '1', 'viscosity', '1', 'vtkValidPointMask', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'Points_Magnitude', '1']
#plotOverLine1Display.SeriesLineThickness = ['arc_length', '2', 'density', '2', 'nonadiabatic_temperature', '2', 'p', '2', 'T', '2', 'viscosity', '2', 'vtkValidPointMask', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'Points_Magnitude', '2']
#plotOverLine1Display.SeriesMarkerStyle = ['arc_length', '0', 'density', '0', 'nonadiabatic_temperature', '0', 'p', '0', 'T', '0', 'viscosity', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']

# destroy lineChartView1
Delete(lineChartView1)
del lineChartView1

# close an empty frame
viewLayout1.Collapse(2)

# set active view
SetActiveView(renderView1)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=slice1)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=slice1)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=plotOverLine1)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=plotOverLine1)

# set active source
SetActiveSource(slice1)

# create a new 'Plot Over Line'
plotOverLine2 = PlotOverLine(Input=slice1,
    Source='High Resolution Line Source')

# init the 'High Resolution Line Source' selected for 'Source'
plotOverLine2.Source.Point1 = [3265797.25, 0.0, -2157930.0]
plotOverLine2.Source.Point2 = [6370500.0, 0.0, 2157930.0]

# set active source
SetActiveSource(plotOverLine1)

# set active source
SetActiveSource(plotOverLine2)

# Properties modified on plotOverLine2.Source
plotOverLine2.Source.Point1 = [3352730.0, 0.0, 898361.0]
plotOverLine2.Source.Point2 = [6153910.0, 0.0, 1648940.0]

# Create a new 'Line Chart View'
lineChartView1 = CreateView('XYChartView')
lineChartView1.ViewSize = [605, 822]

# place view in the layout
viewLayout1.AssignView(2, lineChartView1)

# show data in view
plotOverLine2Display = Show(plotOverLine2, lineChartView1)
# trace defaults for the display properties.
plotOverLine2Display.CompositeDataSetIndex = [0]
plotOverLine2Display.UseIndexForXAxis = 0
plotOverLine2Display.XArrayName = 'arc_length'
plotOverLine2Display.SeriesVisibility = ['density', 'nonadiabatic_temperature', 'p', 'T', 'viscosity', 'thermal_expansivity', 'specific_heat']
#plotOverLine2Display.SeriesLabel = ['arc_length', 'arc_length', 'density', 'density', 'nonadiabatic_temperature', 'nonadiabatic_temperature', 'p', 'p', 'T', 'T', 'viscosity', 'viscosity', 'vtkValidPointMask', 'vtkValidPointMask', 'Points_X', 'Points_X', 'Points_Y', 'Points_Y', 'Points_Z', 'Points_Z', 'Points_Magnitude', 'Points_Magnitude']
#plotOverLine2Display.SeriesColor = ['arc_length', '0', '0', '0', 'density', '0.89', '0.1', '0.11', 'nonadiabatic_temperature', '0.22', '0.49', '0.72', 'p', '0.3', '0.69', '0.29', 'T', '0.6', '0.31', '0.64', 'viscosity', '1', '0.5', '0', 'vtkValidPointMask', '0.65', '0.34', '0.16', 'Points_X', '0', '0', '0', 'Points_Y', '0.89', '0.1', '0.11', 'Points_Z', '0.22', '0.49', '0.72', 'Points_Magnitude', '0.3', '0.69', '0.29']
#plotOverLine2Display.SeriesPlotCorner = ['arc_length', '0', 'density', '0', 'nonadiabatic_temperature', '0', 'p', '0', 'T', '0', 'viscosity', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']
#plotOverLine2Display.SeriesLineStyle = ['arc_length', '1', 'density', '1', 'nonadiabatic_temperature', '1', 'p', '1', 'T', '1', 'viscosity', '1', 'vtkValidPointMask', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'Points_Magnitude', '1']
#plotOverLine2Display.SeriesLineThickness = ['arc_length', '2', 'density', '2', 'nonadiabatic_temperature', '2', 'p', '2', 'T', '2', 'viscosity', '2', 'vtkValidPointMask', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'Points_Magnitude', '2']
#plotOverLine2Display.SeriesMarkerStyle = ['arc_length', '0', 'density', '0', 'nonadiabatic_temperature', '0', 'p', '0', 'T', '0', 'viscosity', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']

# 2021: increase resolution
plotOverLine1.Source.Resolution = 1450

# destroy lineChartView1
Delete(lineChartView1)
del lineChartView1

# close an empty frame
viewLayout1.Collapse(2)

# set active view
SetActiveView(renderView1)

# set active source
SetActiveSource(plotOverLine1)

# save data
SaveData('/Users/visu/Glerum/01092016/aspect/EGU/K08_notomo_2/slices/profile_K08_notomo_2_0_105.csv', proxy=plotOverLine1)

# set active source
SetActiveSource(plotOverLine2)

# save data
SaveData('/Users/visu/Glerum/01092016/aspect/EGU/K08_notomo_2/slices/profile_K08_notomo_2_0_75.csv', proxy=plotOverLine2)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [5913632.76423422, 13087945.195947213, -2627393.624790256]
renderView1.CameraFocalPoint = [4721612.875, 0.0, 0.0]
renderView1.CameraViewUp = [0.8960131604086753, -0.16440943181887602, -0.4124681261666799]
renderView1.CameraParallelScale = 3468739.3273047726
