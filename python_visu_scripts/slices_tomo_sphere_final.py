#### import the simple module from the paraview
from paraview.simple import *
import math
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Partitioned Unstructured Grid Reader'
solution00000pvtu = XMLPartitionedUnstructuredGridReader(FileName=['/Users/visu/Glerum/01092016/aspect/EGU/K08_SEMUCB/solution/solution-00000.pvtu'])
solution00000pvtu.CellArrayStatus = []
solution00000pvtu.PointArrayStatus = ['velocity', 'p', 'T', 'C_1', 'C_2', 'C_3', 'strain_rate', 'nonadiabatic_temperature', 'density', 'thermal_expansivity', 'specific_heat', 'viscosity']

# Properties modified on solution00000pvtu
solution00000pvtu.PointArrayStatus = ['velocity', 'T', 'C_14', 'C_15', 'nonadiabatic_temperature', 'density', 'thermal_expansivity', 'specific_heat', 'viscosity']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# create a new 'Slice'
slice2 = Slice(Input=solution00000pvtu)
slice2.SliceType = 'Sphere'
slice2.Crinkleslice = 0
slice2.Triangulatetheslice = 1
slice2.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice2.SliceType.Center = [4717985.065516494, 0.0, 0.0]
slice2.SliceType.Radius = 1.0

# Properties modified on slice2.SliceType
slice2.SliceType.Center = [0.0, 0.0, 0.0]
slice2.SliceType.Radius = 6366000.0

Hide(slice2, renderView1)

# create a new 'Contour'
contour4 = Contour(Input=slice2)
contour4.ContourBy = ['POINTS', 'C_14']
contour4.Isosurfaces = [0.5]
contour4.PointMergeMethod = 'Uniform Binning'

# properties modified on contour4
contour4.Isosurfaces = [0.1]
# show data in view
contour4Display = Show(contour4, renderView1)
# trace defaults for the display properties.
contour4Display.AmbientColor = [0.0, 0.0, 1.0]
contour4Display.ColorArrayName = ['POINTS', 'velocity']
contour4Display.LookupTable = None
contour4Display.EdgeColor = [0.0, 0.0, 1.0]
contour4Display.GlyphType = 'Arrow'
contour4Display.CubeAxesColor = [0.0, 0.0, 1.0]
contour4Display.SetScaleArray = ['POINTS', 'C_17']
contour4Display.ScaleTransferFunction = 'PiecewiseFunction'
contour4Display.OpacityArray = ['POINTS', 'C_17']
contour4Display.OpacityTransferFunction = 'PiecewiseFunction'

# turn off scalar coloring
ColorBy(contour4Display, None)
# change solid color and line width
contour4Display.DiffuseColor = [1.0, 1.0, 1.0]
contour4Display.LineWidth = 2.0

SetActiveSource(slice2)
# create a new 'Contour'
contour5 = Contour(Input=slice2)
contour5.ContourBy = ['POINTS', 'C_15']
contour5.Isosurfaces = [0.5]
contour5.PointMergeMethod = 'Uniform Binning'

# properties modified on contour5
contour5.Isosurfaces = [0.1]
# show data in view
contour5Display = Show(contour5, renderView1)
# trace defaults for the display properties.
contour5Display.AmbientColor = [0.0, 0.0, 1.0]
contour5Display.ColorArrayName = ['POINTS', 'velocity']
contour5Display.LookupTable = None
contour5Display.EdgeColor = [0.0, 0.0, 1.0]
contour5Display.GlyphType = 'Arrow'
contour5Display.CubeAxesColor = [0.0, 0.0, 1.0]
contour5Display.SetScaleArray = ['POINTS', 'C_18']
contour5Display.ScaleTransferFunction = 'PiecewiseFunction'
contour5Display.OpacityArray = ['POINTS', 'C_18']
contour5Display.OpacityTransferFunction = 'PiecewiseFunction'

# turn off scalar coloring
ColorBy(contour5Display, None)
# change solid color and line width
contour5Display.DiffuseColor = [1.0, 1.0, 1.0]
contour5Display.LineWidth = 2.0
Hide(contour5)

# create a new 'XML Partitioned Unstructured Grid Reader'
spakmanTomo_eqvtk = LegacyVTKReader(FileNames=['/Applications/paraview.app/Contents/data/spakmanTomo_eq.vtk'])

# get color transfer function/color map for 'dv'
dvLUT = GetColorTransferFunction('dv')

# invert the transfer function
dvLUT.InvertTransferFunction()


# show data in view
spakmanTomo_eqvtkDisplay = Show(spakmanTomo_eqvtk, renderView1)
# trace defaults for the display properties.
spakmanTomo_eqvtkDisplay.AmbientColor = [0.0, 0.0, 1.0]
spakmanTomo_eqvtkDisplay.ColorArrayName = ['POINTS', 'dv']
spakmanTomo_eqvtkDisplay.LookupTable = dvLUT
spakmanTomo_eqvtkDisplay.EdgeColor = [0.0, 0.0, 1.0]
spakmanTomo_eqvtkDisplay.GlyphType = 'Arrow'
spakmanTomo_eqvtkDisplay.CubeAxesColor = [0.0, 0.0, 1.0]
spakmanTomo_eqvtkDisplay.ScalarOpacityUnitDistance = 57240.10222431871
spakmanTomo_eqvtkDisplay.SetScaleArray = ['POINTS', 'dv']
spakmanTomo_eqvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
spakmanTomo_eqvtkDisplay.OpacityArray = ['POINTS', 'dv']
spakmanTomo_eqvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# reset view to fit data
renderView1.ResetCamera()

# show color bar/color legend
spakmanTomo_eqvtkDisplay.SetScalarBarVisibility(renderView1, True)

# get opacity transfer function/opacity map for 'dv'
dvPWF = GetOpacityTransferFunction('dv')

#changing interaction mode based on data extents
renderView1.InteractionMode = '3D'

# create a new 'Slice'
slice1 = Slice(Input=spakmanTomo_eqvtk)
slice1.SliceType = 'Sphere'
slice1.Crinkleslice = 0
slice1.Triangulatetheslice = 1
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Center = [0.0, 0.0, 0.0]
slice1.SliceType.Radius = 1.0

# Properties modified on slice1.SliceType
slice1.SliceType.Center = [0.0, 0.0, 0.0]
slice1.SliceType.Radius = 6366000.0

# show data in view
slice1Display = Show(slice1, renderView1)
# trace defaults for the display properties.
slice1Display.AmbientColor = [0.0, 0.0, 1.0]
slice1Display.ColorArrayName = ['POINTS', 'dv']
slice1Display.LookupTable = dvLUT
slice1Display.EdgeColor = [0.0, 0.0, 1.0]
slice1Display.GlyphType = 'Arrow'
slice1Display.CubeAxesColor = [0.0, 0.0, 1.0]
slice1Display.SetScaleArray = ['POINTS', 'dv']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = ['POINTS', 'dv']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'

# hide data in view
Hide(spakmanTomo_eqvtk, renderView1)

# set scalar coloring
ColorBy(slice1Display, ('POINTS', 'dv'))

# rescale color and/or opacity maps used to include current data range
slice1Display.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'velocity'
diffLUT = GetColorTransferFunction('dv')
diffLUT.LockDataRange = 1
diffLUT.InterpretValuesAsCategories = 0
diffLUT.ShowCategoricalColorsinDataRangeOnly = 0
diffLUT.RescaleOnVisibilityChange = 0
diffLUT.EnableOpacityMapping = 0
diffLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 0.05, 0.865003, 0.865003, 0.865003, 0.1, 0.705882, 0.0156863, 0.14902]
diffLUT.UseLogScale = 0
diffLUT.ColorSpace = 'Diverging'
diffLUT.UseBelowRangeColor = 0
diffLUT.BelowRangeColor = [0.0, 0.0, 0.0]
diffLUT.UseAboveRangeColor = 0
diffLUT.AboveRangeColor = [1.0, 1.0, 1.0]
diffLUT.NanColor = [1.0, 1.0, 0.0]
diffLUT.Discretize = 1
diffLUT.NumberOfTableValues = 32
diffLUT.ScalarRangeInitialized = 1.0
diffLUT.HSVWrap = 0
diffLUT.VectorComponent = 0
diffLUT.VectorMode = 'Magnitude'
diffLUT.AllowDuplicateScalars = 1
diffLUT.Annotations = []
diffLUT.ActiveAnnotatedValues = []
diffLUT.IndexedColors = []

# rescale transfer function
diffLUT.RescaleTransferFunction(-2.0, 2.0)
diffLUT.Annotations = ['-2.00', '-2', '-1.0', '-1.0', '0.0', '0', '1.0', '1.0', '2.0', '2']


# invert the transfer function
diffLUT.InvertTransferFunction()

# hide color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# get color legend
diffLUTColorBar = GetScalarBar(diffLUT, renderView1)

## Properties modified on diffLUTColorBar
diffLUTColorBar.Title = ''
diffLUTColorBar.ComponentTitle = ''
diffLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
diffLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
diffLUTColorBar.AutomaticLabelFormat = 0
diffLUTColorBar.LabelFormat = '%-#6.2f'
diffLUTColorBar.RangeLabelFormat = '%6.2f'
diffLUTColorBar.AutoOrient = 0
diffLUTColorBar.Orientation = 'Horizontal'
diffLUTColorBar.Position = [100, 0.85]
diffLUTColorBar.DrawTickMarks = 0
diffLUTColorBar.DrawTickLabels = 0
diffLUTColorBar.AddRangeLabels = 0
diffLUTColorBar.TitleFontSize = 10
diffLUTColorBar.LabelFontSize = 10
diffLUTColorBar.Position = [0.288, 0.05]

# create a new 'Contour'
contour1 = Contour(Input=slice1)
contour1.ContourBy = ['POINTS', 'dv']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# show data in view
contour1Display = Show(contour1, renderView1)
# trace defaults for the display properties.
contour1Display.AmbientColor = [0.0, 0.0, 0.0]
contour1Display.ColorArrayName = ['POINTS', 'dv']
contour1Display.LookupTable = diffLUT
contour1Display.EdgeColor = [0.0, 0.0, 0.0]
contour1Display.GlyphType = 'Arrow'
contour1Display.CubeAxesColor = [0.0, 0.0, 0.0]
contour1Display.SetScaleArray = ['POINTS', 'C_2']
contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
contour1Display.OpacityArray = ['POINTS', 'C_2']
contour1Display.OpacityTransferFunction = 'PiecewiseFunction'

# show color bar/color legend
contour1Display.SetScalarBarVisibility(renderView1, False)
# Properties modified on contour1
contour1.ContourBy = ['POINTS', 'dv']
contour1.Isosurfaces = [-8.0, -7.0, -6.0, 6.0, 7.0, 8.0]

# turn off scalar coloring
ColorBy(contour1Display, None)

# change solid color
contour1Display.DiffuseColor = [0.0, 0.0, 0.0]

diffLUTColorBar = GetScalarBar(diffLUT, renderView1)
diffLUT.RescaleTransferFunction(-2.0, 2.0)
slice1Display.LookupTable = diffLUT
diffLUTColorBar.Position = [0.288, 0.05]

# set active view
SetActiveView(renderView1)

# set camera
# distance from slice
far = 8.0e6
renderView1.CameraFocalPoint = [6371000.0, 0.0, 0.0]
renderView1.CameraPosition = [6371000.0 + far, 0.0, 0.0]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]

# ACG set background
renderView1.Background = [1,1,1] #White

# get Layout
viewLayout1 = GetLayout()

# Add map with indicator line
# create a new 'CSV Reader'
map_lon_latcsv = CSVReader(FileName=['/Users/visu/Glerum/090915/aspect/T_perturb/other_visu_data/map_lon_lat.csv'])

# Properties modified on map_lon_latcsv
map_lon_latcsv.HaveHeaders = 0

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.BlockSize = 1024L

# show data in view
map_lon_latcsvDisplay = Show(map_lon_latcsv, spreadSheetView1)
# trace defaults for the display properties.
map_lon_latcsvDisplay.FieldAssociation = 'Row Data'

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=map_lon_latcsv)
tableToPoints1.XColumn = 'Field 0'
tableToPoints1.YColumn = 'Field 0'
tableToPoints1.ZColumn = 'Field 0'

# Properties modified on tableToPoints1
tableToPoints1.YColumn = 'Field 1'
tableToPoints1.ZColumn = 'Field 2'

# show data in view
tableToPoints1Display = Show(tableToPoints1, spreadSheetView1)

# hide data in view
Hide(map_lon_latcsv, spreadSheetView1)

# destroy spreadSheetView1
Delete(spreadSheetView1)
del spreadSheetView1

# close an empty frame
#viewLayout1.Collapse(3)

# set active view
SetActiveView(renderView1)

# set active source
SetActiveSource(tableToPoints1)

# show data in view
tableToPoints1Display = Show(tableToPoints1, renderView1)
# trace defaults for the display properties.
tableToPoints1Display.AmbientColor = [0.0, 0.0, 0.0]
tableToPoints1Display.ColorArrayName = [None, '']
tableToPoints1Display.EdgeColor = [0.0, 0.0, 0.0]
tableToPoints1Display.GlyphType = 'Arrow'
tableToPoints1Display.CubeAxesColor = [0.0, 0.0, 0.0]
tableToPoints1Display.SetScaleArray = [None, '']
tableToPoints1Display.ScaleTransferFunction = 'PiecewiseFunction'
tableToPoints1Display.OpacityArray = [None, '']
tableToPoints1Display.OpacityTransferFunction = 'PiecewiseFunction'

# change solid color
tableToPoints1Display.DiffuseColor = [0.0, 0.0, 0.0]

# create a new 'Calculator'
calculator1 = Calculator(Input=tableToPoints1)
calculator1.Function = ''

Hide(tableToPoints1)

# Properties modified on calculator1
calculator1.Function = 'coords/6371000.0*6371000.0'

# show data in view
calculator1Display_1 = Show(calculator1, renderView1)
# trace defaults for the display properties.
calculator1Display_1.AmbientColor = [0.0, 0.0, 0.0]
calculator1Display_1.ColorArrayName = [None, '']
calculator1Display_1.EdgeColor = [0.0, 0.0, 0.0]
calculator1Display_1.GlyphType = 'Arrow'
calculator1Display_1.CubeAxesColor = [0.0, 0.0, 0.0]
calculator1Display_1.SetScaleArray = ['POINTS', 'Field 0']
calculator1Display_1.ScaleTransferFunction = 'PiecewiseFunction'
calculator1Display_1.OpacityArray = ['POINTS', 'Field 0']
calculator1Display_1.OpacityTransferFunction = 'PiecewiseFunction'

# change solid color
calculator1Display_1.DiffuseColor = [0.0, 0.0, 0.0]

# Properties modified on calculator1
calculator1.CoordinateResults = 1

# set up 4 other splines as the boundaries of the map
# set spline coordinates
# need to go from latitude to colatitude
lat_bot = -20
lat_top = 20
lat_mid = (lat_bot + lat_top) / 2.0
lon_min = -20
lon_max = 20
lb_x = 6371000.0 * math.cos(math.radians(lon_min)) * math.sin(math.radians(90-lat_bot))
lb_y = 6371000.0 * math.sin(math.radians(lon_min)) * math.sin(math.radians(90-lat_bot))
lb_z = 6371000.0 * math.cos(math.radians(90-lat_bot))
l_x = 6371000.0 * math.cos(math.radians(lon_min)) * math.sin(math.radians(90-lat_mid))
l_y = 6371000.0 * math.sin(math.radians(lon_min)) * math.sin(math.radians(90-lat_mid))
l_z = 6371000.0 * math.cos(math.radians(90-lat_mid))
b_x = 6371000.0 * math.sin(math.radians(90-lat_bot))
b_y = 0.0
b_z = 6371000.0 * math.cos(math.radians(90-lat_bot))
rb_x = 6371000.0 * math.cos(math.radians(lon_max)) * math.sin(math.radians(90-lat_bot))
rb_y = 6371000.0 * math.sin(math.radians(lon_max)) * math.sin(math.radians(90-lat_bot))
rb_z = 6371000.0 * math.cos(math.radians(90-lat_bot))
r_x = 6371000.0 * math.cos(math.radians(lon_max)) * math.sin(math.radians(90-lat_mid))
r_y = 6371000.0 * math.sin(math.radians(lon_max)) * math.sin(math.radians(90-lat_mid))
r_z = 6371000.0 * math.cos(math.radians(90-lat_mid))
lt_x = 6371000.0 * math.cos(math.radians(lon_min)) * math.sin(math.radians(90-lat_top))
lt_y = 6371000.0 * math.sin(math.radians(lon_min)) * math.sin(math.radians(90-lat_top))
lt_z = 6371000.0 * math.cos(math.radians(90-lat_top))
rt_x = 6371000.0 * math.cos(math.radians(lon_max)) * math.sin(math.radians(90-lat_top))
rt_y = 6371000.0 * math.sin(math.radians(lon_max)) * math.sin(math.radians(90-lat_top))
rt_z = 6371000.0 * math.cos(math.radians(90-lat_top))
t_x = 6371000.0 * math.sin(math.radians(90-lat_top))
t_y = 0.0
t_z = 6371000.0 * math.cos(math.radians(90-lat_top))

# create a new 'SplineSource'
splineSource2 = SplineSource()
splineSource2.ParametricFunction = 'Spline'

# Properties modified on splineSource1.ParametricFunction
splineSource2.ParametricFunction.Points = [lb_x, lb_y, lb_z, b_x, b_y, b_z, rb_x, rb_y, rb_z]

# show data in view
splineSource2Display = Show(splineSource2, renderView1)

# change solid color
splineSource2Display.DiffuseColor = [0.0, 0.0, 0.0]

# Properties modified on splineSource7Display
splineSource2Display.LineWidth = 3.0


# create a new 'SplineSource'
splineSource3 = SplineSource()
splineSource3.ParametricFunction = 'Spline'

# Properties modified on splineSource1.ParametricFunction
splineSource3.ParametricFunction.Points = [rb_x, rb_y, rb_z, r_x, r_y, r_z, rt_x, rt_y, rt_z]

# show data in view
splineSource3Display = Show(splineSource3, renderView1)

# change solid color
splineSource3Display.DiffuseColor = [0.0, 0.0, 0.0]

# Properties modified on splineSource7Display
splineSource3Display.LineWidth = 3.0

# create a new 'SplineSource'
splineSource4 = SplineSource()
splineSource4.ParametricFunction = 'Spline'

# Properties modified on splineSource1.ParametricFunction
splineSource4.ParametricFunction.Points = [rt_x, rt_y, rt_z, t_x, t_y, t_z, lt_x, lt_y, lt_z]

# show data in view
splineSource4Display = Show(splineSource4, renderView1)

# change solid color
splineSource4Display.DiffuseColor = [0.0, 0.0, 0.0]

# Properties modified on splineSource7Display
splineSource4Display.LineWidth = 3.0

# create a new 'SplineSource'
splineSource5 = SplineSource()
splineSource5.ParametricFunction = 'Spline'

# Properties modified on splineSource1.ParametricFunction
splineSource5.ParametricFunction.Points = [lt_x, lt_y, lt_z, l_x, l_y, l_z, lb_x, lb_y, lb_z]

# show data in view
splineSource5Display = Show(splineSource5, renderView1)

# change solid color
splineSource5Display.DiffuseColor = [0.0, 0.0, 0.0]

# Properties modified on splineSource7Display
splineSource5Display.LineWidth = 3.0


# create a new 'Line'
line1 = Line()

# show data in view
line1Display = Show(line1, renderView1)
# trace defaults for the display properties.
line1Display.AmbientColor = [0.0, 0.0, 0.0]
line1Display.ColorArrayName = [None, '']
line1Display.EdgeColor = [0.0, 0.0, 0.0]
line1Display.GlyphType = 'Arrow'
line1Display.CubeAxesColor = [0.0, 0.0, 0.0]
line1Display.SetScaleArray = [None, '']
line1Display.ScaleTransferFunction = 'PiecewiseFunction'
line1Display.OpacityArray = [None, '']
line1Display.OpacityTransferFunction = 'PiecewiseFunction'

# Properties modified on line1
line1.Point1 = [lt_x, lt_y, lt_z]
line1.Point2 = [lt_x, lt_y, lt_z]

# change solid color
line1Display.DiffuseColor = [0.0, 0.0, 0.0]

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=line1)

# Properties modified on line1Display
line1Display.LineWidth = 3.0

# create a new 'Line'
line2 = Line()

# show data in view
line2Display = Show(line2, renderView1)
# trace defaults for the display properties.
line2Display.AmbientColor = [0.0, 0.0, 0.0]
line2Display.ColorArrayName = [None, '']
line2Display.EdgeColor = [0.0, 0.0, 0.0]
line2Display.GlyphType = 'Arrow'
line2Display.CubeAxesColor = [0.0, 0.0, 0.0]
line2Display.SetScaleArray = [None, '']
line2Display.ScaleTransferFunction = 'PiecewiseFunction'
line2Display.OpacityArray = [None, '']
line2Display.OpacityTransferFunction = 'PiecewiseFunction'

# Properties modified on line2
line2.Point1 = [rt_x, rt_y, rt_z]
line2.Point2 = [rt_x, rt_y, rt_z]

# change solid color
line2Display.DiffuseColor = [0.0, 0.0, 0.0]

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=line2)

# Properties modified on line2Display
line2Display.LineWidth = 3.0

# create a new 'Line'
line3 = Line()

# show data in view
line3Display = Show(line3, renderView1)
# trace defaults for the display properties.
line3Display.AmbientColor = [0.0, 0.0, 0.0]
line3Display.ColorArrayName = [None, '']
line3Display.EdgeColor = [0.0, 0.0, 0.0]
line3Display.GlyphType = 'Arrow'
line3Display.CubeAxesColor = [0.0, 0.0, 0.0]
line3Display.SetScaleArray = [None, '']
line3Display.ScaleTransferFunction = 'PiecewiseFunction'
line3Display.OpacityArray = [None, '']
line3Display.OpacityTransferFunction = 'PiecewiseFunction'

# Properties modified on line3
line3.Point1 = [lb_x, lb_y, lb_z]
line3.Point2 = [lb_x, lb_y, lb_z]

# change solid color
line3Display.DiffuseColor = [0.0, 0.0, 0.0]

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=line3)

# Properties modified on line3Display
line3Display.LineWidth = 3.0

# create a new 'Line'
line4 = Line()

# show data in view
line4Display = Show(line4, renderView1)
# trace defaults for the display properties.
line4Display.AmbientColor = [0.0, 0.0, 0.0]
line4Display.ColorArrayName = [None, '']
line4Display.EdgeColor = [0.0, 0.0, 0.0]
line4Display.GlyphType = 'Arrow'
line4Display.CubeAxesColor = [0.0, 0.0, 0.0]
line4Display.SetScaleArray = [None, '']
line4Display.ScaleTransferFunction = 'PiecewiseFunction'
line4Display.OpacityArray = [None, '']
line4Display.OpacityTransferFunction = 'PiecewiseFunction'

# Properties modified on line4
line4.Point1 = [rb_x, rb_y, rb_z]
line4.Point2 = [rb_x, rb_y, rb_z]

# change solid color
line4Display.DiffuseColor = [0.0, 0.0, 0.0]

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=line4)

# Properties modified on line4Display
line4Display.LineWidth = 3.0

SetActiveView(renderView1)
#renderView1.ViewSize = [954, 945]
renderView1.ViewSize = [954, 954]

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)
diffLUTColorBar = GetScalarBar(diffLUT, renderView1)
diffLUT.RescaleTransferFunction(-2.0, 2.0)
slice1Display.LookupTable = diffLUT
#diffLUTColorBar.Position = [0.288, 0.07]
diffLUTColorBar.Position = [0.288, 0.11]

# write label
text3 = Text()
text3.Text = 'dV [%]'
text3Display = Show(text3, renderView1)
text3Display.Color = [0.0, 0.0, 0.0]
#text3Display.Position = [0.45, 0.10]
text3Display.Position = [0.45, 0.14]
text3Display.FontSize = 10

# write label
text2 = Text()
text2.Text = '0 km'
text2Display = Show(text2, renderView1)
text2Display.Color = [0.0, 0.0, 0.0]
#text2Display.Position = [0.44, 0.94]
text2Display.Position = [0.47, 0.90]
text2Display.FontSize = 10


SaveScreenshot('/Users/visu/Glerum/01092016/aspect/EGU/K08_SEMUCB/slices/slice_tomo_slab_sphere_K08_SEMUCB__1.png',layout=viewLayout1, magnification=1, quality=100)

filename_base = "/Users/visu/Glerum/01092016/aspect/EGU/K08_SEMUCB/slices/slice_tomo_slab_sphere_K08_SEMUCB__"
count = 1

# Loop over all the other slices, changing only the necessary settings
for i in xrange(1,27):
    count = count + 1
    radius = 6371000.0 - i * 100000.0
    depth = 6371000.0 - radius

    text2.Text = str(int(depth/1000.0)) + ' km'
 
    # set active view
    SetActiveView(renderView1)
    SetActiveSource(slice1)

    # Properties modified on slice1.SliceType
    slice1.SliceType.Center = [0.0, 0.0, 0.0]
    slice1.SliceType.Radius = radius
    slice2.SliceType.Center = [0.0, 0.0, 0.0]
    slice2.SliceType.Radius = radius

    SetActiveSource(line1)
    scale = radius/6371000.0
    line1.Point2 = [lt_x*scale, lt_y*scale, lt_z*scale]
    SetActiveSource(line2)
    line2.Point2 = [rt_x*scale, rt_y*scale, rt_z*scale]
    SetActiveSource(line3)
    line3.Point2 = [lb_x*scale, lb_y*scale, lb_z*scale]
    SetActiveSource(line4)
    line4.Point2 = [rb_x*scale, rb_y*scale, rb_z*scale]

    SetActiveSource(calculator1)
    calculator1.Function = 'coords/6371000.0*'+str(radius)
  
    SetActiveSource(contour1)
    SetActiveSource(contour4)
    #SetActiveSource(contour5)

    # save screenshot
    filename = filename_base + str(count) + ".png"
    SaveScreenshot(filename, layout=viewLayout1, magnification=1, quality=100)


