#### import the simple module from the paraview
from paraview.simple import *
import math
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Partitioned Unstructured Grid Reader'
solution00000pvtu = XMLPartitionedUnstructuredGridReader(FileName=['/Users/visu/Glerum/01092016/aspect/EGU/K08_P06_CSloc/solution/solution-00000.pvtu'])
solution00000pvtu.CellArrayStatus = []
solution00000pvtu.PointArrayStatus = ['velocity', 'p', 'T', 'C_1', 'C_2', 'C_14', 'strain_rate', 'nonadiabatic_temperature', 'density', 'thermal_expansivity', 'specific_heat', 'viscosity']

# Properties modified on solution00000pvtu
# C_14 is the Aegean slab
solution00000pvtu.PointArrayStatus = ['C_14']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
#renderView1.ViewSize = [1908,945]

# show data in view
solution00000pvtuDisplay = Show(solution00000pvtu, renderView1)
# trace defaults for the display properties.
solution00000pvtuDisplay.CubeAxesVisibility = 0
solution00000pvtuDisplay.Representation = 'Surface'
solution00000pvtuDisplay.AmbientColor = [0.0, 0.0, 0.0]
solution00000pvtuDisplay.ColorArrayName = [None, '']
solution00000pvtuDisplay.DiffuseColor = [1.0, 1.0, 1.0]
solution00000pvtuDisplay.LookupTable = None
solution00000pvtuDisplay.MapScalars = 1
solution00000pvtuDisplay.InterpolateScalarsBeforeMapping = 1
solution00000pvtuDisplay.Opacity = 1.0
solution00000pvtuDisplay.PointSize = 2.0
solution00000pvtuDisplay.LineWidth = 1.0
solution00000pvtuDisplay.Interpolation = 'Gouraud'
solution00000pvtuDisplay.Specular = 0.0
solution00000pvtuDisplay.SpecularColor = [1.0, 1.0, 1.0]
solution00000pvtuDisplay.SpecularPower = 100.0
solution00000pvtuDisplay.Ambient = 0.0
solution00000pvtuDisplay.Diffuse = 1.0
solution00000pvtuDisplay.EdgeColor = [0.0, 0.0, 0.0]
solution00000pvtuDisplay.BackfaceRepresentation = 'Follow Frontface'
solution00000pvtuDisplay.BackfaceAmbientColor = [1.0, 1.0, 1.0]
solution00000pvtuDisplay.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
solution00000pvtuDisplay.BackfaceOpacity = 1.0
solution00000pvtuDisplay.Position = [0.0, 0.0, 0.0]
solution00000pvtuDisplay.Scale = [1.0, 1.0, 1.0]
solution00000pvtuDisplay.Orientation = [0.0, 0.0, 0.0]
solution00000pvtuDisplay.Origin = [0.0, 0.0, 0.0]
solution00000pvtuDisplay.Pickable = 1
solution00000pvtuDisplay.Texture = None
solution00000pvtuDisplay.Triangulate = 0
solution00000pvtuDisplay.NonlinearSubdivisionLevel = 1
solution00000pvtuDisplay.GlyphType = 'Arrow'
solution00000pvtuDisplay.CubeAxesColor = [0.0, 0.0, 0.0]
solution00000pvtuDisplay.CubeAxesCornerOffset = 0.0
solution00000pvtuDisplay.CubeAxesFlyMode = 'Closest Triad'
solution00000pvtuDisplay.CubeAxesInertia = 1
solution00000pvtuDisplay.CubeAxesTickLocation = 'Inside'
solution00000pvtuDisplay.CubeAxesXAxisMinorTickVisibility = 1
solution00000pvtuDisplay.CubeAxesXAxisTickVisibility = 1
solution00000pvtuDisplay.CubeAxesXAxisVisibility = 1
solution00000pvtuDisplay.CubeAxesXGridLines = 0
solution00000pvtuDisplay.CubeAxesXTitle = 'X-Axis'
solution00000pvtuDisplay.CubeAxesUseDefaultXTitle = 1
solution00000pvtuDisplay.CubeAxesYAxisMinorTickVisibility = 1
solution00000pvtuDisplay.CubeAxesYAxisTickVisibility = 1
solution00000pvtuDisplay.CubeAxesYAxisVisibility = 1
solution00000pvtuDisplay.CubeAxesYGridLines = 0
solution00000pvtuDisplay.CubeAxesYTitle = 'Y-Axis'
solution00000pvtuDisplay.CubeAxesUseDefaultYTitle = 1
solution00000pvtuDisplay.CubeAxesZAxisMinorTickVisibility = 1
solution00000pvtuDisplay.CubeAxesZAxisTickVisibility = 1
solution00000pvtuDisplay.CubeAxesZAxisVisibility = 1
solution00000pvtuDisplay.CubeAxesZGridLines = 0
solution00000pvtuDisplay.CubeAxesZTitle = 'Z-Axis'
solution00000pvtuDisplay.CubeAxesUseDefaultZTitle = 1
solution00000pvtuDisplay.CubeAxesGridLineLocation = 'All Faces'
solution00000pvtuDisplay.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
solution00000pvtuDisplay.CustomBoundsActive = [0, 0, 0]
solution00000pvtuDisplay.OriginalBoundsRangeActive = [0, 0, 0]
solution00000pvtuDisplay.CustomRange = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
solution00000pvtuDisplay.CustomRangeActive = [0, 0, 0]
solution00000pvtuDisplay.UseAxesOrigin = 0
solution00000pvtuDisplay.AxesOrigin = [0.0, 0.0, 0.0]
solution00000pvtuDisplay.CubeAxesXLabelFormat = '%-#6.3g'
solution00000pvtuDisplay.CubeAxesYLabelFormat = '%-#6.3g'
solution00000pvtuDisplay.CubeAxesZLabelFormat = '%-#6.3g'
solution00000pvtuDisplay.StickyAxes = 0
solution00000pvtuDisplay.CenterStickyAxes = 0
solution00000pvtuDisplay.SelectionCellLabelBold = 0
solution00000pvtuDisplay.SelectionCellLabelColor = [0.0, 1.0, 0.0]
solution00000pvtuDisplay.SelectionCellLabelFontFamily = 'Arial'
solution00000pvtuDisplay.SelectionCellLabelFontSize = 18
solution00000pvtuDisplay.SelectionCellLabelItalic = 0
solution00000pvtuDisplay.SelectionCellLabelJustification = 'Left'
solution00000pvtuDisplay.SelectionCellLabelOpacity = 1.0
solution00000pvtuDisplay.SelectionCellLabelShadow = 0
solution00000pvtuDisplay.SelectionPointLabelBold = 0
solution00000pvtuDisplay.SelectionPointLabelColor = [1.0, 1.0, 0.0]
solution00000pvtuDisplay.SelectionPointLabelFontFamily = 'Arial'
solution00000pvtuDisplay.SelectionPointLabelFontSize = 18
solution00000pvtuDisplay.SelectionPointLabelItalic = 0
solution00000pvtuDisplay.SelectionPointLabelJustification = 'Left'
solution00000pvtuDisplay.SelectionPointLabelOpacity = 1.0
solution00000pvtuDisplay.SelectionPointLabelShadow = 0
solution00000pvtuDisplay.ScalarOpacityUnitDistance = 51152.954947465136
solution00000pvtuDisplay.SelectMapper = 'Projected tetra'
solution00000pvtuDisplay.GaussianRadius = 0.0
solution00000pvtuDisplay.ShaderPreset = 'Sphere'
solution00000pvtuDisplay.Emissive = 0
solution00000pvtuDisplay.ScaleByArray = 0
solution00000pvtuDisplay.SetScaleArray = ['POINTS', 'T']
solution00000pvtuDisplay.ScaleTransferFunction = 'PiecewiseFunction'
solution00000pvtuDisplay.OpacityByArray = 0
solution00000pvtuDisplay.OpacityArray = ['POINTS', 'T']
solution00000pvtuDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'Arrow' selected for 'GlyphType'
solution00000pvtuDisplay.GlyphType.TipResolution = 6
solution00000pvtuDisplay.GlyphType.TipRadius = 0.1
solution00000pvtuDisplay.GlyphType.TipLength = 0.35
solution00000pvtuDisplay.GlyphType.ShaftResolution = 6
solution00000pvtuDisplay.GlyphType.ShaftRadius = 0.03
solution00000pvtuDisplay.GlyphType.Invert = 0

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
solution00000pvtuDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
solution00000pvtuDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.InteractionMode = '3D'

# create a new 'Slice'
slice1 = Slice(Input=solution00000pvtu)
slice1.SliceType = 'Plane'
slice1.Crinkleslice = 0
slice1.Triangulatetheslice = 1
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [4717985.065516494, 0.0, 0.0]
slice1.SliceType.Normal = [1.0, 0.0, 0.0]
slice1.SliceType.Offset = 0.0

# Properties modified on slice1.SliceType
slice1.SliceType.Origin = [5341292.5, 81227.5811174638, -38203.7658714104]
slice1.SliceType.Normal = [0., 0.843311524056685, -0.59]

# show data in view
slice1Display = Show(slice1, renderView1)
# trace defaults for the display properties.
slice1Display.CubeAxesVisibility = 0
slice1Display.Representation = 'Wireframe'
slice1Display.AmbientColor = [0.0, 0.0, 0.0]
slice1Display.ColorArrayName = [None, '']
slice1Display.DiffuseColor = [1.0, 1.0, 1.0]
slice1Display.LookupTable = None
slice1Display.MapScalars = 1
slice1Display.InterpolateScalarsBeforeMapping = 1
slice1Display.Opacity = 1.0
slice1Display.PointSize = 2.0
slice1Display.LineWidth = 1.0
slice1Display.Interpolation = 'Gouraud'
slice1Display.Specular = 0.0
slice1Display.SpecularColor = [1.0, 1.0, 1.0]
slice1Display.SpecularPower = 100.0
slice1Display.Ambient = 0.0
slice1Display.Diffuse = 1.0
slice1Display.EdgeColor = [0.0, 0.0, 0.0]
slice1Display.BackfaceRepresentation = 'Follow Frontface'
slice1Display.BackfaceAmbientColor = [1.0, 1.0, 1.0]
slice1Display.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
slice1Display.BackfaceOpacity = 1.0
slice1Display.Position = [0.0, 0.0, 0.0]
slice1Display.Scale = [1.0, 1.0, 1.0]
slice1Display.Orientation = [0.0, 0.0, 0.0]
slice1Display.Origin = [0.0, 0.0, 0.0]
slice1Display.Pickable = 1
slice1Display.Texture = None
slice1Display.Triangulate = 0
slice1Display.NonlinearSubdivisionLevel = 1
slice1Display.GlyphType = 'Arrow'
slice1Display.CubeAxesColor = [0.0, 0.0, 0.0]
slice1Display.CubeAxesCornerOffset = 0.0
slice1Display.CubeAxesFlyMode = 'Closest Triad'
slice1Display.CubeAxesInertia = 1
slice1Display.CubeAxesTickLocation = 'Inside'
slice1Display.CubeAxesXAxisMinorTickVisibility = 1
slice1Display.CubeAxesXAxisTickVisibility = 1
slice1Display.CubeAxesXAxisVisibility = 1
slice1Display.CubeAxesXGridLines = 0
slice1Display.CubeAxesXTitle = 'X-Axis'
slice1Display.CubeAxesUseDefaultXTitle = 1
slice1Display.CubeAxesYAxisMinorTickVisibility = 1
slice1Display.CubeAxesYAxisTickVisibility = 1
slice1Display.CubeAxesYAxisVisibility = 1
slice1Display.CubeAxesYGridLines = 0
slice1Display.CubeAxesYTitle = 'Y-Axis'
slice1Display.CubeAxesUseDefaultYTitle = 1
slice1Display.CubeAxesZAxisMinorTickVisibility = 1
slice1Display.CubeAxesZAxisTickVisibility = 1
slice1Display.CubeAxesZAxisVisibility = 1
slice1Display.CubeAxesZGridLines = 0
slice1Display.CubeAxesZTitle = 'Z-Axis'
slice1Display.CubeAxesUseDefaultZTitle = 1
slice1Display.CubeAxesGridLineLocation = 'All Faces'
slice1Display.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
slice1Display.CustomBoundsActive = [0, 0, 0]
slice1Display.OriginalBoundsRangeActive = [0, 0, 0]
slice1Display.CustomRange = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
slice1Display.CustomRangeActive = [0, 0, 0]
slice1Display.UseAxesOrigin = 0
slice1Display.AxesOrigin = [0.0, 0.0, 0.0]
slice1Display.CubeAxesXLabelFormat = '%-#6.3g'
slice1Display.CubeAxesYLabelFormat = '%-#6.3g'
slice1Display.CubeAxesZLabelFormat = '%-#6.3g'
slice1Display.StickyAxes = 0
slice1Display.CenterStickyAxes = 0
slice1Display.SelectionCellLabelBold = 0
slice1Display.SelectionCellLabelColor = [0.0, 1.0, 0.0]
slice1Display.SelectionCellLabelFontFamily = 'Arial'
slice1Display.SelectionCellLabelFontSize = 18
slice1Display.SelectionCellLabelItalic = 0
slice1Display.SelectionCellLabelJustification = 'Left'
slice1Display.SelectionCellLabelOpacity = 1.0
slice1Display.SelectionCellLabelShadow = 0
slice1Display.SelectionPointLabelBold = 0
slice1Display.SelectionPointLabelColor = [1.0, 1.0, 0.0]
slice1Display.SelectionPointLabelFontFamily = 'Arial'
slice1Display.SelectionPointLabelFontSize = 18
slice1Display.SelectionPointLabelItalic = 0
slice1Display.SelectionPointLabelJustification = 'Left'
slice1Display.SelectionPointLabelOpacity = 1.0
slice1Display.SelectionPointLabelShadow = 0
slice1Display.GaussianRadius = 0.0
slice1Display.ShaderPreset = 'Sphere'
slice1Display.Emissive = 0
slice1Display.ScaleByArray = 0
slice1Display.SetScaleArray = ['POINTS', 'T']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityByArray = 0
slice1Display.OpacityArray = ['POINTS', 'T']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'Arrow' selected for 'GlyphType'
slice1Display.GlyphType.TipResolution = 6
slice1Display.GlyphType.TipRadius = 0.1
slice1Display.GlyphType.TipLength = 0.35
slice1Display.GlyphType.ShaftResolution = 6
slice1Display.GlyphType.ShaftRadius = 0.03
slice1Display.GlyphType.Invert = 0

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
slice1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
slice1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# hide data in view
Hide(solution00000pvtu, renderView1)

# set scalar coloring
ColorBy(slice1Display, None)
temperatureLUT = GetColorTransferFunction('C_14')

# rescale color and/or opacity maps used to include current data range
slice1Display.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)
slice1Display.AmbientColor = [0.0, 0.0, 1.0]

contour1 = Contour(Input=slice1)
contour1.ContourBy = ['POINTS', 'C_1']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# show data in view
contour1Display = Show(contour1, renderView1)
# trace defaults for the display properties.
contour1Display.AmbientColor = [0.0, 0.0, 0.0]
contour1Display.ColorArrayName = ['POINTS', 'C_14']
contour1Display.LookupTable = temperatureLUT
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
contour1.ContourBy = ['POINTS', 'C_14']
contour1.Isosurfaces = [0.5]

# turn off scalar coloring
ColorBy(contour1Display, None)

# change solid color
contour1Display.DiffuseColor = [1.0, 0.0, 0.0]
contour1Display.LineWidth = 3.0

scale_410 = 5961000./6371000.

# set spline coordinates based on the points with max and min z coordinates (predetermined)
min_x = 5.8276e6 * scale_410
min_y = -1.4018e6 * scale_410
min_z = -2.1579e6 * scale_410
mid_x = 6.3696e6 * scale_410
mid_y = 1.0796e5 * scale_410
mid_z = 0.0
max_x = 5.7714e6 * scale_410
max_y = 1.6177e6 * scale_410
max_z = 2.1579e6 * scale_410

# create a new 'SplineSource' along 410
splineSource6 = SplineSource()
splineSource6.ParametricFunction = 'Spline'

# Properties modified on splineSource6.ParametricFunction
splineSource6.ParametricFunction.Points = [min_x, min_y, min_z, mid_x, mid_y, mid_z, max_x, max_y, max_z]

# show data in view
splineSource6Display = Show(splineSource6, renderView1)
# trace defaults for the display properties.
splineSource6Display.AmbientColor = [0.0, 0.0, 0.0]
splineSource6Display.ColorArrayName = [None, '']
splineSource6Display.EdgeColor = [0.0, 0.0, 0.0]
splineSource6Display.GlyphType = 'Arrow'
splineSource6Display.CubeAxesColor = [0.0, 0.0, 0.0]
splineSource6Display.SetScaleArray = [None, '']
splineSource6Display.ScaleTransferFunction = 'PiecewiseFunction'
splineSource6Display.OpacityArray = [None, '']
splineSource6Display.OpacityTransferFunction = 'PiecewiseFunction'

# change solid color
splineSource6Display.DiffuseColor = [0.0, 0.0, 0.0]

# Properties modified on splineSource6Display
splineSource6Display.LineWidth = 3.0

# set spline coordinates
min_x = 5711000.0 / 5961000.0 * min_x
min_y = 5711000.0 / 5961000.0 * min_y
min_z = 5711000.0 / 5961000.0 * min_z
mid_x = 5711000.0 / 5961000.0 * mid_x
mid_y = 5711000.0 / 5961000.0 * mid_y
mid_z = 5711000.0 / 5961000.0 * mid_z
max_x = 5711000.0 / 5961000.0 * max_x
max_y = 5711000.0 / 5961000.0 * max_y
max_z = 5711000.0 / 5961000.0 * max_z

# create a new 'SplineSource' along 660
splineSource7 = SplineSource()
splineSource7.ParametricFunction = 'Spline'

# Properties modified on splineSource7.ParametricFunction
splineSource7.ParametricFunction.Points = [min_x, min_y, min_z, mid_x, mid_y, mid_z, max_x, max_y, max_z]

# show data in view
splineSource7Display = Show(splineSource7, renderView1)
# trace defaults for the display properties.
splineSource7Display.AmbientColor = [0.0, 0.0, 0.0]
splineSource7Display.ColorArrayName = [None, '']
splineSource7Display.EdgeColor = [0.0, 0.0, 0.0]
splineSource7Display.GlyphType = 'Arrow'
splineSource7Display.CubeAxesColor = [0.0, 0.0, 0.0]
splineSource7Display.SetScaleArray = [None, '']
splineSource7Display.ScaleTransferFunction = 'PiecewiseFunction'
splineSource7Display.OpacityArray = [None, '']
splineSource7Display.OpacityTransferFunction = 'PiecewiseFunction'

# change solid color
splineSource7Display.DiffuseColor = [0.0, 0.0, 0.0]

# Properties modified on splineSource7Display
splineSource7Display.LineWidth = 3.0

# set active view
SetActiveView(renderView1)

# set camera
# distance from slice
far = 4.0e6
origin = slice1.SliceType.Origin
# focal point is center point of slice
scale_slice = 6500
center = [origin[0] / 6371 * scale_slice, origin[1] / 6371 * scale_slice, origin[2] / 6371 * scale_slice ]
renderView1.CameraFocalPoint = [ center[0], center[1], center[2] ]
# position = focal point - X * n_norm
normal = slice1.SliceType.Normal
norm = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]) ** 0.5
normal_norm = [ normal[0] / norm, normal[1] / norm, normal[2] / norm ]
position = [center[0] - far * normal_norm[0], center[1] - far * normal_norm[1], center[2] - far * normal_norm[2]]
renderView1.CameraPosition = [ position[0], position[1], position[2] ]
# View up is the radial vector
renderView1.CameraViewUp = [origin[0], origin[1], origin[2] ]

# ACG set background
renderView1.Background = [1,1,1] #White

# get Layout
viewLayout1 = GetLayout()

# split cell
##viewLayout1.SplitVertical(0, 0.730)
##viewLayout1.SplitVertical(0, 0.5)
#
## set active view
#SetActiveView(None)
#
## get active view
#renderView2 = GetActiveViewOrCreate('RenderView')
##renderView2.ViewSize = [954, 945]
##renderView2.ViewSize = [954, 283]
#renderView2.ViewSize = [954, 473]
#renderView2.AxesGrid = 'GridAxes3DActor'
#renderView2.OrientationAxesOutlineColor = [0.0, 0.0, 0.0]
#renderView2.StereoType = 0
#renderView2.Background = [1.0, 1.0, 1.0]
#
## init the 'GridAxes3DActor'
#renderView2.AxesGrid.GridColor = [0.0, 0.0, 0.0]
#
## place view in the layout
#viewLayout1.AssignView(2, renderView2)
#
## swap view locations (map on top)
#viewLayout1.SwapCells(2, 1)
#
## Add map with indicator line
## create a new 'CSV Reader'
#map_lon_latcsv = CSVReader(FileName=['/Users/visu/Glerum/090915/aspect/T_perturb/other_visu_data/map_lon_lat.csv'])
#
## Properties modified on map_lon_latcsv
#map_lon_latcsv.HaveHeaders = 0
#
## Create a new 'SpreadSheet View'
#spreadSheetView1 = CreateView('SpreadSheetView')
#spreadSheetView1.BlockSize = 1024L
#
## show data in view
#map_lon_latcsvDisplay = Show(map_lon_latcsv, spreadSheetView1)
## trace defaults for the display properties.
#map_lon_latcsvDisplay.FieldAssociation = 'Row Data'
#
## create a new 'Table To Points'
#tableToPoints1 = TableToPoints(Input=map_lon_latcsv)
#tableToPoints1.XColumn = 'Field 0'
#tableToPoints1.YColumn = 'Field 0'
#tableToPoints1.ZColumn = 'Field 0'
#
## Properties modified on tableToPoints1
#tableToPoints1.YColumn = 'Field 1'
#tableToPoints1.ZColumn = 'Field 2'
#
## show data in view
#tableToPoints1Display = Show(tableToPoints1, spreadSheetView1)
#
## hide data in view
#Hide(map_lon_latcsv, spreadSheetView1)
#
## destroy spreadSheetView1
#Delete(spreadSheetView1)
#del spreadSheetView1
#
## set active view
#SetActiveView(renderView2)
#
## set active source
#SetActiveSource(tableToPoints1)
#
## show data in view
#tableToPoints1Display = Show(tableToPoints1, renderView2)
## trace defaults for the display properties.
#tableToPoints1Display.AmbientColor = [0.0, 0.0, 0.0]
#tableToPoints1Display.ColorArrayName = [None, '']
#tableToPoints1Display.EdgeColor = [0.0, 0.0, 0.0]
#tableToPoints1Display.GlyphType = 'Arrow'
#tableToPoints1Display.CubeAxesColor = [0.0, 0.0, 0.0]
#tableToPoints1Display.SetScaleArray = [None, '']
#tableToPoints1Display.ScaleTransferFunction = 'PiecewiseFunction'
#tableToPoints1Display.OpacityArray = [None, '']
#tableToPoints1Display.OpacityTransferFunction = 'PiecewiseFunction'
#
## change solid color
#tableToPoints1Display.DiffuseColor = [0.0, 0.0, 0.0]
#
## create a new 'SplineSource'
#splineSource1 = SplineSource()
#splineSource1.ParametricFunction = 'Spline'
#
## set spline coordinates
#min_x = 6371000.0 / 5711000.0 * min_x
#min_y = 6371000.0 / 5711000.0 * min_y
#min_z = 6371000.0 / 5711000.0 * min_z
#mid_x = 6371000.0 / 5711000.0 * mid_x
#mid_y = 6371000.0 / 5711000.0 * mid_y
#mid_z = 6371000.0 / 5711000.0 * mid_z
#max_x = 6371000.0 / 5711000.0 * max_x
#max_y = 6371000.0 / 5711000.0 * max_y
#max_z = 6371000.0 / 5711000.0 * max_z
#
#
## Properties modified on splineSource1.ParametricFunction
#splineSource1.ParametricFunction.Points = [min_x, min_y, min_z, mid_x, mid_y, mid_z, max_x, max_y, max_z]
#
## show data in view
#splineSource1Display = Show(splineSource1, renderView2)
## trace defaults for the display properties.
#splineSource1Display.AmbientColor = [0.0, 0.0, 0.0]
#splineSource1Display.ColorArrayName = [None, '']
#splineSource1Display.EdgeColor = [0.0, 0.0, 0.0]
#splineSource1Display.GlyphType = 'Arrow'
#splineSource1Display.CubeAxesColor = [0.0, 0.0, 0.0]
#splineSource1Display.SetScaleArray = [None, '']
#splineSource1Display.ScaleTransferFunction = 'PiecewiseFunction'
#splineSource1Display.OpacityArray = [None, '']
#splineSource1Display.OpacityTransferFunction = 'PiecewiseFunction'
#
## change solid color
#splineSource1Display.DiffuseColor = [0.0, 0.0, 1.0]
#
## Properties modified on splineSource1Display
#splineSource1Display.LineWidth = 3.0
#
#
## set up 4 other splines as the boundaries of the map
## set spline coordinates
## need to go from latitude to colatitude
##lon_min = max(-20-5,-20)
##lon_max = min(-20+5,20)
#lon_min = -14
#lon_max = 16
##lat_mid = (lat_bot + lat_top) / 2.0
#lon_mid = 1
#lat_bot = -20
#lat_mid = 0
#lat_top = 20
#lb_x = 6371000.0 * math.cos(math.radians(lon_min)) * math.sin(math.radians(90-lat_bot))
#lb_y = 6371000.0 * math.sin(math.radians(lon_min)) * math.sin(math.radians(90-lat_bot))
#lb_z = 6371000.0 * math.cos(math.radians(90-lat_bot))
#l_x = 6371000.0 * math.cos(math.radians(lon_min)) * math.sin(math.radians(90-lat_mid))
#l_y = 6371000.0 * math.sin(math.radians(lon_min)) * math.sin(math.radians(90-lat_mid))
#l_z = 6371000.0 * math.cos(math.radians(90-lat_mid))
#b_x = 6371000.0 * math.cos(math.radians(lon_mid)) * math.sin(math.radians(90-lat_bot))
#b_y = 6371000.0 * math.sin(math.radians(lon_mid)) * math.sin(math.radians(90-lat_bot))
#b_z = 6371000.0 * math.cos(math.radians(90-lat_bot))
#rb_x = 6371000.0 * math.cos(math.radians(lon_max)) * math.sin(math.radians(90-lat_bot))
#rb_y = 6371000.0 * math.sin(math.radians(lon_max)) * math.sin(math.radians(90-lat_bot))
#rb_z = 6371000.0 * math.cos(math.radians(90-lat_bot))
#r_x = 6371000.0 * math.cos(math.radians(lon_max)) * math.sin(math.radians(90-lat_mid))
#r_y = 6371000.0 * math.sin(math.radians(lon_max)) * math.sin(math.radians(90-lat_mid))
#r_z = 6371000.0 * math.cos(math.radians(90-lat_mid))
#lt_x = 6371000.0 * math.cos(math.radians(lon_min)) * math.sin(math.radians(90-lat_top))
#lt_y = 6371000.0 * math.sin(math.radians(lon_min)) * math.sin(math.radians(90-lat_top))
#lt_z = 6371000.0 * math.cos(math.radians(90-lat_top))
#rt_x = 6371000.0 * math.cos(math.radians(lon_max)) * math.sin(math.radians(90-lat_top))
#rt_y = 6371000.0 * math.sin(math.radians(lon_max)) * math.sin(math.radians(90-lat_top))
#rt_z = 6371000.0 * math.cos(math.radians(90-lat_top))
#t_x = 6371000.0 * math.cos(math.radians(lon_mid)) * math.sin(math.radians(90-lat_top))
#t_y = 6371000.0 * math.sin(math.radians(lon_mid)) * math.sin(math.radians(90-lat_top))
#t_z = 6371000.0 * math.cos(math.radians(90-lat_top))
#
## create a new 'SplineSource'
#splineSource2 = SplineSource()
#splineSource2.ParametricFunction = 'Spline'
#
## Properties modified on splineSource1.ParametricFunction
#splineSource2.ParametricFunction.Points = [lb_x, lb_y, lb_z, b_x, b_y, b_z, rb_x, rb_y, rb_z]
#
## show data in view
#splineSource2Display = Show(splineSource2, renderView2)
#
## change solid color
#splineSource2Display.DiffuseColor = [0.0, 0.0, 0.0]
#
## Properties modified on splineSource7Display
#splineSource2Display.LineWidth = 3.0
#
#
## create a new 'SplineSource'
#splineSource3 = SplineSource()
#splineSource3.ParametricFunction = 'Spline'
#
## Properties modified on splineSource1.ParametricFunction
#splineSource3.ParametricFunction.Points = [rb_x, rb_y, rb_z, r_x, r_y, r_z, rt_x, rt_y, rt_z]
#
## show data in view
#splineSource3Display = Show(splineSource3, renderView2)
#
## change solid color
#splineSource3Display.DiffuseColor = [0.0, 0.0, 0.0]
#
## Properties modified on splineSource7Display
#splineSource3Display.LineWidth = 3.0
#
## create a new 'SplineSource'
#splineSource4 = SplineSource()
#splineSource4.ParametricFunction = 'Spline'
#
## Properties modified on splineSource1.ParametricFunction
#splineSource4.ParametricFunction.Points = [rt_x, rt_y, rt_z, t_x, t_y, t_z, lt_x, lt_y, lt_z]
#
## show data in view
#splineSource4Display = Show(splineSource4, renderView2)
#
## change solid color
#splineSource4Display.DiffuseColor = [0.0, 0.0, 0.0]
#
## Properties modified on splineSource7Display
#splineSource4Display.LineWidth = 3.0
#
## create a new 'SplineSource'
#splineSource5 = SplineSource()
#splineSource5.ParametricFunction = 'Spline'
#
## Properties modified on splineSource1.ParametricFunction
#splineSource5.ParametricFunction.Points = [lt_x, lt_y, lt_z, l_x, l_y, l_z, lb_x, lb_y, lb_z]
#
## show data in view
#splineSource5Display = Show(splineSource5, renderView2)
#
## change solid color
#splineSource5Display.DiffuseColor = [0.0, 0.0, 0.0]
#
## Properties modified on splineSource7Display
#splineSource5Display.LineWidth = 3.0
#
## create a new 'Clip'
#clip1 = Clip(Input=tableToPoints1)
#clip1.ClipType = 'Plane'
#clip1.Scalars = [None, '']
#
## Properties modified on clip1.ClipType
#clip1.ClipType.Origin = [l_x, l_y, l_z]
#clip1.ClipType.Normal = [-l_y, l_x, t_z]
#clip1.InsideOut = 0
#
## show data in view
#clip1Display = Show(clip1, renderView2)
#
## Properties modified on clip1Display
#clip1Display.DiffuseColor = [0.0, 0.0, 0.0]
#
## hide data in view
#Hide(tableToPoints1, renderView2)
#
## create a new 'Clip'
#clip2 = Clip(Input=clip1)
#clip2.ClipType = 'Plane'
#clip2.Scalars = [None, '']
#
## Properties modified on clip2.ClipType
#clip2.ClipType.Origin = [r_x, r_y, r_z]
#clip2.ClipType.Normal = [-r_y, r_x, r_z]
#clip2.InsideOut = 1
#
## show data in view
#clip2Display = Show(clip2, renderView2)
#
## Properties modified on clip2Display
#clip2Display.DiffuseColor = [0.0, 0.0, 0.0]
#
## hide data in view
#Hide(clip1, renderView2)
#
# current camera placement for renderView2
#renderView2.CameraFocalPoint = [mid_x, mid_y, mid_z]
#renderView2.CameraViewUp = [-mid_y, mid_x, mid_z]
#far2 = 6.0e6
#radial = [mid_x, mid_y, mid_z]
#norm_r = (radial[0] * radial[0] + radial[1] * radial[1] + radial[2] * radial[2]) ** 0.5
#radial_norm = [ radial[0] / norm_r, radial[1] / norm_r, radial[2] / norm_r ]
#position = [radial[0] + far2 * radial_norm[0], radial[1] + far2 * radial_norm[1], radial[2] + far2 * radial_norm[2]]
#renderView2.CameraPosition = [ position[0], position[1], position[2] ]

SetActiveView(renderView1)
renderView1.ViewSize = [1054, 473]

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, False)
slice1Display.LookupTable = temperatureLUT

# write label
text3 = Text()
text3.Text = 'FE Mesh'
text3Display = Show(text3, renderView1)
text3Display.Color = [0.0, 0.0, 0.0]
text3Display.Position = [0.42, 0.13]
text3Display.FontSize = 10

filename_base = "/Users/visu/Glerum/01092016/aspect/EGU/K08_P06_CSloc/slices/slice_mesh_K08_P06_CSloc__"
count = 1
filename = filename_base + str(count) + ".png"

# save screenshot
SaveScreenshot(filename,layout=viewLayout1, magnification=2, quality=100)
