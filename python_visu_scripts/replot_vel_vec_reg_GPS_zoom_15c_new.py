#### import the simple module from the paraview
from paraview.simple import *
import math
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Partitioned Unstructured Grid Reader'
solution00000pvtu = XMLPartitionedUnstructuredGridReader(FileName=['/Users/visu/Glerum/01092016/aspect/EGU/K08_P06_CSloc_new/solution/solution-00000.pvtu'])
solution00000pvtu.PointArrayStatus = ['velocity', 'p', 'T', 'C_1', 'C_2', 'C_3', 'C_4', 'C_5', 'C_6', 'C_7', 'C_8', 'C_9', 'C_10', 'C_11', 'C_12', 'C_13', 'C_14', 'C_15', 'C_16', 'C_17', 'C_18', 'strain_rate', 'density', 'thermal_expansivity', 'specific_heat', 'viscosity', 'viscosity_ratio', 'nonadiabatic_temperature', 'stress_xx', 'stress_yy', 'stress_zz', 'stress_xy', 'stress_xz', 'stress_yz', 'trace_strain_rate', 'dynamic_topography']

# Properties modified on solution00000pvtu
solution00000pvtu.PointArrayStatus = ['velocity','C_14','C_12','C_10','C_8','C_6', 'C_1', 'C_3']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1028, 638]

# Add GPS
# create a new 'Legacy VTK Reader'
points_vel_cart_newvtk = LegacyVTKReader(FileNames=['/Users/visu/Glerum/01092016/aspect/modelK/other_visu_data/pred_vel_GPS_Nocq_Doubr.vtk'])
# create a new 'Calculator'
calculator4 = Calculator(Input=points_vel_cart_newvtk)
calculator4.Function = ''
 
# Properties modified on calculator1
calculator4.CoordinateResults = 1
calculator4.Function = 'iHat*6366/6371*coordsX+jHat*6366/6371*coordsY+kHat*6366/6371*coordsZ'
 
Hide(points_vel_cart_newvtk, renderView1)



# create a new 'Glyph'
glyph2 = Glyph(Input=calculator4,
    GlyphType='Arrow')
glyph2.Scalars = ['POINTS', 'vel_mag']
glyph2.Vectors = ['POINTS', 'velocity']
glyph2.ScaleFactor = 1134255.0
glyph2.GlyphTransform = 'Transform2'

# Properties modified on glyph2
glyph2.ScaleMode = 'vector'
glyph2.ScaleFactor = 7000000.0

# show data in view
glyph2Display = Show(glyph2, renderView1)
# trace defaults for the display properties.
glyph2Display.AmbientColor = [0.0, 0.0, 1.0]
glyph2Display.ColorArrayName = [None, '']
glyph2Display.EdgeColor = [0.0, 0.0, 1.0]
glyph2Display.GlyphType = 'Arrow'
glyph2Display.CubeAxesColor = [0.0, 0.0, 1.0]
glyph2Display.SetScaleArray = [None, '']
glyph2Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph2Display.OpacityArray = [None, '']
glyph2Display.OpacityTransferFunction = 'PiecewiseFunction'

# hide data in view
Hide(points_vel_cart_newvtk, renderView1)
# Properties modified on glyph2
glyph2.GlyphMode = 'Every Nth Point'

# Properties modified on glyph2Display
glyph2Display.SetScaleArray = [None, 'vel_mag']

# Properties modified on glyph2Display
glyph2Display.OpacityArray = [None, 'vel_mag']

# Properties modified on glyph2
#glyph2.Stride = 5
glyph2.Stride = 2
#glyph2Display.SetRepresentationType('Wireframe')

# change solid color
#glyph2Display.DiffuseColor = [1.0, 0.0, 1.0]
#glyph2Display.AmbientColor = [1.0, 0.0, 1.0]
#glyph2Display.DiffuseColor = [0.545, 0.902, 0.706]
#glyph2Display.AmbientColor = [0.545, 0.902, 0.706]
glyph2Display.DiffuseColor = [0.369, 0.153, 0.157] #Brown
glyph2Display.AmbientColor = [0.369, 0.153, 0.157]
glyph2Display.AmbientColor = [0.545, 0.902, 0.706] #Green
glyph2Display.DiffuseColor = [0.545, 0.902, 0.706]

# show data in view
solution00000pvtuDisplay = Show(solution00000pvtu, renderView1)
# trace defaults for the display properties.
solution00000pvtuDisplay.AmbientColor = [0.0, 0.0, 1.0]
solution00000pvtuDisplay.ColorArrayName = [None, '']
solution00000pvtuDisplay.EdgeColor = [0.0, 0.0, 1.0]
solution00000pvtuDisplay.GlyphType = 'Arrow'
solution00000pvtuDisplay.CubeAxesColor = [0.0, 0.0, 1.0]
solution00000pvtuDisplay.ScalarOpacityUnitDistance = 52856.9360269387
solution00000pvtuDisplay.SetScaleArray = [None, '']
solution00000pvtuDisplay.ScaleTransferFunction = 'PiecewiseFunction'
solution00000pvtuDisplay.OpacityArray = [None, '']
solution00000pvtuDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# reset view to fit data
renderView1.ResetCamera()

# create a new 'Clip'
clip1 = Clip(Input=solution00000pvtu)
clip1.ClipType = 'Plane'
clip1.Scalars = ['FIELD', 'TIME']

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Origin = [4721612.875, 0.0, 0.0]

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip1)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=clip1)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip1)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=clip1)

# Properties modified on clip1
clip1.ClipType = 'Sphere'

# Properties modified on clip1.ClipType
clip1.ClipType.Center = [0.0, 0.0, 0.0]
clip1.ClipType.Radius = 6351000.0

# show data in view
clip1Display = Show(clip1, renderView1)
# trace defaults for the display properties.
clip1Display.AmbientColor = [0.0, 0.0, 1.0]
clip1Display.ColorArrayName = [None, '']
clip1Display.EdgeColor = [0.0, 0.0, 1.0]
clip1Display.GlyphType = 'Arrow'
clip1Display.CubeAxesColor = [0.0, 0.0, 1.0]
clip1Display.ScalarOpacityUnitDistance = 57752.9367665238
clip1Display.SetScaleArray = [None, '']
clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
clip1Display.OpacityArray = [None, '']
clip1Display.OpacityTransferFunction = 'PiecewiseFunction'

# hide data in view
Hide(solution00000pvtu, renderView1)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip1)

# hide data in view
Hide(clip1, renderView1)

# set active source
SetActiveSource(solution00000pvtu)

# create a new 'Slice'
slice1 = Slice(Input=solution00000pvtu)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [4721612.875, 0.0, 0.0]

# Properties modified on slice1
slice1.SliceType = 'Sphere'

# Properties modified on slice1.SliceType
slice1.SliceType.Center = [0.0, 0.0, 0.0]
slice1.SliceType.Radius = 6366000.0

# show data in view
slice1Display = Show(slice1, renderView1)
# trace defaults for the display properties.
slice1Display.AmbientColor = [0.0, 0.0, 1.0]
slice1Display.ColorArrayName = [None, '']
slice1Display.EdgeColor = [0.0, 0.0, 1.0]
slice1Display.GlyphType = 'Arrow'
slice1Display.CubeAxesColor = [0.0, 0.0, 1.0]
slice1Display.SetScaleArray = [None, '']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = [None, '']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'

# hide data in view
Hide(solution00000pvtu, renderView1)

# set scalar coloring
ColorBy(slice1Display, ('POINTS', 'velocity'))

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'velocity'
velocityLUT = GetColorTransferFunction('velocity')

# get opacity transfer function/opacity map for 'velocity'
velocityPWF = GetOpacityTransferFunction('velocity')

# Rescale transfer function
velocityLUT.RescaleTransferFunction(0.0, 0.03)

# Rescale transfer function
velocityPWF.RescaleTransferFunction(0.0, 0.03)

# get color legend
velocityLUTColorBar = GetScalarBar(velocityLUT, renderView1)

## Properties modified on velocityLUTColorBar
velocityLUTColorBar.Title = ''
velocityLUTColorBar.ComponentTitle = ''
velocityLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
velocityLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
velocityLUTColorBar.AutomaticLabelFormat = 0
velocityLUTColorBar.LabelFormat = '%-#6.2f'
velocityLUTColorBar.RangeLabelFormat = '%6.2f'
velocityLUTColorBar.AutoOrient = 0
velocityLUTColorBar.Orientation = 'Horizontal'
velocityLUTColorBar.Position = [100, 0.85]
velocityLUTColorBar.DrawTickMarks = 0
velocityLUTColorBar.DrawTickLabels = 0
velocityLUTColorBar.AddRangeLabels = 0
velocityLUTColorBar.TitleFontSize = 10
velocityLUTColorBar.LabelFontSize = 10

velocityLUT.Annotations = ['0.00', '0', '0.01', '1.0', '0.02', '2.0', '0.03', '3.0']

# get opacity transfer function/opacity map for 'velocity'
velocityPWF = GetOpacityTransferFunction('velocity')
velocityPWF.Points = [0.0, 0.699999988079071, 0.5, 0.0, 0.1, 1.0, 0.5, 0.0]
velocityPWF.AllowDuplicateScalars = 1
velocityPWF.ScalarRangeInitialized = 1

# create a new 'Calculator'
calculator1 = Calculator(Input=slice1)
calculator1.Function = ''

# Properties modified on calculator1
calculator1.ResultArrayName = 'Vel mag'
calculator1.Function = 'sqrt(velocity_X*velocity_X+velocity_Y*velocity_Y+velocity_Z*velocity_Z)'

# create a new 'Contour'
contour1 = Contour(Input=calculator1)
contour1.ContourBy = ['POINTS', 'Vel mag']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# show data in view
contour1Display = Show(contour1, renderView1)
# trace defaults for the display properties.
contour1Display.AmbientColor = [0.0, 0.0, 0.0]
contour1Display.ColorArrayName = ['POINTS', 'diff']
contour1Display.LookupTable = velocityLUT
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
contour1.ContourBy = ['POINTS', 'Vel mag']
contour1.Isosurfaces = [0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28]

# turn off scalar coloring
ColorBy(contour1Display, None)

# change solid color
contour1Display.DiffuseColor = [0.0, 0.0, 0.0]

# create a new 'Contour'
contour2 = Contour(Input=calculator1)
contour2.ContourBy = ['POINTS', 'Vel mag']
contour2.Isosurfaces = [0.5]
contour2.PointMergeMethod = 'Uniform Binning'

# show data in view
contour2Display = Show(contour2, renderView1)
# trace defaults for the display properties.
contour2Display.AmbientColor = [0.0, 0.0, 0.0]
contour2Display.ColorArrayName = ['POINTS', 'C_15']
contour2Display.LookupTable = velocityLUT
contour2Display.EdgeColor = [0.0, 0.0, 0.0]
contour2Display.GlyphType = 'Arrow'
contour2Display.CubeAxesColor = [0.0, 0.0, 0.0]
contour2Display.SetScaleArray = ['POINTS', 'C_15']
contour2Display.ScaleTransferFunction = 'PiecewiseFunction'
contour2Display.OpacityArray = ['POINTS', 'C_15']
contour2Display.OpacityTransferFunction = 'PiecewiseFunction'

# show color bar/color legend
contour2Display.SetScalarBarVisibility(renderView1, False)
# Properties modified on contour2
contour2.ContourBy = ['POINTS', 'C_6']
contour2.Isosurfaces = [0.5]

# turn off scalar coloring
ColorBy(contour2Display, None)

# change solid color
contour2Display.DiffuseColor = [1.0, 1.0, 1.0]
contour2Display.LineWidth = 2.0

Hide(contour2,renderView1)

# create a new 'Contour'
contour3 = Contour(Input=calculator1)
contour3.ContourBy = ['POINTS', 'Vel mag']
contour3.Isosurfaces = [0.5]
contour3.PointMergeMethod = 'Uniform Binning'

# show data in view
contour3Display = Show(contour3, renderView1)
# trace defaults for the display properties.
contour3Display.AmbientColor = [0.0, 0.0, 0.0]
contour3Display.ColorArrayName = ['POINTS', 'C_15']
contour3Display.LookupTable = velocityLUT
contour3Display.EdgeColor = [0.0, 0.0, 0.0]
contour3Display.GlyphType = 'Arrow'
contour3Display.CubeAxesColor = [0.0, 0.0, 0.0]
contour3Display.SetScaleArray = ['POINTS', 'C_15']
contour3Display.ScaleTransferFunction = 'PiecewiseFunction'
contour3Display.OpacityArray = ['POINTS', 'C_15']
contour3Display.OpacityTransferFunction = 'PiecewiseFunction'

# show color bar/color legend
contour3Display.SetScalarBarVisibility(renderView1, False)
# Properties modified on contour3
contour3.ContourBy = ['POINTS', 'C_12']
contour3.Isosurfaces = [0.5]

# turn off scalar coloring
ColorBy(contour3Display, None)

# change solid color
contour3Display.DiffuseColor = [1.0, 1.0, 1.0]
contour3Display.LineWidth = 2.0

Hide(contour3,renderView1)

# create a new 'Contour'
contour4 = Contour(Input=calculator1)
contour4.ContourBy = ['POINTS', 'Vel mag']
contour4.Isosurfaces = [0.5]
contour4.PointMergeMethod = 'Uniform Binning'

# show data in view
contour4Display = Show(contour4, renderView1)
# trace defaults for the display properties.
contour4Display.AmbientColor = [0.0, 0.0, 0.0]
contour4Display.ColorArrayName = ['POINTS', 'C_15']
contour4Display.LookupTable = velocityLUT
contour4Display.EdgeColor = [0.0, 0.0, 0.0]
contour4Display.GlyphType = 'Arrow'
contour4Display.CubeAxesColor = [0.0, 0.0, 0.0]
contour4Display.SetScaleArray = ['POINTS', 'C_15']
contour4Display.ScaleTransferFunction = 'PiecewiseFunction'
contour4Display.OpacityArray = ['POINTS', 'C_15']
contour4Display.OpacityTransferFunction = 'PiecewiseFunction'

# show color bar/color legend
contour4Display.SetScalarBarVisibility(renderView1, False)
# Properties modified on contour4
contour4.ContourBy = ['POINTS', 'C_10']
contour4.Isosurfaces = [0.5]

# turn off scalar coloring
ColorBy(contour4Display, None)

# change solid color
contour4Display.DiffuseColor = [1.0, 1.0, 1.0]
contour4Display.LineWidth = 2.0

Hide(contour4,renderView1)

# create a new 'Contour'
contour5 = Contour(Input=calculator1)
contour5.ContourBy = ['POINTS', 'Vel mag']
contour5.Isosurfaces = [0.5]
contour5.PointMergeMethod = 'Uniform Binning'

# show data in view
contour5Display = Show(contour5, renderView1)
# trace defaults for the display properties.
contour5Display.AmbientColor = [0.0, 0.0, 0.0]
contour5Display.ColorArrayName = ['POINTS', 'C_15']
contour5Display.LookupTable = velocityLUT
contour5Display.EdgeColor = [0.0, 0.0, 0.0]
contour5Display.GlyphType = 'Arrow'
contour5Display.CubeAxesColor = [0.0, 0.0, 0.0]
contour5Display.SetScaleArray = ['POINTS', 'C_15']
contour5Display.ScaleTransferFunction = 'PiecewiseFunction'
contour5Display.OpacityArray = ['POINTS', 'C_15']
contour5Display.OpacityTransferFunction = 'PiecewiseFunction'

# show color bar/color legend
contour5Display.SetScalarBarVisibility(renderView1, False)
# Properties modified on contour5
contour5.ContourBy = ['POINTS', 'C_8']
contour5.Isosurfaces = [0.5]

# turn off scalar coloring
ColorBy(contour5Display, None)

# change solid color
contour5Display.DiffuseColor = [1.0, 1.0, 1.0]
contour5Display.LineWidth = 2.0

Hide(contour5,renderView1)

calculator5 = Calculator(Input=calculator1)
calculator5.Function = ''

# Properties modified on calculator5
# Scale from 5mm to 2cm per year
calculator5.Function = 'C_12+C_10+C_8+C_6'

# show data in view
#calculator5Display = Show(calculator5, spreadSheetView1)

# create a new 'Contour'
contour6 = Contour(Input=calculator5)
contour6.ContourBy = ['POINTS', 'Vel mag']
contour6.Isosurfaces = [0.5]
contour6.PointMergeMethod = 'Uniform Binning'

# show data in view
contour6Display = Show(contour6, renderView1)
# trace defaults for the display properties.
contour6Display.AmbientColor = [0.0, 0.0, 0.0]
contour6Display.ColorArrayName = ['POINTS', 'C_15']
contour6Display.LookupTable = velocityLUT
contour6Display.EdgeColor = [0.0, 0.0, 0.0]
contour6Display.GlyphType = 'Arrow'
contour6Display.CubeAxesColor = [0.0, 0.0, 0.0]
contour6Display.SetScaleArray = ['POINTS', 'C_15']
contour6Display.ScaleTransferFunction = 'PiecewiseFunction'
contour6Display.OpacityArray = ['POINTS', 'C_15']
contour6Display.OpacityTransferFunction = 'PiecewiseFunction'

# show color bar/color legend
contour6Display.SetScalarBarVisibility(renderView1, False)
# Properties modified on contour6
contour6.ContourBy = ['POINTS', 'Result']
contour6.Isosurfaces = [0.5]

# turn off scalar coloring
ColorBy(contour6Display, None)

# change solid color
contour6Display.DiffuseColor = [1.0, 1.0, 1.0]
contour6Display.LineWidth = 2.0

Hide(contour5)
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
#Delete(spreadSheetView1)
#del spreadSheetView1

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

# Hide
Hide(tableToPoints1)

# create a new 'Calculator'
calculator2 = Calculator(Input=tableToPoints1)
calculator2.Function = ''

# Properties modified on calculator1
calculator2.Function = 'coords/6371000.0*6366000.0'

# show data in view
calculator2Display_1 = Show(calculator2, renderView1)
# trace defaults for the display properties.
calculator2Display_1.AmbientColor = [0.0, 0.0, 0.0]
calculator2Display_1.ColorArrayName = [None, '']
calculator2Display_1.EdgeColor = [0.0, 0.0, 0.0]
calculator2Display_1.GlyphType = 'Arrow'
calculator2Display_1.CubeAxesColor = [0.0, 0.0, 0.0]
calculator2Display_1.SetScaleArray = ['POINTS', 'Field 0']
calculator2Display_1.ScaleTransferFunction = 'PiecewiseFunction'
calculator2Display_1.OpacityArray = ['POINTS', 'Field 0']
calculator2Display_1.OpacityTransferFunction = 'PiecewiseFunction'

# change solid color
calculator2Display_1.DiffuseColor = [0.0, 0.0, 0.0]

# Properties modified on calculator1
calculator2.CoordinateResults = 1
# Properties modified on calculator1Display
calculator2Display_1.PointSize = 3.0



# reset view to fit data
renderView1.ResetCamera()

# create a new 'Resample With Dataset'
resampleWithDataset1 = ResampleWithDataset(Input=clip1,
    Source=calculator4)

# Properties modified on resampleWithDataset1
resampleWithDataset1.Tolerance = 2.22044604925031e-16

# show data in view
resampleWithDataset1Display = Show(resampleWithDataset1, renderView1)
# trace defaults for the display properties.
resampleWithDataset1Display.AmbientColor = [0.0, 0.0, 1.0]
resampleWithDataset1Display.ColorArrayName = [None, '']
resampleWithDataset1Display.EdgeColor = [0.0, 0.0, 1.0]
resampleWithDataset1Display.GlyphType = 'Arrow'
resampleWithDataset1Display.CubeAxesColor = [0.0, 0.0, 1.0]
resampleWithDataset1Display.SetScaleArray = ['POINTS', 'vtkValidPointMask']
resampleWithDataset1Display.ScaleTransferFunction = 'PiecewiseFunction'
resampleWithDataset1Display.OpacityArray = ['POINTS', 'vtkValidPointMask']
resampleWithDataset1Display.OpacityTransferFunction = 'PiecewiseFunction'
# hide data in view
Hide(resampleWithDataset1, renderView1)
# set active source
SetActiveSource(resampleWithDataset1)

calculator6 = Calculator(Input=resampleWithDataset1)
calculator6.Function = ''

# Properties modified on calculator5
# Scale from 5mm to 2cm per year
calculator6.Function = 'C_12+C_10+C_8+C_6'

# create a new 'Threshold'
threshold1 = Threshold(Input=calculator6)
threshold1.Scalars = ['POINTS', 'Result']
threshold1.ThresholdRange = [-0.01628371886909008, 0.1265261024236679]

# Properties modified on threshold1
threshold1.ThresholdRange = [-2.5, 0.5]

# show data in view
threshold1Display = Show(threshold1, renderView1)
# trace defaults for the display properties.
threshold1Display.AmbientColor = [0.0, 0.0, 0.0]
threshold1Display.ColorArrayName = [None, '']
threshold1Display.EdgeColor = [0.0, 0.0, 0.0]
threshold1Display.GlyphType = 'Arrow'
threshold1Display.CubeAxesColor = [0.0, 0.0, 0.0]
threshold1Display.ScalarOpacityUnitDistance = 177559.8661261182
threshold1Display.SetScaleArray = ['POINTS', 'result']
threshold1Display.ScaleTransferFunction = 'PiecewiseFunction'
threshold1Display.OpacityArray = ['POINTS', 'result']
threshold1Display.OpacityTransferFunction = 'PiecewiseFunction'



# create a new 'Glyph'
glyph1 = Glyph(Input=threshold1,
    GlyphType='Arrow')
glyph1.Scalars = ['POINTS', 'None']
glyph1.Vectors = ['POINTS', 'None']
glyph1.ScaleFactor = 435460.05000000005
glyph1.GlyphTransform = 'Transform2'

# Properties modified on glyph1
glyph1.Vectors = ['POINTS', 'velocity']
glyph1.ScaleMode = 'vector'
glyph1.ScaleFactor = 7000000.0
glyph1.GlyphMode = 'Every Nth Point'

# show data in view
glyph1Display = Show(glyph1, renderView1)
# trace defaults for the display properties.
glyph1Display.AmbientColor = [0.0, 0.0, 1.0]
glyph1Display.ColorArrayName = ['POINTS', 'velocity']
glyph1Display.EdgeColor = [0.0, 0.0, 1.0]
glyph1Display.LookupTable = velocityLUT
glyph1Display.GlyphType = 'Arrow'
glyph1Display.CubeAxesColor = [0.0, 0.0, 1.0]
glyph1Display.SetScaleArray = ['POINTS', 'vtkValidPointMask']
glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph1Display.OpacityArray = ['POINTS', 'vtkValidPointMask']
glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# turn off scalar coloring
ColorBy(glyph1Display, None)

# change solid color
#glyph1Display.DiffuseColor = [1.0, 1.0, 0.0]
#glyph1Display.DiffuseColor = [0.369, 0.153, 0.157]
glyph1Display.DiffuseColor = [0.545, 0.902, 0.706]

# Add a scale vector
# create a new 'CSV Reader'
scale_vectorcsv = CSVReader(FileName=['/Users/visu/Glerum/01092016/aspect/modelM/python_visu_scripts/scale_vector_zoom.csv'])

# show data in view
scale_vectorcsvDisplay = Show(scale_vectorcsv, spreadSheetView1)
# trace defaults for the display properties.
scale_vectorcsvDisplay.FieldAssociation = 'Row Data'

# create a new 'Table To Points'
tableToPoints3 = TableToPoints(Input=scale_vectorcsv)
tableToPoints3.XColumn = '#x'
tableToPoints3.YColumn = '#x'
tableToPoints3.ZColumn = '#x'

# Properties modified on tableToPoints3
tableToPoints3.YColumn = 'y'
tableToPoints3.ZColumn = 'z'

# show data in view
tableToPoints3Display = Show(tableToPoints3, spreadSheetView1)

# hide data in view
Hide(scale_vectorcsv, spreadSheetView1)

# create a new 'Calculator'
calculator3 = Calculator(Input=tableToPoints3)
calculator3.Function = ''

# Properties modified on calculator3
# Scale from 5mm to 2cm per year
calculator3.Function = 'iHat*vx/0.005*0.02+jHat*vy/0.005*0.02+kHat/0.005*0.02*vz (point_at_-4_-8_6366000_vel 5mm/yr)'

# show data in view
calculator3Display = Show(calculator3, spreadSheetView1)

# hide data in view
Hide(tableToPoints3, spreadSheetView1)

# destroy spreadSheetView1
Delete(spreadSheetView1)
del spreadSheetView1


# set active view
SetActiveView(renderView1)

# create a new 'Glyph'
glyph4 = Glyph(Input=calculator3,
    GlyphType='Arrow')
glyph4.Scalars = ['POINTS', 'None']
glyph4.Vectors = ['POINTS', 'Result']
glyph4.ScaleFactor = 0.1
glyph4.GlyphTransform = 'Transform2'

# Properties modified on glyph4
glyph4.ScaleMode = 'vector'
glyph4.ScaleFactor = 7000000.0
glyph4.GlyphMode = 'All Points'

# show data in view
glyph4Display = Show(glyph4, renderView1)
# trace defaults for the display properties.
glyph4Display.AmbientColor = [0.0, 0.0, 1.0]
glyph4Display.ColorArrayName = [None, '']
glyph4Display.EdgeColor = [0.0, 0.0, 1.0]
glyph4Display.GlyphType = 'Arrow'
glyph4Display.CubeAxesColor = [0.0, 0.0, 1.0]
glyph4Display.SetScaleArray = ['POINTS', 'vx']
glyph4Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph4Display.OpacityArray = ['POINTS', 'vx']
glyph4Display.OpacityTransferFunction = 'PiecewiseFunction'

glyph4Display.DiffuseColor = [0.0, 0.0, 0.0]

# set active view
SetActiveView(renderView1)

# set up 4 other splines as the boundaries of the map
# set spline coordinates
# need to go from latitude to colatitude
#lat_bot = -20
#lat_top = 20
lat_bot = -10
lat_top = 6
lat_mid = (lat_bot + lat_top) / 2.0
#lon_min = -20
#lon_max = 20
lon_min = -7
lon_max = 13
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

SetActiveView(renderView1)
#renderView1.ViewSize = [954, 954]
renderView1.ViewSize = [954, 763]

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)
#velocityLUTColorBar = GetScalarBar(velocityLUT, renderView1)
#velocityLUT.RescaleTransferFunction(0.0, 0.03)
slice1Display.LookupTable = velocityLUT
#velocityLUTColorBar.Position = [0.288, 0.07]
velocityLUTColorBar.Position = [0.288, 0.11]

# write label
text2 = Text()
text2.Text = '5 km'
text2Display = Show(text2, renderView1)
text2Display.Color = [0.0, 0.0, 0.0]
#text2Display.Position = [0.46, 0.94]
text2Display.Position = [0.49, 0.90]
text2Display.FontSize = 10

# write label
text3 = Text()
text3.Text = 'Velocity [cm/yr]'
text3Display = Show(text3, renderView1)
text3Display.Color = [0.0, 0.0, 0.0]
#text3Display.Position = [0.42, 0.10]
text3Display.Position = [0.42, 0.14]
text3Display.FontSize = 10

# write label
text4 = Text()
text4.Text = '[2.0 cm/yr]'
text4Display = Show(text4, renderView1)
text4Display.Color = [0.0, 0.0, 0.0]
text4Display.Position = [0.11, 0.14]
text4Display.FontSize = 10


# set camera
# distance from slice
far = 1.55
domain_center_lon = (lon_min+lon_max)/2.0
domain_center_lat = lat_mid
domain_center_radius = 6371000.0
domain_center_x = domain_center_radius * math.cos(math.radians(domain_center_lon)) * math.sin(math.radians(90-domain_center_lat))
domain_center_y = domain_center_radius * math.sin(math.radians(domain_center_lon)) * math.sin(math.radians(90-domain_center_lat))
domain_center_z = domain_center_radius * math.cos(math.radians(90-domain_center_lat))
renderView1.CameraFocalPoint = [domain_center_x, domain_center_y, domain_center_z]
renderView1.CameraPosition = [domain_center_x * far, domain_center_y * far, domain_center_z * far]
renderView1.CameraViewUp = [0.0, 0.0, domain_center_x]
renderView1.Background = [1,1,1] #White
# get Layout
viewLayout1 = GetLayout()

SaveScreenshot('/Users/visu/Glerum/01092016/aspect/EGU/K08_P06_CSloc_new/slices/slice_vel_reg_sphere_zoom_K08_P06_CSloc_new__1_new.png',layout=viewLayout1, magnification=1, quality=100)

