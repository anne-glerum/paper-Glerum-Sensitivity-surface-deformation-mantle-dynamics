#### import the simple module from the paraview
from paraview.simple import *
import math
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

print 'Loading first data set'

# create a new 'XML Partitioned Unstructured Grid Reader'
solution00000pvtu = XMLPartitionedUnstructuredGridReader(FileName=['/Users/visu/Glerum/01092016/aspect/EGU/K08_P06_CSloc_freeslip_notomo/solution/solution-00000.pvtu'])
solution00000pvtu.PointArrayStatus = ['velocity', 'C_1', 'C_3', 'C_6', 'C_8', 'C_10', 'C_12', 'C_14'] 

# Properties modified on solution00000pvtu
solution00000pvtu.PointArrayStatus = ['velocity','C_14','C_12','C_10','C_8','C_6', 'C_1', 'C_3']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1028, 638]

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

print 'Loading second data set'
# create a new 'XML Partitioned Unstructured Grid Reader'
solution00001pvtu = XMLPartitionedUnstructuredGridReader(FileName=['/Users/visu/Glerum/01092016/aspect/EGU/K08_P06_CSloc/solution/solution-00000.pvtu'])
solution00001pvtu.PointArrayStatus = ['velocity','C_14','C_12','C_10','C_8','C_6', 'C_1', 'C_3']

# Properties modified on solution00001pvtu
solution00001pvtu.PointArrayStatus = ['velocity']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
#first renderView1.ViewSize = [1028, 638]

# show data in view
solution00001pvtuDisplay = Show(solution00001pvtu, renderView1)
# trace defaults for the display properties.
solution00001pvtuDisplay.AmbientColor = [0.0, 0.0, 1.0]
solution00001pvtuDisplay.ColorArrayName = [None, '']
solution00001pvtuDisplay.EdgeColor = [0.0, 0.0, 1.0]
solution00001pvtuDisplay.GlyphType = 'Arrow'
solution00001pvtuDisplay.CubeAxesColor = [0.0, 0.0, 1.0]
solution00001pvtuDisplay.ScalarOpacityUnitDistance = 52856.9360269387
solution00001pvtuDisplay.SetScaleArray = [None, '']
solution00001pvtuDisplay.ScaleTransferFunction = 'PiecewiseFunction'
solution00001pvtuDisplay.OpacityArray = [None, '']
solution00001pvtuDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# reset view to fit data
renderView1.ResetCamera()

# set active source
SetActiveSource(solution00000pvtu)

print 'Adding Euler pole to model field data'
calculator9 = Calculator(Input=solution00000pvtu)
calculator9.Function = ''

# Properties modified on calculator1
calculator9.ResultArrayName = 'Abs_vel'
# add absolute Nb motion Doubrovine 2012
# -44.9 3.72 0.1608
calculator9.Function = 'iHat*(velocity_X+(-1.97685e-09*coordsZ-1.82087e-10*coordsY))+jHat*(velocity_Y+(1.82087e-10*coordsX-1.98376e-09*coordsZ))+kHat*(velocity_Z+(1.98376e-09*coordsY+1.97685e-09*coordsX))'

Hide(calculator9, renderView1)

print 'Subtracting data set 2 from set 1'
programmableFilter1 = ProgrammableFilter(Input=[calculator9, solution00001pvtu])
programmableFilter1.Script = ''
programmableFilter1.RequestInformationScript = ''
programmableFilter1.RequestUpdateExtentScript = ''
programmableFilter1.PythonPath = ''

# Properties modified on programmableFilter1
# Subtract second pvtu velocity from first
programmableFilter1.Script = "vel0=inputs[0].PointData['Abs_vel']\nvel1=inputs[1].PointData['velocity']\noutput.PointData.append(vel0-vel1,'diff')"
programmableFilter1.RequestInformationScript = ''
programmableFilter1.RequestUpdateExtentScript = ''
programmableFilter1.PythonPath = ''

# show data in view
programmableFilter1Display = Show(programmableFilter1, renderView1)

print 'Creating slice from velocity difference'
# create a new 'Slice'
slice1 = Slice(Input=programmableFilter1)
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
Hide(solution00001pvtu, renderView1)
Hide(programmableFilter1, renderView1)

# set scalar coloring
ColorBy(slice1Display, ('POINTS', 'diff'))

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'velocity'
velocityLUT = GetColorTransferFunction('diff')

# get opacity transfer function/opacity map for 'velocity'
velocityPWF = GetOpacityTransferFunction('diff')

# Rescale transfer function
#velocityLUT.RescaleTransferFunction(0.0, 0.005)
velocityLUT.RescaleTransferFunction(0.0, 0.01)
#velocityLUT.RescaleTransferFunction(0.0, 0.0025)

# Rescale transfer function
#velocityPWF.RescaleTransferFunction(0.0, 0.005)
velocityPWF.RescaleTransferFunction(0.0, 0.01)
#velocityPWF.RescaleTransferFunction(0.0, 0.0025)

# get color legend
velocityLUTColorBar = GetScalarBar(velocityLUT, renderView1)

## Properties modified on velocityLUTColorBar
velocityLUTColorBar.Title = ''
velocityLUTColorBar.ComponentTitle = ''
velocityLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
velocityLUTColorBar.LabelColor = [1.0, 1.0, 1.0]
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

#velocityLUT.Annotations = ['0.00', '0.0', '0.0025', '0.25', '0.005', '0.50', '0.0075', '0.75', '0.01', '1.0']
velocityLUT.Annotations = ['0.00', '0.0', '0.005', '0.5', '0.01', '1.0']
#velocityLUT.Annotations = ['0.00', '0.00', '0.00125', '1.25', '0.0025', '2.50']
#velocityLUT.Annotations = ['0.00', '0.0', '0.0025', '0.25', '0.005', '0.50']

# get opacity transfer function/opacity map for 'velocity'
velocityPWF = GetOpacityTransferFunction('diff')
velocityPWF.Points = [0.0, 0.699999988079071, 0.5, 0.0, 0.1, 1.0, 0.5, 0.0]
velocityPWF.AllowDuplicateScalars = 1
velocityPWF.ScalarRangeInitialized = 1

print 'Calculating velocity diff magnitude'
# create a new 'Calculator'
calculator1 = Calculator(Input=slice1)
calculator1.Function = ''

# Properties modified on calculator1
calculator1.ResultArrayName = 'Vel diff mag'
calculator1.Function = 'sqrt(diff_X*diff_X+diff_Y*diff_Y+diff_Z*diff_Z)'

Hide(calculator1, renderView1)

print 'Contouring velocity diff magnitude'
# create a new 'Contour'
contour1 = Contour(Input=calculator1)
contour1.ContourBy = ['POINTS', 'Vel diff mag']
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
contour1.ContourBy = ['POINTS', 'Vel diff mag']
contour1.Isosurfaces = [0.011,0.012, 0.013,0.014,0.015]

# turn off scalar coloring
ColorBy(contour1Display, None)

# change solid color
contour1Display.DiffuseColor = [0.0, 0.0, 0.0]

print 'Summing weak zone fields'
# create a new 'Slice'
slice2 = Slice(Input=solution00000pvtu)
slice2.SliceType = 'Plane'
slice2.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice2.SliceType.Origin = [4721612.875, 0.0, 0.0]

# Properties modified on slice2
slice2.SliceType = 'Sphere'

# Properties modified on slice2.SliceType
slice2.SliceType.Center = [0.0, 0.0, 0.0]
slice2.SliceType.Radius = 6366000.0

# show data in view
slice2Display = Show(slice2, renderView1)
# trace defaults for the display properties.
slice2Display.AmbientColor = [0.0, 0.0, 1.0]
slice2Display.ColorArrayName = [None, '']
slice2Display.EdgeColor = [0.0, 0.0, 1.0]
slice2Display.GlyphType = 'Arrow'
slice2Display.CubeAxesColor = [0.0, 0.0, 1.0]
slice2Display.SetScaleArray = [None, '']
slice2Display.ScaleTransferFunction = 'PiecewiseFunction'
slice2Display.OpacityArray = [None, '']
slice2Display.OpacityTransferFunction = 'PiecewiseFunction'
# create new calculator
calculator5 = Calculator(Input=slice2)
calculator5.Function = ''

# Properties modified on calculator5
calculator5.Function = 'C_12+C_10+C_8+C_6'
#calculator5.ResultArrayName = 'WZ'

Hide(slice2,renderView1)
Hide(calculator5, renderView1)

print 'Contouring weak zones'
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
contour6Display.LineWidth = 4.0
## create a new 'Contour'
#contour6 = Contour(Input=calculator5)
#contour6.ContourBy = ['POINTS', 'WZ']
#contour6.Isosurfaces = [0.5]
#contour6.PointMergeMethod = 'Uniform Binning'
#print 'show contours'
## show data in view
#contour6Display = Show(contour6, renderView1)
## trace defaults for the display properties.
#contour6Display.AmbientColor = [0.0, 0.0, 0.0]
#contour6Display.ColorArrayName = ['POINTS', 'C_15']
#contour6Display.LookupTable = velocityLUT
#contour6Display.EdgeColor = [0.0, 0.0, 0.0]
#contour6Display.GlyphType = 'Arrow'
#contour6Display.CubeAxesColor = [0.0, 0.0, 0.0]
#contour6Display.SetScaleArray = ['POINTS', 'C_15']
#contour6Display.ScaleTransferFunction = 'PiecewiseFunction'
#contour6Display.OpacityArray = ['POINTS', 'C_15']
#contour6Display.OpacityTransferFunction = 'PiecewiseFunction'
#print 'Change contours'
## show color bar/color legend
#contour6Display.SetScalarBarVisibility(renderView1, False)
## Properties modified on contour6
#contour6.ContourBy = ['POINTS', 'WZ']
#contour6.Isosurfaces = [0.5]
#
## turn off scalar coloring
#ColorBy(contour6Display, None)
#
## change solid color
#contour6Display.DiffuseColor = [1.0, 1.0, 1.0]
#contour6Display.LineWidth = 2.0

print 'Read in map points'
# Add map with indicator line
# create a new 'CSV Reader'
#map_lon_latcsv = CSVReader(FileName=['/Users/visu/Glerum/090915/aspect/T_perturb/other_visu_data/map_lon_lat.csv'])

# Properties modified on map_lon_latcsv
#map_lon_latcsv.HaveHeaders = 0

# Create a new 'SpreadSheet View'
#spreadSheetView1 = CreateView('SpreadSheetView')
#spreadSheetView1.BlockSize = 1024L

# show data in view
#map_lon_latcsvDisplay = Show(map_lon_latcsv, spreadSheetView1)
# trace defaults for the display properties.
#map_lon_latcsvDisplay.FieldAssociation = 'Row Data'

# create a new 'Table To Points'
#tableToPoints1 = TableToPoints(Input=map_lon_latcsv)
#tableToPoints1.XColumn = 'Field 0'
#tableToPoints1.YColumn = 'Field 0'
#tableToPoints1.ZColumn = 'Field 0'

# Properties modified on tableToPoints1
#tableToPoints1.YColumn = 'Field 1'
#tableToPoints1.ZColumn = 'Field 2'

# show data in view
#tableToPoints1Display = Show(tableToPoints1, spreadSheetView1)

# hide data in view
#Hide(map_lon_latcsv, spreadSheetView1)

# set active view
#SetActiveView(renderView1)

# set active source
#SetActiveSource(tableToPoints1)

# show data in view
#tableToPoints1Display = Show(tableToPoints1, renderView1)
# trace defaults for the display properties.
#tableToPoints1Display.AmbientColor = [0.0, 0.0, 0.0]
#tableToPoints1Display.ColorArrayName = [None, '']

# Hide
#Hide(tableToPoints1)
map3vtk = LegacyVTKReader(FileNames=['/Users/visu/Glerum/01092016/aspect/EGU/python_visu_scripts/map3.vtk'])
# get color transfer function/color map for 'cell_scalars'
cell_scalarsLUT = GetColorTransferFunction('cell_scalars')
# show data in view
map3vtkDisplay = Show(map3vtk, renderView1)

# create a new 'Calculator'
#calculator2 = Calculator(Input=tableToPoints1)
calculator2 = Calculator(Input=map3vtk)
calculator2.Function = ''

# Properties modified on calculator1
calculator2.CoordinateResults = 1
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

# Properties modified on calculator1Display
calculator2Display_1.PointSize = 3.0
calculator2Display_1.LineWidth = 2.5

Hide(map3vtk, renderView1)

print 'Reading first point values'
# Add map with indicator line
# create a new 'CSV Reader'
point_values0csv = CSVReader(FileName=['/Users/visu/Glerum/01092016/aspect/EGU/K08_P06_CSloc_freeslip_notomo/point_values.txt'])

# Properties modified on point_values0csv
#point_values0csv.HaveHeaders = 1
point_values0csv.FieldDelimiterCharacters = ' '

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.BlockSize = 1024L

# show data in view
point_values0csvDisplay = Show(point_values0csv, spreadSheetView1)
# trace defaults for the display properties.
point_values0csvDisplay.FieldAssociation = 'Row Data'

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=point_values0csv)
tableToPoints1.XColumn = '#<time>'
tableToPoints1.YColumn = '#<time>'
tableToPoints1.ZColumn = '#<time>'

# Properties modified on tableToPoints1
tableToPoints1.XColumn = '<evaluation_point_x>'
tableToPoints1.YColumn = '<evaluation_point_y>'
tableToPoints1.ZColumn = '<evaluation_point_z>'

# show data in view
tableToPoints1Display = Show(tableToPoints1, spreadSheetView1)

# hide data in view
Hide(point_values0csv, spreadSheetView1)

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
calculator4 = Calculator(Input=tableToPoints1)
calculator4.Function = ''

# Properties modified on calculator1
calculator4.ResultArrayName = 'Vel_1'
calculator4.Function = 'iHat*<velocity_x>+jHat*<velocity_y>+kHat*<velocity_z>'

Hide(calculator4, renderView1)

print 'Adding Euler pole to model data'
calculator7 = Calculator(Input=calculator4)
calculator7.Function = ''

# Properties modified on calculator1
calculator7.ResultArrayName = 'Abs_vel'
# add absolute Nb motion Doubrovine 2012
# -44.9 3.72 0.1608
calculator7.Function = 'iHat*(Vel_1_X+(-1.97685e-09*coordsZ-1.82087e-10*coordsY))+jHat*(Vel_1_Y+(1.82087e-10*coordsX-1.98376e-09*coordsZ))+kHat*(Vel_1_Z+(1.98376e-09*coordsY+1.97685e-09*coordsX))'

Hide(calculator7, renderView1)
print 'Reading second point values'
# Add map with indicator line
# create a new 'CSV Reader'
point_values1csv = CSVReader(FileName=['/Users/visu/Glerum/01092016/aspect/EGU/K08_P06_CSloc/point_values.txt'])

# Properties modified on point_values0csv
#point_values1csv.HaveHeaders = 1
point_values1csv.FieldDelimiterCharacters = ' '

# Create a new 'SpreadSheet View'
spreadSheetView2 = CreateView('SpreadSheetView')
spreadSheetView2.BlockSize = 1024L

# show data in view
point_values1csvDisplay = Show(point_values1csv, spreadSheetView2)
# trace defaults for the display properties.
point_values1csvDisplay.FieldAssociation = 'Row Data'

# create a new 'Table To Points'
tableToPoints2 = TableToPoints(Input=point_values1csv)
tableToPoints2.XColumn = '#<time>'
tableToPoints2.YColumn = '#<time>'
tableToPoints2.ZColumn = '#<time>'

# Properties modified on tableToPoints1
tableToPoints2.XColumn = '<evaluation_point_x>'
tableToPoints2.YColumn = '<evaluation_point_y>'
tableToPoints2.ZColumn = '<evaluation_point_z>'

# show data in view
tableToPoints2Display = Show(tableToPoints2, spreadSheetView2)

# hide data in view
Hide(point_values1csv, spreadSheetView2)

# set active view
SetActiveView(renderView1)

# set active source
SetActiveSource(tableToPoints2)

# show data in view
tableToPoints2Display = Show(tableToPoints2, renderView1)
# trace defaults for the display properties.
tableToPoints2Display.AmbientColor = [0.0, 0.0, 0.0]
tableToPoints2Display.ColorArrayName = [None, '']

# Hide
Hide(tableToPoints2)

# create a new 'Calculator'
calculator8 = Calculator(Input=tableToPoints2)
calculator8.Function = ''

# Properties modified on calculator1
calculator8.ResultArrayName = 'Vel 2'
calculator8.Function = 'iHat*<velocity_x>+jHat*<velocity_y>+kHat*<velocity_z>'

Hide(calculator8, renderView1)

print 'Computing velocity vector diff'
programmableFilter2 = ProgrammableFilter(Input=[calculator7, calculator8])
programmableFilter2.Script = ''
programmableFilter2.RequestInformationScript = ''
programmableFilter2.RequestUpdateExtentScript = ''
programmableFilter2.PythonPath = ''

# Properties modified on programmableFilter2
# Subtract second pvtu velocity from first
programmableFilter2.Script = "vel0=inputs[0].PointData['Abs_vel']\nvel1=inputs[1].PointData['Vel 2']\noutput.PointData.append(vel0-vel1,'GPS diff')"
#programmableFilter2.Script = "vel0=inputs[0].PointData['Vel 1']\nvel1=inputs[1].PointData['Vel 2']\noutput.PointData.append(vel0-vel1,'GPS diff')\nprint vel0\nprint vel1"
programmableFilter2.RequestInformationScript = ''
programmableFilter2.RequestUpdateExtentScript = ''
programmableFilter2.PythonPath = ''

# show data in view
programmableFilter2Display = Show(programmableFilter2, renderView1)

Hide(programmableFilter2, renderView1)
Delete(spreadSheetView2)
del spreadSheetView2

print 'Creating vel diff vector glyph'
# create a new 'Glyph'
glyph2 = Glyph(Input=programmableFilter2,
    GlyphType='Arrow')
glyph2.Scalars = ['POINTS', 'vel_mag']
glyph2.Vectors = ['POINTS', 'GPS diff'] ##GPS diff
glyph2.ScaleFactor = 1134255.0
glyph2.GlyphTransform = 'Transform2'

# Properties modified on glyph2
glyph2.ScaleMode = 'vector'
glyph2.ScaleFactor = 25000000.0

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

# Properties modified on glyph2
glyph2.GlyphMode = 'Every Nth Point'

# Properties modified on glyph2Display
glyph2Display.SetScaleArray = [None, 'vel_mag']

# Properties modified on glyph2Display
glyph2Display.OpacityArray = [None, 'vel_mag']

# Properties modified on glyph2
#glyph2.Stride = 5
glyph2.Stride = 1
glyph2Display.SetRepresentationType('Surface')

# change solid color
#glyph2Display.DiffuseColor = [21.0, 119.0, 40.0]
#glyph2Display.AmbientColor = [21.0, 119.0, 40.0]
#glyph2Display.DiffuseColor = [1., 0.0, 0.0]
glyph2Display.DiffuseColor = [0.953, 0.251, 0.576]
#glyph2Display.AmbientColor = [0.2, 1.0, 0.3]


# reset view to fit data
renderView1.ResetCamera()

print 'Reading scale vector'
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
# Scale from 5mm to 2mm per year
#calculator3.Function = 'iHat*vx/0.005*0.001+jHat*vy/0.005*0.001+kHat/0.005*0.001*vz (point_at_-4_-8_6366000_vel 5mm/yr)'
#calculator3.Function = 'iHat*vx/0.005*0.005+jHat*vy/0.005*0.005+kHat/0.005*0.005*vz (point_at_-4_-8_6366000_vel 5mm/yr)'
calculator3.Function = 'iHat*vx/0.005*0.01+jHat*vy/0.005*0.01+kHat/0.005*0.01*vz (point_at_-4_-8_6366000_vel 5mm/yr)'

# show data in view
calculator3Display = Show(calculator3, spreadSheetView1)

# hide data in view
Hide(tableToPoints3, spreadSheetView1)

# destroy spreadSheetView1
Delete(spreadSheetView1)
del spreadSheetView1

Hide(calculator3, renderView1)

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
glyph4.ScaleFactor = glyph2.ScaleFactor
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

#glyph4Display.DiffuseColor = [1.0, 1.0, 1.0]
glyph4Display.DiffuseColor = [0.953, 0.251, 0.576]

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
renderView1.ViewSize = [954, 763]
#renderView1.ViewSize = [1908, 1526]

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
#text2Display.Color = [0.0, 0.0, 0.0]
text2Display.Color = [1.0, 1.0, 1.0]
#text2Display.Position = [0.46, 0.94]
text2Display.Position = [0.475, 0.90]
text2Display.FontSize = 10

# write label
text3 = Text()
text3.Text = 'Velocity [cm/yr]'
text3Display = Show(text3, renderView1)
text3Display.Color = [1.0, 1.0, 1.0]
text3Display.Position = [0.425, 0.14]
text3Display.FontSize = 10

# write label
text4 = Text()
text4.Text = '[1 cm/yr]'
text4Display = Show(text4, renderView1)
#text4Display.Color = [0.0, 0.0, 0.0]
text4Display.Color = [1.0, 1.0, 1.0]
text4Display.Position = [0.13, 0.14]
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

SaveScreenshot('/Users/visu/Glerum/01092016/aspect/EGU/K08_P06_CSloc_freeslip_notomo/slices/slice_diffvel_reg_sphere_zoom_K08_P06_CSloc_freeslip_notomo__1_new.png',
               layout=viewLayout1, magnification=2, quality=100)
#               viewOrLayout=viewLayout1, 
#               ImageResolution=[3816, 3052],
#               FontScaling='Do not scale fonts',
#               TransparentBackground=1,
#               ImageQuality=100)

