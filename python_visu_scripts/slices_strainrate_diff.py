#### import the simple module from the paraview
from paraview.simple import *
import math
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Partitioned Unstructured Grid Reader'
solution00000pvtu = XMLPartitionedUnstructuredGridReader(FileName=['/Users/visu/Glerum/01092016/aspect/EGU/K08_Nbx4/solution/solution-00000.pvtu'])
solution00000pvtu.CellArrayStatus = []
solution00000pvtu.PointArrayStatus = ['velocity','strain_rate', 'C_14','C_12','C_10','C_8','C_6','C_15']

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

# create a new 'XML Partitioned Unstructured Grid Reader'
solution00001pvtu = XMLPartitionedUnstructuredGridReader(FileName=['/Users/visu/Glerum/01092016/aspect/EGU/K08_P06_CSloc/solution/solution-00000.pvtu'])
solution00001pvtu.CellArrayStatus = []
solution00001pvtu.PointArrayStatus = ['velocity','strain_rate','C_14','C_12','C_10','C_8','C_6', 'C_1', 'C_3']
solution00001pvtuDisplay = Show(solution00001pvtu, renderView1)

programmableFilter1 = ProgrammableFilter(Input=[solution00000pvtu, solution00001pvtu])
programmableFilter1.Script = ''
programmableFilter1.RequestInformationScript = ''
programmableFilter1.RequestUpdateExtentScript = ''
programmableFilter1.PythonPath = ''

# Properties modified on programmableFilter1
# Subtract second pvtu velocity from first
programmableFilter1.Script = "e0=inputs[0].PointData['strain_rate']\ne1=inputs[1].PointData['strain_rate']\noutput.PointData.append(abs(e0-e1),'diff')"
programmableFilter1.RequestInformationScript = ''
programmableFilter1.RequestUpdateExtentScript = ''
programmableFilter1.PythonPath = ''

# show data in view
programmableFilter1Display = Show(programmableFilter1, renderView1)

# hide data in view
Hide(solution00000pvtu, renderView1)
Hide(solution00001pvtu, renderView1)
Hide(programmableFilter1, renderView1)

diff_range = programmableFilter1.PointData.GetArray("diff").GetRange()
print "strain rate diff range ", diff_range

# create a new 'Slice'
#slice1 = Slice(Input=solution00000pvtu)
slice1 = Slice(Input=programmableFilter1)
slice1.SliceType = 'Sphere'
slice1.Crinkleslice = 0
slice1.Triangulatetheslice = 1
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Center = [4717985.065516494, 0.0, 0.0]
slice1.SliceType.Radius = 1.0

# Properties modified on slice1.SliceType
slice1.SliceType.Center = [0.0, 0.0, 0.0]
slice1.SliceType.Radius = 6366000.0

# show data in view
slice1Display = Show(slice1, renderView1)
# trace defaults for the display properties.
slice1Display.CubeAxesVisibility = 0
slice1Display.Representation = 'Surface'
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
ColorBy(slice1Display, ('POINTS', 'diff'))
# get color transfer function/color map for 'velocity'
diffLUT = GetColorTransferFunction('diff')
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
diffLUT.RescaleTransferFunction(5e-18, 5e-14)
# due to bug in log scale bar and annotations, calculate non-log position:
point1=5e-18
point2=5e-16
#point3=1e-7
point4=4.95e-14
xmin=5e-18
xmax=5e-14
annot1=xmin + (xmax-xmin)*(math.log10(point1)-math.log10(xmin))/(math.log10(xmax)-math.log10(xmin))
annot2=xmin + (xmax-xmin)*(math.log10(point2)-math.log10(xmin))/(math.log10(xmax)-math.log10(xmin))
##annot3=xmin + (xmax-xmin)*(math.log10(point3)-math.log10(xmin))/(math.log10(xmax)-math.log10(xmin))
annot4=xmin + (xmax-xmin)*(math.log10(point4)-math.log10(xmin))/(math.log10(xmax)-math.log10(xmin))
##diffLUT.Annotations = [str(annot1), str(point1), str(annot2), str(point2), str(annot3), str(point3), str(annot4), str(point4)]
diffLUT.Annotations = [str(annot1), str(point1), str(annot2), str(point2), str(annot4), "5e-14"]
#diffLUT.Annotations = [str(annot1), str(point1), str(annot2), str(point2), str(annot4), str(point4)]
diffLUT.MapControlPointsToLogSpace()
diffLUT.UseLogScale = 1

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
contour1.ContourBy = ['POINTS', 'diff']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# show data in view
contour1Display = Show(contour1, renderView1)
# trace defaults for the display properties.
contour1Display.AmbientColor = [0.0, 0.0, 0.0]
contour1Display.ColorArrayName = ['POINTS', 'strain_rate']

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
contour1.ContourBy = ['POINTS', 'strain_rate']
contour1.Isosurfaces = [5e-20, 5e-19]

# turn off scalar coloring
ColorBy(contour1Display, None)

# change solid color
contour1Display.DiffuseColor = [0.0, 0.0, 0.0]

diffLUTColorBar = GetScalarBar(diffLUT, renderView1)
diffLUT.RescaleTransferFunction(5e-18, 5e-14)
#diffLUT.RescaleTransferFunction(-1e-14, 1e-14)
#diffLUT.RescaleTransferFunction(-18, -14)
slice1Display.LookupTable = diffLUT
diffLUTColorBar.Position = [0.288, 0.05]

# rescale color and/or opacity maps used to include current data range
slice1Display.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# hide color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# set active view
SetActiveView(renderView1)

# create a new 'Slice'
slice2 = Slice(Input=solution00000pvtu)
slice2.SliceType = 'Sphere'
slice2.Crinkleslice = 0
slice2.Triangulatetheslice = 1
slice2.SliceOffsetValues = [0.0]
slice2.SliceType.Center = [0.0, 0.0, 0.0]
slice2.SliceType.Radius = 6366000.0

calculator5 = Calculator(Input=slice2)
calculator5.Function = 'C_12+C_10+C_8+C_6'

Hide(slice2,renderView1)
Hide(calculator5, renderView1)

print 'Contouring weak zones'
# create a new 'Contour'
contour6 = Contour(Input=calculator5)
contour6.ContourBy = ['POINTS', 'Result']
contour6.Isosurfaces = [0.5]
contour6.PointMergeMethod = 'Uniform Binning'

# show data in view
contour6Display = Show(contour6, renderView1)
# trace defaults for the display properties.
contour6Display.AmbientColor = [0.0, 0.0, 0.0]
contour6Display.ColorArrayName = ['POINTS', 'C_15']
contour6Display.LookupTable = diffLUT
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
contour6.Isosurfaces = [0.5]

# turn off scalar coloring
ColorBy(contour6Display, None)

# change solid color
contour6Display.DiffuseColor = [1.0, 1.0, 1.0]
contour6Display.LineWidth = 4.0

contour4 = Contour(Input=slice2)
contour4.ContourBy = ['POINTS', 'C_14']
contour4.Isosurfaces = [0.5]
contour4.PointMergeMethod = 'Uniform Binning'

# show data in view
contour4Display = Show(contour4, renderView1)
contour4Display.AmbientColor = [0.0, 0.0, 0.0]
contour4Display.ColorArrayName = ['POINTS', 'C_15']
contour4Display.LookupTable = diffLUT
contour4Display.EdgeColor = [0.0, 0.0, 0.0]
contour4Display.GlyphType = 'Arrow'
contour4Display.CubeAxesColor = [0.0, 0.0, 0.0]
contour4Display.SetScaleArray = ['POINTS', 'C_15']
contour4Display.ScaleTransferFunction = 'PiecewiseFunction'
contour4Display.OpacityArray = ['POINTS', 'C_15']
contour4Display.OpacityTransferFunction = 'PiecewiseFunction'
contour4Display.SetScalarBarVisibility(renderView1, False)

# turn off scalar coloring
ColorBy(contour4Display, None)
# change solid color and line width
contour4Display.DiffuseColor = [1.0, 1.0, 1.0]
contour4Display.LineWidth = 4.0

# Hide slab contour in lithosphere
Hide(contour4,renderView1)
SetActiveSource(slice2)
# create a new 'Contour'
contour5 = Contour(Input=slice2)
contour5.ContourBy = ['POINTS', 'C_15']
contour5.Isosurfaces = [0.5]
contour5.PointMergeMethod = 'Uniform Binning'

# show data in view
contour5Display = Show(contour5, renderView1)
# trace defaults for the display properties.
contour5Display.AmbientColor = [0.0, 0.0, 1.0]
contour5Display.ColorArrayName = ['POINTS', 'velocity']
contour5Display.LookupTable = None
contour5Display.EdgeColor = [0.0, 0.0, 1.0]
contour5Display.GlyphType = 'Arrow'
contour5Display.CubeAxesColor = [0.0, 0.0, 1.0]
contour5Display.SetScaleArray = ['POINTS', 'C_15']
contour5Display.ScaleTransferFunction = 'PiecewiseFunction'
contour5Display.OpacityArray = ['POINTS', 'C_15']
contour5Display.OpacityTransferFunction = 'PiecewiseFunction'

# turn off scalar coloring
ColorBy(contour5Display, None)
# change solid color and line width
contour5Display.DiffuseColor = [1.0, 1.0, 1.0]
contour5Display.LineWidth = 4.0

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

map3vtk = LegacyVTKReader(FileNames=['/Users/visu/Glerum/01092016/aspect/EGU/python_visu_scripts/map3.vtk'])
# get color transfer function/color map for 'cell_scalars'
cell_scalarsLUT = GetColorTransferFunction('cell_scalars')
# show data in view
map3vtkDisplay = Show(map3vtk, renderView1)

# create a new 'Calculator'
calculator4 = Calculator(Input=map3vtk)
calculator4.Function = ''

# Properties modified on calculator1
calculator4.CoordinateResults = 1
calculator4.Function = 'coords/6371000.0*6366000.0'

# show data in view
calculator4Display_1 = Show(calculator4, renderView1)
# trace defaults for the display properties.
calculator4Display_1.AmbientColor = [0.0, 0.0, 0.0]
calculator4Display_1.ColorArrayName = [None, '']
calculator4Display_1.EdgeColor = [0.0, 0.0, 0.0]
calculator4Display_1.GlyphType = 'Arrow'
calculator4Display_1.CubeAxesColor = [0.0, 0.0, 0.0]
calculator4Display_1.SetScaleArray = ['POINTS', 'Field 0']
calculator4Display_1.ScaleTransferFunction = 'PiecewiseFunction'
calculator4Display_1.OpacityArray = ['POINTS', 'Field 0']
calculator4Display_1.OpacityTransferFunction = 'PiecewiseFunction'

# change solid color
calculator4Display_1.DiffuseColor = [0.0, 0.0, 0.0]

# Properties modified on calculator1Display
calculator4Display_1.PointSize = 3.0
calculator4Display_1.LineWidth = 2.5

Hide(map3vtk, renderView1)
# set active view
SetActiveView(renderView1)

scale = 6172000.


# set spline coordinates based on the points with max and min z coordinates (predetermined)
lon_min = -10.43 #-7.93
lon_max = 12.43 #10.93
lon_mid = 1 #1.5
lat_min = -20.
lat_max = 20.
lat_mid = 0.
min_x = scale * math.cos(math.radians(lon_min)) * math.sin(math.radians(90-lat_min))
min_y = scale * math.sin(math.radians(lon_min)) * math.sin(math.radians(90-lat_min))
min_z = scale * math.cos(math.radians(90-lat_min))
mid_x = scale * math.cos(math.radians(lon_mid)) * math.sin(math.radians(90-lat_mid))
mid_y = scale * math.sin(math.radians(lon_mid)) * math.sin(math.radians(90-lat_mid))
mid_z = scale * math.cos(math.radians(90-lat_mid))
max_x = scale * math.cos(math.radians(lon_max)) * math.sin(math.radians(90-lat_max))
max_y = scale * math.sin(math.radians(lon_max)) * math.sin(math.radians(90-lat_max))
max_z = scale * math.cos(math.radians(90-lat_max))
print min_x, min_y, min_z, mid_x, mid_y, mid_z, max_x, max_y, max_z
a=max_x-min_x
b=max_y-min_y
c=max_z-min_z
d=mid_x
e=mid_y
f=mid_z
g=b*f-e*c
h=-a*f+d*c
j=a*e-d*b
R=math.sqrt(g*g+h*h+j*j)
print g/R,h/R,j/R

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
Hide(splineSource6,renderView1)

SetActiveView(renderView1)
renderView1.ViewSize = [954, 763]

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)
diffLUTColorBar = GetScalarBar(diffLUT, renderView1)
slice1Display.LookupTable = diffLUT
diffLUTColorBar.Position = [0.288, 0.07]

# write label
text2 = Text()
text2.Text = '5 km'
text2Display = Show(text2, renderView1)
text2Display.Color = [1.0, 1.0, 1.0]
text2Display.Color = [0.0, 0.0, 0.0]
text2Display.Position = [0.47, 0.94]
text2Display.Position = [0.475, 0.90]
text2Display.FontSize = 10

# write label
text3 = Text()
text3.Text = 'Strain rate difference [1/s]'
text3Display = Show(text3, renderView1)
text3Display.Color = [0.0, 0.0, 0.0]
text3Display.Color = [1.0, 1.0, 1.0]
text3Display.Position = [0.36, 0.10]
text3Display.FontSize = 10

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
splineSource5Display.LineWidth = 3.0
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

diffLUT.RescaleTransferFunction(5e-18, 5e-14)
#diffLUT.RescaleTransferFunction(-1e-14, 1e-14)
#diffLUT.RescaleTransferFunction(-18, -14)
# save screenshot
SaveScreenshot('/Users/visu/Glerum/01092016/aspect/EGU/K08_Nbx4/slices/slice_diff_strainrate_sphere_K08_Nbx4_ref_K08_P06_CSloc_zoom__1.png',layout=viewLayout1, magnification=2, quality=100)

# Loop over all the other slices, changing only the necessary settings

filename_base = "/Users/visu/Glerum/01092016/aspect/EGU/K08_Nbx4/slices/slice_diff_strainrate_sphere_K08_Nbx4_ref_K08_P06_CSloc_zoom__"
count = 1
#count = 4

for i in xrange(1,5):
    count = count + 1
#    count = count + 4
    radius = 6371000.0 - i * 50000.0
    depth = 6371000.0 - radius
    loc = max(0.1,depth*2./2900000.0)
    print loc
    text3Display.Position = [0.36, 0.10]
    text3Display.Color = [0.0, 0.0, 0.0]

  #  if (depth > 55000.0):

    if (depth > 105000.0):
        Show(contour4,renderView1)
      #velocityLUT.Annotations = ['0.00', '0.0', '0.025', '2.5', '0.05', '5.0', '0.075', '7.5', '0.10', '10']
      #text4Display.Color = [0.0, 0.0, 0.0]
      #text3Display.Color = [0.0, 0.0, 0.0]
      #text2Display.Color = [0.0, 0.0, 0.0]
      #velocityLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
    
#    if (depth > 505000.0):
      # hide data in view
      #calculator8.Function = 'iHat*radius/6371000*coordsX+jHat*radius/6371000*coordsY+kHat*radius/6371000*coordsZ'
      #calculator3.Function = 'iHat*vx/0.005*0.03+jHat*vy/0.005*0.03+kHat/0.005*0.03*vz (point_at_-5.25_-9.75_6366000_vel 5mm/yr)'
      #text4.Text = '[3.0 cm/yr]'
      #velocityLUT.RescaleTransferFunction(0.0, 0.03)
      #velocityLUT.Annotations = ['0.00', '0', '0.01', '1.0', '0.02', '2.0', '0.03', '3.0', '0.04', '4', '0.05', '5']
      #velocityLUT.Annotations = ['0.00', '0', '0.01', '1.0', '0.02', '2.0', '0.03', '3.0']
    
    

    text2.Text = str(int(depth/1000.0)) + ' km'
 
    # set active view
    SetActiveView(renderView1)
    SetActiveSource(slice1)

    # Properties modified on slice1.SliceType
    slice1.SliceType.Center = [0.0, 0.0, 0.0]
    slice1.SliceType.Radius = radius
    slice2.SliceType.Center = [0.0, 0.0, 0.0]
    slice2.SliceType.Radius = radius

    scale = radius/6371000.0

    SetActiveSource(calculator4)
    calculator4.Function = 'coords/6371000.0*'+str(radius)

    # save screenshot
    filename = filename_base + str(count) + ".png"
    SaveScreenshot(filename, layout=viewLayout1, magnification=2, quality=100)



