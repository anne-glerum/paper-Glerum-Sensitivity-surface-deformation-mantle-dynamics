#### import the simple module from the paraview
from paraview.simple import *
import math
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Partitioned Unstructured Grid Reader'
solution00000pvtu = XMLPartitionedUnstructuredGridReader(FileName=['/Users/visu/Glerum/01092016/aspect/EGU/K08_SEMUCB/solution/solution-00000.pvtu'])
solution00000pvtu.CellArrayStatus = []
solution00000pvtu.PointArrayStatus = ['velocity', 'p', 'T', 'C_1', 'C_2', 'C_14', 'strain_rate', 'nonadiabatic_temperature', 'density', 'thermal_expansivity', 'specific_heat', 'viscosity']

# Properties modified on solution00000pvtu
# C_14 is the Aegean slab
solution00000pvtu.PointArrayStatus = ['velocity', 'C_14', 'viscosity']

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
# 6371029.5 111206.7 3.9e-10
slice1.SliceType.Origin = [6371029.5, 111206.7, 3.9e-10]
slice1.SliceType.Normal = [-0.015327735, 0.87812535, -0.4781850365]
#slice1.SliceType.Origin = [5341292.5, 81227.5811174638, -38203.7658714104]
#slice1.SliceType.Normal = [0., 0.843311524056685, -0.59]

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
ColorBy(slice1Display, ('POINTS', 'viscosity'))

# rescale color and/or opacity maps used to include current data range
slice1Display.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

viscosityLUT = GetColorTransferFunction('viscosity')
viscosityLUT.LockDataRange = 1
viscosityLUT.InterpretValuesAsCategories = 0
viscosityLUT.ShowCategoricalColorsinDataRangeOnly = 0
viscosityLUT.RescaleOnVisibilityChange = 0
viscosityLUT.EnableOpacityMapping = 0
viscosityLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 0.05, 0.865003, 0.865003, 0.865003, 0.1, 0.705882, 0.0156863, 0.14902]
viscosityLUT.UseLogScale = 0
viscosityLUT.ColorSpace = 'Diverging'
viscosityLUT.UseBelowRangeColor = 0
viscosityLUT.BelowRangeColor = [0.0, 0.0, 0.0]
viscosityLUT.UseAboveRangeColor = 0
viscosityLUT.AboveRangeColor = [1.0, 1.0, 1.0]
viscosityLUT.NanColor = [1.0, 1.0, 0.0]
viscosityLUT.Discretize = 1
viscosityLUT.NumberOfTableValues = 32
viscosityLUT.ScalarRangeInitialized = 1.0
viscosityLUT.HSVWrap = 0
viscosityLUT.VectorComponent = 0
viscosityLUT.VectorMode = 'Magnitude'
viscosityLUT.AllowDuplicateScalars = 1
viscosityLUT.Annotations = []
viscosityLUT.ActiveAnnotatedValues = []
viscosityLUT.IndexedColors = []
# get color transfer function/color map for 'velocity'
velocityLUT = GetColorTransferFunction('velocity')
velocityLUT.LockDataRange = 1
velocityLUT.InterpretValuesAsCategories = 0
velocityLUT.ShowCategoricalColorsinDataRangeOnly = 0
velocityLUT.RescaleOnVisibilityChange = 0
velocityLUT.EnableOpacityMapping = 1
velocityLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 0.05, 0.865003, 0.865003, 0.865003, 0.1, 0.705882, 0.0156863, 0.14902]
velocityLUT.UseLogScale = 0
velocityLUT.ColorSpace = 'Diverging'
velocityLUT.UseBelowRangeColor = 0
velocityLUT.BelowRangeColor = [0.0, 0.0, 0.0]
velocityLUT.UseAboveRangeColor = 0
velocityLUT.AboveRangeColor = [1.0, 1.0, 1.0]
velocityLUT.NanColor = [1.0, 1.0, 0.0]
velocityLUT.Discretize = 1
velocityLUT.NumberOfTableValues = 32
velocityLUT.ScalarRangeInitialized = 1.0
velocityLUT.HSVWrap = 0
velocityLUT.VectorComponent = 0
velocityLUT.VectorMode = 'Magnitude'
velocityLUT.AllowDuplicateScalars = 1
velocityLUT.Annotations = []
velocityLUT.ActiveAnnotatedValues = []
velocityLUT.IndexedColors = []

# set viscosity color function
viscosityLUT.ApplyPreset('Plasma', True)
viscosityPWF = GetOpacityTransferFunction('viscosity')
# rescale transfer function
viscosityLUT.RescaleTransferFunction(1e19, 1e24)
# due to bug in log scale bar and annotations, calculate non-log position:
point1=1e19
point2=5e20
point3=5e22
point4=1e24
xmin=1e19
xmax=1e24
annot1=xmin + (xmax-xmin)*(math.log10(point1)-math.log10(xmin))/(math.log10(xmax)-math.log10(xmin))
annot2=xmin + (xmax-xmin)*(math.log10(point2)-math.log10(xmin))/(math.log10(xmax)-math.log10(xmin))
annot3=xmin + (xmax-xmin)*(math.log10(point3)-math.log10(xmin))/(math.log10(xmax)-math.log10(xmin))
annot4=xmin + (xmax-xmin)*(math.log10(point4)-math.log10(xmin))/(math.log10(xmax)-math.log10(xmin))
viscosityLUT.Annotations = [str(annot1), str(point1), str(annot2), str(point2), str(annot3), str(point3), str(annot4), str(point4)]
viscosityLUT.MapControlPointsToLogSpace()
viscosityLUT.UseLogScale = 1

# get opacity transfer function/opacity map for 'velocity'
velocityPWF = GetOpacityTransferFunction('velocity')
velocityPWF.Points = [0.0, 0.699999988079071, 0.5, 0.0, 0.1, 1.0, 0.5, 0.0]
velocityPWF.AllowDuplicateScalars = 1
velocityPWF.ScalarRangeInitialized = 1

# hide color bar/color legend
#slice1Display.SetScalarBarVisibility(renderView1, True)

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

# get color legend
viscosityLUTColorBar = GetScalarBar(viscosityLUT, renderView1)

## Properties modified on viscosityLUTColorBar
viscosityLUTColorBar.Title = ''
viscosityLUTColorBar.ComponentTitle = ''
viscosityLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
viscosityLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
viscosityLUTColorBar.AutomaticLabelFormat = 0
viscosityLUTColorBar.LabelFormat = '%-#1.0e'
viscosityLUTColorBar.RangeLabelFormat = '%1.0e'
viscosityLUTColorBar.AutoOrient = 0
viscosityLUTColorBar.Orientation = 'Horizontal'
viscosityLUTColorBar.Position = [100, 0.85]
viscosityLUTColorBar.DrawTickMarks = 0
viscosityLUTColorBar.DrawTickLabels = 0
viscosityLUTColorBar.AddRangeLabels = 0
viscosityLUTColorBar.TitleFontSize = 10
viscosityLUTColorBar.LabelFontSize = 10

# create a new 'Glyph'
glyph1 = Glyph(Input=slice1,
    GlyphType='Arrow')
glyph1.Scalars = ['POINTS', 'None']
glyph1.Vectors = ['POINTS', 'None']
glyph1.Orient = 1
glyph1.ScaleMode = 'off'
glyph1.ScaleFactor = 412349.4755170329
glyph1.GlyphMode = 'Uniform Spatial Distribution'
glyph1.MaximumNumberOfSamplePoints = 5000
glyph1.Seed = 10339
glyph1.Stride = 1
glyph1.GlyphTransform = 'Transform2'

# init the 'Arrow' selected for 'GlyphType'
glyph1.GlyphType.TipResolution = 6
glyph1.GlyphType.TipRadius = 0.1
glyph1.GlyphType.TipLength = 0.35
glyph1.GlyphType.ShaftResolution = 6
glyph1.GlyphType.ShaftRadius = 0.03
glyph1.GlyphType.Invert = 0

# init the 'Transform2' selected for 'GlyphTransform'
glyph1.GlyphTransform.Translate = [0.0, 0.0, 0.0]
glyph1.GlyphTransform.Rotate = [0.0, 0.0, 0.0]
glyph1.GlyphTransform.Scale = [1.0, 1.0, 1.0]

# Properties modified on glyph1
glyph1.Vectors = ['POINTS', 'velocity']
glyph1.ScaleMode = 'vector'
#CHANGE
glyph1.ScaleFactor = 4000000.0
#glyph1.ScaleFactor = 2000000.0
glyph1.GlyphMode = 'Every Nth Point'
glyph1.Stride = 30

# show data in view
glyph1Display = Show(glyph1, renderView1)
# trace defaults for the display properties.
glyph1Display.CubeAxesVisibility = 0
glyph1Display.Representation = 'Surface'
glyph1Display.AmbientColor = [0.0, 0.0, 0.0]
glyph1Display.ColorArrayName = ['POINTS', 'velocity']
glyph1Display.DiffuseColor = [1.0, 1.0, 1.0]
glyph1Display.LookupTable = velocityLUT
glyph1Display.MapScalars = 1
glyph1Display.InterpolateScalarsBeforeMapping = 1
glyph1Display.Opacity = 1.0
glyph1Display.PointSize = 2.0
glyph1Display.LineWidth = 1.0
glyph1Display.Interpolation = 'Gouraud'
glyph1Display.Specular = 0.0
glyph1Display.SpecularColor = [1.0, 1.0, 1.0]
glyph1Display.SpecularPower = 100.0
glyph1Display.Ambient = 0.0
glyph1Display.Diffuse = 1.0
glyph1Display.EdgeColor = [0.0, 0.0, 0.0]
glyph1Display.BackfaceRepresentation = 'Follow Frontface'
glyph1Display.BackfaceAmbientColor = [1.0, 1.0, 1.0]
glyph1Display.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
glyph1Display.BackfaceOpacity = 1.0
glyph1Display.Position = [0.0, 0.0, 0.0]
glyph1Display.Scale = [1.0, 1.0, 1.0]
glyph1Display.Orientation = [0.0, 0.0, 0.0]
glyph1Display.Origin = [0.0, 0.0, 0.0]
glyph1Display.Pickable = 1
glyph1Display.Texture = None
glyph1Display.Triangulate = 0
glyph1Display.NonlinearSubdivisionLevel = 1
glyph1Display.GlyphType = 'Arrow'
glyph1Display.CubeAxesColor = [0.0, 0.0, 0.0]
glyph1Display.CubeAxesCornerOffset = 0.0
glyph1Display.CubeAxesFlyMode = 'Closest Triad'
glyph1Display.CubeAxesInertia = 1
glyph1Display.CubeAxesTickLocation = 'Inside'
glyph1Display.CubeAxesXAxisMinorTickVisibility = 1
glyph1Display.CubeAxesXAxisTickVisibility = 1
glyph1Display.CubeAxesXAxisVisibility = 1
glyph1Display.CubeAxesXGridLines = 0
glyph1Display.CubeAxesXTitle = 'X-Axis'
glyph1Display.CubeAxesUseDefaultXTitle = 1
glyph1Display.CubeAxesYAxisMinorTickVisibility = 1
glyph1Display.CubeAxesYAxisTickVisibility = 1
glyph1Display.CubeAxesYAxisVisibility = 1
glyph1Display.CubeAxesYGridLines = 0
glyph1Display.CubeAxesYTitle = 'Y-Axis'
glyph1Display.CubeAxesUseDefaultYTitle = 1
glyph1Display.CubeAxesZAxisMinorTickVisibility = 1
glyph1Display.CubeAxesZAxisTickVisibility = 1
glyph1Display.CubeAxesZAxisVisibility = 1
glyph1Display.CubeAxesZGridLines = 0
glyph1Display.CubeAxesZTitle = 'Z-Axis'
glyph1Display.CubeAxesUseDefaultZTitle = 1
glyph1Display.CubeAxesGridLineLocation = 'All Faces'
glyph1Display.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
glyph1Display.CustomBoundsActive = [0, 0, 0]
glyph1Display.OriginalBoundsRangeActive = [0, 0, 0]
glyph1Display.CustomRange = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
glyph1Display.CustomRangeActive = [0, 0, 0]
glyph1Display.UseAxesOrigin = 0
glyph1Display.AxesOrigin = [0.0, 0.0, 0.0]
glyph1Display.CubeAxesXLabelFormat = '%-#6.3g'
glyph1Display.CubeAxesYLabelFormat = '%-#6.3g'
glyph1Display.CubeAxesZLabelFormat = '%-#6.3g'
glyph1Display.StickyAxes = 0
glyph1Display.CenterStickyAxes = 0
glyph1Display.SelectionCellLabelBold = 0
glyph1Display.SelectionCellLabelColor = [0.0, 1.0, 0.0]
glyph1Display.SelectionCellLabelFontFamily = 'Arial'
glyph1Display.SelectionCellLabelFontSize = 18
glyph1Display.SelectionCellLabelItalic = 0
glyph1Display.SelectionCellLabelJustification = 'Left'
glyph1Display.SelectionCellLabelOpacity = 1.0
glyph1Display.SelectionCellLabelShadow = 0
glyph1Display.SelectionPointLabelBold = 0
glyph1Display.SelectionPointLabelColor = [1.0, 1.0, 0.0]
glyph1Display.SelectionPointLabelFontFamily = 'Arial'
glyph1Display.SelectionPointLabelFontSize = 18
glyph1Display.SelectionPointLabelItalic = 0
glyph1Display.SelectionPointLabelJustification = 'Left'
glyph1Display.SelectionPointLabelOpacity = 1.0
glyph1Display.SelectionPointLabelShadow = 0
glyph1Display.GaussianRadius = 0.0
glyph1Display.ShaderPreset = 'Sphere'
glyph1Display.Emissive = 0
glyph1Display.ScaleByArray = 0
glyph1Display.SetScaleArray = ['POINTS', 'T']
glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph1Display.OpacityByArray = 0
glyph1Display.OpacityArray = ['POINTS', 'T']
glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'Arrow' selected for 'GlyphType'
glyph1Display.GlyphType.TipResolution = 6
glyph1Display.GlyphType.TipRadius = 0.1
glyph1Display.GlyphType.TipLength = 0.35
glyph1Display.GlyphType.ShaftResolution = 6
glyph1Display.GlyphType.ShaftRadius = 0.03
glyph1Display.GlyphType.Invert = 0

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
glyph1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
glyph1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, False)

velocityLUTColorBar.Position = [0.288, 0.05]

# hide color bar/color legend
#glyph1Display.SetScalarBarVisibility(renderView1, False)

# create a new 'Calculator'
calculator1 = Calculator(Input=slice1)
calculator1.Function = ''

# Properties modified on calculator1
calculator1.ResultArrayName = 'Vel mag'
calculator1.Function = 'sqrt(velocity_X*velocity_X+velocity_Y*velocity_Y+velocity_Z*velocity_Z)'

# create a new 'Contour'
contour1 = Contour(Input=calculator1)
contour1.ContourBy = ['POINTS', 'C_1']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# show data in view
contour1Display = Show(contour1, renderView1)
# trace defaults for the display properties.
contour1Display.AmbientColor = [0.0, 0.0, 0.0]
contour1Display.ColorArrayName = ['POINTS', 'velocity']
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
#contour1.ContourBy = ['POINTS', 'Vel mag']
contour1.ContourBy = ['POINTS', 'C_14']
#contour1.Isosurfaces = [0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22]
contour1.Isosurfaces = [0.5]

# turn off scalar coloring
ColorBy(contour1Display, None)

# change solid color
contour1Display.DiffuseColor = [0.0, 0.0, 0.0]

Hide(calculator1, renderView1)
#Hide(contour1, renderView1)

# create a new 'Contour'
contour2 = Contour(Input=slice1)
contour2.ContourBy = ['POINTS', 'viscosity']
contour2.PointMergeMethod = 'Uniform Binning'
contour2.Isosurfaces = [5e19,5e20,1e22, 1e23, 5e23]

# show data in view
contour2Display = Show(contour2, renderView1)
# trace defaults for the display properties.
contour2Display.AmbientColor = [0.0, 0.0, 0.0]
contour2Display.ColorArrayName = ['POINTS', 'viscosity']
contour2Display.LookupTable = viscosityLUT
contour2Display.EdgeColor = [0.0, 0.0, 0.0]
contour2Display.GlyphType = 'Arrow'
contour2Display.CubeAxesColor = [0.0, 0.0, 0.0]
contour2Display.SetScaleArray = ['POINTS', 'C_2']
contour2Display.ScaleTransferFunction = 'PiecewiseFunction'
contour2Display.OpacityArray = ['POINTS', 'C_2']
contour2Display.OpacityTransferFunction = 'PiecewiseFunction'

# show color bar/color legend
contour2Display.SetScalarBarVisibility(renderView1, False)

# turn off scalar coloring
ColorBy(contour2Display, None)

# change solid color
contour2Display.DiffuseColor = [0.0, 0.0, 0.0]

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, False)
viscosityLUTColorBar = GetScalarBar(viscosityLUT, renderView1)
viscosityLUT.RescaleTransferFunction(1e19, 1e24)
glyph1Display.LookupTable = viscosityLUT
viscosityLUTColorBar.Position = [0.288, 0.06]

print 'reading disk grid'
vertical_slice_gridcsv = CSVReader(FileName=['/Users/visu/Glerum/01092016/aspect/EGU/python_visu_scripts/disk_grid_eq_6370_100km_05deg.csv'])

spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.BlockSize = 1024L

# show data in view
vertical_slice_gridcsvDisplay = Show(vertical_slice_gridcsv, spreadSheetView1)
# trace defaults for the display properties.
vertical_slice_gridcsvDisplay.FieldAssociation = 'Row Data'
print 'Table to points1'
# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=vertical_slice_gridcsv)
tableToPoints1.XColumn = '#x'
tableToPoints1.YColumn = 'y'
tableToPoints1.ZColumn = 'z'

# show data in view
tableToPoints1Display = Show(tableToPoints1, spreadSheetView1)

# hide data in view
Hide(vertical_slice_gridcsv, spreadSheetView1)
Hide(tableToPoints1, spreadSheetView1)
#SetActiveSource(vertical_slice_gridcsv)

print 'converting grid'
calculator2 = Calculator(Input=tableToPoints1)
#calculator2.Function = 'coordsX*iHat+jHat*(coordsY*0.4781850365-coordsZ*0.878259113)+kHat*(0.878259113*coordsY+coordsZ*0.4781850365)'
calculator2.Function = 'iHat*(coordsX*(0.4781850365+6370030/6371000*6370030/6371000*(1-0.4781850365))+coordsY*(111189/6371000*6370030/6371000*(1-0.4781850365)+0.00130671/6371000*0.878259113)+coordsZ*(6370030/6371000*-0.00130671/6371000*(1-0.4781850365)+111189/6371000*0.878259113))+jHat*(coordsX*(111189/6371000*6370030/6371000*(1-0.4781850365)-0.00130671/6371000*0.878259113)+coordsY*(0.4781850365+111189/6371000*111189/6371000*(1-0.4781850365))+coordsZ*(111189/6371000*-0.00130671/6371000*(1-0.4781850365)-6370030/6371000*0.878259113))+kHat*(coordsX*(6370030/6371000*-0.00130671/6371000*(1-0.4781850365)-111189/6371000*0.878259113)+coordsY*(111189/6371000*-0.00130671/6371000*(1-0.4781850365)+6370030/6371000*0.878259113)+coordsZ*(0.4781850365+0.00130671/6371000*0.0013067/6371000*(1.0-0.4781850365)))'
calculator2.CoordinateResults = 1
# show data in view
calculator2Display = Show(calculator2, spreadSheetView1)

Hide(calculator2, spreadSheetView1)

Delete(spreadSheetView1)
del spreadSheetView1
SetActiveView(renderView1)
calculator2Display = Show(calculator2, renderView1)
Hide(calculator2, renderView1)

print 'Done converting grid'
# create a new 'Generate Surface Normals'
generateSurfaceNormals1 = GenerateSurfaceNormals(Input=slice1)

# show data in view
generateSurfaceNormals1Display = Show(generateSurfaceNormals1, renderView1)
# trace defaults for the display properties.
generateSurfaceNormals1Display.AmbientColor = [0.0, 0.0, 0.0]
generateSurfaceNormals1Display.ColorArrayName = [None, '']
generateSurfaceNormals1Display.EdgeColor = [0.0, 0.0, 0.0]
generateSurfaceNormals1Display.GlyphType = 'Arrow'
generateSurfaceNormals1Display.CubeAxesColor = [0.0, 0.0, 0.0]
generateSurfaceNormals1Display.SetScaleArray = [None, '']
generateSurfaceNormals1Display.ScaleTransferFunction = 'PiecewiseFunction'
generateSurfaceNormals1Display.OpacityArray = [None, '']
generateSurfaceNormals1Display.OpacityTransferFunction = 'PiecewiseFunction'

# create a new 'Python Calculator'
pythonCalculator1 = PythonCalculator(Input=generateSurfaceNormals1)
pythonCalculator1.Expression = ''

# Properties modified on pythonCalculator1
pythonCalculator1.Expression = 'dot(Normals,norm(velocity))'

# show data in view
pythonCalculator1Display = Show(pythonCalculator1, renderView1)
# trace defaults for the display properties.
pythonCalculator1Display.AmbientColor = [0.0, 0.0, 0.0]
pythonCalculator1Display.ColorArrayName = [None, '']
pythonCalculator1Display.EdgeColor = [0.0, 0.0, 0.0]
pythonCalculator1Display.GlyphType = 'Arrow'
pythonCalculator1Display.CubeAxesColor = [0.0, 0.0, 0.0]
pythonCalculator1Display.SetScaleArray = ['POINTS', 'result']
pythonCalculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
pythonCalculator1Display.OpacityArray = ['POINTS', 'result']
pythonCalculator1Display.OpacityTransferFunction = 'PiecewiseFunction'

# hide data in view
Hide(generateSurfaceNormals1, renderView1)

# create a new 'Threshold'
threshold1 = Threshold(Input=pythonCalculator1)
threshold1.Scalars = ['POINTS', 'result']
threshold1.ThresholdRange = [-0.01628371886909008, 0.1265261024236679]

# Properties modified on threshold1
threshold1.ThresholdRange = [-10.0, -0.0]

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

# create a new 'Threshold'
#threshold3 = Threshold(Input=threshold1)
threshold3 = Threshold(Input=slice1)
threshold3.Scalars = ['POINTS', 'C_14']
threshold3.ThresholdRange = [-0.01628371886909008, 0.1265261024236679]

# Properties modified on threshold3
threshold3.ThresholdRange = [0.5, 1.5]

# show data in view
threshold3Display = Show(threshold3, renderView1)
# trace defaults for the display properties.
threshold3Display.AmbientColor = [0.0, 0.0, 0.0]
threshold3Display.ColorArrayName = [None, '']
threshold3Display.EdgeColor = [0.0, 0.0, 0.0]
threshold3Display.GlyphType = 'Arrow'
threshold3Display.CubeAxesColor = [0.0, 0.0, 0.0]
threshold3Display.ScalarOpacityUnitDistance = 177559.8661261182
threshold3Display.SetScaleArray = ['POINTS', 'result']
threshold3Display.ScaleTransferFunction = 'PiecewiseFunction'
threshold3Display.OpacityArray = ['POINTS', 'result']
threshold3Display.OpacityTransferFunction = 'PiecewiseFunction'

print 'Resampling'
resampleWithDataset1 = ResampleWithDataset(Input=threshold3,
    Source=calculator2)

## Properties modified on resampleWithDataset1
resampleWithDataset1.Tolerance = 2.22044604925031e-16

## show data in view
resampleWithDataset1Display = Show(resampleWithDataset1, renderView1)
## trace defaults for the display properties.
resampleWithDataset1Display.AmbientColor = [0.0, 0.0, 1.0]
resampleWithDataset1Display.ColorArrayName = [None, '']
resampleWithDataset1Display.EdgeColor = [0.0, 0.0, 1.0]
resampleWithDataset1Display.GlyphType = 'Arrow'
resampleWithDataset1Display.CubeAxesColor = [0.0, 0.0, 1.0]
resampleWithDataset1Display.SetScaleArray = ['POINTS', 'vtkValidPointMask']
resampleWithDataset1Display.ScaleTransferFunction = 'PiecewiseFunction'
resampleWithDataset1Display.OpacityArray = ['POINTS', 'vtkValidPointMask']
resampleWithDataset1Display.OpacityTransferFunction = 'PiecewiseFunction'
## hide data in view
Hide(resampleWithDataset1, renderView1)
# hide data in view
Hide(pythonCalculator1, renderView1)
# create a new 'Glyph'
#glyph3 = Glyph(Input=threshold3,
glyph3 = Glyph(Input=resampleWithDataset1,
    GlyphType='Arrow')
glyph3.Scalars = ['POINTS', 'None']
glyph3.Vectors = ['POINTS', 'None']
glyph3.ScaleFactor = 435802.05000000005
glyph3.GlyphTransform = 'Transform2'

# Properties modified on glyph1
glyph3.Vectors = ['POINTS', 'velocity']
glyph3.ScaleMode = 'vector'
glyph3.ScaleFactor = glyph1.ScaleFactor
glyph3.GlyphMode = 'Every Nth Point'
glyph3.Stride = 1

# show data in view
glyph3Display = Show(glyph3, renderView1)
# trace defaults for the display properties.
glyph3Display.AmbientColor = [0.0, 0.0, 0.0]
glyph3Display.ColorArrayName = [None, '']
glyph3Display.EdgeColor = [0.0, 0.0, 0.0]
glyph3Display.GlyphType = 'Arrow'
glyph3Display.CubeAxesColor = [0.0, 0.0, 0.0]
glyph3Display.SetScaleArray = ['POINTS', 'result']
glyph3Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph3Display.OpacityArray = ['POINTS', 'result']
glyph3Display.OpacityTransferFunction = 'PiecewiseFunction'

# hide data in view
Hide(threshold1, renderView1)

# set active source
SetActiveSource(pythonCalculator1)

# create a new 'Threshold'
threshold2 = Threshold(Input=pythonCalculator1)
threshold2.Scalars = ['POINTS', 'result']
threshold2.ThresholdRange = [-0.01628371886909008, 0.1265261024236679]

# Properties modified on threshold2
threshold2.ThresholdRange = [0.0, 10.0]

# show data in view
threshold2Display = Show(threshold2, renderView1)
# trace defaults for the display properties.
threshold2Display.AmbientColor = [0.0, 0.0, 0.0]
threshold2Display.ColorArrayName = [None, '']
threshold2Display.EdgeColor = [0.0, 0.0, 0.0]
threshold2Display.GlyphType = 'Arrow'
threshold2Display.CubeAxesColor = [0.0, 0.0, 0.0]
threshold2Display.ScalarOpacityUnitDistance = 206567.89145085454
threshold2Display.SetScaleArray = ['POINTS', 'result']
threshold2Display.ScaleTransferFunction = 'PiecewiseFunction'
threshold2Display.OpacityArray = ['POINTS', 'result']
threshold2Display.OpacityTransferFunction = 'PiecewiseFunction'

# create a new 'Threshold'
threshold4 = Threshold(Input=threshold2)
threshold4.Scalars = ['POINTS', 'C_14']
threshold4.ThresholdRange = [-0.01628371886909008, 0.1265261024236679]

# Properties modified on threshold4
threshold4.ThresholdRange = [0.5, 1.5]

# show data in view
threshold4Display = Show(threshold4, renderView1)
# trace defaults for the display properties.
threshold4Display.AmbientColor = [0.0, 0.0, 0.0]
threshold4Display.ColorArrayName = [None, '']
threshold4Display.EdgeColor = [0.0, 0.0, 0.0]
threshold4Display.GlyphType = 'Arrow'
threshold4Display.CubeAxesColor = [0.0, 0.0, 0.0]
threshold4Display.ScalarOpacityUnitDistance = 177559.8661261182
threshold4Display.SetScaleArray = ['POINTS', 'result']
threshold4Display.ScaleTransferFunction = 'PiecewiseFunction'
threshold4Display.OpacityArray = ['POINTS', 'result']
threshold4Display.OpacityTransferFunction = 'PiecewiseFunction'

# hide data in view
Hide(pythonCalculator1, renderView1)
Hide(threshold3, renderView1)
Hide(threshold4, renderView1)

# create a new 'Glyph'
glyph2 = Glyph(Input=threshold4,
    GlyphType='Arrow')
glyph2.Scalars = ['POINTS', 'None']
glyph2.Vectors = ['POINTS', 'None']
glyph2.ScaleFactor = 380784.92500000005
glyph2.GlyphTransform = 'Transform2'

# Properties modified on glyph2
glyph2.Vectors = ['POINTS', 'velocity']
glyph2.ScaleMode = 'vector'
glyph2.ScaleFactor = glyph1.ScaleFactor
glyph2.GlyphMode = 'Every Nth Point'
glyph2.Stride = 10

# show data in view
glyph2Display = Show(glyph2, renderView1)
# trace defaults for the display properties.
glyph2Display.AmbientColor = [0.0, 0.0, 0.0]
glyph2Display.ColorArrayName = [None, '']
glyph2Display.EdgeColor = [0.0, 0.0, 0.0]
glyph2Display.GlyphType = 'Arrow'
glyph2Display.CubeAxesColor = [0.0, 0.0, 0.0]
glyph2Display.SetScaleArray = ['POINTS', 'result']
glyph2Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph2Display.OpacityArray = ['POINTS', 'result']
glyph2Display.OpacityTransferFunction = 'PiecewiseFunction'

# hide data in view
Hide(threshold2, renderView1)
Hide(glyph2, renderView1)

# change solid color
#glyph2Display.DiffuseColor = [0.0, 1.0, 0.0]
# for presi make everything pink
#glyph2Display.DiffuseColor = [1.0, 0.0, 1.0]
# for presi make everything pink
glyph2Display.DiffuseColor = [0.545, 0.902, 0.706]

# set active source
SetActiveSource(glyph1)

# change solid color
#glyph3Display.DiffuseColor = [1.0, 0.0, 1.0]
glyph3Display.DiffuseColor = [0.545, 0.902, 0.706] # green
glyph3.GlyphType.TipRadius = 0.2
glyph3.GlyphType.TipResolution = 12
glyph3.GlyphType.ShaftRadius = 0.09
glyph3.GlyphType.ShaftResolution = 12 

# hide data in view
Hide(glyph1, renderView1)

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

min_x = 5888785.15988 * scale_410
min_y = -1083981.8546 * scale_410
min_z = -2179352.35327 * scale_410
mid_x = 6371029.51354 * scale_410
mid_y = 111206.733818 * scale_410
mid_z = 3.9e-10
max_x = 5847367.44991 * scale_410
max_y = 1288837.16041 * scale_410
max_z = 2179352.35327 * scale_410

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
scale_slice = 5500
center = [origin[0] / 6371 * scale_slice, origin[1] / 6371 * scale_slice, origin[2] / 6371 * scale_slice ]
renderView1.CameraFocalPoint = [ center[0], center[1], center[2] ]
# position = focal point - X * n_norm
normal = slice1.SliceType.Normal
norm = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]) ** 0.5
normal_norm = [ normal[0] / norm, normal[1] / norm, normal[2] / norm ]
position = [center[0] - far * normal_norm[0], center[1] - far * normal_norm[1], center[2] - far * normal_norm[2]]
renderView1.CameraPosition = [ position[0], position[1], position[2] ]
# View up is the radial vector
#renderView1.CameraViewUp = [origin[0], origin[1], origin[2] ]
renderView1.CameraViewUp = [origin[0], min_y+(max_y-min_y)/2, min_z+(max_z-min_z)/2 ]

# ACG set background
renderView1.Background = [1,1,1] #White

# get Layout
viewLayout1 = GetLayout()

# split cell
#viewLayout1.SplitVertical(0, 0.730)
#viewLayout1.SplitVertical(0, 0.5)
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
## current camera placement for renderView2
#renderView2.CameraFocalPoint = [mid_x, mid_y, mid_z]
#renderView2.CameraViewUp = [-mid_y, mid_x, mid_z]
#far2 = 6.0e6
#radial = [mid_x, mid_y, mid_z]
#norm_r = (radial[0] * radial[0] + radial[1] * radial[1] + radial[2] * radial[2]) ** 0.5
#radial_norm = [ radial[0] / norm_r, radial[1] / norm_r, radial[2] / norm_r ]
#position = [radial[0] + far2 * radial_norm[0], radial[1] + far2 * radial_norm[1], radial[2] + far2 * radial_norm[2]]
#renderView2.CameraPosition = [ position[0], position[1], position[2] ]
#
SetActiveView(renderView1)
#renderView1.ViewSize = [954, 945]
#renderView1.ViewSize = [954, 764]
#renderView1.ViewSize = [954, 473]
##renderView1.ViewSize = [1054, 473]
renderView1.ViewSize = [954, 477]

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, False)
slice1Display.SetScalarBarVisibility(renderView1, True)
viscosityLUTColorBar = GetScalarBar(viscosityLUT, renderView1)
viscosityLUT.RescaleTransferFunction(1e19,1e24)
slice1Display.LookupTable = viscosityLUT
viscosityLUTColorBar.Position = [0.288, 0.08]
viscosityLUT.Annotations = [str(annot1), str(point1), str(annot2), str(point2), str(annot3), str(point3), str(annot4), str(point4)]

# write label
text3 = Text()
text3.Text = 'Viscosity [Pa s]'
text3Display = Show(text3, renderView1)
text3Display.Color = [0.0, 0.0, 0.0]
text3Display.Position = [0.42, 0.13]
text3Display.FontSize = 10

# write label
text4 = Text()
#CHANGE
#text4.Text = '[10.0 cm/yr]'
text4.Text = '[5.0 cm/yr]'
text4Display = Show(text4, renderView1)
text4Display.Color = [1.0, 1.0, 1.0]
text4Display.Color = [0.0, 0.0, 0.0]
text4Display.Position = [0.11, 0.13]
text4Display.FontSize = 10

scale_vectorcsv = CSVReader(FileName=['/Users/visu/Glerum/01092016/aspect/EGU/python_visu_scripts/scale_vector_vertical_slice.csv'])

spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.BlockSize = 1024L

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
#CHANGE
calculator3.Function = 'iHat*vx/0.005*0.05+jHat*vy/0.005*0.05+kHat/0.005*0.05*vz (vel 5mm/yr)'
#calculator3.Function = 'iHat*vx/0.005*0.10+jHat*vy/0.005*0.10+kHat/0.005*0.10*vz (vel 5mm/yr)'

# show data in view
calculator3Display = Show(calculator3, spreadSheetView1)

# hide data in view
Hide(tableToPoints3, spreadSheetView1)

# destroy spreadSheetView1
Delete(spreadSheetView1)
del spreadSheetView1
# create a new 'Glyph'
glyph4 = Glyph(Input=calculator3,
    GlyphType='Arrow')
glyph4.Scalars = ['POINTS', 'None']
glyph4.Vectors = ['POINTS', 'Result']
glyph4.ScaleFactor = 0.1
glyph4.GlyphTransform = 'Transform2'

# Properties modified on glyph4
glyph4.ScaleMode = 'vector'
glyph4.ScaleFactor = glyph1.ScaleFactor
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

glyph4Display.DiffuseColor = [1.0, 1.0, 1.0]
glyph4.GlyphType.TipRadius = 0.2
glyph4.GlyphType.TipResolution = 12
glyph4.GlyphType.ShaftRadius = 0.09
glyph4.GlyphType.ShaftResolution = 12 
viscosityLUT.ApplyPreset('Plasma', True)

Hide(contour2, renderView1)
filename_base = "/Users/visu/Glerum/01092016/aspect/EGU/K08_SEMUCB/slices/slice_vel_slab_K08_SEMUCB__"
count = 1
filename = filename_base + str(count) + ".png"

# save screenshot
SaveScreenshot(filename,layout=viewLayout1, magnification=2, quality=100)
