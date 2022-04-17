#### import the simple module from the paraview
from paraview.simple import *
import math
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Partitioned Unstructured Grid Reader'
solution00000pvtu = XMLPartitionedUnstructuredGridReader(FileName=['/Users/visu/Glerum/01092016/aspect/modelM_paper2/K08_Nb_fixed/solution/solution-00000.pvtu'])
solution00000pvtu.CellArrayStatus = []
solution00000pvtu.PointArrayStatus = ['velocity', 'p', 'T', 'C_14', 'C_12', 'C_13', 'C_10', 'C_11', 'C_8', 'C_9', 'C_6', 'C_7', 'strain_rate', 'nonadiabatic_temperature', 'density', 'thermal_expansivity', 'specific_heat', 'viscosity']

# Properties modified on solution00000pvtu
solution00000pvtu.PointArrayStatus = ['velocity', 'T', 'C_14', 'C_12', 'C_13', 'C_10', 'C_11', 'C_8', 'C_9', 'C_6', 'C_7', 'nonadiabatic_temperature', 'density', 'thermal_expansivity', 'specific_heat', 'viscosity']

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
ColorBy(slice1Display, ('POINTS', 'velocity'))

# rescale color and/or opacity maps used to include current data range
slice1Display.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

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

# rescale transfer function
velocityLUT.RescaleTransferFunction(0.0, 0.06)
velocityLUT.Annotations = ['0.00', '0', '0.02', '2.0', '0.04', '4.0', '0.06', '6.0', '0.08', '8', '0.10', '10']

# get opacity transfer function/opacity map for 'velocity'
velocityPWF = GetOpacityTransferFunction('velocity')
velocityPWF.Points = [0.0, 0.699999988079071, 0.5, 0.0, 0.1, 1.0, 0.5, 0.0]
velocityPWF.AllowDuplicateScalars = 1
velocityPWF.ScalarRangeInitialized = 1

# hide color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

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

####################
Hide(slice1,renderView1)

# create a new 'Calculator'
calculator3 = Calculator(Input=slice1)
calculator3.Function = ''

# Properties modified on calculator1
calculator3.ResultArrayName = 'Vel vec'
#calculator1.Function = 'sqrt(velocity_X*velocity_X+velocity_Y*velocity_Y+velocity_Z*velocity_Z)'
# add absolute Nb motion Doubrovine 2012
# -44.9 3.72 0.1608
calculator3.Function = 'iHat*(velocity_X+(-1.97685e-09*coordsZ-1.82087e-10*coordsY))+jHat*(velocity_Y+(1.82087e-10*coordsX-1.98376e-09*coordsZ))+kHat*(velocity_Z+(1.98376e-09*coordsY+1.97685e-09*coordsX))'
# remove twice the absolute Nb motion Doubrovine 2012
# -44.9 3.72 -0.3216
# Ex = -0.3216 * 1e-6 * PI / 180.0 * cos(3.72) * cos(-44.9) =-3.96752e-09
# Ey = -0.3216 * 1e-6 * PI / 180.0 * cos(3.72) * sin(-44.9) = 3.9537e-09
# Ez = -0.3216 * 1e-6 * PI / 180.0 * sin(3.72)              =-3.64174e-10
#calculator3.Function = 'iHat*(velocity_X+(3.9537e-09*coordsZ+3.64174e-10*coordsY))+jHat*(velocity_Y+(-3.64174e-10*coordsX+3.96752e-09*coordsZ))+kHat*(velocity_Z+(-3.96752e-09*coordsY-3.9537e-09*coordsX))'
# add absolute Eu motion based on Kreemer 14 + Doubrovine 12
# Eu : -38.63 21.61 0.1232
# Ex = 0.1232 * 1e-6 * PI / 180.0 * cos(21.61) * cos(-38.63) =  1.56179982e-9
# Ey = 0.1232 * 1e-6 * PI / 180.0 * cos(21.61) * sin(-38.63) = -1.24810786e-9
# Ez = 0.1232 * 1e-6 * PI / 180.0 * sin(21.61)               =  7.9190714e-10
#calculator3.Function = 'iHat*(velocity_X+(-1.24810786e-09*coordsZ-7.9190714e-10*coordsY))+jHat*(velocity_Y+(7.9190714e-10*coordsX-1.56179982e-09*coordsZ))+kHat*(velocity_Z+(1.56179982e-09*coordsY+1.24810786*coordsX))'

# create a new 'Glyph'
glyph1 = Glyph(Input=calculator3,
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
glyph1.Vectors = ['POINTS', 'Vel vec']
glyph1.ScaleMode = 'vector'
glyph1.ScaleFactor = 6000000.0
glyph1.GlyphMode = 'Every Nth Point'
glyph1.Stride = 50
###################

## create a new 'Glyph'
#glyph1 = Glyph(Input=slice1,
#    GlyphType='Arrow')
#glyph1.Scalars = ['POINTS', 'None']
#glyph1.Vectors = ['POINTS', 'None']
#glyph1.Orient = 1
#glyph1.ScaleMode = 'off'
#glyph1.ScaleFactor = 412349.4755170329
#glyph1.GlyphMode = 'Uniform Spatial Distribution'
#glyph1.MaximumNumberOfSamplePoints = 5000
#glyph1.Seed = 10339
#glyph1.Stride = 1
#glyph1.GlyphTransform = 'Transform2'
#
## init the 'Arrow' selected for 'GlyphType'
#glyph1.GlyphType.TipResolution = 6
#glyph1.GlyphType.TipRadius = 0.1
#glyph1.GlyphType.TipLength = 0.35
#glyph1.GlyphType.ShaftResolution = 6
#glyph1.GlyphType.ShaftRadius = 0.03
#glyph1.GlyphType.Invert = 0
#
## init the 'Transform2' selected for 'GlyphTransform'
#glyph1.GlyphTransform.Translate = [0.0, 0.0, 0.0]
#glyph1.GlyphTransform.Rotate = [0.0, 0.0, 0.0]
#glyph1.GlyphTransform.Scale = [1.0, 1.0, 1.0]
#
## Properties modified on glyph1
#glyph1.Vectors = ['POINTS', 'velocity']
#glyph1.ScaleMode = 'vector'
#glyph1.ScaleFactor = 6000000.0
#glyph1.GlyphMode = 'Every Nth Point'
#glyph1.Stride = 50


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
#glyph1Display.SetScalarBarVisibility(renderView1, True)

ColorBy(glyph1Display, None)
# change solid color
#glyph1Display.DiffuseColor = [0.0, 0.0, 0.0]
glyph1Display.DiffuseColor = [1.0, 1.0, 0.0]

velocityLUTColorBar.Position = [0.288, 0.05]

# hide color bar/color legend
#glyph1Display.SetScalarBarVisibility(renderView1, False)

# create a new 'Calculator'
calculator1 = Calculator(Input=slice1)
calculator1.Function = ''

# Properties modified on calculator1
calculator1.ResultArrayName = 'Vel mag'
#calculator1.Function = 'sqrt(velocity_X*velocity_X+velocity_Y*velocity_Y+velocity_Z*velocity_Z)'
# add absolute Nb motion Doubrovine 2012
# -44.9 3.72 0.1608
calculator1.Function = 'sqrt((velocity_X+(-1.97685e-09*coordsZ-1.82087e-10*coordsY))^2+(velocity_Y+(1.82087e-10*coordsX-1.98376e-09*coordsZ))^2+(velocity_Z+(1.98376e-09*coordsY+1.97685e-09*coordsX))^2)'
# remove double abs Nb motion D12
#calculator1.Function = 'sqrt((velocity_X+(3.9537e-09*coordsZ+3.64174e-10*coordsY))^2+(velocity_Y+(-3.64174e-10*coordsX+3.96752e-09*coordsZ))^2+(velocity_Z+(-3.96752e-09*coordsY-3.9537e-09*coordsX))^2)'

vel_magLUT = GetColorTransferFunction('Vel mag')
calculator1Display = Show(calculator1, renderView1)
ColorBy(calculator1Display, 'Vel mag')
# rescale transfer function
vel_magLUT.RescaleTransferFunction(0.0, 0.06)
vel_magLUT.Annotations = ['0.00', '0', '0.02', '2.0', '0.04', '4.0', '0.06', '6.0', '0.08', '8', '0.10', '10']

# get color legend
vel_magLUTColorBar = GetScalarBar(vel_magLUT, renderView1)

## Properties modified on velocityLUTColorBar
vel_magLUTColorBar.Title = ''
vel_magLUTColorBar.ComponentTitle = ''
vel_magLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
vel_magLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
vel_magLUTColorBar.AutomaticLabelFormat = 0
vel_magLUTColorBar.LabelFormat = '%-#6.2f'
vel_magLUTColorBar.RangeLabelFormat = '%6.2f'
vel_magLUTColorBar.AutoOrient = 0
vel_magLUTColorBar.Orientation = 'Horizontal'
vel_magLUTColorBar.Position = [100, 0.85]
vel_magLUTColorBar.DrawTickMarks = 0
vel_magLUTColorBar.DrawTickLabels = 0
vel_magLUTColorBar.AddRangeLabels = 0
vel_magLUTColorBar.TitleFontSize = 10
vel_magLUTColorBar.LabelFontSize = 10

vel_magLUT.EnableOpacityMapping = 1
vel_magPWF = GetOpacityTransferFunction('Vel mag')
vel_magPWF.Points = [0.0, 0.699999988079071, 0.5, 0.0, 0.1, 1.0, 0.5, 0.0]
vel_magPWF.AllowDuplicateScalars = 1
vel_magPWF.ScalarRangeInitialized = 1

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
contour1.ContourBy = ['POINTS', 'Vel mag']
contour1.Isosurfaces = [0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18]
#contour1.Isosurfaces = [0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22]

# turn off scalar coloring
ColorBy(contour1Display, None)

# change solid color
contour1Display.DiffuseColor = [0.0, 0.0, 0.0]

#Hide(calculator1, renderView1)

# create a new 'Contour'
contour3 = Contour(Input=slice1)
contour3.ContourBy = ['POINTS', 'C_1']
contour3.Isosurfaces = [0.5]
contour3.PointMergeMethod = 'Uniform Binning'

# properties modified on contour3
contour3.ContourBy = ['POINTS', 'C_14']
contour3.Isosurfaces = [0.1]
# show data in view
contour3Display = Show(contour3, renderView1)
# trace defaults for the display properties.
contour3Display.AmbientColor = [0.0, 0.0, 1.0]
contour3Display.ColorArrayName = ['POINTS', 'velocity']
contour3Display.LookupTable = velocityLUT
contour3Display.EdgeColor = [0.0, 0.0, 1.0]
contour3Display.GlyphType = 'Arrow'
contour3Display.CubeAxesColor = [0.0, 0.0, 1.0]
contour3Display.SetScaleArray = ['POINTS', 'C_1']
contour3Display.ScaleTransferFunction = 'PiecewiseFunction'
contour3Display.OpacityArray = ['POINTS', 'C_1']
contour3Display.OpacityTransferFunction = 'PiecewiseFunction'

# turn off scalar coloring
ColorBy(contour3Display, None)
# change solid color and line width
contour3Display.DiffuseColor = [1.0, 1.0, 1.0]
contour3Display.LineWidth = 2.0

calculator5 = Calculator(Input=slice1)
calculator5.Function = ''

# Properties modified on calculator5
# Scale from 5mm to 2cm per year
calculator5.Function = 'C_12+C_13+C_10+C_11+C_8+C_9+C_6+C_7'

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

#Hide(contour3,renderView1)

#glyph1Display.SetScalarBarVisibility(renderView1, True)
#velocityLUTColorBar = GetScalarBar(velocityLUT, renderView1)
#velocityLUT.RescaleTransferFunction(0.0, 0.06)
#glyph1Display.LookupTable = velocityLUT
#velocityLUTColorBar.Position = [0.288, 0.05]
vel_magLUTColorBar = GetScalarBar(vel_magLUT, renderView1)
vel_magLUT.RescaleTransferFunction(0.0, 0.06)
glyph1Display.LookupTable = velocityLUT
vel_magLUTColorBar.Position = [0.288, 0.05]

# set active source
SetActiveSource(glyph1)

# hide data in view
#Hide(glyph1, renderView1)

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
calculator2 = Calculator(Input=tableToPoints1)
calculator2.Function = ''

Hide(tableToPoints1)

# Properties modified on calculator1
calculator2.Function = 'coords/6371000.0*6371000.0'

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

# set up 4 other splines as the boundaries of the map
# set spline coordinates
# need to go from latitude to colatitude
lat_bot = -20
lat_top = 20
#lat_bot = -13
#lat_top = 7
lat_mid = (lat_bot + lat_top) / 2.0
lon_min = -20
lon_max = 20
#lon_min = -10
#lon_max = 13
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

# hide data in view
#Hide(tableToPoints1, renderView2)

SetActiveView(renderView1)
renderView1.ViewSize = [954, 954]

# show color bar/color legend
#slice1Display.SetScalarBarVisibility(renderView1, True)
#velocityLUTColorBar = GetScalarBar(velocityLUT, renderView1)
#velocityLUT.RescaleTransferFunction(0.0, 0.06)
#slice1Display.LookupTable = velocityLUT
#velocityLUTColorBar.Position = [0.288, 0.07]
calculator1Display.SetScalarBarVisibility(renderView1, True)
vel_magLUTColorBar = GetScalarBar(vel_magLUT, renderView1)
vel_magLUT.RescaleTransferFunction(0.0, 0.06)
calculator1Display.LookupTable = vel_magLUT
vel_magLUTColorBar.Position = [0.288, 0.07]

# write label
text2 = Text()
text2.Text = '0 km'
text2Display = Show(text2, renderView1)
text2Display.Color = [0.0, 0.0, 0.0]
text2Display.Position = [0.46, 0.94]
text2Display.FontSize = 10

# write label
text3 = Text()
text3.Text = 'Velocity [cm/yr]'
text3Display = Show(text3, renderView1)
text3Display.Color = [0.0, 0.0, 0.0]
text3Display.Position = [0.42, 0.10]
text3Display.FontSize = 10

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

# Add GPS
# create a new 'Legacy VTK Reader'
points_vel_cart_newvtk = LegacyVTKReader(FileNames=['/Users/visu/Glerum/01092016/aspect/aegean_plate_motion/other_visu_data/GPS_Noc_Doubr_Nb_frame.vtk'])
#points_vel_cart_newvtk = LegacyVTKReader(FileNames=['/Users/visu/Glerum/01092016/aspect/aegean_plate_motion/other_visu_data/GPS_Noc_Nb_fixed.vtk'])
# create a new 'Glyph'
glyph2 = Glyph(Input=points_vel_cart_newvtk,
    GlyphType='Arrow')
glyph2.Scalars = ['POINTS', 'vel_mag']
glyph2.Vectors = ['POINTS', 'velocity']
glyph2.ScaleFactor = 1134255.0
glyph2.GlyphTransform = 'Transform2'

# Properties modified on glyph2
glyph2.ScaleMode = 'vector'
glyph2.ScaleFactor = 6000000.0

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
glyph2.Stride = 10

# Properties modified on glyph2
glyph2.Stride = 5

# change solid color
#glyph2Display.DiffuseColor = [1.0, 0.0, 0.0]

glyph2Display.DiffuseColor = [1.0, 0.0, 1.0]










# save screenshot
SaveScreenshot('/Users/visu/Glerum/01092016/aspect/modelM_paper2/K08_Nb_fixed/slices/slice_vel_sphere_absNb_K08_Nb_fixed__1.png',layout=viewLayout1, magnification=1, quality=100)
#SaveScreenshot('/Applications/paraview.app/Contents/data/K08_P07/slices/slice_sphere_K08_P07__4.png',layout=viewLayout1, magnification=1, quality=100)

# Loop over all the other slices, changing only the necessary settings

filename_base = "/Users/visu/Glerum/01092016/aspect/modelM_paper2/K08_Nb_fixed/slices/slice_vel_sphere_absNb_K08_Nb_fixed__"
#filename_base = "/Applications/paraview.app/Contents/data/K08_P07/slices/slice_sphere_K08_P07__"
count = 1
#count = 4

for i in xrange(1,27):
    count = count + 1
#    count = count + 4
    radius = 6371000.0 - i * 50000.0
    depth = 6371000.0 - radius

    if (depth > 10000.0):
      # hide data in view
      Hide(glyph2, renderView1)
#    if (depth > 12000.0):
#      Show(contour3, renderView1)    
    

    text2.Text = str(int(depth/1000.0)) + ' km'
 
    # set active view
    SetActiveView(renderView1)
    SetActiveSource(slice1)

    # Properties modified on slice1.SliceType
    slice1.SliceType.Center = [0.0, 0.0, 0.0]
    slice1.SliceType.Radius = radius

    SetActiveSource(line1)
    scale = radius/6371000.0
    line1.Point2 = [lt_x*scale, lt_y*scale, lt_z*scale]
    SetActiveSource(line2)
    line2.Point2 = [rt_x*scale, rt_y*scale, rt_z*scale]
    SetActiveSource(line3)
    line3.Point2 = [lb_x*scale, lb_y*scale, lb_z*scale]
    SetActiveSource(line4)
    line4.Point2 = [rb_x*scale, rb_y*scale, rb_z*scale]

    SetActiveSource(calculator2)
    calculator2.Function = 'coords/6371000.0*'+str(radius)

    # save screenshot
    filename = filename_base + str(count) + ".png"
    SaveScreenshot(filename, layout=viewLayout1, magnification=1, quality=100)



