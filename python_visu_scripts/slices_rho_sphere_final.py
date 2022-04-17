#### import the simple module from the paraview
from paraview.simple import *
import math
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Partitioned Unstructured Grid Reader'
solution00000pvtu = XMLPartitionedUnstructuredGridReader(FileName=['/Users/visu/Glerum/01092016/aspect/EGU/K08_P06_CSloc/solution/solution-00000.pvtu'])
solution00000pvtu.CellArrayStatus = []
solution00000pvtu.PointArrayStatus = ['velocity', 'p', 'T', 'C_1', 'C_2', 'C_3', 'strain_rate', 'nonadiabatic_temperature', 'density', 'thermal_expansivity', 'specific_heat', 'viscosity', 'shear_heating', 'adiabatic_heating']

# Properties modified on solution00000pvtu
solution00000pvtu.PointArrayStatus = ['density']

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
solution00000pvtu_1 = XMLPartitionedUnstructuredGridReader(FileName=['/Users/visu/Glerum/01092016/aspect/EGU/adiab_tbl_noslab3/solution/solution-00000.pvtu'])
solution00000pvtu_1.CellArrayStatus = []
solution00000pvtu_1.PointArrayStatus = ['velocity', 'p', 'T', 'C_1', 'C_2', 'C_3', 'strain_rate', 'nonadiabatic_temperature', 'density', 'thermal_expansivity', 'specific_heat', 'viscosity']

# Properties modified on solution00000pvtu_1
solution00000pvtu_1.PointArrayStatus = ['density']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
solution00000pvtu_1Display = Show(solution00000pvtu_1, renderView1)
# trace defaults for the display properties.
solution00000pvtu_1Display.CubeAxesVisibility = 0
solution00000pvtu_1Display.Representation = 'Surface'
solution00000pvtu_1Display.AmbientColor = [0.0, 0.0, 0.0]
solution00000pvtu_1Display.ColorArrayName = [None, '']
solution00000pvtu_1Display.DiffuseColor = [1.0, 1.0, 1.0]
solution00000pvtu_1Display.LookupTable = None
solution00000pvtu_1Display.MapScalars = 1
solution00000pvtu_1Display.InterpolateScalarsBeforeMapping = 1
solution00000pvtu_1Display.Opacity = 1.0
solution00000pvtu_1Display.PointSize = 2.0
solution00000pvtu_1Display.LineWidth = 1.0
solution00000pvtu_1Display.Interpolation = 'Gouraud'
solution00000pvtu_1Display.Specular = 0.0
solution00000pvtu_1Display.SpecularColor = [1.0, 1.0, 1.0]
solution00000pvtu_1Display.SpecularPower = 100.0
solution00000pvtu_1Display.Ambient = 0.0
solution00000pvtu_1Display.Diffuse = 1.0
solution00000pvtu_1Display.EdgeColor = [0.0, 0.0, 0.0]
solution00000pvtu_1Display.BackfaceRepresentation = 'Follow Frontface'
solution00000pvtu_1Display.BackfaceAmbientColor = [1.0, 1.0, 1.0]
solution00000pvtu_1Display.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
solution00000pvtu_1Display.BackfaceOpacity = 1.0
solution00000pvtu_1Display.Position = [0.0, 0.0, 0.0]
solution00000pvtu_1Display.Scale = [1.0, 1.0, 1.0]
solution00000pvtu_1Display.Orientation = [0.0, 0.0, 0.0]
solution00000pvtu_1Display.Origin = [0.0, 0.0, 0.0]
solution00000pvtu_1Display.Pickable = 1
solution00000pvtu_1Display.Texture = None
solution00000pvtu_1Display.Triangulate = 0
solution00000pvtu_1Display.NonlinearSubdivisionLevel = 1
solution00000pvtu_1Display.GlyphType = 'Arrow'
solution00000pvtu_1Display.CubeAxesColor = [0.0, 0.0, 0.0]
solution00000pvtu_1Display.CubeAxesCornerOffset = 0.0
solution00000pvtu_1Display.CubeAxesFlyMode = 'Closest Triad'
solution00000pvtu_1Display.CubeAxesInertia = 1
solution00000pvtu_1Display.CubeAxesTickLocation = 'Inside'
solution00000pvtu_1Display.CubeAxesXAxisMinorTickVisibility = 1
solution00000pvtu_1Display.CubeAxesXAxisTickVisibility = 1
solution00000pvtu_1Display.CubeAxesXAxisVisibility = 1
solution00000pvtu_1Display.CubeAxesXGridLines = 0
solution00000pvtu_1Display.CubeAxesXTitle = 'X-Axis'
solution00000pvtu_1Display.CubeAxesUseDefaultXTitle = 1
solution00000pvtu_1Display.CubeAxesYAxisMinorTickVisibility = 1
solution00000pvtu_1Display.CubeAxesYAxisTickVisibility = 1
solution00000pvtu_1Display.CubeAxesYAxisVisibility = 1
solution00000pvtu_1Display.CubeAxesYGridLines = 0
solution00000pvtu_1Display.CubeAxesYTitle = 'Y-Axis'
solution00000pvtu_1Display.CubeAxesUseDefaultYTitle = 1
solution00000pvtu_1Display.CubeAxesZAxisMinorTickVisibility = 1
solution00000pvtu_1Display.CubeAxesZAxisTickVisibility = 1
solution00000pvtu_1Display.CubeAxesZAxisVisibility = 1
solution00000pvtu_1Display.CubeAxesZGridLines = 0
solution00000pvtu_1Display.CubeAxesZTitle = 'Z-Axis'
solution00000pvtu_1Display.CubeAxesUseDefaultZTitle = 1
solution00000pvtu_1Display.CubeAxesGridLineLocation = 'All Faces'
solution00000pvtu_1Display.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
solution00000pvtu_1Display.CustomBoundsActive = [0, 0, 0]
solution00000pvtu_1Display.OriginalBoundsRangeActive = [0, 0, 0]
solution00000pvtu_1Display.CustomRange = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
solution00000pvtu_1Display.CustomRangeActive = [0, 0, 0]
solution00000pvtu_1Display.UseAxesOrigin = 0
solution00000pvtu_1Display.AxesOrigin = [0.0, 0.0, 0.0]
solution00000pvtu_1Display.CubeAxesXLabelFormat = '%-#6.3g'
solution00000pvtu_1Display.CubeAxesYLabelFormat = '%-#6.3g'
solution00000pvtu_1Display.CubeAxesZLabelFormat = '%-#6.3g'
solution00000pvtu_1Display.StickyAxes = 0
solution00000pvtu_1Display.CenterStickyAxes = 0
solution00000pvtu_1Display.SelectionCellLabelBold = 0
solution00000pvtu_1Display.SelectionCellLabelColor = [0.0, 1.0, 0.0]
solution00000pvtu_1Display.SelectionCellLabelFontFamily = 'Arial'
solution00000pvtu_1Display.SelectionCellLabelFontSize = 18
solution00000pvtu_1Display.SelectionCellLabelItalic = 0
solution00000pvtu_1Display.SelectionCellLabelJustification = 'Left'
solution00000pvtu_1Display.SelectionCellLabelOpacity = 1.0
solution00000pvtu_1Display.SelectionCellLabelShadow = 0
solution00000pvtu_1Display.SelectionPointLabelBold = 0
solution00000pvtu_1Display.SelectionPointLabelColor = [1.0, 1.0, 0.0]
solution00000pvtu_1Display.SelectionPointLabelFontFamily = 'Arial'
solution00000pvtu_1Display.SelectionPointLabelFontSize = 18
solution00000pvtu_1Display.SelectionPointLabelItalic = 0
solution00000pvtu_1Display.SelectionPointLabelJustification = 'Left'
solution00000pvtu_1Display.SelectionPointLabelOpacity = 1.0
solution00000pvtu_1Display.SelectionPointLabelShadow = 0
solution00000pvtu_1Display.ScalarOpacityUnitDistance = 51152.954947465136
solution00000pvtu_1Display.SelectMapper = 'Projected tetra'
solution00000pvtu_1Display.GaussianRadius = 0.0
solution00000pvtu_1Display.ShaderPreset = 'Sphere'
solution00000pvtu_1Display.Emissive = 0
solution00000pvtu_1Display.ScaleByArray = 0
solution00000pvtu_1Display.SetScaleArray = ['POINTS', 'T']
solution00000pvtu_1Display.ScaleTransferFunction = 'PiecewiseFunction'
solution00000pvtu_1Display.OpacityByArray = 0
solution00000pvtu_1Display.OpacityArray = ['POINTS', 'T']
solution00000pvtu_1Display.OpacityTransferFunction = 'PiecewiseFunction'

## init the 'Arrow' selected for 'GlyphType'
#solution00000pvtu_1Display.GlyphType.TipResolution = 6
#solution00000pvtu_1Display.GlyphType.TipRadius = 0.1
#solution00000pvtu_1Display.GlyphType.TipLength = 0.35
#solution00000pvtu_1Display.GlyphType.ShaftResolution = 6
#solution00000pvtu_1Display.GlyphType.ShaftRadius = 0.03
#solution00000pvtu_1Display.GlyphType.Invert = 0
#
## init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
#solution00000pvtu_1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
#
## init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
#solution00000pvtu_1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
#
## reset view to fit data
#renderView1.ResetCamera()
#
##changing interaction mode based on data extents
#renderView1.InteractionMode = '3D'
#
# create a new 'Programmable Filter'
programmableFilter1 = ProgrammableFilter(Input=[solution00000pvtu, solution00000pvtu_1])
programmableFilter1.Script = ''
programmableFilter1.RequestInformationScript = ''
programmableFilter1.RequestUpdateExtentScript = ''
programmableFilter1.PythonPath = ''

# Properties modified on programmableFilter1
# Subtract second pvtu velocity from first
programmableFilter1.Script = "rho0=inputs[0].PointData['density']\nrho1=inputs[1].PointData['density']\noutput.PointData.append(rho0-rho1,'diff')"
programmableFilter1.RequestInformationScript = ''
programmableFilter1.RequestUpdateExtentScript = ''
programmableFilter1.PythonPath = ''

# show data in view
programmableFilter1Display = Show(programmableFilter1, renderView1)

# hide data in view
Hide(solution00000pvtu, renderView1)
Hide(solution00000pvtu_1, renderView1)
Hide(programmableFilter1, renderView1)

# create a new 'Slice'
slice1 = Slice(Input=programmableFilter1)
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
#Hide(solution00000pvtu, renderView1)

# set scalar coloring
ColorBy(slice1Display, ('POINTS', 'diff'))

# rescale color and/or opacity maps used to include current data range
slice1Display.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

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
diffLUT.RescaleTransferFunction(-80.0, 80.0)
diffLUT.Annotations = ['-80.0', '-80', '-40.0', '-40', '0.0', '0', '40.0', '40', '80.0', '80']

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
contour1Display.ColorArrayName = ['POINTS', 'density']
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
contour1.ContourBy = ['POINTS', 'diff']
contour1.Isosurfaces = [-240.0, -200.0, -160.0, -120., -80.0, -40.0, 40.0, 80.0, 120., 160., 200.0, 240.0]

# turn off scalar coloring
ColorBy(contour1Display, None)

# change solid color
contour1Display.DiffuseColor = [0.0, 0.0, 0.0]

diffLUTColorBar = GetScalarBar(diffLUT, renderView1)
diffLUT.RescaleTransferFunction(-80.0, 80.0)
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
calculator1.Function = 'coords/6371000.0*6366000.0'

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
diffLUT.RescaleTransferFunction(-80.0, 80.0)
slice1Display.LookupTable = diffLUT
diffLUTColorBar.Position = [0.288, 0.11]

# write label
text3 = Text()
text3.Text = 'dDensity [kg/m3]'
text3Display = Show(text3, renderView1)
text3Display.Color = [0.0, 0.0, 0.0]
text3Display.Position = [0.43, 0.14]
text3Display.FontSize = 10

# write label
text2 = Text()
text2.Text = '0 km'
text2Display = Show(text2, renderView1)
text2Display.Color = [0.0, 0.0, 0.0]
text2Display.Position = [0.47, 0.90]
text2Display.FontSize = 10

# save screenshot
SaveScreenshot('/Users/visu/Glerum/01092016/aspect/EGU/K08_P06_CSloc/slices/slice_rho_sphere_K08_P06_CSloc__1.png',layout=viewLayout1, magnification=1, quality=100)

# Loop over all the other slices, changing only the necessary settings

filename_base = "/Users/visu/Glerum/01092016/aspect/EGU/K08_P06_CSloc/slices/slice_rho_sphere_K08_P06_CSloc__"
count = 1

for i in xrange(1,30):
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

    # save screenshot
    filename = filename_base + str(count) + ".png"
    SaveScreenshot(filename, layout=viewLayout1, magnification=1, quality=100)


