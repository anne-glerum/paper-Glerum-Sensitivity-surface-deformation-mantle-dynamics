#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
solution00000pvtu = XMLPartitionedUnstructuredGridReader(FileName=['/Users/visu/Glerum/01092016/aspect/EGU/K08_SEMUCB/solution/solution-00000.pvtu'])
solution00000pvtu.CellArrayStatus = []
solution00000pvtu.PointArrayStatus = ['velocity', 'C_14', 'viscosity']
renderView1 = GetActiveViewOrCreate('RenderView')

# get color transfer function/color map for 'viscosity'
viscosityLUT = GetColorTransferFunction('viscosity')

# get opacity transfer function/opacity map for 'viscosity'
viscosityPWF = GetOpacityTransferFunction('viscosity')

# create a new 'Contour'
contour2 = Contour(Input=solution00000pvtu)
contour2.ContourBy = ['POINTS', 'C_1']
contour2.Isosurfaces = [0.502566933631897]
contour2.PointMergeMethod = 'Uniform Binning'

# Properties modified on contour2
contour2.ContourBy = ['POINTS', 'C_14']

# uncomment following to set a specific view size
# renderView1.ViewSize = [1033, 638]

# show data in view
contour2Display = Show(contour2, renderView1)
# trace defaults for the display properties.
contour2Display.AmbientColor = [0.0, 0.0, 1.0]
contour2Display.ColorArrayName = ['POINTS', 'viscosity']
contour2Display.LookupTable = viscosityLUT
contour2Display.EdgeColor = [0.0, 0.0, 1.0]
contour2Display.GlyphType = 'Arrow'
contour2Display.CubeAxesColor = [0.0, 0.0, 1.0]
contour2Display.SetScaleArray = ['POINTS', 'C_1']
contour2Display.ScaleTransferFunction = 'PiecewiseFunction'
contour2Display.OpacityArray = ['POINTS', 'C_1']
contour2Display.OpacityTransferFunction = 'PiecewiseFunction'

# hide data in view
Hide(solution00000pvtu, renderView1)

# show color bar/color legend
contour2Display.SetScalarBarVisibility(renderView1, True)

# Rescale transfer function
viscosityLUT.RescaleTransferFunction(9.99999998051e+18, 8.93726292845e+23)

# Properties modified on contour2
contour2.Isosurfaces = [0.5]

# Rescale transfer function
viscosityLUT.RescaleTransferFunction(9.99999998051e+18, 8.93726292845e+23)

# create a new 'Calculator'
calculator1 = Calculator(Input=contour2)
calculator1.Function = ''

# Properties modified on calculator1
calculator1.Function = '(6371000-sqrt(coordsX*coordsX+coordsY*coordsY+coordsZ*coordsZ))/1000'

# get color transfer function/color map for 'Result'
resultLUT = GetColorTransferFunction('Result')

# show data in view
calculator1Display = Show(calculator1, renderView1)
# trace defaults for the display properties.
calculator1Display.AmbientColor = [0.0, 0.0, 1.0]
calculator1Display.ColorArrayName = ['POINTS', 'Result']
calculator1Display.LookupTable = resultLUT
calculator1Display.EdgeColor = [0.0, 0.0, 1.0]
calculator1Display.GlyphType = 'Arrow'
calculator1Display.CubeAxesColor = [0.0, 0.0, 1.0]
calculator1Display.SetScaleArray = ['POINTS', 'Result']
calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
calculator1Display.OpacityArray = ['POINTS', 'Result']
calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'

# hide data in view
Hide(contour2, renderView1)

# show color bar/color legend
calculator1Display.SetScalarBarVisibility(renderView1, True)

# get opacity transfer function/opacity map for 'Result'
resultPWF = GetOpacityTransferFunction('Result')

# get color legend/bar for resultLUT in view renderView1
resultLUTColorBar = GetScalarBar(resultLUT, renderView1)

# Properties modified on resultLUTColorBar
resultLUTColorBar.Title = 'Depth [km]'
resultLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
resultLUTColorBar.TitleFontSize = 10
resultLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
resultLUTColorBar.LabelFontSize = 10

# Properties modified on calculator1

# Rescale transfer function
resultLUT.RescaleTransferFunction(0.499871158286, 718400.544602)

# Rescale transfer function
resultPWF.RescaleTransferFunction(0.499871158286, 718400.544602)

# rescale color and/or opacity maps used to exactly fit the current data range
calculator1Display.RescaleTransferFunctionToDataRange(False)

# Properties modified on resultLUTColorBar
resultLUTColorBar.AutomaticLabelFormat = 0
resultLUTColorBar.LabelFormat = '%-#3.0f'
resultLUTColorBar.RangeLabelFormat = '%3.0f'

# Rescale transfer function
resultLUT.RescaleTransferFunction(0.0, 720.0)

# Rescale transfer function
resultPWF.RescaleTransferFunction(0.0, 720.0)

# Properties modified on resultLUTColorBar
resultLUTColorBar.LabelFormat = '%#3.0f'

# Properties modified on resultLUTColorBar
resultLUTColorBar.LabelFormat = '%3.0f'

# reset view to fit data
renderView1.ResetCamera()

# Properties modified on resultLUT
resultLUT.Annotations = ['', '']
resultLUT.IndexedColors = [0.0, 0.0, 0.0]

# Properties modified on resultLUT
resultLUT.Annotations = ['0', '']

# Properties modified on resultLUT
resultLUT.Annotations = ['0', '0']

# Properties modified on resultLUT
resultLUT.Annotations = ['0', '0', '', '']
resultLUT.IndexedColors = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# Properties modified on resultLUT
resultLUT.Annotations = ['0', '0', '180', '']

# Properties modified on resultLUT
resultLUT.Annotations = ['0', '0', '180', '180']

# Properties modified on resultLUT
resultLUT.Annotations = ['0', '0', '180', '180', '', '']
resultLUT.IndexedColors = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# Properties modified on resultLUT
resultLUT.Annotations = ['0', '0', '180', '180', '360', '']

# Properties modified on resultLUT
resultLUT.Annotations = ['0', '0', '180', '180', '360', '360']

# Properties modified on resultLUT
resultLUT.Annotations = ['0', '0', '180', '180', '360', '360', '', '']
resultLUT.IndexedColors = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# Properties modified on resultLUT
resultLUT.Annotations = ['0', '0', '180', '180', '360', '360', '540', '']

# Properties modified on resultLUT
resultLUT.Annotations = ['0', '0', '180', '180', '360', '360', '540', '540']

# Properties modified on resultLUT
resultLUT.Annotations = ['0', '0', '180', '180', '360', '360', '540', '540', '', '']
resultLUT.IndexedColors = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# Properties modified on resultLUT
resultLUT.Annotations = ['0', '0', '180', '180', '360', '360', '540', '540', '720', '']

# Properties modified on resultLUT
resultLUT.Annotations = ['0', '0', '180', '180', '360', '360', '540', '540', '720', '720']

# Properties modified on resultLUTColorBar
resultLUTColorBar.AddRangeLabels = 0

# Properties modified on resultLUTColorBar
resultLUTColorBar.DrawTickLabels = 0

resultLUTColorBar.Orientation = 'Horizontal'
resultLUTColorBar.Position = [0.2, 0.06]
# current camera placement for renderView1
renderView1.CameraPosition = [7745692.735918371, 1676557.1694394003, -78439.76626821003]
renderView1.CameraFocalPoint = [5998906.499999999, -77042.54687500252, -273030.70312500047]
renderView1.CameraViewUp = [0.7029961754214453, -0.6756164793843724, -0.2221232768691438]
renderView1.CameraParallelScale = 777537.5320806042

renderView1.Background = [1,1,1] #White
viewLayout1 = GetLayout()
# save screenshot
filename_base = "/Users/visu/Glerum/01092016/aspect/EGU/K08_SEMUCB/slices/contour_slab_K08_SEMUCB__"
count = 1
filename = filename_base + str(count) + ".png"

# save screenshot
SaveScreenshot(filename,layout=viewLayout1, magnification=1, quality=100)

#### saving camera placements for all active views

