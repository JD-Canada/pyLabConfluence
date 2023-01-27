# trace generated using paraview version 5.10.0
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
paraviewfoam = FindSource('paraview.foam')

# create a new 'Calculator'
calculator1 = Calculator(registrationName='Calculator1', Input=paraviewfoam)
calculator1.Function = ''

# Properties modified on paraviewfoam
paraviewfoam.CellArrays = ['U', 'UMean', 'UPrime2Mean', 'alpha.dense', 'alpha.denseMean', 'nut', 'p', 'p_rgh', 'turbulencePropertiesR', 'turbulenceProperties:R', 'turbulenceProperties:RMean', 'yPlus']

# Properties modified on calculator1
calculator1.ResultArrayName = 'LESindex'
calculator1.Function = '(0.5*(UPrime2Mean_XX+UPrime2Mean_YY+UPrime2Mean_ZZ ))/(0.5*(UPrime2Mean_XX+UPrime2Mean_YY+UPrime2Mean_ZZ+"turbulenceProperties:RMean_XX"+"turbulenceProperties:RMean_YY"+"turbulenceProperties:RMean_ZZ"))'

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
calculator1Display = Show(calculator1, renderView1, 'UnstructuredGridRepresentation')

# get color transfer function/color map for 'LESindex'
lESindexLUT = GetColorTransferFunction('LESindex')

# get opacity transfer function/opacity map for 'LESindex'
lESindexPWF = GetOpacityTransferFunction('LESindex')

# trace defaults for the display properties.
calculator1Display.Representation = 'Surface'
calculator1Display.ColorArrayName = ['POINTS', 'LESindex']
calculator1Display.LookupTable = lESindexLUT
calculator1Display.SelectTCoordArray = 'None'
calculator1Display.SelectNormalArray = 'None'
calculator1Display.SelectTangentArray = 'None'
calculator1Display.OSPRayScaleArray = 'LESindex'
calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
calculator1Display.SelectOrientationVectors = 'U'
calculator1Display.ScaleFactor = 0.2046410232782364
calculator1Display.SelectScaleArray = 'LESindex'
calculator1Display.GlyphType = 'Arrow'
calculator1Display.GlyphTableIndexArray = 'LESindex'
calculator1Display.GaussianRadius = 0.010232051163911819
calculator1Display.SetScaleArray = ['POINTS', 'LESindex']
calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
calculator1Display.OpacityArray = ['POINTS', 'LESindex']
calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
calculator1Display.DataAxesGrid = 'GridAxesRepresentation'
calculator1Display.PolarAxes = 'PolarAxesRepresentation'
calculator1Display.ScalarOpacityFunction = lESindexPWF
calculator1Display.ScalarOpacityUnitDistance = 0.011713658603110485
calculator1Display.OpacityArrayName = ['POINTS', 'LESindex']

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
calculator1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0000000407565455, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
calculator1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0000000407565455, 1.0, 0.5, 0.0]

# hide data in view
Hide(paraviewfoam, renderView1)

# show color bar/color legend
calculator1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

# get layout
layout1 = GetLayout()

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(2034, 1156)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.CameraPosition = [0.34514957667120005, 0.610089528552282, 1.088605601858783]
renderView1.CameraFocalPoint = [0.28000000119209284, -6.447283160294883e-18, 0.035000000149011584]
renderView1.CameraViewUp = [-0.24173856360007737, -0.8331651362855087, 0.49739151837050993]
renderView1.CameraParallelScale = 0.3818296737353555

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).