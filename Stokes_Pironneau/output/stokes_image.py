# trace generated using paraview version 5.6.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

import os
path_ = os.getcwd()
# create a new 'PVD Reader'
u_mesh0pvd = PVDReader(FileName=path_ +'/u_mesh0.pvd')
u_mesh0pvd.PointArrays = ['x']

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1547, 713]

# show data in view
u_mesh0pvdDisplay = Show(u_mesh0pvd, renderView1)

# trace defaults for the display properties.
u_mesh0pvdDisplay.Representation = 'Surface'
u_mesh0pvdDisplay.AmbientColor = [0.0, 0.0, 0.0]
u_mesh0pvdDisplay.ColorArrayName = [None, '']
u_mesh0pvdDisplay.DiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_mesh0pvdDisplay.BackfaceDiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_mesh0pvdDisplay.OSPRayScaleArray = 'x'
u_mesh0pvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
u_mesh0pvdDisplay.SelectOrientationVectors = 'x'
u_mesh0pvdDisplay.ScaleFactor = 0.1
u_mesh0pvdDisplay.SelectScaleArray = 'None'
u_mesh0pvdDisplay.GlyphType = 'Arrow'
u_mesh0pvdDisplay.GlyphTableIndexArray = 'None'
u_mesh0pvdDisplay.GaussianRadius = 0.005
u_mesh0pvdDisplay.SetScaleArray = ['POINTS', 'x']
u_mesh0pvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
u_mesh0pvdDisplay.OpacityArray = ['POINTS', 'x']
u_mesh0pvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
u_mesh0pvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
u_mesh0pvdDisplay.SelectionCellLabelFontFile = ''
u_mesh0pvdDisplay.SelectionPointLabelFontFile = ''
u_mesh0pvdDisplay.PolarAxes = 'PolarAxesRepresentation'
u_mesh0pvdDisplay.ScalarOpacityUnitDistance = 0.07032022880262821

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
u_mesh0pvdDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
u_mesh0pvdDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
u_mesh0pvdDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
u_mesh0pvdDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
u_mesh0pvdDisplay.DataAxesGrid.XTitleFontFile = ''
u_mesh0pvdDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
u_mesh0pvdDisplay.DataAxesGrid.YTitleFontFile = ''
u_mesh0pvdDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
u_mesh0pvdDisplay.DataAxesGrid.ZTitleFontFile = ''
u_mesh0pvdDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
u_mesh0pvdDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
u_mesh0pvdDisplay.DataAxesGrid.XLabelFontFile = ''
u_mesh0pvdDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
u_mesh0pvdDisplay.DataAxesGrid.YLabelFontFile = ''
u_mesh0pvdDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
u_mesh0pvdDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
u_mesh0pvdDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
u_mesh0pvdDisplay.PolarAxes.PolarAxisTitleFontFile = ''
u_mesh0pvdDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
u_mesh0pvdDisplay.PolarAxes.PolarAxisLabelFontFile = ''
u_mesh0pvdDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
u_mesh0pvdDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
u_mesh0pvdDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
u_mesh0pvdDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 0.5, 10000.0]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# create a new 'PVD Reader'
u_mesh1pvd = PVDReader(FileName=path_ +'/u_mesh1.pvd')
u_mesh1pvd.PointArrays = ['x']

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

# show data in view
u_mesh1pvdDisplay = Show(u_mesh1pvd, renderView1)

# trace defaults for the display properties.
u_mesh1pvdDisplay.Representation = 'Surface'
u_mesh1pvdDisplay.AmbientColor = [0.0, 0.0, 0.0]
u_mesh1pvdDisplay.ColorArrayName = [None, '']
u_mesh1pvdDisplay.DiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_mesh1pvdDisplay.BackfaceDiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_mesh1pvdDisplay.OSPRayScaleArray = 'x'
u_mesh1pvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
u_mesh1pvdDisplay.SelectOrientationVectors = 'x'
u_mesh1pvdDisplay.ScaleFactor = 0.03723139999999999
u_mesh1pvdDisplay.SelectScaleArray = 'None'
u_mesh1pvdDisplay.GlyphType = 'Arrow'
u_mesh1pvdDisplay.GlyphTableIndexArray = 'None'
u_mesh1pvdDisplay.GaussianRadius = 0.0018615699999999997
u_mesh1pvdDisplay.SetScaleArray = ['POINTS', 'x']
u_mesh1pvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
u_mesh1pvdDisplay.OpacityArray = ['POINTS', 'x']
u_mesh1pvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
u_mesh1pvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
u_mesh1pvdDisplay.SelectionCellLabelFontFile = ''
u_mesh1pvdDisplay.SelectionPointLabelFontFile = ''
u_mesh1pvdDisplay.PolarAxes = 'PolarAxesRepresentation'
u_mesh1pvdDisplay.ScalarOpacityUnitDistance = 0.039251038035751876

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
u_mesh1pvdDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
u_mesh1pvdDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
u_mesh1pvdDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
u_mesh1pvdDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
u_mesh1pvdDisplay.DataAxesGrid.XTitleFontFile = ''
u_mesh1pvdDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
u_mesh1pvdDisplay.DataAxesGrid.YTitleFontFile = ''
u_mesh1pvdDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
u_mesh1pvdDisplay.DataAxesGrid.ZTitleFontFile = ''
u_mesh1pvdDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
u_mesh1pvdDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
u_mesh1pvdDisplay.DataAxesGrid.XLabelFontFile = ''
u_mesh1pvdDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
u_mesh1pvdDisplay.DataAxesGrid.YLabelFontFile = ''
u_mesh1pvdDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
u_mesh1pvdDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
u_mesh1pvdDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
u_mesh1pvdDisplay.PolarAxes.PolarAxisTitleFontFile = ''
u_mesh1pvdDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
u_mesh1pvdDisplay.PolarAxes.PolarAxisLabelFontFile = ''
u_mesh1pvdDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
u_mesh1pvdDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
u_mesh1pvdDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
u_mesh1pvdDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(u_mesh1pvdDisplay, ('POINTS', 'x', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
u_mesh1pvdDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
u_mesh1pvdDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'x'
xLUT = GetColorTransferFunction('x')

# get opacity transfer function/opacity map for 'x'
xPWF = GetOpacityTransferFunction('x')

# set active source
SetActiveSource(u_mesh0pvd)

# set scalar coloring
ColorBy(u_mesh0pvdDisplay, ('POINTS', 'x', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
u_mesh0pvdDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
u_mesh0pvdDisplay.SetScalarBarVisibility(renderView1, True)

# change representation type
u_mesh0pvdDisplay.SetRepresentationType('Surface With Edges')

# set active source
SetActiveSource(u_mesh1pvd)

# change representation type
u_mesh1pvdDisplay.SetRepresentationType('Surface With Edges')

# Properties modified on u_mesh1pvdDisplay
u_mesh1pvdDisplay.EdgeColor = [1.0, 0.0, 0.0]

# set active source
SetActiveSource(u_mesh0pvd)

# Properties modified on u_mesh0pvdDisplay
u_mesh0pvdDisplay.EdgeColor = [0.0, 0.0, 0.0]

# get color legend/bar for xLUT in view renderView1
xLUTColorBar = GetScalarBar(xLUT, renderView1)


# change scalar bar placement
xLUTColorBar.Orientation = 'Vertical'
xLUTColorBar.Position = [0.785126050420168, 0.08550573514077156]
xLUTColorBar.ScalarBarLength = 0.8028290854667779
xLUTColorBar.ScalarBarThickness = 5
xLUTColorBar.RangeLabelFormat = '%-#.2f'


# Properties modified on xLUTColorBar
xLUTColorBar.TitleFontSize = 5
xLUTColorBar.LabelFontSize = 3
xLUTColorBar.UseCustomLabels = 1
xLUTColorBar.CustomLabels = [0.5, 1.0, 1.5]


# show color bar/color legend
u_mesh0pvdDisplay.SetScalarBarVisibility(renderView1, False)


# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 0.5, 10000.0]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
renderView1.CameraParallelScale = 0.584385769575659

# save screenshot
SaveScreenshot(path_ + '/initial_stokes.png', renderView1, ImageResolution=[1547, 959])

animationScene1.GoToLast()

# show color bar/color legend
u_mesh0pvdDisplay.SetScalarBarVisibility(renderView1, True)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 0.5, 10000.0]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
renderView1.CameraParallelScale = 0.584385769575659

# Properties modified on xLUTColorBar
xLUTColorBar.Title = 'u'


# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 0.5, 10000.0]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
renderView1.CameraParallelScale = 0.584385769575659

# save screenshot
SaveScreenshot(path_ + '/final_stokes.png', renderView1, ImageResolution=[1547, 959])

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 0.5, 10000.0]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
renderView1.CameraParallelScale = 0.584385769575659

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).


os.system("convert {0}/final_stokes.png -trim {0}/final_stokes.png".format(path_))
os.system("convert {0}/initial_stokes.png -trim {0}/initial_stokes.png".format(path_))
