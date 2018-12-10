### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

import os;
pwd = os.getcwd()

# create a new 'PVD Reader'
u_0pvd = PVDReader(FileName=pwd + '/u_0.pvd')
u_0pvd.PointArrays = ['x']

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1547, 807]

# show data in view
u_0pvdDisplay = Show(u_0pvd, renderView1)

# trace defaults for the display properties.
u_0pvdDisplay.Representation = 'Surface'
u_0pvdDisplay.AmbientColor = [0.0, 0.0, 0.0]
u_0pvdDisplay.ColorArrayName = [None, '']
u_0pvdDisplay.DiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_0pvdDisplay.BackfaceDiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_0pvdDisplay.OSPRayScaleArray = 'x'
u_0pvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
u_0pvdDisplay.SelectOrientationVectors = 'x'
u_0pvdDisplay.ScaleFactor = 0.1
u_0pvdDisplay.SelectScaleArray = 'None'
u_0pvdDisplay.GlyphType = 'Arrow'
u_0pvdDisplay.GlyphTableIndexArray = 'None'
u_0pvdDisplay.GaussianRadius = 0.005
u_0pvdDisplay.SetScaleArray = ['POINTS', 'x']
u_0pvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
u_0pvdDisplay.OpacityArray = ['POINTS', 'x']
u_0pvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
u_0pvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
u_0pvdDisplay.SelectionCellLabelFontFile = ''
u_0pvdDisplay.SelectionPointLabelFontFile = ''
u_0pvdDisplay.PolarAxes = 'PolarAxesRepresentation'
u_0pvdDisplay.ScalarOpacityUnitDistance = 0.04396504373333011

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
u_0pvdDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
u_0pvdDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
u_0pvdDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
u_0pvdDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
u_0pvdDisplay.DataAxesGrid.XTitleFontFile = ''
u_0pvdDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
u_0pvdDisplay.DataAxesGrid.YTitleFontFile = ''
u_0pvdDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
u_0pvdDisplay.DataAxesGrid.ZTitleFontFile = ''
u_0pvdDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
u_0pvdDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
u_0pvdDisplay.DataAxesGrid.XLabelFontFile = ''
u_0pvdDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
u_0pvdDisplay.DataAxesGrid.YLabelFontFile = ''
u_0pvdDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
u_0pvdDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
u_0pvdDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
u_0pvdDisplay.PolarAxes.PolarAxisTitleFontFile = ''
u_0pvdDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
u_0pvdDisplay.PolarAxes.PolarAxisLabelFontFile = ''
u_0pvdDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
u_0pvdDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
u_0pvdDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
u_0pvdDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

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

# set scalar coloring
ColorBy(u_0pvdDisplay, ('POINTS', 'x', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
u_0pvdDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
u_0pvdDisplay.SetScalarBarVisibility(renderView1, False)

# get color transfer function/color map for 'x'
xLUT = GetColorTransferFunction('x')

# get opacity transfer function/opacity map for 'x'
xPWF = GetOpacityTransferFunction('x')

u_pvds = []
for i in range(1,10):
    u_pvds.append(PVDReader(FileName=pwd + '/u_{}.pvd'.format(i)))
    u_pvds[i-1].PointArrays = ['x']

for i in range(1,10):
    # show data in view
    u_ipvdDisplay = Show(u_pvds[i-1], renderView1)

    # trace defaults for the display properties.
    u_ipvdDisplay.Representation = 'Surface'
    u_ipvdDisplay.AmbientColor = [0.0, 0.0, 0.0]
    u_ipvdDisplay.ColorArrayName = [None, '']
    u_ipvdDisplay.DiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
    u_ipvdDisplay.BackfaceDiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
    u_ipvdDisplay.OSPRayScaleArray = 'x'
    u_ipvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    u_ipvdDisplay.SelectOrientationVectors = 'x'
    u_ipvdDisplay.ScaleFactor = 0.0214
    u_ipvdDisplay.SelectScaleArray = 'None'
    u_ipvdDisplay.GlyphType = 'Arrow'
    u_ipvdDisplay.GlyphTableIndexArray = 'None'
    u_ipvdDisplay.GaussianRadius = 0.0010699999999999998
    u_ipvdDisplay.SetScaleArray = ['POINTS', 'x']
    u_ipvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    u_ipvdDisplay.OpacityArray = ['POINTS', 'x']
    u_ipvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    u_ipvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
    u_ipvdDisplay.SelectionCellLabelFontFile = ''
    u_ipvdDisplay.SelectionPointLabelFontFile = ''
    u_ipvdDisplay.PolarAxes = 'PolarAxesRepresentation'
    u_ipvdDisplay.ScalarOpacityUnitDistance = 0.020385860120059825
    
    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    u_ipvdDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]
    
    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    u_ipvdDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]
    
    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    u_ipvdDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]
    
    # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
    u_ipvdDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
    u_ipvdDisplay.DataAxesGrid.XTitleFontFile = ''
    u_ipvdDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
    u_ipvdDisplay.DataAxesGrid.YTitleFontFile = ''
    u_ipvdDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
    u_ipvdDisplay.DataAxesGrid.ZTitleFontFile = ''
    u_ipvdDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
    u_ipvdDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
    u_ipvdDisplay.DataAxesGrid.XLabelFontFile = ''
    u_ipvdDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
    u_ipvdDisplay.DataAxesGrid.YLabelFontFile = ''
    u_ipvdDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
    u_ipvdDisplay.DataAxesGrid.ZLabelFontFile = ''
    
    # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
    u_ipvdDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
    u_ipvdDisplay.PolarAxes.PolarAxisTitleFontFile = ''
    u_ipvdDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
    u_ipvdDisplay.PolarAxes.PolarAxisLabelFontFile = ''
    u_ipvdDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
    u_ipvdDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
    u_ipvdDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
    u_ipvdDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''


    # update the view to ensure updated data information
    #renderView1.Update()

    # set active source
    SetActiveSource(u_pvds[i-1])


    # change representation type
    u_ipvdDisplay.SetRepresentationType('Wireframe')

    # change solid color
    u_ipvdDisplay.AmbientColor = [1.0, 0.0, 0.0]


    # rescale color and/or opacity maps used to include current data range
    #u_ipvdDisplay.RescaleTransferFunctionToDataRange(True, False)

    # show color bar/color legend
    # u_ipvdDisplay.SetScalarBarVisibility(renderView1, True)
    # u_ipvdDisplay.RescaleTransferFunctionToDataRange(False, True)

    # turn off scalar coloring
    ColorBy(u_ipvdDisplay, None)
renderView1.Update()

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(xLUT, renderView1)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 0.5, 10000.0]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
renderView1.CameraParallelScale = 0.7071067811865476


# save screenshot
SaveScreenshot(pwd + '/initial_rotation.png', renderView1, ImageResolution=[1547, 807])

animationScene1.GoToLast()

# set active source
SetActiveSource(u_0pvd)

# get color legend/bar for xLUT in view renderView1
xLUTColorBar = GetScalarBar(xLUT, renderView1)

# Properties modified on xLUTColorBar
xLUTColorBar.AutoOrient = 0
xLUTColorBar.Orientation = 'Vertical'
xLUTColorBar.ComponentTitle = 'Magnitude'
# Properties modified on xLUTColorBar
xLUTColorBar.Position = [0.68, 0.15]
xLUTColorBar.ScalarBarLength = 0.67
xLUTColorBar.ScalarBarThickness = 5
xLUTColorBar.UseCustomLabels = 1
xLUTColorBar.CustomLabels = [0,0.2,0.4,0.625]
xLUTColorBar.RangeLabelFormat = '%-#.2f'
xLUTColorBar.TitleFontSize = 5
xLUTColorBar.LabelFontSize = 3
xLUTColorBar.Title = 'u'
for prop in xLUTColorBar:
    print type(prop), prop.GetXMLLabel()
# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5245340642790872, 0.5052572994883759, 10000.0]
renderView1.CameraFocalPoint = [0.5245340642790872, 0.5052572994883759, 0.0]
renderView1.CameraParallelScale = 0.7071067811865476

SaveScreenshot(pwd + '/final_rotation.png', renderView1, ImageResolution=[1547, 807])

os.system("convert {0} -trim {0}".format(pwd+'/final_rotation.png'))
os.system("convert {0} -trim {0}".format(pwd+'/initial_rotation.png'))
          
#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
