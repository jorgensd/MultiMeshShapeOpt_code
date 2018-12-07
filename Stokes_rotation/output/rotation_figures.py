# trace generated using paraview version 5.6.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'PVD Reader'
u_0pvd = PVDReader(FileName='/home/dokken/Documents/MultiMeshShapeOpt_code/Stokes_rotation/output/u_0.pvd')
u_0pvd.PointArrays = ['x']

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1547, 835]

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
renderView1.CameraPosition = [0.5, 0.5, 10000.0]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# create a new 'PVD Reader'
u_1pvd = PVDReader(FileName='/home/dokken/Documents/MultiMeshShapeOpt_code/Stokes_rotation/output/u_1.pvd')
u_1pvd.PointArrays = ['x']

# create a new 'PVD Reader'
u_2pvd = PVDReader(FileName='/home/dokken/Documents/MultiMeshShapeOpt_code/Stokes_rotation/output/u_2.pvd')
u_2pvd.PointArrays = ['x']

# create a new 'PVD Reader'
u_3pvd = PVDReader(FileName='/home/dokken/Documents/MultiMeshShapeOpt_code/Stokes_rotation/output/u_3.pvd')
u_3pvd.PointArrays = ['x']

# create a new 'PVD Reader'
u_4pvd = PVDReader(FileName='/home/dokken/Documents/MultiMeshShapeOpt_code/Stokes_rotation/output/u_4.pvd')
u_4pvd.PointArrays = ['x']

# create a new 'PVD Reader'
u_5pvd = PVDReader(FileName='/home/dokken/Documents/MultiMeshShapeOpt_code/Stokes_rotation/output/u_5.pvd')
u_5pvd.PointArrays = ['x']

# create a new 'PVD Reader'
u_6pvd = PVDReader(FileName='/home/dokken/Documents/MultiMeshShapeOpt_code/Stokes_rotation/output/u_6.pvd')
u_6pvd.PointArrays = ['x']

# create a new 'PVD Reader'
u_7pvd = PVDReader(FileName='/home/dokken/Documents/MultiMeshShapeOpt_code/Stokes_rotation/output/u_7.pvd')
u_7pvd.PointArrays = ['x']

# create a new 'PVD Reader'
u_8pvd = PVDReader(FileName='/home/dokken/Documents/MultiMeshShapeOpt_code/Stokes_rotation/output/u_8.pvd')
u_8pvd.PointArrays = ['x']

# create a new 'PVD Reader'
u_9pvd = PVDReader(FileName='/home/dokken/Documents/MultiMeshShapeOpt_code/Stokes_rotation/output/u_9.pvd')
u_9pvd.PointArrays = ['x']

# show data in view
u_1pvdDisplay = Show(u_1pvd, renderView1)

# trace defaults for the display properties.
u_1pvdDisplay.Representation = 'Surface'
u_1pvdDisplay.AmbientColor = [0.0, 0.0, 0.0]
u_1pvdDisplay.ColorArrayName = [None, '']
u_1pvdDisplay.DiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_1pvdDisplay.BackfaceDiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_1pvdDisplay.OSPRayScaleArray = 'x'
u_1pvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
u_1pvdDisplay.SelectOrientationVectors = 'x'
u_1pvdDisplay.ScaleFactor = 0.021400000000000002
u_1pvdDisplay.SelectScaleArray = 'None'
u_1pvdDisplay.GlyphType = 'Arrow'
u_1pvdDisplay.GlyphTableIndexArray = 'None'
u_1pvdDisplay.GaussianRadius = 0.00107
u_1pvdDisplay.SetScaleArray = ['POINTS', 'x']
u_1pvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
u_1pvdDisplay.OpacityArray = ['POINTS', 'x']
u_1pvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
u_1pvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
u_1pvdDisplay.SelectionCellLabelFontFile = ''
u_1pvdDisplay.SelectionPointLabelFontFile = ''
u_1pvdDisplay.PolarAxes = 'PolarAxesRepresentation'
u_1pvdDisplay.ScalarOpacityUnitDistance = 0.02038586012005983

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
u_1pvdDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
u_1pvdDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
u_1pvdDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
u_1pvdDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
u_1pvdDisplay.DataAxesGrid.XTitleFontFile = ''
u_1pvdDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
u_1pvdDisplay.DataAxesGrid.YTitleFontFile = ''
u_1pvdDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
u_1pvdDisplay.DataAxesGrid.ZTitleFontFile = ''
u_1pvdDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
u_1pvdDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
u_1pvdDisplay.DataAxesGrid.XLabelFontFile = ''
u_1pvdDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
u_1pvdDisplay.DataAxesGrid.YLabelFontFile = ''
u_1pvdDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
u_1pvdDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
u_1pvdDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
u_1pvdDisplay.PolarAxes.PolarAxisTitleFontFile = ''
u_1pvdDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
u_1pvdDisplay.PolarAxes.PolarAxisLabelFontFile = ''
u_1pvdDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
u_1pvdDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
u_1pvdDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
u_1pvdDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data in view
u_3pvdDisplay = Show(u_3pvd, renderView1)

# trace defaults for the display properties.
u_3pvdDisplay.Representation = 'Surface'
u_3pvdDisplay.AmbientColor = [0.0, 0.0, 0.0]
u_3pvdDisplay.ColorArrayName = [None, '']
u_3pvdDisplay.DiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_3pvdDisplay.BackfaceDiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_3pvdDisplay.OSPRayScaleArray = 'x'
u_3pvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
u_3pvdDisplay.SelectOrientationVectors = 'x'
u_3pvdDisplay.ScaleFactor = 0.0214
u_3pvdDisplay.SelectScaleArray = 'None'
u_3pvdDisplay.GlyphType = 'Arrow'
u_3pvdDisplay.GlyphTableIndexArray = 'None'
u_3pvdDisplay.GaussianRadius = 0.0010699999999999998
u_3pvdDisplay.SetScaleArray = ['POINTS', 'x']
u_3pvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
u_3pvdDisplay.OpacityArray = ['POINTS', 'x']
u_3pvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
u_3pvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
u_3pvdDisplay.SelectionCellLabelFontFile = ''
u_3pvdDisplay.SelectionPointLabelFontFile = ''
u_3pvdDisplay.PolarAxes = 'PolarAxesRepresentation'
u_3pvdDisplay.ScalarOpacityUnitDistance = 0.020385860120059825

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
u_3pvdDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
u_3pvdDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
u_3pvdDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
u_3pvdDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
u_3pvdDisplay.DataAxesGrid.XTitleFontFile = ''
u_3pvdDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
u_3pvdDisplay.DataAxesGrid.YTitleFontFile = ''
u_3pvdDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
u_3pvdDisplay.DataAxesGrid.ZTitleFontFile = ''
u_3pvdDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
u_3pvdDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
u_3pvdDisplay.DataAxesGrid.XLabelFontFile = ''
u_3pvdDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
u_3pvdDisplay.DataAxesGrid.YLabelFontFile = ''
u_3pvdDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
u_3pvdDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
u_3pvdDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
u_3pvdDisplay.PolarAxes.PolarAxisTitleFontFile = ''
u_3pvdDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
u_3pvdDisplay.PolarAxes.PolarAxisLabelFontFile = ''
u_3pvdDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
u_3pvdDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
u_3pvdDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
u_3pvdDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data in view
u_7pvdDisplay = Show(u_7pvd, renderView1)

# trace defaults for the display properties.
u_7pvdDisplay.Representation = 'Surface'
u_7pvdDisplay.AmbientColor = [0.0, 0.0, 0.0]
u_7pvdDisplay.ColorArrayName = [None, '']
u_7pvdDisplay.DiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_7pvdDisplay.BackfaceDiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_7pvdDisplay.OSPRayScaleArray = 'x'
u_7pvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
u_7pvdDisplay.SelectOrientationVectors = 'x'
u_7pvdDisplay.ScaleFactor = 0.021400000000000002
u_7pvdDisplay.SelectScaleArray = 'None'
u_7pvdDisplay.GlyphType = 'Arrow'
u_7pvdDisplay.GlyphTableIndexArray = 'None'
u_7pvdDisplay.GaussianRadius = 0.00107
u_7pvdDisplay.SetScaleArray = ['POINTS', 'x']
u_7pvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
u_7pvdDisplay.OpacityArray = ['POINTS', 'x']
u_7pvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
u_7pvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
u_7pvdDisplay.SelectionCellLabelFontFile = ''
u_7pvdDisplay.SelectionPointLabelFontFile = ''
u_7pvdDisplay.PolarAxes = 'PolarAxesRepresentation'
u_7pvdDisplay.ScalarOpacityUnitDistance = 0.02038586012005983

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
u_7pvdDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
u_7pvdDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
u_7pvdDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
u_7pvdDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
u_7pvdDisplay.DataAxesGrid.XTitleFontFile = ''
u_7pvdDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
u_7pvdDisplay.DataAxesGrid.YTitleFontFile = ''
u_7pvdDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
u_7pvdDisplay.DataAxesGrid.ZTitleFontFile = ''
u_7pvdDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
u_7pvdDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
u_7pvdDisplay.DataAxesGrid.XLabelFontFile = ''
u_7pvdDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
u_7pvdDisplay.DataAxesGrid.YLabelFontFile = ''
u_7pvdDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
u_7pvdDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
u_7pvdDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
u_7pvdDisplay.PolarAxes.PolarAxisTitleFontFile = ''
u_7pvdDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
u_7pvdDisplay.PolarAxes.PolarAxisLabelFontFile = ''
u_7pvdDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
u_7pvdDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
u_7pvdDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
u_7pvdDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data in view
u_9pvdDisplay = Show(u_9pvd, renderView1)

# trace defaults for the display properties.
u_9pvdDisplay.Representation = 'Surface'
u_9pvdDisplay.AmbientColor = [0.0, 0.0, 0.0]
u_9pvdDisplay.ColorArrayName = [None, '']
u_9pvdDisplay.DiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_9pvdDisplay.BackfaceDiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_9pvdDisplay.OSPRayScaleArray = 'x'
u_9pvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
u_9pvdDisplay.SelectOrientationVectors = 'x'
u_9pvdDisplay.ScaleFactor = 0.0214
u_9pvdDisplay.SelectScaleArray = 'None'
u_9pvdDisplay.GlyphType = 'Arrow'
u_9pvdDisplay.GlyphTableIndexArray = 'None'
u_9pvdDisplay.GaussianRadius = 0.0010699999999999998
u_9pvdDisplay.SetScaleArray = ['POINTS', 'x']
u_9pvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
u_9pvdDisplay.OpacityArray = ['POINTS', 'x']
u_9pvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
u_9pvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
u_9pvdDisplay.SelectionCellLabelFontFile = ''
u_9pvdDisplay.SelectionPointLabelFontFile = ''
u_9pvdDisplay.PolarAxes = 'PolarAxesRepresentation'
u_9pvdDisplay.ScalarOpacityUnitDistance = 0.020385860120059825

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
u_9pvdDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
u_9pvdDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
u_9pvdDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
u_9pvdDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
u_9pvdDisplay.DataAxesGrid.XTitleFontFile = ''
u_9pvdDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
u_9pvdDisplay.DataAxesGrid.YTitleFontFile = ''
u_9pvdDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
u_9pvdDisplay.DataAxesGrid.ZTitleFontFile = ''
u_9pvdDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
u_9pvdDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
u_9pvdDisplay.DataAxesGrid.XLabelFontFile = ''
u_9pvdDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
u_9pvdDisplay.DataAxesGrid.YLabelFontFile = ''
u_9pvdDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
u_9pvdDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
u_9pvdDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
u_9pvdDisplay.PolarAxes.PolarAxisTitleFontFile = ''
u_9pvdDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
u_9pvdDisplay.PolarAxes.PolarAxisLabelFontFile = ''
u_9pvdDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
u_9pvdDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
u_9pvdDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
u_9pvdDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data in view
u_2pvdDisplay = Show(u_2pvd, renderView1)

# trace defaults for the display properties.
u_2pvdDisplay.Representation = 'Surface'
u_2pvdDisplay.AmbientColor = [0.0, 0.0, 0.0]
u_2pvdDisplay.ColorArrayName = [None, '']
u_2pvdDisplay.DiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_2pvdDisplay.BackfaceDiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_2pvdDisplay.OSPRayScaleArray = 'x'
u_2pvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
u_2pvdDisplay.SelectOrientationVectors = 'x'
u_2pvdDisplay.ScaleFactor = 0.0214
u_2pvdDisplay.SelectScaleArray = 'None'
u_2pvdDisplay.GlyphType = 'Arrow'
u_2pvdDisplay.GlyphTableIndexArray = 'None'
u_2pvdDisplay.GaussianRadius = 0.0010699999999999998
u_2pvdDisplay.SetScaleArray = ['POINTS', 'x']
u_2pvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
u_2pvdDisplay.OpacityArray = ['POINTS', 'x']
u_2pvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
u_2pvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
u_2pvdDisplay.SelectionCellLabelFontFile = ''
u_2pvdDisplay.SelectionPointLabelFontFile = ''
u_2pvdDisplay.PolarAxes = 'PolarAxesRepresentation'
u_2pvdDisplay.ScalarOpacityUnitDistance = 0.020385860120059825

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
u_2pvdDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
u_2pvdDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
u_2pvdDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
u_2pvdDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
u_2pvdDisplay.DataAxesGrid.XTitleFontFile = ''
u_2pvdDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
u_2pvdDisplay.DataAxesGrid.YTitleFontFile = ''
u_2pvdDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
u_2pvdDisplay.DataAxesGrid.ZTitleFontFile = ''
u_2pvdDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
u_2pvdDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
u_2pvdDisplay.DataAxesGrid.XLabelFontFile = ''
u_2pvdDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
u_2pvdDisplay.DataAxesGrid.YLabelFontFile = ''
u_2pvdDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
u_2pvdDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
u_2pvdDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
u_2pvdDisplay.PolarAxes.PolarAxisTitleFontFile = ''
u_2pvdDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
u_2pvdDisplay.PolarAxes.PolarAxisLabelFontFile = ''
u_2pvdDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
u_2pvdDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
u_2pvdDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
u_2pvdDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data in view
u_4pvdDisplay = Show(u_4pvd, renderView1)

# trace defaults for the display properties.
u_4pvdDisplay.Representation = 'Surface'
u_4pvdDisplay.AmbientColor = [0.0, 0.0, 0.0]
u_4pvdDisplay.ColorArrayName = [None, '']
u_4pvdDisplay.DiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_4pvdDisplay.BackfaceDiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_4pvdDisplay.OSPRayScaleArray = 'x'
u_4pvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
u_4pvdDisplay.SelectOrientationVectors = 'x'
u_4pvdDisplay.ScaleFactor = 0.021400000000000002
u_4pvdDisplay.SelectScaleArray = 'None'
u_4pvdDisplay.GlyphType = 'Arrow'
u_4pvdDisplay.GlyphTableIndexArray = 'None'
u_4pvdDisplay.GaussianRadius = 0.00107
u_4pvdDisplay.SetScaleArray = ['POINTS', 'x']
u_4pvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
u_4pvdDisplay.OpacityArray = ['POINTS', 'x']
u_4pvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
u_4pvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
u_4pvdDisplay.SelectionCellLabelFontFile = ''
u_4pvdDisplay.SelectionPointLabelFontFile = ''
u_4pvdDisplay.PolarAxes = 'PolarAxesRepresentation'
u_4pvdDisplay.ScalarOpacityUnitDistance = 0.020385860120059825

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
u_4pvdDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
u_4pvdDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
u_4pvdDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
u_4pvdDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
u_4pvdDisplay.DataAxesGrid.XTitleFontFile = ''
u_4pvdDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
u_4pvdDisplay.DataAxesGrid.YTitleFontFile = ''
u_4pvdDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
u_4pvdDisplay.DataAxesGrid.ZTitleFontFile = ''
u_4pvdDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
u_4pvdDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
u_4pvdDisplay.DataAxesGrid.XLabelFontFile = ''
u_4pvdDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
u_4pvdDisplay.DataAxesGrid.YLabelFontFile = ''
u_4pvdDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
u_4pvdDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
u_4pvdDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
u_4pvdDisplay.PolarAxes.PolarAxisTitleFontFile = ''
u_4pvdDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
u_4pvdDisplay.PolarAxes.PolarAxisLabelFontFile = ''
u_4pvdDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
u_4pvdDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
u_4pvdDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
u_4pvdDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data in view
u_6pvdDisplay = Show(u_6pvd, renderView1)

# trace defaults for the display properties.
u_6pvdDisplay.Representation = 'Surface'
u_6pvdDisplay.AmbientColor = [0.0, 0.0, 0.0]
u_6pvdDisplay.ColorArrayName = [None, '']
u_6pvdDisplay.DiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_6pvdDisplay.BackfaceDiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_6pvdDisplay.OSPRayScaleArray = 'x'
u_6pvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
u_6pvdDisplay.SelectOrientationVectors = 'x'
u_6pvdDisplay.ScaleFactor = 0.0214
u_6pvdDisplay.SelectScaleArray = 'None'
u_6pvdDisplay.GlyphType = 'Arrow'
u_6pvdDisplay.GlyphTableIndexArray = 'None'
u_6pvdDisplay.GaussianRadius = 0.0010699999999999998
u_6pvdDisplay.SetScaleArray = ['POINTS', 'x']
u_6pvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
u_6pvdDisplay.OpacityArray = ['POINTS', 'x']
u_6pvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
u_6pvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
u_6pvdDisplay.SelectionCellLabelFontFile = ''
u_6pvdDisplay.SelectionPointLabelFontFile = ''
u_6pvdDisplay.PolarAxes = 'PolarAxesRepresentation'
u_6pvdDisplay.ScalarOpacityUnitDistance = 0.020385860120059825

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
u_6pvdDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
u_6pvdDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
u_6pvdDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
u_6pvdDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
u_6pvdDisplay.DataAxesGrid.XTitleFontFile = ''
u_6pvdDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
u_6pvdDisplay.DataAxesGrid.YTitleFontFile = ''
u_6pvdDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
u_6pvdDisplay.DataAxesGrid.ZTitleFontFile = ''
u_6pvdDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
u_6pvdDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
u_6pvdDisplay.DataAxesGrid.XLabelFontFile = ''
u_6pvdDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
u_6pvdDisplay.DataAxesGrid.YLabelFontFile = ''
u_6pvdDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
u_6pvdDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
u_6pvdDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
u_6pvdDisplay.PolarAxes.PolarAxisTitleFontFile = ''
u_6pvdDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
u_6pvdDisplay.PolarAxes.PolarAxisLabelFontFile = ''
u_6pvdDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
u_6pvdDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
u_6pvdDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
u_6pvdDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data in view
u_8pvdDisplay = Show(u_8pvd, renderView1)

# trace defaults for the display properties.
u_8pvdDisplay.Representation = 'Surface'
u_8pvdDisplay.AmbientColor = [0.0, 0.0, 0.0]
u_8pvdDisplay.ColorArrayName = [None, '']
u_8pvdDisplay.DiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_8pvdDisplay.BackfaceDiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_8pvdDisplay.OSPRayScaleArray = 'x'
u_8pvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
u_8pvdDisplay.SelectOrientationVectors = 'x'
u_8pvdDisplay.ScaleFactor = 0.0214
u_8pvdDisplay.SelectScaleArray = 'None'
u_8pvdDisplay.GlyphType = 'Arrow'
u_8pvdDisplay.GlyphTableIndexArray = 'None'
u_8pvdDisplay.GaussianRadius = 0.0010699999999999998
u_8pvdDisplay.SetScaleArray = ['POINTS', 'x']
u_8pvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
u_8pvdDisplay.OpacityArray = ['POINTS', 'x']
u_8pvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
u_8pvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
u_8pvdDisplay.SelectionCellLabelFontFile = ''
u_8pvdDisplay.SelectionPointLabelFontFile = ''
u_8pvdDisplay.PolarAxes = 'PolarAxesRepresentation'
u_8pvdDisplay.ScalarOpacityUnitDistance = 0.020385860120059825

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
u_8pvdDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
u_8pvdDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
u_8pvdDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
u_8pvdDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
u_8pvdDisplay.DataAxesGrid.XTitleFontFile = ''
u_8pvdDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
u_8pvdDisplay.DataAxesGrid.YTitleFontFile = ''
u_8pvdDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
u_8pvdDisplay.DataAxesGrid.ZTitleFontFile = ''
u_8pvdDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
u_8pvdDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
u_8pvdDisplay.DataAxesGrid.XLabelFontFile = ''
u_8pvdDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
u_8pvdDisplay.DataAxesGrid.YLabelFontFile = ''
u_8pvdDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
u_8pvdDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
u_8pvdDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
u_8pvdDisplay.PolarAxes.PolarAxisTitleFontFile = ''
u_8pvdDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
u_8pvdDisplay.PolarAxes.PolarAxisLabelFontFile = ''
u_8pvdDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
u_8pvdDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
u_8pvdDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
u_8pvdDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data in view
u_5pvdDisplay = Show(u_5pvd, renderView1)

# trace defaults for the display properties.
u_5pvdDisplay.Representation = 'Surface'
u_5pvdDisplay.AmbientColor = [0.0, 0.0, 0.0]
u_5pvdDisplay.ColorArrayName = [None, '']
u_5pvdDisplay.DiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_5pvdDisplay.BackfaceDiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_5pvdDisplay.OSPRayScaleArray = 'x'
u_5pvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
u_5pvdDisplay.SelectOrientationVectors = 'x'
u_5pvdDisplay.ScaleFactor = 0.0214
u_5pvdDisplay.SelectScaleArray = 'None'
u_5pvdDisplay.GlyphType = 'Arrow'
u_5pvdDisplay.GlyphTableIndexArray = 'None'
u_5pvdDisplay.GaussianRadius = 0.0010699999999999998
u_5pvdDisplay.SetScaleArray = ['POINTS', 'x']
u_5pvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
u_5pvdDisplay.OpacityArray = ['POINTS', 'x']
u_5pvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
u_5pvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
u_5pvdDisplay.SelectionCellLabelFontFile = ''
u_5pvdDisplay.SelectionPointLabelFontFile = ''
u_5pvdDisplay.PolarAxes = 'PolarAxesRepresentation'
u_5pvdDisplay.ScalarOpacityUnitDistance = 0.020385860120059825

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
u_5pvdDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
u_5pvdDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
u_5pvdDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
u_5pvdDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
u_5pvdDisplay.DataAxesGrid.XTitleFontFile = ''
u_5pvdDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
u_5pvdDisplay.DataAxesGrid.YTitleFontFile = ''
u_5pvdDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
u_5pvdDisplay.DataAxesGrid.ZTitleFontFile = ''
u_5pvdDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
u_5pvdDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
u_5pvdDisplay.DataAxesGrid.XLabelFontFile = ''
u_5pvdDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
u_5pvdDisplay.DataAxesGrid.YLabelFontFile = ''
u_5pvdDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
u_5pvdDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
u_5pvdDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
u_5pvdDisplay.PolarAxes.PolarAxisTitleFontFile = ''
u_5pvdDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
u_5pvdDisplay.PolarAxes.PolarAxisLabelFontFile = ''
u_5pvdDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
u_5pvdDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
u_5pvdDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
u_5pvdDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(u_1pvd)

# set active source
SetActiveSource(u_0pvd)

# set scalar coloring
ColorBy(u_0pvdDisplay, ('POINTS', 'x', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
u_0pvdDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
u_0pvdDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'x'
xLUT = GetColorTransferFunction('x')
xLUT.RGBPoints = [0.0, 0.267004, 0.004874, 0.329415, 0.0024512499999999977, 0.26851, 0.009605, 0.335427, 0.0049018749999999904, 0.269944, 0.014625, 0.341379, 0.007353124999999988, 0.271305, 0.019942, 0.347269, 0.009803749999999981, 0.272594, 0.025563, 0.353093, 0.01225499999999998, 0.273809, 0.031497, 0.358853, 0.014705624999999974, 0.274952, 0.037752, 0.364543, 0.01715687499999997, 0.276022, 0.044167, 0.370164, 0.019608124999999966, 0.277018, 0.050344, 0.375715, 0.022058749999999964, 0.277941, 0.056324, 0.381191, 0.02450999999999996, 0.278791, 0.062145, 0.386592, 0.026960624999999988, 0.279566, 0.067836, 0.391917, 0.029411874999999952, 0.280267, 0.073417, 0.397163, 0.03186249999999998, 0.280894, 0.078907, 0.402329, 0.034313749999999976, 0.281446, 0.08432, 0.407414, 0.03676499999999998, 0.281924, 0.089666, 0.412415, 0.03921562500000001, 0.282327, 0.094955, 0.417331, 0.041666874999999895, 0.282656, 0.100196, 0.42216, 0.044117499999999775, 0.28291, 0.105393, 0.426902, 0.04656875000000004, 0.283091, 0.110553, 0.431554, 0.049019374999999914, 0.283197, 0.11568, 0.436115, 0.0514706249999998, 0.283229, 0.120777, 0.440584, 0.053921875000000064, 0.283187, 0.125848, 0.44496, 0.056372499999999943, 0.283072, 0.130895, 0.449241, 0.05882374999999983, 0.282884, 0.13592, 0.453427, 0.061274375000000075, 0.282623, 0.140926, 0.457517, 0.06372562499999998, 0.28229, 0.145912, 0.46151, 0.06617624999999985, 0.281887, 0.150881, 0.465405, 0.0686275000000001, 0.281412, 0.155834, 0.469201, 0.07107812499999999, 0.280868, 0.160771, 0.472899, 0.07352937499999987, 0.280255, 0.165693, 0.476498, 0.07598062499999976, 0.279574, 0.170599, 0.479997, 0.07843125000000002, 0.278826, 0.17549, 0.483397, 0.08088249999999989, 0.278012, 0.180367, 0.486697, 0.08333312499999979, 0.277134, 0.185228, 0.489898, 0.08578437500000004, 0.276194, 0.190074, 0.493001, 0.08823499999999992, 0.275191, 0.194905, 0.496005, 0.09068624999999982, 0.274128, 0.199721, 0.498911, 0.09313750000000008, 0.273006, 0.20452, 0.501721, 0.09558812499999995, 0.271828, 0.209303, 0.504434, 0.09803937499999983, 0.270595, 0.214069, 0.507052, 0.10048999999999972, 0.269308, 0.218818, 0.509577, 0.10294124999999998, 0.267968, 0.223549, 0.512008, 0.10539187499999986, 0.26658, 0.228262, 0.514349, 0.10784312499999975, 0.265145, 0.232956, 0.516599, 0.110294375, 0.263663, 0.237631, 0.518762, 0.11274499999999989, 0.262138, 0.242286, 0.520837, 0.11519624999999976, 0.260571, 0.246922, 0.522828, 0.11764687500000001, 0.258965, 0.251537, 0.524736, 0.12009812499999992, 0.257322, 0.25613, 0.526563, 0.12254874999999979, 0.255645, 0.260703, 0.528312, 0.12500000000000003, 0.253935, 0.265254, 0.529983, 0.12745124999999996, 0.252194, 0.269783, 0.531579, 0.12990187499999983, 0.250425, 0.27429, 0.533103, 0.1323531249999997, 0.248629, 0.278775, 0.534556, 0.13480374999999994, 0.246811, 0.283237, 0.535941, 0.13725499999999985, 0.244972, 0.287675, 0.53726, 0.13970562499999972, 0.243113, 0.292092, 0.538516, 0.142156875, 0.241237, 0.296485, 0.539709, 0.14460812499999987, 0.239346, 0.300855, 0.540844, 0.14705874999999974, 0.237441, 0.305202, 0.541921, 0.14951, 0.235526, 0.309527, 0.542944, 0.1519606249999999, 0.233603, 0.313828, 0.543914, 0.15441187499999978, 0.231674, 0.318106, 0.544834, 0.15686249999999965, 0.229739, 0.322361, 0.545706, 0.15931374999999992, 0.227802, 0.326594, 0.546532, 0.16176499999999977, 0.225863, 0.330805, 0.547314, 0.16421562499999967, 0.223925, 0.334994, 0.548053, 0.16666687499999994, 0.221989, 0.339161, 0.548752, 0.1691174999999998, 0.220057, 0.343307, 0.549413, 0.17156874999999974, 0.21813, 0.347432, 0.550038, 0.17401937499999995, 0.21621, 0.351535, 0.550627, 0.17647062499999983, 0.214298, 0.355619, 0.551184, 0.17892187499999973, 0.212395, 0.359683, 0.55171, 0.1813725, 0.210503, 0.363727, 0.552206, 0.1838237499999999, 0.208623, 0.367752, 0.552675, 0.18627437499999974, 0.206756, 0.371758, 0.553117, 0.18872562499999965, 0.204903, 0.375746, 0.553533, 0.1911762499999999, 0.203063, 0.379716, 0.553925, 0.19362749999999976, 0.201239, 0.38367, 0.554294, 0.19607812499999963, 0.19943, 0.387607, 0.554642, 0.19852937499999995, 0.197636, 0.391528, 0.554969, 0.2009806249999998, 0.19586, 0.395433, 0.555276, 0.2034312499999997, 0.1941, 0.399323, 0.555565, 0.20588249999999997, 0.192357, 0.403199, 0.555836, 0.2083331249999998, 0.190631, 0.407061, 0.556089, 0.21078437499999972, 0.188923, 0.41091, 0.556326, 0.21323499999999995, 0.187231, 0.414746, 0.556547, 0.21568624999999989, 0.185556, 0.41857, 0.556753, 0.21813749999999973, 0.183898, 0.422383, 0.556944, 0.2205881249999996, 0.182256, 0.426184, 0.55712, 0.2230393749999999, 0.180629, 0.429975, 0.557282, 0.22548999999999977, 0.179019, 0.433756, 0.55743, 0.22794124999999965, 0.177423, 0.437527, 0.557565, 0.23039187499999988, 0.175841, 0.44129, 0.557685, 0.2328431249999998, 0.174274, 0.445044, 0.557792, 0.23529437499999967, 0.172719, 0.448791, 0.557885, 0.23774499999999993, 0.171176, 0.45253, 0.557965, 0.24019624999999983, 0.169646, 0.456262, 0.55803, 0.2426468749999997, 0.168126, 0.459988, 0.558082, 0.24509812499999956, 0.166617, 0.463708, 0.558119, 0.24754874999999985, 0.165117, 0.467423, 0.558141, 0.24999999999999972, 0.163625, 0.471133, 0.558148, 0.2524512499999997, 0.162142, 0.474838, 0.55814, 0.25490187499999983, 0.160665, 0.47854, 0.558115, 0.2573531249999997, 0.159194, 0.482237, 0.558073, 0.25980374999999967, 0.157729, 0.485932, 0.558013, 0.26225499999999985, 0.15627, 0.489624, 0.557936, 0.2647056249999998, 0.154815, 0.493313, 0.55784, 0.2671568749999997, 0.153364, 0.497, 0.557724, 0.26960812499999987, 0.151918, 0.500685, 0.557587, 0.2720587499999998, 0.150476, 0.504369, 0.55743, 0.2745099999999997, 0.149039, 0.508051, 0.55725, 0.27696062499999957, 0.147607, 0.511733, 0.557049, 0.2794118749999998, 0.14618, 0.515413, 0.556823, 0.2818624999999997, 0.144759, 0.519093, 0.556572, 0.2843137499999996, 0.143343, 0.522773, 0.556295, 0.2867649999999998, 0.141935, 0.526453, 0.555991, 0.28921562499999975, 0.140536, 0.530132, 0.555659, 0.2916668749999996, 0.139147, 0.533812, 0.555298, 0.2941174999999998, 0.13777, 0.537492, 0.554906, 0.29656874999999977, 0.136408, 0.541173, 0.554483, 0.29901937499999964, 0.135066, 0.544853, 0.554029, 0.30147062499999955, 0.133743, 0.548535, 0.553541, 0.3039218749999998, 0.132444, 0.552216, 0.553018, 0.30637249999999966, 0.131172, 0.555899, 0.552459, 0.30882374999999956, 0.129933, 0.559582, 0.551864, 0.3112743749999998, 0.128729, 0.563265, 0.551229, 0.3137256249999997, 0.127568, 0.566949, 0.550556, 0.3161762499999996, 0.126453, 0.570633, 0.549841, 0.31862749999999984, 0.125394, 0.574318, 0.549086, 0.3210781249999997, 0.124395, 0.578002, 0.548287, 0.3235293749999996, 0.123463, 0.581687, 0.547445, 0.32598062499999947, 0.122606, 0.585371, 0.546557, 0.32843124999999973, 0.121831, 0.589055, 0.545623, 0.3308824999999997, 0.121148, 0.592739, 0.544641, 0.3333331249999995, 0.120565, 0.596422, 0.543611, 0.33578437499999975, 0.120092, 0.600104, 0.54253, 0.3382349999999996, 0.119738, 0.603785, 0.5414, 0.3406862499999995, 0.119512, 0.607464, 0.540218, 0.3431374999999995, 0.119423, 0.611141, 0.538982, 0.34558812499999964, 0.119483, 0.614817, 0.537692, 0.34803937499999954, 0.119699, 0.61849, 0.536347, 0.35048999999999947, 0.120081, 0.622161, 0.534946, 0.35294124999999965, 0.120638, 0.625828, 0.533488, 0.3553918749999996, 0.12138, 0.629492, 0.531973, 0.3578431249999995, 0.122312, 0.633153, 0.530398, 0.36029437499999967, 0.123444, 0.636809, 0.528763, 0.3627449999999996, 0.12478, 0.640461, 0.527068, 0.3651962499999995, 0.126326, 0.644107, 0.525311, 0.36764687499999943, 0.128087, 0.647749, 0.523491, 0.37009812499999817, 0.130067, 0.651384, 0.521608, 0.3725487499999995, 0.132268, 0.655014, 0.519661, 0.3750000000000001, 0.134692, 0.658636, 0.517649, 0.37745125000000074, 0.137339, 0.662252, 0.515571, 0.3799018749999984, 0.14021, 0.665859, 0.513427, 0.38235312499999907, 0.143303, 0.669459, 0.511215, 0.3848037500000004, 0.146616, 0.67305, 0.508936, 0.38725500000000107, 0.150148, 0.676631, 0.506589, 0.3897056249999987, 0.153894, 0.680203, 0.504172, 0.3921568749999993, 0.157851, 0.683765, 0.501686, 0.394608125, 0.162016, 0.687316, 0.499129, 0.39705875000000135, 0.166383, 0.690856, 0.496502, 0.3995099999999982, 0.170948, 0.694384, 0.493803, 0.4019606249999996, 0.175707, 0.6979, 0.491033, 0.40441187500000025, 0.180653, 0.701402, 0.488189, 0.4068624999999979, 0.185783, 0.704891, 0.485273, 0.40931374999999853, 0.19109, 0.708366, 0.482284, 0.41176499999999916, 0.196571, 0.711827, 0.479221, 0.41421562500000053, 0.202219, 0.715272, 0.476084, 0.4166668750000011, 0.20803, 0.718701, 0.472873, 0.4191174999999988, 0.214, 0.722114, 0.469588, 0.42156874999999944, 0.220124, 0.725509, 0.466226, 0.42401937500000075, 0.226397, 0.728888, 0.462789, 0.4264706249999977, 0.232815, 0.732247, 0.459277, 0.42892187499999834, 0.239374, 0.735588, 0.455688, 0.43137249999999977, 0.24607, 0.73891, 0.452024, 0.4338237500000003, 0.252899, 0.742211, 0.448284, 0.436274374999998, 0.259857, 0.745492, 0.444467, 0.4387256249999987, 0.266941, 0.748751, 0.440573, 0.44117624999999994, 0.274149, 0.751988, 0.436601, 0.44362750000000056, 0.281477, 0.755203, 0.432552, 0.4460781249999983, 0.288921, 0.758394, 0.428426, 0.44852937499999884, 0.296479, 0.761561, 0.424223, 0.4509806249999995, 0.304148, 0.764704, 0.419943, 0.4534312500000009, 0.311925, 0.767822, 0.415586, 0.45588249999999786, 0.319809, 0.770914, 0.411152, 0.4583331249999991, 0.327796, 0.77398, 0.40664, 0.4607843749999998, 0.335885, 0.777018, 0.402049, 0.4632350000000012, 0.344074, 0.780029, 0.397381, 0.465686249999998, 0.35236, 0.783011, 0.392636, 0.4681374999999987, 0.360741, 0.785964, 0.387814, 0.4705881250000001, 0.369214, 0.788888, 0.382914, 0.4730393750000007, 0.377779, 0.791781, 0.377939, 0.47548999999999836, 0.386433, 0.794644, 0.372886, 0.477941249999999, 0.395174, 0.797475, 0.367757, 0.4803918750000004, 0.404001, 0.800275, 0.362552, 0.482843125000001, 0.412913, 0.803041, 0.357269, 0.4852943749999979, 0.421908, 0.805774, 0.35191, 0.4877449999999993, 0.430983, 0.808473, 0.346476, 0.4901962499999999, 0.440137, 0.811138, 0.340967, 0.4926468750000012, 0.449368, 0.813768, 0.335384, 0.4950981249999982, 0.458674, 0.816363, 0.329727, 0.4975487499999996, 0.468053, 0.818921, 0.323998, 0.5000000000000001, 0.477504, 0.821444, 0.318195, 0.5024512500000009, 0.487026, 0.823929, 0.312321, 0.5049018749999984, 0.496615, 0.826376, 0.306377, 0.5073531249999991, 0.506271, 0.828786, 0.300362, 0.5098037500000006, 0.515992, 0.831158, 0.294279, 0.5122550000000011, 0.525776, 0.833491, 0.288127, 0.5147056249999987, 0.535621, 0.835785, 0.281908, 0.5171568749999994, 0.545524, 0.838039, 0.275626, 0.519608125, 0.555484, 0.840254, 0.269281, 0.5220587499999978, 0.565498, 0.84243, 0.262877, 0.5245099999999984, 0.575563, 0.844566, 0.256415, 0.5269606249999996, 0.585678, 0.846661, 0.249897, 0.5294118750000003, 0.595839, 0.848717, 0.243329, 0.531862499999998, 0.606045, 0.850733, 0.236712, 0.5343137499999986, 0.616293, 0.852709, 0.230052, 0.5367649999999993, 0.626579, 0.854645, 0.223353, 0.5392156250000005, 0.636902, 0.856542, 0.21662, 0.5416668750000012, 0.647257, 0.8584, 0.209861, 0.5441174999999989, 0.657642, 0.860219, 0.203082, 0.5465687499999994, 0.668054, 0.861999, 0.196293, 0.5490193750000009, 0.678489, 0.863742, 0.189503, 0.5514706249999979, 0.688944, 0.865448, 0.182725, 0.5539218749999985, 0.699415, 0.867117, 0.175971, 0.5563724999999997, 0.709898, 0.868751, 0.169257, 0.5588237500000004, 0.720391, 0.87035, 0.162603, 0.5612743749999981, 0.730889, 0.871916, 0.156029, 0.5637256249999987, 0.741388, 0.873449, 0.149561, 0.56617625, 0.751884, 0.874951, 0.143228, 0.5686275000000007, 0.762373, 0.876424, 0.137064, 0.5710781249999983, 0.772852, 0.877868, 0.131109, 0.573529374999999, 0.783315, 0.879285, 0.125405, 0.5759806249999996, 0.79376, 0.880678, 0.120005, 0.578431250000001, 0.804182, 0.882046, 0.114965, 0.5808824999999979, 0.814576, 0.883393, 0.110347, 0.5833331249999992, 0.82494, 0.88472, 0.106217, 0.5857843749999999, 0.83527, 0.886029, 0.102646, 0.5882350000000012, 0.845561, 0.887322, 0.099702, 0.5906862499999981, 0.85581, 0.888601, 0.097452, 0.5931374999999988, 0.866013, 0.889868, 0.095953, 0.5955881250000001, 0.876168, 0.891125, 0.09525, 0.5980393750000007, 0.886271, 0.892374, 0.095374, 0.6004899999999984, 0.89632, 0.893616, 0.096335, 0.6029412499999991, 0.906311, 0.894855, 0.098125, 0.6053918750000004, 0.916242, 0.896091, 0.100717, 0.607843125000001, 0.926106, 0.89733, 0.104071, 0.610294374999998, 0.935904, 0.89857, 0.108131, 0.6127449999999993, 0.945636, 0.899815, 0.112838, 0.61519625, 0.9553, 0.901065, 0.118128, 0.6176468749999976, 0.964894, 0.902323, 0.123941, 0.6200981249999983, 0.974417, 0.90359, 0.130215, 0.6225487499999997, 0.983868, 0.904867, 0.136897, 0.6250000000000002, 0.993248, 0.906157, 0.143936]
xLUT.ColorSpace = 'RGB'
xLUT.AboveRangeColor = [1.0, 1.0, 1.0]
xLUT.NanColor = [1.0, 0.0, 0.0]
xLUT.NumberOfTableValues = 512
xLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'x'
xPWF = GetOpacityTransferFunction('x')
xPWF.Points = [0.0, 0.0, 0.5, 0.0, 0.6250000000000002, 1.0, 0.5, 0.0]
xPWF.ScalarRangeInitialized = 1

# set active source
SetActiveSource(u_1pvd)

# set scalar coloring
ColorBy(u_1pvdDisplay, ('POINTS', 'x', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
u_1pvdDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
u_1pvdDisplay.SetScalarBarVisibility(renderView1, True)

# set active source
SetActiveSource(u_2pvd)

# set scalar coloring
ColorBy(u_2pvdDisplay, ('POINTS', 'x', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
u_2pvdDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
u_2pvdDisplay.SetScalarBarVisibility(renderView1, True)

# set active source
SetActiveSource(u_3pvd)

# set scalar coloring
ColorBy(u_3pvdDisplay, ('POINTS', 'x', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
u_3pvdDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
u_3pvdDisplay.SetScalarBarVisibility(renderView1, True)

# set active source
SetActiveSource(u_4pvd)

# set scalar coloring
ColorBy(u_4pvdDisplay, ('POINTS', 'x', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
u_4pvdDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
u_4pvdDisplay.SetScalarBarVisibility(renderView1, True)

# set active source
SetActiveSource(u_5pvd)

# set scalar coloring
ColorBy(u_5pvdDisplay, ('POINTS', 'x', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
u_5pvdDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
u_5pvdDisplay.SetScalarBarVisibility(renderView1, True)

# set active source
SetActiveSource(u_6pvd)

# set scalar coloring
ColorBy(u_6pvdDisplay, ('POINTS', 'x', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
u_6pvdDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
u_6pvdDisplay.SetScalarBarVisibility(renderView1, True)

# set active source
SetActiveSource(u_7pvd)

# set scalar coloring
ColorBy(u_7pvdDisplay, ('POINTS', 'x', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
u_7pvdDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
u_7pvdDisplay.SetScalarBarVisibility(renderView1, True)

# set active source
SetActiveSource(u_8pvd)

# set scalar coloring
ColorBy(u_8pvdDisplay, ('POINTS', 'x', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
u_8pvdDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
u_8pvdDisplay.SetScalarBarVisibility(renderView1, True)

# set active source
SetActiveSource(u_9pvd)

# set scalar coloring
ColorBy(u_9pvdDisplay, ('POINTS', 'x', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
u_9pvdDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
u_9pvdDisplay.SetScalarBarVisibility(renderView1, True)

# create a new 'PVD Reader'
u_1pvd_1 = PVDReader(FileName='/home/dokken/Documents/MultiMeshShapeOpt_code/Stokes_rotation/output/u_1.pvd')
u_1pvd_1.PointArrays = ['x']

# create a new 'PVD Reader'
u_2pvd_1 = PVDReader(FileName='/home/dokken/Documents/MultiMeshShapeOpt_code/Stokes_rotation/output/u_2.pvd')
u_2pvd_1.PointArrays = ['x']

# create a new 'PVD Reader'
u_3pvd_1 = PVDReader(FileName='/home/dokken/Documents/MultiMeshShapeOpt_code/Stokes_rotation/output/u_3.pvd')
u_3pvd_1.PointArrays = ['x']

# create a new 'PVD Reader'
u_4pvd_1 = PVDReader(FileName='/home/dokken/Documents/MultiMeshShapeOpt_code/Stokes_rotation/output/u_4.pvd')
u_4pvd_1.PointArrays = ['x']

# create a new 'PVD Reader'
u_5pvd_1 = PVDReader(FileName='/home/dokken/Documents/MultiMeshShapeOpt_code/Stokes_rotation/output/u_5.pvd')
u_5pvd_1.PointArrays = ['x']

# create a new 'PVD Reader'
u_6pvd_1 = PVDReader(FileName='/home/dokken/Documents/MultiMeshShapeOpt_code/Stokes_rotation/output/u_6.pvd')
u_6pvd_1.PointArrays = ['x']

# create a new 'PVD Reader'
u_7pvd_1 = PVDReader(FileName='/home/dokken/Documents/MultiMeshShapeOpt_code/Stokes_rotation/output/u_7.pvd')
u_7pvd_1.PointArrays = ['x']

# create a new 'PVD Reader'
u_8pvd_1 = PVDReader(FileName='/home/dokken/Documents/MultiMeshShapeOpt_code/Stokes_rotation/output/u_8.pvd')
u_8pvd_1.PointArrays = ['x']

# create a new 'PVD Reader'
u_9pvd_1 = PVDReader(FileName='/home/dokken/Documents/MultiMeshShapeOpt_code/Stokes_rotation/output/u_9.pvd')
u_9pvd_1.PointArrays = ['x']

# show data in view
u_1pvd_1Display = Show(u_1pvd_1, renderView1)

# trace defaults for the display properties.
u_1pvd_1Display.Representation = 'Surface'
u_1pvd_1Display.AmbientColor = [0.0, 0.0, 0.0]
u_1pvd_1Display.ColorArrayName = [None, '']
u_1pvd_1Display.DiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_1pvd_1Display.BackfaceDiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_1pvd_1Display.OSPRayScaleArray = 'x'
u_1pvd_1Display.OSPRayScaleFunction = 'PiecewiseFunction'
u_1pvd_1Display.SelectOrientationVectors = 'x'
u_1pvd_1Display.ScaleFactor = 0.021400000000000002
u_1pvd_1Display.SelectScaleArray = 'None'
u_1pvd_1Display.GlyphType = 'Arrow'
u_1pvd_1Display.GlyphTableIndexArray = 'None'
u_1pvd_1Display.GaussianRadius = 0.00107
u_1pvd_1Display.SetScaleArray = ['POINTS', 'x']
u_1pvd_1Display.ScaleTransferFunction = 'PiecewiseFunction'
u_1pvd_1Display.OpacityArray = ['POINTS', 'x']
u_1pvd_1Display.OpacityTransferFunction = 'PiecewiseFunction'
u_1pvd_1Display.DataAxesGrid = 'GridAxesRepresentation'
u_1pvd_1Display.SelectionCellLabelFontFile = ''
u_1pvd_1Display.SelectionPointLabelFontFile = ''
u_1pvd_1Display.PolarAxes = 'PolarAxesRepresentation'
u_1pvd_1Display.ScalarOpacityUnitDistance = 0.02038586012005983

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
u_1pvd_1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
u_1pvd_1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
u_1pvd_1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
u_1pvd_1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
u_1pvd_1Display.DataAxesGrid.XTitleFontFile = ''
u_1pvd_1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
u_1pvd_1Display.DataAxesGrid.YTitleFontFile = ''
u_1pvd_1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
u_1pvd_1Display.DataAxesGrid.ZTitleFontFile = ''
u_1pvd_1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
u_1pvd_1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
u_1pvd_1Display.DataAxesGrid.XLabelFontFile = ''
u_1pvd_1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
u_1pvd_1Display.DataAxesGrid.YLabelFontFile = ''
u_1pvd_1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
u_1pvd_1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
u_1pvd_1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
u_1pvd_1Display.PolarAxes.PolarAxisTitleFontFile = ''
u_1pvd_1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
u_1pvd_1Display.PolarAxes.PolarAxisLabelFontFile = ''
u_1pvd_1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
u_1pvd_1Display.PolarAxes.LastRadialAxisTextFontFile = ''
u_1pvd_1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
u_1pvd_1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data in view
u_2pvd_1Display = Show(u_2pvd_1, renderView1)

# trace defaults for the display properties.
u_2pvd_1Display.Representation = 'Surface'
u_2pvd_1Display.AmbientColor = [0.0, 0.0, 0.0]
u_2pvd_1Display.ColorArrayName = [None, '']
u_2pvd_1Display.DiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_2pvd_1Display.BackfaceDiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_2pvd_1Display.OSPRayScaleArray = 'x'
u_2pvd_1Display.OSPRayScaleFunction = 'PiecewiseFunction'
u_2pvd_1Display.SelectOrientationVectors = 'x'
u_2pvd_1Display.ScaleFactor = 0.0214
u_2pvd_1Display.SelectScaleArray = 'None'
u_2pvd_1Display.GlyphType = 'Arrow'
u_2pvd_1Display.GlyphTableIndexArray = 'None'
u_2pvd_1Display.GaussianRadius = 0.0010699999999999998
u_2pvd_1Display.SetScaleArray = ['POINTS', 'x']
u_2pvd_1Display.ScaleTransferFunction = 'PiecewiseFunction'
u_2pvd_1Display.OpacityArray = ['POINTS', 'x']
u_2pvd_1Display.OpacityTransferFunction = 'PiecewiseFunction'
u_2pvd_1Display.DataAxesGrid = 'GridAxesRepresentation'
u_2pvd_1Display.SelectionCellLabelFontFile = ''
u_2pvd_1Display.SelectionPointLabelFontFile = ''
u_2pvd_1Display.PolarAxes = 'PolarAxesRepresentation'
u_2pvd_1Display.ScalarOpacityUnitDistance = 0.020385860120059825

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
u_2pvd_1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
u_2pvd_1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
u_2pvd_1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
u_2pvd_1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
u_2pvd_1Display.DataAxesGrid.XTitleFontFile = ''
u_2pvd_1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
u_2pvd_1Display.DataAxesGrid.YTitleFontFile = ''
u_2pvd_1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
u_2pvd_1Display.DataAxesGrid.ZTitleFontFile = ''
u_2pvd_1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
u_2pvd_1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
u_2pvd_1Display.DataAxesGrid.XLabelFontFile = ''
u_2pvd_1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
u_2pvd_1Display.DataAxesGrid.YLabelFontFile = ''
u_2pvd_1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
u_2pvd_1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
u_2pvd_1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
u_2pvd_1Display.PolarAxes.PolarAxisTitleFontFile = ''
u_2pvd_1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
u_2pvd_1Display.PolarAxes.PolarAxisLabelFontFile = ''
u_2pvd_1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
u_2pvd_1Display.PolarAxes.LastRadialAxisTextFontFile = ''
u_2pvd_1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
u_2pvd_1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data in view
u_7pvd_1Display = Show(u_7pvd_1, renderView1)

# trace defaults for the display properties.
u_7pvd_1Display.Representation = 'Surface'
u_7pvd_1Display.AmbientColor = [0.0, 0.0, 0.0]
u_7pvd_1Display.ColorArrayName = [None, '']
u_7pvd_1Display.DiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_7pvd_1Display.BackfaceDiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_7pvd_1Display.OSPRayScaleArray = 'x'
u_7pvd_1Display.OSPRayScaleFunction = 'PiecewiseFunction'
u_7pvd_1Display.SelectOrientationVectors = 'x'
u_7pvd_1Display.ScaleFactor = 0.021400000000000002
u_7pvd_1Display.SelectScaleArray = 'None'
u_7pvd_1Display.GlyphType = 'Arrow'
u_7pvd_1Display.GlyphTableIndexArray = 'None'
u_7pvd_1Display.GaussianRadius = 0.00107
u_7pvd_1Display.SetScaleArray = ['POINTS', 'x']
u_7pvd_1Display.ScaleTransferFunction = 'PiecewiseFunction'
u_7pvd_1Display.OpacityArray = ['POINTS', 'x']
u_7pvd_1Display.OpacityTransferFunction = 'PiecewiseFunction'
u_7pvd_1Display.DataAxesGrid = 'GridAxesRepresentation'
u_7pvd_1Display.SelectionCellLabelFontFile = ''
u_7pvd_1Display.SelectionPointLabelFontFile = ''
u_7pvd_1Display.PolarAxes = 'PolarAxesRepresentation'
u_7pvd_1Display.ScalarOpacityUnitDistance = 0.02038586012005983

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
u_7pvd_1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
u_7pvd_1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
u_7pvd_1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
u_7pvd_1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
u_7pvd_1Display.DataAxesGrid.XTitleFontFile = ''
u_7pvd_1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
u_7pvd_1Display.DataAxesGrid.YTitleFontFile = ''
u_7pvd_1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
u_7pvd_1Display.DataAxesGrid.ZTitleFontFile = ''
u_7pvd_1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
u_7pvd_1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
u_7pvd_1Display.DataAxesGrid.XLabelFontFile = ''
u_7pvd_1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
u_7pvd_1Display.DataAxesGrid.YLabelFontFile = ''
u_7pvd_1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
u_7pvd_1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
u_7pvd_1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
u_7pvd_1Display.PolarAxes.PolarAxisTitleFontFile = ''
u_7pvd_1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
u_7pvd_1Display.PolarAxes.PolarAxisLabelFontFile = ''
u_7pvd_1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
u_7pvd_1Display.PolarAxes.LastRadialAxisTextFontFile = ''
u_7pvd_1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
u_7pvd_1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data in view
u_4pvd_1Display = Show(u_4pvd_1, renderView1)

# trace defaults for the display properties.
u_4pvd_1Display.Representation = 'Surface'
u_4pvd_1Display.AmbientColor = [0.0, 0.0, 0.0]
u_4pvd_1Display.ColorArrayName = [None, '']
u_4pvd_1Display.DiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_4pvd_1Display.BackfaceDiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_4pvd_1Display.OSPRayScaleArray = 'x'
u_4pvd_1Display.OSPRayScaleFunction = 'PiecewiseFunction'
u_4pvd_1Display.SelectOrientationVectors = 'x'
u_4pvd_1Display.ScaleFactor = 0.021400000000000002
u_4pvd_1Display.SelectScaleArray = 'None'
u_4pvd_1Display.GlyphType = 'Arrow'
u_4pvd_1Display.GlyphTableIndexArray = 'None'
u_4pvd_1Display.GaussianRadius = 0.00107
u_4pvd_1Display.SetScaleArray = ['POINTS', 'x']
u_4pvd_1Display.ScaleTransferFunction = 'PiecewiseFunction'
u_4pvd_1Display.OpacityArray = ['POINTS', 'x']
u_4pvd_1Display.OpacityTransferFunction = 'PiecewiseFunction'
u_4pvd_1Display.DataAxesGrid = 'GridAxesRepresentation'
u_4pvd_1Display.SelectionCellLabelFontFile = ''
u_4pvd_1Display.SelectionPointLabelFontFile = ''
u_4pvd_1Display.PolarAxes = 'PolarAxesRepresentation'
u_4pvd_1Display.ScalarOpacityUnitDistance = 0.020385860120059825

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
u_4pvd_1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
u_4pvd_1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
u_4pvd_1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
u_4pvd_1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
u_4pvd_1Display.DataAxesGrid.XTitleFontFile = ''
u_4pvd_1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
u_4pvd_1Display.DataAxesGrid.YTitleFontFile = ''
u_4pvd_1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
u_4pvd_1Display.DataAxesGrid.ZTitleFontFile = ''
u_4pvd_1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
u_4pvd_1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
u_4pvd_1Display.DataAxesGrid.XLabelFontFile = ''
u_4pvd_1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
u_4pvd_1Display.DataAxesGrid.YLabelFontFile = ''
u_4pvd_1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
u_4pvd_1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
u_4pvd_1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
u_4pvd_1Display.PolarAxes.PolarAxisTitleFontFile = ''
u_4pvd_1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
u_4pvd_1Display.PolarAxes.PolarAxisLabelFontFile = ''
u_4pvd_1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
u_4pvd_1Display.PolarAxes.LastRadialAxisTextFontFile = ''
u_4pvd_1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
u_4pvd_1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data in view
u_5pvd_1Display = Show(u_5pvd_1, renderView1)

# trace defaults for the display properties.
u_5pvd_1Display.Representation = 'Surface'
u_5pvd_1Display.AmbientColor = [0.0, 0.0, 0.0]
u_5pvd_1Display.ColorArrayName = [None, '']
u_5pvd_1Display.DiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_5pvd_1Display.BackfaceDiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_5pvd_1Display.OSPRayScaleArray = 'x'
u_5pvd_1Display.OSPRayScaleFunction = 'PiecewiseFunction'
u_5pvd_1Display.SelectOrientationVectors = 'x'
u_5pvd_1Display.ScaleFactor = 0.0214
u_5pvd_1Display.SelectScaleArray = 'None'
u_5pvd_1Display.GlyphType = 'Arrow'
u_5pvd_1Display.GlyphTableIndexArray = 'None'
u_5pvd_1Display.GaussianRadius = 0.0010699999999999998
u_5pvd_1Display.SetScaleArray = ['POINTS', 'x']
u_5pvd_1Display.ScaleTransferFunction = 'PiecewiseFunction'
u_5pvd_1Display.OpacityArray = ['POINTS', 'x']
u_5pvd_1Display.OpacityTransferFunction = 'PiecewiseFunction'
u_5pvd_1Display.DataAxesGrid = 'GridAxesRepresentation'
u_5pvd_1Display.SelectionCellLabelFontFile = ''
u_5pvd_1Display.SelectionPointLabelFontFile = ''
u_5pvd_1Display.PolarAxes = 'PolarAxesRepresentation'
u_5pvd_1Display.ScalarOpacityUnitDistance = 0.020385860120059825

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
u_5pvd_1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
u_5pvd_1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
u_5pvd_1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
u_5pvd_1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
u_5pvd_1Display.DataAxesGrid.XTitleFontFile = ''
u_5pvd_1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
u_5pvd_1Display.DataAxesGrid.YTitleFontFile = ''
u_5pvd_1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
u_5pvd_1Display.DataAxesGrid.ZTitleFontFile = ''
u_5pvd_1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
u_5pvd_1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
u_5pvd_1Display.DataAxesGrid.XLabelFontFile = ''
u_5pvd_1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
u_5pvd_1Display.DataAxesGrid.YLabelFontFile = ''
u_5pvd_1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
u_5pvd_1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
u_5pvd_1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
u_5pvd_1Display.PolarAxes.PolarAxisTitleFontFile = ''
u_5pvd_1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
u_5pvd_1Display.PolarAxes.PolarAxisLabelFontFile = ''
u_5pvd_1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
u_5pvd_1Display.PolarAxes.LastRadialAxisTextFontFile = ''
u_5pvd_1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
u_5pvd_1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data in view
u_3pvd_1Display = Show(u_3pvd_1, renderView1)

# trace defaults for the display properties.
u_3pvd_1Display.Representation = 'Surface'
u_3pvd_1Display.AmbientColor = [0.0, 0.0, 0.0]
u_3pvd_1Display.ColorArrayName = [None, '']
u_3pvd_1Display.DiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_3pvd_1Display.BackfaceDiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_3pvd_1Display.OSPRayScaleArray = 'x'
u_3pvd_1Display.OSPRayScaleFunction = 'PiecewiseFunction'
u_3pvd_1Display.SelectOrientationVectors = 'x'
u_3pvd_1Display.ScaleFactor = 0.0214
u_3pvd_1Display.SelectScaleArray = 'None'
u_3pvd_1Display.GlyphType = 'Arrow'
u_3pvd_1Display.GlyphTableIndexArray = 'None'
u_3pvd_1Display.GaussianRadius = 0.0010699999999999998
u_3pvd_1Display.SetScaleArray = ['POINTS', 'x']
u_3pvd_1Display.ScaleTransferFunction = 'PiecewiseFunction'
u_3pvd_1Display.OpacityArray = ['POINTS', 'x']
u_3pvd_1Display.OpacityTransferFunction = 'PiecewiseFunction'
u_3pvd_1Display.DataAxesGrid = 'GridAxesRepresentation'
u_3pvd_1Display.SelectionCellLabelFontFile = ''
u_3pvd_1Display.SelectionPointLabelFontFile = ''
u_3pvd_1Display.PolarAxes = 'PolarAxesRepresentation'
u_3pvd_1Display.ScalarOpacityUnitDistance = 0.020385860120059825

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
u_3pvd_1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
u_3pvd_1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
u_3pvd_1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
u_3pvd_1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
u_3pvd_1Display.DataAxesGrid.XTitleFontFile = ''
u_3pvd_1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
u_3pvd_1Display.DataAxesGrid.YTitleFontFile = ''
u_3pvd_1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
u_3pvd_1Display.DataAxesGrid.ZTitleFontFile = ''
u_3pvd_1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
u_3pvd_1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
u_3pvd_1Display.DataAxesGrid.XLabelFontFile = ''
u_3pvd_1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
u_3pvd_1Display.DataAxesGrid.YLabelFontFile = ''
u_3pvd_1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
u_3pvd_1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
u_3pvd_1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
u_3pvd_1Display.PolarAxes.PolarAxisTitleFontFile = ''
u_3pvd_1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
u_3pvd_1Display.PolarAxes.PolarAxisLabelFontFile = ''
u_3pvd_1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
u_3pvd_1Display.PolarAxes.LastRadialAxisTextFontFile = ''
u_3pvd_1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
u_3pvd_1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data in view
u_8pvd_1Display = Show(u_8pvd_1, renderView1)

# trace defaults for the display properties.
u_8pvd_1Display.Representation = 'Surface'
u_8pvd_1Display.AmbientColor = [0.0, 0.0, 0.0]
u_8pvd_1Display.ColorArrayName = [None, '']
u_8pvd_1Display.DiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_8pvd_1Display.BackfaceDiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_8pvd_1Display.OSPRayScaleArray = 'x'
u_8pvd_1Display.OSPRayScaleFunction = 'PiecewiseFunction'
u_8pvd_1Display.SelectOrientationVectors = 'x'
u_8pvd_1Display.ScaleFactor = 0.0214
u_8pvd_1Display.SelectScaleArray = 'None'
u_8pvd_1Display.GlyphType = 'Arrow'
u_8pvd_1Display.GlyphTableIndexArray = 'None'
u_8pvd_1Display.GaussianRadius = 0.0010699999999999998
u_8pvd_1Display.SetScaleArray = ['POINTS', 'x']
u_8pvd_1Display.ScaleTransferFunction = 'PiecewiseFunction'
u_8pvd_1Display.OpacityArray = ['POINTS', 'x']
u_8pvd_1Display.OpacityTransferFunction = 'PiecewiseFunction'
u_8pvd_1Display.DataAxesGrid = 'GridAxesRepresentation'
u_8pvd_1Display.SelectionCellLabelFontFile = ''
u_8pvd_1Display.SelectionPointLabelFontFile = ''
u_8pvd_1Display.PolarAxes = 'PolarAxesRepresentation'
u_8pvd_1Display.ScalarOpacityUnitDistance = 0.020385860120059825

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
u_8pvd_1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
u_8pvd_1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
u_8pvd_1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
u_8pvd_1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
u_8pvd_1Display.DataAxesGrid.XTitleFontFile = ''
u_8pvd_1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
u_8pvd_1Display.DataAxesGrid.YTitleFontFile = ''
u_8pvd_1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
u_8pvd_1Display.DataAxesGrid.ZTitleFontFile = ''
u_8pvd_1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
u_8pvd_1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
u_8pvd_1Display.DataAxesGrid.XLabelFontFile = ''
u_8pvd_1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
u_8pvd_1Display.DataAxesGrid.YLabelFontFile = ''
u_8pvd_1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
u_8pvd_1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
u_8pvd_1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
u_8pvd_1Display.PolarAxes.PolarAxisTitleFontFile = ''
u_8pvd_1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
u_8pvd_1Display.PolarAxes.PolarAxisLabelFontFile = ''
u_8pvd_1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
u_8pvd_1Display.PolarAxes.LastRadialAxisTextFontFile = ''
u_8pvd_1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
u_8pvd_1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data in view
u_9pvd_1Display = Show(u_9pvd_1, renderView1)

# trace defaults for the display properties.
u_9pvd_1Display.Representation = 'Surface'
u_9pvd_1Display.AmbientColor = [0.0, 0.0, 0.0]
u_9pvd_1Display.ColorArrayName = [None, '']
u_9pvd_1Display.DiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_9pvd_1Display.BackfaceDiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_9pvd_1Display.OSPRayScaleArray = 'x'
u_9pvd_1Display.OSPRayScaleFunction = 'PiecewiseFunction'
u_9pvd_1Display.SelectOrientationVectors = 'x'
u_9pvd_1Display.ScaleFactor = 0.0214
u_9pvd_1Display.SelectScaleArray = 'None'
u_9pvd_1Display.GlyphType = 'Arrow'
u_9pvd_1Display.GlyphTableIndexArray = 'None'
u_9pvd_1Display.GaussianRadius = 0.0010699999999999998
u_9pvd_1Display.SetScaleArray = ['POINTS', 'x']
u_9pvd_1Display.ScaleTransferFunction = 'PiecewiseFunction'
u_9pvd_1Display.OpacityArray = ['POINTS', 'x']
u_9pvd_1Display.OpacityTransferFunction = 'PiecewiseFunction'
u_9pvd_1Display.DataAxesGrid = 'GridAxesRepresentation'
u_9pvd_1Display.SelectionCellLabelFontFile = ''
u_9pvd_1Display.SelectionPointLabelFontFile = ''
u_9pvd_1Display.PolarAxes = 'PolarAxesRepresentation'
u_9pvd_1Display.ScalarOpacityUnitDistance = 0.020385860120059825

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
u_9pvd_1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
u_9pvd_1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
u_9pvd_1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
u_9pvd_1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
u_9pvd_1Display.DataAxesGrid.XTitleFontFile = ''
u_9pvd_1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
u_9pvd_1Display.DataAxesGrid.YTitleFontFile = ''
u_9pvd_1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
u_9pvd_1Display.DataAxesGrid.ZTitleFontFile = ''
u_9pvd_1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
u_9pvd_1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
u_9pvd_1Display.DataAxesGrid.XLabelFontFile = ''
u_9pvd_1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
u_9pvd_1Display.DataAxesGrid.YLabelFontFile = ''
u_9pvd_1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
u_9pvd_1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
u_9pvd_1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
u_9pvd_1Display.PolarAxes.PolarAxisTitleFontFile = ''
u_9pvd_1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
u_9pvd_1Display.PolarAxes.PolarAxisLabelFontFile = ''
u_9pvd_1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
u_9pvd_1Display.PolarAxes.LastRadialAxisTextFontFile = ''
u_9pvd_1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
u_9pvd_1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data in view
u_6pvd_1Display = Show(u_6pvd_1, renderView1)

# trace defaults for the display properties.
u_6pvd_1Display.Representation = 'Surface'
u_6pvd_1Display.AmbientColor = [0.0, 0.0, 0.0]
u_6pvd_1Display.ColorArrayName = [None, '']
u_6pvd_1Display.DiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_6pvd_1Display.BackfaceDiffuseColor = [0.6039215686274509, 0.6039215686274509, 0.6039215686274509]
u_6pvd_1Display.OSPRayScaleArray = 'x'
u_6pvd_1Display.OSPRayScaleFunction = 'PiecewiseFunction'
u_6pvd_1Display.SelectOrientationVectors = 'x'
u_6pvd_1Display.ScaleFactor = 0.0214
u_6pvd_1Display.SelectScaleArray = 'None'
u_6pvd_1Display.GlyphType = 'Arrow'
u_6pvd_1Display.GlyphTableIndexArray = 'None'
u_6pvd_1Display.GaussianRadius = 0.0010699999999999998
u_6pvd_1Display.SetScaleArray = ['POINTS', 'x']
u_6pvd_1Display.ScaleTransferFunction = 'PiecewiseFunction'
u_6pvd_1Display.OpacityArray = ['POINTS', 'x']
u_6pvd_1Display.OpacityTransferFunction = 'PiecewiseFunction'
u_6pvd_1Display.DataAxesGrid = 'GridAxesRepresentation'
u_6pvd_1Display.SelectionCellLabelFontFile = ''
u_6pvd_1Display.SelectionPointLabelFontFile = ''
u_6pvd_1Display.PolarAxes = 'PolarAxesRepresentation'
u_6pvd_1Display.ScalarOpacityUnitDistance = 0.020385860120059825

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
u_6pvd_1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
u_6pvd_1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
u_6pvd_1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.3286008936764144, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
u_6pvd_1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
u_6pvd_1Display.DataAxesGrid.XTitleFontFile = ''
u_6pvd_1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
u_6pvd_1Display.DataAxesGrid.YTitleFontFile = ''
u_6pvd_1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
u_6pvd_1Display.DataAxesGrid.ZTitleFontFile = ''
u_6pvd_1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
u_6pvd_1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
u_6pvd_1Display.DataAxesGrid.XLabelFontFile = ''
u_6pvd_1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
u_6pvd_1Display.DataAxesGrid.YLabelFontFile = ''
u_6pvd_1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
u_6pvd_1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
u_6pvd_1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
u_6pvd_1Display.PolarAxes.PolarAxisTitleFontFile = ''
u_6pvd_1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
u_6pvd_1Display.PolarAxes.PolarAxisLabelFontFile = ''
u_6pvd_1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
u_6pvd_1Display.PolarAxes.LastRadialAxisTextFontFile = ''
u_6pvd_1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
u_6pvd_1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(u_1pvd_1)

# change representation type
u_1pvd_1Display.SetRepresentationType('Surface With Edges')

# Properties modified on u_1pvd_1Display
u_1pvd_1Display.EdgeColor = [1.0, 0.0, 0.0]

# change representation type
u_1pvd_1Display.SetRepresentationType('Wireframe')

# change solid color
u_1pvd_1Display.AmbientColor = [1.0, 0.0, 0.0]

# set active source
SetActiveSource(u_2pvd_1)

# change representation type
u_2pvd_1Display.SetRepresentationType('Wireframe')

# change solid color
u_2pvd_1Display.AmbientColor = [1.0, 0.0, 0.0]

# set active source
SetActiveSource(u_3pvd_1)

# change representation type
u_3pvd_1Display.SetRepresentationType('Wireframe')

# change solid color
u_3pvd_1Display.AmbientColor = [1.0, 0.0, 0.0]

# set active source
SetActiveSource(u_4pvd_1)

# change representation type
u_4pvd_1Display.SetRepresentationType('Wireframe')

# change solid color
u_4pvd_1Display.AmbientColor = [1.0, 0.0, 0.0]

# set active source
SetActiveSource(u_5pvd_1)

# change representation type
u_5pvd_1Display.SetRepresentationType('Wireframe')

# change solid color
u_5pvd_1Display.AmbientColor = [1.0, 0.0, 0.0]

# set active source
SetActiveSource(u_6pvd_1)

# change representation type
u_6pvd_1Display.SetRepresentationType('Wireframe')

# change solid color
u_6pvd_1Display.AmbientColor = [1.0, 0.0, 0.0]

# set active source
SetActiveSource(u_7pvd_1)

# change representation type
u_7pvd_1Display.SetRepresentationType('Wireframe')

# change solid color
u_7pvd_1Display.AmbientColor = [1.0, 0.0, 0.0]

# set active source
SetActiveSource(u_8pvd_1)

# change representation type
u_8pvd_1Display.SetRepresentationType('Wireframe')

# change solid color
u_8pvd_1Display.AmbientColor = [1.0, 0.0, 0.0]

# set active source
SetActiveSource(u_9pvd_1)

# change representation type
u_9pvd_1Display.SetRepresentationType('Wireframe')

# change solid color
u_9pvd_1Display.AmbientColor = [1.0, 0.0, 0.0]

# get color legend/bar for xLUT in view renderView1
xLUTColorBar = GetScalarBar(xLUT, renderView1)
xLUTColorBar.Orientation = 'Horizontal'
xLUTColorBar.WindowLocation = 'AnyLocation'
xLUTColorBar.Position = [0.8466290182450043, 0.1298071219396452]
xLUTColorBar.Title = 'x'
xLUTColorBar.ComponentTitle = 'Magnitude'
xLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
xLUTColorBar.TitleFontFile = ''
xLUTColorBar.TitleBold = 1
xLUTColorBar.TitleFontSize = 25
xLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
xLUTColorBar.LabelFontFile = ''
xLUTColorBar.LabelBold = 1
xLUTColorBar.LabelFontSize = 20
xLUTColorBar.RangeLabelFormat = '%-#6.3e'
xLUTColorBar.ScalarBarLength = 0.802829085466778

# change scalar bar placement
xLUTColorBar.Position = [0.7302747842437114, 0.08070532553245954]
xLUTColorBar.ScalarBarLength = 0.8028290854667782

# set active source
SetActiveSource(u_1pvd)

# set active source
SetActiveSource(u_0pvd)

# Properties modified on xLUTColorBar
xLUTColorBar.Title = 'u'

# change scalar bar placement
xLUTColorBar.Position = [0.7296283718325931, 0.08070532553245954]

# change scalar bar placement
xLUTColorBar.Position = [0.7296283718325931, 0.08190293032287871]
xLUTColorBar.ScalarBarLength = 0.8028290854667783

# change scalar bar placement
xLUTColorBar.ScalarBarLength = 0.6375596243889319

# change scalar bar placement
xLUTColorBar.ScalarBarLength = 0.6663021393589916

# change scalar bar placement
xLUTColorBar.Position = [0.7121752367323991, 0.11663346924503443]

# change scalar bar placement
xLUTColorBar.Position = [0.7031254629767429, 0.1322023315204835]

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

# hide color bar/color legend
u_0pvdDisplay.SetScalarBarVisibility(renderView1, False)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.506774675747895, 0.5355670476764489, 10000.0]
renderView1.CameraFocalPoint = [0.506774675747895, 0.5355670476764489, 0.0]
renderView1.CameraParallelScale = 0.7071067811865476

# save screenshot
SaveScreenshot('/home/dokken/Documents/MultiMeshShapeOpt_code/Stokes_rotation/output/initial_rotation.png', renderView1, ImageResolution=[1547, 835])

animationScene1.GoToLast()

# show color bar/color legend
u_0pvdDisplay.SetScalarBarVisibility(renderView1, True)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.506774675747895, 0.5355670476764489, 10000.0]
renderView1.CameraFocalPoint = [0.506774675747895, 0.5355670476764489, 0.0]
renderView1.CameraParallelScale = 0.7071067811865476

# save screenshot
SaveScreenshot('/home/dokken/Documents/MultiMeshShapeOpt_code/Stokes_rotation/output/final_rotation.png', renderView1, ImageResolution=[1547, 835])

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.506774675747895, 0.5355670476764489, 10000.0]
renderView1.CameraFocalPoint = [0.506774675747895, 0.5355670476764489, 0.0]
renderView1.CameraParallelScale = 0.7071067811865476

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).