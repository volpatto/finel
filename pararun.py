#!/usr/bin/env python
# -*- coding: utf-8 -*-

import paraview.simple as pvs
import subprocess as sp
import shlex as sh
import sys
import os
import glob

# Verify if input case name argument is provided
try:
    casename = sys.argv[1]
except:
    sys.exit('*** Input casename error ***')

# Setting bash command
if casename[-1] != '/': casename = casename + '/'
cmd_line = '../finel.exe input.dat'
args = sh.split(cmd_line)

# Setting case dir
target_dir = sp.check_output('pwd').rstrip() + '/' + casename

# Executing the bash command
try:
    p1 = sp.check_call(args,shell=True,cwd=target_dir)
except:
    sys.exit('*** Simulator run error ***')

# Reading .vtk solution file
vtkfiles = glob.glob(target_dir + '/solution*.vtk')
try:
    solution_vtk = []
    for target_file in vtkfiles:
        solution_vtk.append(pvs.LegacyVTKReader(FileNames=[target_file]))
except:
    sys.exit('*** Paraview vtk reading error ***')

# Plotting .vtk solution file
try:
    # Getting the Data Points from solution file
    solution_data = solution_vtk[0].PointData
    solution_key = solution_data.keys()
    # Applying WarpByScalar filter in the retrieved scalar data
    warp = pvs.WarpByScalar(Input=solution_vtk[0],Scalars=solution_key[0],XYPlane=1)
    # Setting the view
    view = pvs.GetActiveViewOrCreate('RenderView')
    # Setting the background color (white)
    view.Background = [1,1,1]
    # Creating a display object
    display = pvs.Show(warp)
    # Setting the colorbar associated to the scalar values
    pvs.ColorBy(display, ('POINTS', solution_key[0]))
    # Rescale to data range
    display.RescaleTransferFunctionToDataRange(True)
    # Display the colorbar
    display.SetScalarBarVisibility(view, True)
    # Setting the view size
    view.ViewSize = [800, 650]
    view.CameraPosition = [0.0, 0.0, 1]
    # Getting color transfer function/color map
    warp_LUT = pvs.GetColorTransferFunction(solution_key[0])
    # Getting color legend/bar in the view
    warp_LUT_ColorBar = pvs.GetScalarBar(warp_LUT, view)
    # Setting colors to Title and Label of ScalarBar
    warp_LUT_ColorBar.TitleColor = [0.0, 0.0, 0.0]
    warp_LUT_ColorBar.LabelColor = [0.0, 0.0, 0.0]
    # Setting font size to label and title of ScalarBar
    warp_LUT_ColorBar.LabelFontSize = 9
    warp_LUT_ColorBar.TitleFontSize = 11
    # Getting render view to change axes label color
    rv = pvs.GetRenderView()
    # Changing axes label color to black
    rv.OrientationAxesLabelColor = [0, 0, 0]
    # Activating a Render object from Active View
    render = pvs.Render()
    # Screenshot 
    pvs.WriteImage('solution.png')
    # Display result in user's screen
    pvs.Interact()
except:
    sys.exit('*** Paraview plotting solution error ***')
