#!/usr/bin/env python

import numpy as np
import os.path
import os
import pandas as pd

from landlab.components import SoilInfiltrationGreenAmpt
from landlab.components import SinkFiller 
from landlab.io import read_esri_ascii
from landlab import imshow_grid
from landlab.components.overland_flow import OverlandFlow
from landlab import RasterModelGrid, CLOSED_BOUNDARY, FIXED_VALUE_BOUNDARY
from matplotlib import pyplot as plt

__all__ = ['loadGrid','runModel','listFiles','getLLFilename','decodeLLFilename',
            'writeLLData','readLLData','filesToDf']



# SoilInfiltrationGreenAmpt
#
# defaults
#    hydraulic_conductivity=0.005 [The soil effective hydraulic conductivity m/s]
#    soil_bulk_density=1590. [The dry bulk density of the soil kg/m**3]
#    rock_density=2650. [The density of the soil constituent material (i.e., lacking porosity) kg/m**3]
#    initial_soil_moisture_content=0.15 [The fraction of the initial pore space filled with water m**3/m**3, 0. to 1.] 
#    soil_type='sandy loam',
#    volume_fraction_coarse_fragments=0.2 [The fraction of the soil made up of rocky fragments with 
#                                          very little porosity, with diameter > 2 mm m**3/m**3, 0. to 1.]
#    coarse_sed_flag=False [If this flag is set to true, the fraction of coarse material in the
#                           soil column with be used as a correction for phi, the porosity factor.]
#    surface_water_minimum_depth=1.e-8 [A minimum water depth to stabilize the solutions for surface flood
#                                       modelling. Leave as the default in most normal use cases m]
#    soil_pore_size_distribution_index [index describing the distribution of pore sizes in the soil,
#            and controlling effective hydraulic conductivity at varying water contents, 
#            following Brooks and Corey (1964). Can be set by soil_type. Typically denoted "lambda"]
#    soil_bubbling_pressure [The bubbling capillary pressure of the soil, controlling effective
#            hydraulic conductivity at varying water contents, following Brooks
#            and Corey (1964). Can be set by soil_type. Typically denoted "h_b"]
#    wetting_front_capillary_pressure_head [The effective head at the wetting front in the soil driven by
#            capillary pressure in the soil pores. If not set, will be calculated by the component
#            from the pore size distribution and bubbling pressure, following Brooks and Corey]
         
# A soil type to automatically set soil_pore_size_distribution_index
#            and soil_bubbling_pressure, using mean values from Rawls et al.,1992.        
#    SOIL_PROPS = {
#        'sand': (0.694, 0.0726),
#        'loamy sand': (0.553, 0.0869),
#        'sandy loam': (0.378, 0.1466),
#        'loam': (0.252, 0.1115),
#        'silt loam': (0.234, 0.2076),
#        'sandy clay loam': (0.319, 0.2808),
#        'clay loam': (0.242, 0.2589),
#        'silty clay loam': (0.177, 0.3256),
#        'sandy clay': (0.223, 0.2917),
#        'silty clay': (0.150, 0.3419),
#        'clay': (0.165, 0.3730),
#    }

# boundary conditions [right, left, top, bottom]
boundaryDict = {
'rightOpen': [FIXED_VALUE_BOUNDARY, CLOSED_BOUNDARY, CLOSED_BOUNDARY, CLOSED_BOUNDARY],
'leftOpen': [CLOSED_BOUNDARY, FIXED_VALUE_BOUNDARY, CLOSED_BOUNDARY, CLOSED_BOUNDARY],
'topOpen': [CLOSED_BOUNDARY, CLOSED_BOUNDARY, FIXED_VALUE_BOUNDARY, CLOSED_BOUNDARY],
'bottomOpen': [CLOSED_BOUNDARY, CLOSED_BOUNDARY, CLOSED_BOUNDARY, FIXED_VALUE_BOUNDARY],
'allOpen': [FIXED_VALUE_BOUNDARY, FIXED_VALUE_BOUNDARY, FIXED_VALUE_BOUNDARY, FIXED_VALUE_BOUNDARY],
'allClosed': [CLOSED_BOUNDARY, CLOSED_BOUNDARY, CLOSED_BOUNDARY, CLOSED_BOUNDARY],
}

###########################################################################
def loadGrid(topopath, dem):
        
    # load topography to grid

    (rmg, z) = read_esri_ascii(os.path.join(topopath,dem), name='topographic__elevation')
 
    return rmg, z

###########################################################################
def runModel(topopath, dem, 
            soil_type, min_h, max_h, 
            starting_precip_mmhr, storm_duration,
            soilWaterInfiltrationDepth, surfaceWaterDepth,
            mannings_n, 
            boundary, monitorLink, 
            runHours, reportInterval, 
            showProgress = False, showPlots = False, 
            soilInf = True, sinkFiller = False, 
            noDataValue=0,outlet_node=0):

    rmg, z =  loadGrid(topopath, dem)

    if 'noData' in boundary:
        if noDataValue != 0:
            no_data = noDataValue
        else:
            no_data=min(z)
        rmg.set_watershed_boundary_condition(z,nodata_value=no_data)
        outlet_node = np.where(rmg.status_at_node==1)
    elif 'setSelf' in boundary:
        rmg.set_watershed_boundary_condition_outlet_id(outlet_node,z)
    else:
        for i, edge in enumerate([rmg.nodes_at_right_edge, rmg.nodes_at_left_edge, rmg.nodes_at_top_edge, rmg.nodes_at_bottom_edge]):
            rmg.status_at_node[edge] = boundaryDict[boundary][i]

    
    # Initial values in grid:

    # Depth of water above the surface
    h = rmg.add_ones('node', 'surface_water__depth')
    h *= surfaceWaterDepth

 
    # Water column height above the surface previously absorbed into the soil. 
    # Note that this is NOT the actual depth of the wetted front, which also depends on the porosity.
    d = rmg.add_ones('node', 'soil_water_infiltration__depth', dtype=float)
    d *= soilWaterInfiltrationDepth

    # Hydraulic Conductivity
    # the ease with which a fluid (usually water) can move through pore spaces or fractures.
    # typical values for conductivity: http://www.aqtesolv.com/aquifer-tests/aquifer_properties.htm
    # randomize hydraulic conductivity
    hydraulic_conductivity = rmg.ones('node')*min_h
    hydraulic_conductivity += np.random.rand(z.size)*(max_h - min_h)  
        
    # sinkfiller
    if sinkFiller:
        sf = SinkFiller(rmg, routing='D4', apply_slope=True, fill_slope=1.e-5)
        sf.fill_pits()
    
    # OverlandFlow Component
    # The OverlandFlow component is built off the de Almeida et al., (2012) algorithm for urban flood inundation, 
    # and is most stable in flat environments. Because of this, instabilities can arise when trying to apply the 
    # algorithm to steep landscapes. To adapt this model for use across a variety of terrains, stability criteria 
    # (following Coulthard et al., 2013) is implemented to using the steep_slopes flag.
    # steep_slopes = True: This method reduces flow discharge to keep flow subcritical according to the Froude number less than or equal to 1.0. 
    #                      For more information, see Adams et al., (in prep for Geoscientific Model Development).
    # of = OverlandFlow(rmg, alpha = alpha, mannings_n = mannings_n, steep_slopes = True, default_fixed_links = True)
    of = OverlandFlow(rmg, mannings_n = mannings_n, steep_slopes = True)

    # SoilInfiltrationGreenAmpt component
    if soilInf:
        SI = SoilInfiltrationGreenAmpt(rmg, soil_type = soil_type, hydraulic_conductivity=hydraulic_conductivity)

    # storm [m/s]
    # converted from mm hr-1 to m s-1
    starting_precip_ms = starting_precip_mmhr * (2.77778 * 10 ** -7)

    # Simulation run time and reporting
    model_run_time = 3600.0 * runHours   
    report_interval = 3600. * reportInterval

    if model_run_time < storm_duration:
        model_run_time = storm_duration * 1.5             

    ## Run simulation:
    # for model stability we do not start at time 0
    elapsed_time = 1.
    next_report = 0.  
    fig = 1

    rmg['link']['water__discharge'] = np.zeros(rmg.number_of_links) 
    rmg['node']['water__depth'] = np.zeros(rmg.number_of_nodes)
    rmg['node']['surface_water__discharge'] = np.zeros(rmg.number_of_nodes)

    # Lists for saving data
    discharge_at_outlet = []
    hydrograph_time = []

    while elapsed_time < model_run_time:
        # Setting the adaptive time step
        of.dt = of.calc_time_step()

        # The storm starts when the model starts. While the elapsed time is less
        # than the storm duration, we add water to the system as rainfall.
        if elapsed_time < (storm_duration):
            of.rainfall_intensity =  starting_precip_ms   
        else: # elapsed time exceeds the storm duration, rainfall ceases.
            of.rainfall_intensity = 0.0

        of.overland_flow() # Generating overland flow based on the deAlmeida solution.

        # Variables listed here are updated by the component at the grid locations listed.
        # surface_water__discharge, link, [m^2 s^-1] : At each link in grid, surface_water__discharge is calculated using the de Almeida et al., (2012) equation. Discharge is a function of the water depth, adaptive time step, surface water slope and Manning’s roughness coefficient.
        # surface_water__depth, node, [m] : At each node in the grid, surface_water__depth is updated using the surface_water__discharge on links connected to a given node.    
        
        # To translate the discharge values calculated on Landlab links to nodes, 
        # values on links (of.q) are summed and mapped to their node neighbors using 
        # the method of.discharge_mapper. Using the convert_to_volume flag, 
        # these discharge values are converted from units of [m2 s-1] to [m3 s-1].
        #rmg.at_node['surface_water__discharge'] = of.discharge_mapper(of.q, convert_to_volume=True)

        # Append time and discharge to their lists to save data and for plotting.
        hydrograph_time.append(elapsed_time/3600.)
        discharge_at_outlet.append(np.abs(of.q[monitorLink]) * rmg.dx) # append discharge in m^3/s
         
        if soilInf:
            SI.run_one_step(of.dt)
        
        # Report progress
        
        if elapsed_time >= next_report:
            if showProgress:
                print('elapsed time = {:.2f} s'.format(elapsed_time))
            next_report += report_interval
        
            if showPlots:
                plt.figure(fig)
                fig = fig + 1
                imshow_grid(rmg, 'surface_water__depth', cmap='Blues')
                plt.title('soil_water_depth (m) at {:.2f} s'.format(elapsed_time))
                plt.figure(fig)
                fig = fig + 1
                imshow_grid(rmg, 'soil_water_infiltration__depth', cmap='Blues')
                plt.title('soil_water_infiltration__depth (m) at {:.2f} s'.format(elapsed_time))

        # Updating elapsed_time
        elapsed_time += of.dt
        
    return rmg, hydraulic_conductivity, discharge_at_outlet, hydrograph_time, h, d
    
    
def listFiles(root, patterns='*', recurse=1, return_folders=0, useRegex=False):
    """Lists the files/directories meeting specific requirement
        Returns a list of file paths to files in a file system, searching a 
        directory structure along the specified path, looking for files 
        that matches the glob pattern. If specified, the search will continue 
        into sub-directories.  A list of matching names is returned. The 
        function supports a local or network reachable filesystem, but not URLs.
    Args:
        | root (string): directory root from where the search must take place
        | patterns (string): glob/regex pattern for filename matching. Multiple pattens 
          may be present, each one separated by ;
        | recurse (unt): flag to indicate if subdirectories must also be searched (optional)
        | return_folders (int): flag to indicate if folder names must also be returned (optional)
        | useRegex (bool): flag to indicate if patterns areregular expression strings (optional)
    Returns:
        | A list with matching file/directory names
    Raises:
        | No exception is raised.
    """
    import sys
    import fnmatch    
    if useRegex:
        import re

    # Expand patterns from semicolon-separated string to list
    pattern_list = patterns.split(';')
    filenames = []
    filertn = []

    if sys.version_info[0] < 3:

        # Collect input and output arguments into one bunch
        class Bunch(object):
            def __init__(self, **kwds): self.__dict__.update(kwds)
        arg = Bunch(recurse=recurse, pattern_list=pattern_list,
                                return_folders=return_folders, results=[])

        def visit(arg, dirname, files):
            # Append to arg.results all relevant files (and perhaps folders)
            for name in files:
                fullname = os.path.normpath(os.path.join(dirname, name))
                if arg.return_folders or os.path.isfile(fullname):
                    for pattern in arg.pattern_list:
                        if useRegex:
                            #search returns None is pattern not found
                            regex = re.compile(pattern)
                            if regex.search(name):
                                arg.results.append(fullname)
                                break
                        else:
                            if fnmatch.fnmatch(name, pattern):
                                arg.results.append(fullname)
                                break
            # Block recursion if recursion was disallowed
            if not arg.recurse: files[:]=[]
        os.path.walk(root, visit, arg)
        return arg.results

    else: #python 3
        for dirpath,dirnames,files in os.walk(root):
            if dirpath==root or recurse:
                for filen in files:
                    # filenames.append(os.path.abspath(os.path.join(os.getcwd(),dirpath,filen)))
                    filenames.append(os.path.relpath(os.path.join(dirpath,filen)))
                if return_folders:
                    for dirn in dirnames:
                        # filenames.append(os.path.abspath(os.path.join(os.getcwd(),dirpath,dirn)))
                        filenames.append(os.path.relpath(os.path.join(dirpath,dirn)))

        for name in filenames:
            if return_folders or os.path.isfile(name):
                for pattern in pattern_list:
                    if useRegex:
                        #search returns None is pattern not found
                        regex = re.compile(pattern)
                        if regex.search(name):
                            filertn.append(name)
                            break
                    else:
                        # split only the filename to compare, discard folder path
                        if fnmatch.fnmatch(os.path.basename(name), pattern):
                            filertn.append(name)
                            break
    return filertn



def getLLFilename(Slope,ManningN,HCond,InfilDepth):
    """Create filename from input parameters
    """
    return '{}_{:.3e}_{:.3e}_{:.3e}.npz'.format(Slope,ManningN,HCond,InfilDepth)

def decodeLLFilename(filename):
    """Decode parameters from filename
    """
    filename = filename.split('\\')[-1]
    filename = filename.replace('.npz','')
    lstf = filename.split('_')
    Slope = lstf[0]
    ManningN = float(lstf[1])
    HCond = float(lstf[2])
    InfilDepth =float(lstf[3])
    return Slope,ManningN,HCond,InfilDepth

def writeLLData(filename,hydrograph_time,discharge_at_outlet,h,d):
    """Write four arrays to npz file
    """
    np.savez(filename,
             hydrograph_time = hydrograph_time,
             discharge_at_outlet = discharge_at_outlet,
             h = h,
             d = d)
    
def readLLData(filename):
    """Read four parameters from npz file
    """
    ndic = np.load(filename)
    return ndic['hydrograph_time'],ndic['discharge_at_outlet'],ndic['h'],ndic['d']

def filesToDf(root):
    df = pd.DataFrame()
    filenames = listFiles(root, patterns='*.npz', recurse=0, return_folders=0, useRegex=False)
    for filename in filenames:
        Slope,ManningN,HCond,InfilDepth = decodeLLFilename(filename)
        df = df.append(pd.DataFrame({
            'Slope':Slope,
            'ManningN':ManningN,
            'HCond':HCond,
            'InfilDepth':InfilDepth,
            'filename':filename}, index=[0]))
        
    return df
        
