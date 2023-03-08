## 2D Extruded Geometry Draping Code
# Developed by: Rex Sherratt
# email: asherra@uwo.ca

# Import packages
import pandas as p
import numpy as np
import matplotlib.pyplot as plt
import math
import time

# Import functions
from draping_calcs import *

# Specify draping data
drape_sequence = [45.]
domain_thickness = 0.003 # Remember to match this with domain units

# Specify domain attributes
extrude_direction = [0, 0, 1]
domain_type = '2D' # Currently can choose from 'square' or '2D'

if domain_type == 'square':
    print('Domain is set to square, specifying constant thickness direction')
    thickness_direction = [0., 0., 1.]

# Read in mesh details
internal_mesh = 'hatChannel_input/cellCenters.csv'
upper_mesh = 'hatChannel_input/upperSurfaceNormals.csv'
lower_mesh = 'hatChannel_input/lowerSurfaceNormals.csv'

# Template file
template_drape = 'template_files/drapeDirectionTemplate'
template_thickness = 'template_files/thicknessDirectionTemplate'

# Output location:
write_file_drape = 'hatChannel_output/drapeDirection'
write_file_thickness = 'hatChannel_output/thicknessDir'

# Boundary Conditions
boundaries = ['uppermold', 'lowermold', 'wall']
conditions = ['zeroGradient','zeroGradient','zeroGradient']

# if 2D geometry
outline_2D = 'hatChannel_input/hatChannel_2D_profile.csv'

# Execute draping code
if domain_type == 'square':
    
    # Calculate normal distance from upper surface
    print('Calculating normal distance from upper surface')
    normal_distance, layer_location = calc_normal_distance_square(internal_mesh, upper_mesh, domain_thickness, drape_sequence)
    
    # Specify draping angle
    print('Calculating drape angle')
    drape_angles = calc_drape_angle(thickness_direction, layer_location, drape_sequence, extrude_direction, domain_type)
    
    # Write to file
    print('Writing to file Drape Dir')
    write_OF_file(drape_angles, template_drape, write_file_drape, boundaries, conditions)
    
    # Write to file
    print('Writing to file thicknessDir')
    write_OF_file(thickness_direction, template_thickness, write_file_thickness, boundaries, conditions)
    
elif domain_type == '2D':
    print('Calculating thickness direction')
    thickness_direction = calc_thickness_2D(internal_mesh, upper_mesh, lower_mesh, outline_2D)
    
    print('Calculating normal distance from upper surface')
    normal_distance, layer_location = calc_normal_distance_2D(internal_mesh, upper_mesh, outline_2D, domain_thickness, drape_sequence)
    
    # Specify draping angle
    print('Calculating drape angle')
    drape_angles = calc_drape_angle(thickness_direction, layer_location, drape_sequence, extrude_direction, domain_type)
    
    # Write to file
    print('Writing to file Drape Dir')
    write_OF_file(drape_angles, template_drape, write_file_drape, boundaries, conditions)
    
    # Write to file
    print('Writing to file thicknessDir')
    write_OF_file(thickness_direction, template_thickness, write_file_thickness, boundaries, conditions)




