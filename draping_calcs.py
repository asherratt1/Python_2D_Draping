import pandas as p
import numpy as np
import math
    
def read_cell_centers_2D(df):
    df = p.read_csv(df, skiprows = 1, header = None, names = ['X', 'Y'] )
    cellcenter = df.to_numpy()
    return cellcenter
        
def read_cell_centers_3D(df):
    df = p.read_csv(df, skiprows = 1, header = None, names = ['X', 'Y', 'Z'] )
    cellcenter = df.to_numpy()
    return cellcenter

def read_normal_vectors(df):
    df = p.read_csv(df, skiprows = 1, header = None, names = ['N1', 'N2', 'N3', 'X', 'Y', 'Z'] )
    cellcenter = df.loc[:,'X' : 'Z'].to_numpy()
    vecnorm = df.loc[:,'N1' : 'N3'].to_numpy()
    return cellcenter, vecnorm

def write_OF_file(drape_dir, template, file, boundaries, condition):
    # Convert from float to string for file writing
    drape_string = drape_dir.astype(str)
    drape_str = p.DataFrame(drape_string,dtype = "string")

    # Copy template over to new file
    with open(template) as f:
        lines = f.readlines()
        with open(file,"w") as f1:
            f1.writelines(lines)

    # Write OpenFOAM file and drape data:
    with open(file,"a") as f1:
        
        f1.write('\n')
        f1.write('(\n')
        for i in range(len(drape_dir)):
            f1.write('(' + drape_str[0][i] + ' ' + drape_str[1][i] + ' ' + drape_str[2][i] + ')\n')
        
    # Write boundary conditions:
    with open(file,"a") as f2:
        f2.write(');\n')
        f2.write('boundaryField\n')
        f2.write('{\n')
        for i,bnd in enumerate(boundaries):
            f2.write(bnd+'\n')
            f2.write('{\n')
            f2.write('type '+condition[i]+';\n')
            f2.write('}\n')
        
    with open(file,"a") as f3:
        f3.write('}\n')
    
def min_dist(v1, v2):
    dist = np.subtract(v1, v2)
    vec_norm_min = np.amin(np.linalg.norm(dist, axis = 1))
    index = np.argmin(np.linalg.norm(dist, axis = 1))
    return vec_norm_min, index
        
def vec_min_index(v1, v2):
    dist = np.subtract(v1, v2)
    vec_norm = np.linalg.norm(dist, axis = 1)
    vec_min_location = np.argpartition(vec_norm, 2)
    return vec_min_location[0:2]
        
def calc_layer_thickness(domain_thickness, drape_sequence):
    layer_thickness = domain_thickness / len(drape_sequence)
    return layer_thickness
        

def calc_thickness_2D(cc_file, um_file, lm_file, outline):
        
    # Read file data
    cell_centers = read_cell_centers_3D(cc_file)
    upper_centers, upper_normal = read_normal_vectors(um_file)
    lower_centers, lower_normal = read_normal_vectors(lm_file)
        
    thickness_direction = np.zeros((len(cell_centers), 3))
        
    for num, cells in enumerate(cell_centers, start = 0):
        min_dist_upper, index_upper = min_dist(upper_centers, cells)
        min_dist_lower, index_lower = min_dist(lower_centers, cells)
            
        # Calculate weighting factor for distance:
        dist_total = min_dist_upper + min_dist_lower
            
        w1 = min_dist_upper / dist_total
        w2 = min_dist_lower / dist_total
            
        upper_vec = upper_normal[index_upper]
        lower_vec = lower_normal[index_lower] * -1
            
        norm_add = np.add((1-w1) * upper_vec, (1-w2) * lower_vec)
            
        thickness_direction[num] = norm_add / np.linalg.norm(norm_add)
            
    return thickness_direction
        
    
def calc_normal_distance_square(cc_file, um_file, domain_thickness, drape_sequence):

    # Read file data
    cell_centers = read_cell_centers_3D(cc_file)

    upper_centers, upper_normal = read_normal_vectors(um_file)
        
    layer_thickness = calc_layer_thickness(domain_thickness, drape_sequence)
        
    layer_devisions = np.linspace(1, len(drape_sequence), len(drape_sequence)) * layer_thickness
        
    min_distance = np.zeros((len(cell_centers), 1))
    layer_location = np.zeros((len(cell_centers), 1))
        
    # Get smallest distance:
    for num, cells in enumerate(cell_centers, start = 0):
        #vec_min_dist = domain_thickness - cells[2]
        vec_min_dist = cells[2]
        min_distance[num] = vec_min_dist
            
        #for index, layers in enumerate(layer_devisions):
        #if vec_min_dist <= layer_devisions[0]:
            #layer_location[num] = 0
        if vec_min_dist > layer_devisions[0] and vec_min_dist <= layer_devisions[1]:
            layer_location[num] = 1
            print('True')
        #elif vec_min_dist > layer_devisions[1] and vec_min_dist <= layer_devisions[2]:
            #layer_location[num] = 2
    return min_distance, layer_location
        
def calc_normal_distance_2D(cc_file, um_file, outline_2D, domain_thickness, drape_sequence):
        
    # Read file data
    cell_centers = read_cell_centers_3D(cc_file)
    upper_centers, upper_normal = read_normal_vectors(um_file)
    upper2D = read_cell_centers_2D(outline_2D)
        
    layer_thickness = calc_layer_thickness(domain_thickness, drape_sequence)
        
    layer_devisions = np.linspace(1, len(drape_sequence), len(drape_sequence)) * layer_thickness
        
    min_distance = np.zeros((len(cell_centers), 1))
    layer_location = np.zeros((len(cell_centers), 1))
        
    for num, cells in enumerate(cell_centers, start = 0):
        point_index = vec_min_index(upper2D, cell_centers[[num],:][:,[0,1]])
            
        xy_point1 = upper2D[point_index[0]]
        xy_point2 = upper2D[point_index[1]]
            
        slope1 = (xy_point2[1] - xy_point1[1]) / (xy_point2[0] - xy_point1[0])
            
        b1 = xy_point2[1] - slope1 * xy_point2[0]
            
        if slope1 == 0: #two points on surface have same y coordinate
            dist = xy_point2[0] - cell_centers[num][0]
            if dist < 0:
                dist = dist * -1
        elif np.isnan(slope1) == True:
            dist = xy_point2[1] - cell_centers[num][1]
            if dist < 0:
                dist = dist * -1
        else:
            slope2 = -1 / slope1
            b2 = cell_centers[num][1] - slope2 * cell_centers[num][0]
                
            x_intercept = (b2 - b1) / (slope1 - slope2)
            y_intercept = slope2 * x_intercept + b2
            sol = [x_intercept, y_intercept]
                
            dist = np.linalg.norm(np.subtract(sol, cell_centers[[num],:][:,[0,1]]))
                
            if dist < 0:
                dist = dist * -1
            
        min_distance[num] = dist
                    
        for index, layers in enumerate(layer_devisions):
            if dist <= layers:
                layer_location[num] = index
                break;
                   
                    
    return min_distance, layer_location
         
        
def calc_drape_angle(thickness_direction, layer_location, drape_sequence, extrude_direction, domain_type):
    
    drape_angles = np.zeros((len(layer_location),3))
        
    surface_parallel = np.cross(thickness_direction, extrude_direction)
        
    for index, x in enumerate(layer_location):
            
        angle_deg = drape_sequence[int(x.item())]
        angle_rad = math.radians(angle_deg)
        
        if domain_type == 'square':
            y_comp = math.sin(angle_rad) * np.asarray(surface_parallel)
            x_comp = math.cos(angle_rad) * np.asarray(extrude_direction)
        elif domain_type == '2D':
            y_comp = math.sin(angle_rad) * np.asarray(surface_parallel[index])
            x_comp = math.cos(angle_rad) * np.asarray(extrude_direction)
            
        drape_angles[index] = np.add(x_comp, y_comp) / np.linalg.norm(np.add(x_comp, y_comp))
            
    return drape_angles
        
        
        
