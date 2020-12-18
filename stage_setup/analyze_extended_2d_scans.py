"""
filename: analyze_extended_2d_scans.py

This script is meant to analyze data originating from four 2d scans of the magnetic field,
where the current configurations are (111), (100), (010) and (001) up to an overall factor.
Before each scan, the coils have to be demagnetized, such that the resulting field lies on 
a virgin hysteresis curve. 


Author: Nicholas Meinhardt (QZabre)
        nmeinhar@student.ethz.ch

Date: 23.11.2020
"""

#%%
# standard library imports
import numpy as np
import serial
from time import sleep, time
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
from scipy.optimize import leastsq, curve_fit

# local imports
# local imports
try:
    from modules.general_functions import transform_between_sensor_stage_coordinates
except ModuleNotFoundError:
    import sys
    sys.path.insert(1, os.path.join(sys.path[0], '..'))
finally:
    from modules.general_functions import transform_between_sensor_stage_coordinates
    from modules.analysis_tools import extract_raw_data_from_2d_scan, plot_2d_scan, get_relative_inplane_angles, get_direction_vector


#%%
# retrieve data from previous scans
directory = '../test_data/2d_scans_different_fields/set1'

# filenames = ['20_11_23_18-32-01_2d_scan_(1_0_0).csv',   # (100)   
#             '20_11_23_18-43-38_2d_scan_(0_1_0).csv',    # (010) 
#             '20_11_24_08-35-21_2d_scan_(0_0_1).csv',    # (001)   
#             '20_11_24_08-47-10_2d_scan_(1_1_1).csv']    # (111)  
filenames = ['20_11_23_17-10-07_2d_scan_(1_0_0).csv',   # (100)   
            '20_11_23_17-37-25_2d_scan_(0_1_0).csv',    # (010) 
            '20_11_23_17-56-46_2d_scan_(0_0_1).csv',    # (001)   
            '20_11_23_18-16-10_2d_scan_(1_1_1).csv']    # (111)  

coils = []
positions_all = []
B_field_all = []
filepath_all = []

for data_filename in filenames:
    data_filepath = os.path.join(directory, data_filename)

    positions, B_field = extract_raw_data_from_2d_scan(data_filepath)

    positions_all.append(positions)
    B_field_all.append(B_field)
    filepath_all.append(data_filepath)

positions_all = np.array(positions_all)
B_field_all = np.array(B_field_all)
filepath_all = np.array(filepath_all)


# %%
# using scans with only a single coil switched on, estimate angles and plot them wrt to index
N_pts = positions_all.shape[1]

angles = np.zeros((N_pts,3))
for i in range(N_pts):
    directions = np.array([[np.zeros(3), B_field_all[0, i]], 
                        [np.zeros(3), B_field_all[1, i]], 
                        [np.zeros(3), B_field_all[2, i]]])

    
    angles[i] = get_relative_inplane_angles(directions)

# estimate sensor position with smallest distance between estimated angles and 120°
i_chosen = np.argmin(np.linalg.norm(np.abs(angles-120), axis=1))
print('angles closest to 120°: {}'.format(np.round(angles[i_chosen], 3)))

stage_pos = transform_between_sensor_stage_coordinates(positions_all[0, i_chosen])
print('Corresponding stage position: {}'.format(np.round(stage_pos, 3)))

# plot angles wrt to iteration 
fig, ax = plt.subplots()
ax.axhline(120, color='k', alpha=0.5)

ax.plot(np.arange(N_pts), angles[:,0], label='angle btw. coils 1+2')
ax.plot(np.arange(N_pts), angles[:,1], label='angle btw. coils 1+3')
ax.plot(np.arange(N_pts), angles[:,2], label='angle btw. coils 2+3')

ax.axvline(i_chosen, color='r', linestyle='--')

ax.legend()
ax.set_xlabel('number iterations')
ax.set_ylabel('angles between coils, $\\alpha$ [°]')

plt.tight_layout()

# save image
image_path = os.path.join(directory, 'relative_angles.png')
fig.savefig(image_path, dpi=300)

plt.show()



#%%
# plot measured fields

# choose which component should be plotted
plot_component = 'z'

# set the center_position to None or [x,y] in stage coordinate system
# center_position = [5.0, 15.9]
center_position = stage_pos[:2]

for i in range(4):
    # generate figure
    fig, ax = plot_2d_scan(positions_all[i], B_field_all[i],  Cont=True, Scat_Mag=False, levels=None,  
                        plot_component=plot_component, center_position=center_position)

    # show legend 
    ax.legend(loc='lower right')

    # save figure
    image_path = '{}_2dplot_{}.png'.format(os.path.splitext(filepath_all[i])[0], plot_component)
    fig.savefig(image_path, dpi=300)

    plt.show()




# %%
