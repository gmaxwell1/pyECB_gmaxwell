#%%
import numpy as np
import os

from modules.analysis_tools import extract_raw_data_from_file, generate_plots


#%%
data_directory = './data_sets/z_field_meas_set_2'
data_filename = '20_10_06_17-10-30_B_vs_I_along_z.csv'
data_filepath = os.path.join(data_directory, data_filename)

I, mean_data_specific_sensor, std_data_specific_sensor, expected_fields = extract_raw_data_from_file(data_filepath)

generate_plots(I, mean_data_specific_sensor, std_data_specific_sensor, expected_fields, 
                flag_xaxis='I', flags_yaxis='zmpa', plot_delta_sim=False, save_image=False,
                directory=data_directory, height_per_plot=2)


