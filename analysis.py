#%%
import numpy as np
import os
from scipy.optimize import curve_fit
kj
from modules.analysis_tools import *


#%%
data_directory = './test_data/xy_field_meas_set_2'
data_filename = '20_10_08_14-18-22_B_field_vs_I_new.csv'
data_filepath = os.path.join(data_directory, data_filename)

I, mean_data_specific_sensor, std_data_specific_sensor, expected_fields = extract_raw_data_from_file(data_filepath)


_ = generate_plots(I, mean_data_specific_sensor, std_data_specific_sensor, expected_fields, 
                flag_xaxis='I1', flags_yaxis='zmptf', plot_delta_sim=False, save_image=False,
                directory=data_directory, height_per_plot=2, show_labels=True, distance=1.7,
                xlim=None, ylim_field_abs=(0,80), ylim_field_z=(-80,80))

plot_I_vs_B(I, mean_data_specific_sensor, std_data_specific_sensor, expected_fields, directory=data_directory,
                save_image=False, xlim=(0,90), ylim=None)


#%%

def fit_B_vs_single_current(I_single_coil, B_magnitudes):

    # fit functions
    fit_functions_collection = [abs_sigmoid, abs_brillouin_fct, abs_lin_and_const]
    fit_labels = ['sigmoid', 'brillouin', 'linconst']

    params = []
    params_err = []

    for k in range(len(fit_functions_collection)):
        p, p_cov = curve_fit(fit_functions_collection[k], I_single_coil, B_magnitudes, p0=[1,50])
        p_err = np.sqrt(np.diag(p_cov))

        params.append(p)
        params_err.append(p_err)

    predictions = np.array([params[0][1], params[1][1], params[2][0]*params[2][1]])
    predictions_errors = np.array([params_err[0][1], params_err[1][1], 
                                np.sqrt((params[2][1]*params_err[2][0])**2 + (params[2][0]*params_err[2][1])**2)])

    print('p = {} +- {}'.format(np.round(predictions,2), np.round(predictions_errors,2)))
    selected_fit = np.argmin(predictions_errors/predictions)

    B_fit = fit_functions_collection[selected_fit](I_single_coil, *params[selected_fit])

    B_saturation = predictions[selected_fit]
    B_saturation_error = predictions_errors[selected_fit]

    return abs(B_saturation), B_saturation_error, B_fit, fit_labels[selected_fit]


mean_magnitudes = np.linalg.norm(mean_data_specific_sensor, axis=1)

for i in range(3):
    print('I_{}'.format(i+1))
    try:
        B_sat, B_sat_error, B_fit, which_fit = fit_B_vs_single_current(I[:,i], mean_magnitudes)
        print('{:.2f} mT +- {:.2f} mT ({})'.format(B_sat, B_sat_error, which_fit))
    except RuntimeError:
        print('no saturation could be estimated, take maximum value of B_mag = {:.2f}'.format(np.max(mean_magnitudes)))
    




#%%
fig, ax = plt.subplots()

x = np.linspace(-5,5,num=100)
y = x**2
ax.plot(x,y)
ax.set_xlim(None)
fig.show()
