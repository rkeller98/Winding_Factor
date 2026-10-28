close all; 
clc; 

eval_set = load_data(); 
meas = eval_set.measurement; 

[kde, meta] = k_factor_contour_calculation(meas);
fig = k_factor_contour_analysis_plot(meas, meta);

