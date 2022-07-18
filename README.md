# Virotherapy_in_collagen
These are the scripts to run virotherapy PDE model 

Discription of files:

Virus_diff_model2.m 
This is the main code that contains the model and the pdepe solver used to esimate the solutions to the model. If you do not want to input values, you can run the code with an arbitrary value for P, eg. Virus_diff_model2(1). If you want to input a value or an array, then assing this to P.

All the initial conditions used in the paper are included in the initial conditions section of this code.


To fit the model to data, we run Mckee_nonlin_fitter.m. This code calls Total_pop_cal_fitting which then calls Virus_diff_model2.m.
Note that Total_pop_cal_fitting requires the data Mckee_treatment_data.mat. This code estimates the tumour fold change in time. Please also ensure that you are using the correct initial conditions in Virus_diff_model2.m. The initial conditions here are Experiment 2.  
