% Matlab build-in solver lsqnonlin 
% Use this to fit Virotherapy model to Mckee data
% This is a simple code that calls Total_pop_cal_fitting

%InitialGuess = [2.4614 0.2916]; % initial guess for squeeze prob and scalling factor.
%InitialGuess = [0.2127    0.2732];  % initial guess for tumour growth in pbs and collagenase 
InitialGuess = [3.5*10^(-8) 0.079];% initial guess fbor infection rate beta and virus decay delta_v

lb = [0 0]; % restricting lower boundary to 0
ub = [ ];
options=optimoptions('lsqnonlin','display','iter-detailed','FiniteDifferenceStepSize',1e-10);
[aa,resnorm,residual,exitflag,output,lambda,jacobian]=lsqnonlin(@Total_pop_cal_fitting,InitialGuess,lb,ub);

conf = nlparci(aa,residual,'jacobian',jacobian) % confidence intervals

