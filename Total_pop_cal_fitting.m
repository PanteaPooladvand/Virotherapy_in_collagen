function [Diff] = Total_pop_cal_fitting

% Use one the the following to fit Experiment 2 to data. 
%p = [aa 0.85 25; aa 1/3 20]; % Use this to fit PDB and Collagenase data [r comp collagen_density t_end]
%p = [aa 0.85 30 ;aa 1/3 52]; % Use this to fit beta and deltaV simulataneously to both high and low collagen data [beta deltaV collagen_density t_end]



Diff = [];

UTotal = [];
VTotal = [];
ITotal = [];
ETotal = [];

load('Mckee_treatment_data.mat')


for i=1:1
    


P = p1(i,:); % Also use this to fit beta and deltaV simultaneously to high and low collagen data*******


[x, t, sol] = Virus_diff_model2(P); %PDE with modified diffusion

% Calculating tumour population

trapintegrationvector = ones(length(x), 1);
trapintegrationvector(1) = 0.5;
trapintegrationvector(end) = 0.5;
xmatrix = diag(x);

% Assume equal spacing of r
dr = (x(end) - x(1))/(length(x) - 1);

% This is like integral of N r dr dtheta
intmatrix = xmatrix.^2*trapintegrationvector*dr*4*pi;

U=sol(:,:,1);           % tumour cell solutions
V=sol(:,:,2);           % virus
I=sol(:,:,3);           % infected cells
E=sol(:,:,4);           % ECM


UT = U*intmatrix;       % calculate totals
VT = V*intmatrix;
IT = I*intmatrix;
ET = E*intmatrix;


UT(UT<0)=0;             % eliminate negative solutions due to numerics
VT(VT<0)=0;
IT(IT<0)=0;


% use the next 2 lines to fit model to PBS and collagen data simultaneously************
%untreated_growth_data = {PBS ; Collagenase}; 
%A = repmat(untreated_growth_data{i}(:,1),[1,length(t)]); 


% use the next 2 lines to fit model to virotherapy data simultaneously********
treated_growth_data = {Mckee_virus ; Mckee_virus_collag}; 
A = repmat(treated_growth_data{i}(:,1),[1,length(t)]); 

[~,closestIndex] = min(abs(A'-t'));
Tvolnew = UT([closestIndex]); % this array contains the output values from the model that match those from the virus data at the given timepoints.
TumourVol = Tvolnew/Tvolnew(1);

TumourVolDiff = (TumourVol - treated_growth_data{i}(:,2));


Diff = [Diff ; TumourVolDiff];


  figure
  plot(t,(UT)/(UT(1)),'Linewidth',4)
  hold on
  scatter(treated_growth_data{i}(:,1),treated_growth_data{i}(:,2),80)
  set(gca,'FontSize',40,'fontweight','bold')
  xlabel('Time (days)')
  ylabel('Tumour fold change')
  axis square
  box off
  hold off

end


end
