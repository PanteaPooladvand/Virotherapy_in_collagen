% This is the code that runs the PDE virotherapy model with both Fickian and
% non-Fickian diffusion terms


function [x, t, sol] = Virus_diff_model2(P) % model with explicit collagen 

    m = 2;  % Geometry of the system, 0=Cartisian, 2=Spherical
    L=15;   % Domain is 60 for experiment 1 and 15 for experiment 2

 
    dx = 0.1; % spatial increments
    dt =0.1; % time increments
    
    x = 0:dx:L;  % Set domeain use x=-60:dx:L for Experiment 1
    t = 0:dt:20; % Set time, use 0:dt:1800 for Experiment 1 and 24,30,52 for Exeriment 2
  
    options = odeset('RelTol',1e-6,'AbsTol',1e-12);
    
    sol = pdepe(m,@viruspdes,@pdeic,@pdebc,x,t,options,P); % PDE solver
    
    
    % --------------------------------------------------------------
    % Set parameter values
    % Experiment 1, rates are in /sec and micrometers
    % Experiment 2, rates are in /day and milimeters
    
function [c,f,s] = viruspdes(x,t,u,DuDx,P)
    
    
    % P = [2.4614 0.2916] % fit of ind P(2) for virus diffusion in collagen alone (Mckee_data.mat)
    % P = [0.0116 0.7355];%[0.6951  0.0049]; % fit of [du r] to Mckee_collagenase data from Mckee_treatment_dat.mat
    %aa = [0.2169 0.3101 3.6937e-08 0.0896];
    
    du = 0.0360;    % diffusion of tumour cells mm^2/day 
    dc = 0;         % collagen diffusion is no longer consider in the model, however, at times we need a small diffusion here to avoid computational problems
    dv = 0.8640;    %10; % max virus diffusion rate 10 micro^2/s or 0.8640 mm^2/day
    a = 150;        % reproduction of virus / cell 
    r = 0.2169;     % tumour growth 
    rc = 0.1;       % collagen growth rate         
    g = 0.0896;     % death of virus /day 0.04 
    h = 1;          % death of infected cells 
    b = 3.6937e-08.*(1-u(4));   % internalization (b = 4/3499 where systme changes stability)  
    k = 10^6;       % Tumour cell carrying capacity 
    kappa = 1;      % max virus diffusion rate
    ind = 0.2916;   % n squeeze prob exponent 
    p = kappa*(1-u(4)^(ind)); % squeeze probability
    m =  1;         % hill-type exponent 
    Kc = k;         % hill-type coefficient 
    comp = 0.3101;  %0.2440; 
  
    c = [1; 1; 1; 1];
    
    % Use extra doses of virus for Experiment 2
    
   % if t>=2 && t<=3 
   %    else
            Sup=0;
   %    end

    % Diffusion and supply function
    
    f = [du 0 0 0; 0 dv*p 0 dv*kappa*ind*u(4)^(ind-1)*u(2); 0 0 du 0; 0 0 0 dc]* DuDx;  %Diffusion when particle struggles to jump in.
    %f = [du 0 0 0 ; 0 dv*p 0 -kappa*u(2)*dv ; 0 0 du 0 ; 0 0 0 dc]* DuDx;              % Diffusion when a particle struggles to jump out
    %f = [du; p*dv; du; dc].*DuDx;                                                      % Fickian diffusion
    
    s = [(r*u(1)*(1-((u(1)+u(3))/k+comp.*u(4)))-b*u(1)*u(2)) ; (-g*u(2)-b*(u(1)+u(3))*u(2) + a*h*u(3))+Sup; b*u(1)*u(2)-h*u(3); rc*u(1)^m./(Kc^m+u(1)^m).*u(4).*(1-(u(1)/k+comp.*u(4)))]; %(1-u(1)/k-u(4))(1-((u(1)+u(3))/k+comp.*u(4)))
    % --------------------------------------------------------------
function u0 = pdeic(x,P) % inital conditions
  
    
    % Initial conditions
    
    % Experiment 1: 
    %Important, change coordinates to Cartisian, set m=0 above
    % conditions for Mckee figure 1D
    
    %a = -0.5574; b = 0.3198; c = 0.7780; coll = a./(1+exp(-b*x))+c; 
    %uu = 0; vv = 1*abs(x)<10; II =0;
    
    
    
    % Experiment 2: Fitting model to data. Run Mckee_nonlin_fitter. 
    % We also use Total_pop_cal. Make sure to include the P values in this
    % code as scaling factors for collagen. Collagen needs to be high and low.
    % Important, coordiantes in spherical geometry, set m=2 above
    % Initial conditions for virus treament with collagenase
    
    %col = (0.15*sin(5*x)+0.75);
    %uu = 10^6*(1-coll);%.*(heaviside(3-x)) 
    %vv = 2.4*10^5*heaviside(1-x);
    %coll = col*P(3);  We multiply this by P(3) when fitting to data.
    %II = 0; 
     
    % Intital condition for global and local sensitivity analysis is the same as
    % above with only one exception. 
    % For local sensitivity we use 0.83 and 1/3 for high and low collagen respectively. 
    % For r_u in high collagen, we allow for very small diffusion of collagen dc=0.001 to avoid numerical problems.
    % When simulating the exponential collagen in global sensitivity
    % analysis use
    
    % coll = exp(-0.4568.*x); % change scaling factor in front inline the methods section 
    % uu = 10^6*(1-coll).*(heaviside(3-x));
    % col = 0.85.*coll;
    % vv = 2.4*10^5*heaviside(1-x);
    %II = 0; 
     
     % simulate desmoplastic reaction for ueno simulations
    
     %coll = -1*(square(1*x+2,5));   % immature collagen
     %coll = -1*(square(3*x,12)); % Intermediate collagen
     %coll = -1*(square(10*x)); % Mature collagen
     %col = (coll+1.1)/2.5; % rescale the collagen for the different cases
     %II = 0;
     %uu = 10^6*(1-col);
     %vv = 2.4*10^5*heaviside(1-x);
     %II = 0; 
     
     % simulate desmoplastic reaction for varying collagen density matureand immature
     % for these conditions, P is inputed from Total_pop_cal
     
     %col1 = -1*(square(10*x)); 
     %coll = col/P(1)+P(2);
     %col1 = -1*(square(1*x+2,5));
     %col = col1/P(1)+P(2);
     %uu = 10^6*(1-col);
     %vv = 2.4*10^5*heaviside(1-x);
     %II = 0; 
 
     
   u0 = [ uu ; vv ; II ; col]; % Initial conditions 
   
% --------------------------------------------------------------
function [pl,ql,pr,qr] = pdebc(xl,ul,xr,ur,~,t) % No flux boundary conditions
    pl = [0; 0; 0; 0];
    ql = [1; 1; 1; 1];
    pr = [0; 0; 0; 0];
    qr = [1; 1; 1; 1];
    
   