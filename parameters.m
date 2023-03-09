%% parameters
% realated to dimentions and mass of BV

% outer lid
param.r_outerLid = 0.073; % [m]
param.h_outerLid = 0.295; % [m]

% vehicle house
param.r_house = 0.065; %[m]
% param.y = 0.102;
param.y = 0.097;    %0.080;%% Initial guess, adjustable legth between wehicle house bottom, and outer lid. range = [0.010, 0.135]m (fresh water 0.090 is found to be a god value)
param.Areal_vehicle_bottom = pi*param.r_house^2;
param.Areal_vehicle_top = pi*param.r_outerLid^2;
% piston
param.h_pistonMax = 0.062; % [m]
param.h_pistonBlock = 0.011; % [m]
param.h_piston = 0.051;% (max - block)
param.r_piston = 0.02; % [m]

param.Areal_piston = pi*param.r_piston^2; %piston areal

% total mass
param.mass = 6.225;  %NEW MEASURMENT 3 + 1.54 + 1.75 + 0.07; %[kg] measured mass (6.36kg)
% param.mass = 6.300; %[kg]measured 02.03.2023 at kybernetisk værksted 
% stepper motor


%other stuff...?
