%% parameters
% realated to dimentions and mass of BV

% outer lid
param.r_outerLid = 0.073; % [m]
param.h_outerLid = 0.295; % [m]

% vehicle house
param.r_house = 0.065; %[m]
% param.h_bot = 0.102;

% h_bot is adjustable legth between wehicle house bottom, and outer lid. range = [0.010, 0.135]m (fresh water 0.0972 is found to be a god value)
param.h_bot = 0.0972; % good value for fresh water tank (rho = 999)    
% param.h_bot seems to have an deviation of about 1.5mm (sim 96.5 == real 95mm)
%param.h_bot = 0.089; %experimental for sea water
%param.h_bot = 0.0924; % experimental with new safety rig 

param.Areal_vehicle_bottom = pi*param.r_house^2;
param.Areal_vehicle_top = pi*param.r_outerLid^2;
% piston
param.h_pistonMax = 0.062; % [m]
param.h_pistonBlock = 0.011; % [m]
param.h_piston = 0.051;% (max - block)
param.r_piston_block = 0.015;
param.r_piston = 0.02; % [m]

param.Areal_piston = pi*param.r_piston^2; %piston areal

% total mass
% param.mass = 6.225;  %NEW MEASURMENT 3 + 1.54 + 1.75 + 0.07; %[kg] measured mass (6.36kg)
% param.mass = 6.300; %[kg]measured 02.03.2023 at kybernetisk værksted 
% stepper motor
param.mass = 0.075 + 4.4470 + 1.746; % [6.268 kg] measured 13.03.20023 at kybernetisk verksted led + house + outer lid 

if flag.useSafetyRig
    param.safetyRig_mass =  (0.196); % Rig and screw/nut
    param.mass = param.mass + param.safetyRig_mass*flag.n_safetyRig;
    
    param.safetyRig_volume = (110987e-9 + (0.030*pi*0.008^2 + 0.008*pi*0.011^2*2) ) * flag.n_safetyRig; %volume of rig (found with solidwork) + volume of bolt and nut;
end

if flag.fishTag
    param.tagHolder_mass = 0.022; % measured mass of fish tag holder + tag
    param.mass = param.mass + param.tagHolder_mass;
    param.tagHolder_volume = 36273.00e-9; % volume of fish tag holder (including tag) found with solidwork
end

%param.mass = 6.492; 
%other stuff...?
