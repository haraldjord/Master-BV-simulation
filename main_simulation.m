%% constants
clear all; close all; clc;
%%
addpath(fullfile(pwd,"src"));
fid = fopen('BV-info.txt', 'wt');
fprintf(fid, 'BV related information from last simulation\n');



%%%%%%%%%%% Flags
flag.searchPresetVolume = true;  % search for preset volume, false use param.h_bot 
flag.useCTDprofile = true;       % boolean variable, (true/false) idicate wheter ctd profile is available, false use constant density of fresh water. 
flag.useSafetyRig = false;         % boolean variable, indicate wheter holder for safety line beeing used.
flag.n_safetyRig = 0;             % number of holder for safety line (1 or 2)
flag.fishTag = true;               % boolean variable, indicate wheter fish tag holder is monted.
fprintf(fid, 'Flags: \nsearchPresetVolume: %d \nuseCTDprofile: %d\n', flag.searchPresetVolume, flag.useCTDprofile);
fprintf(fid, 'useSafetyRig: %d \t number of safety holders %d\nuse fish tag holder: %d\n\n', flag.useSafetyRig, flag.n_safetyRig, flag.fishTag);
%%%%%%%%%%
CTD = load("CTD-probe/borsa14.04.2023.mat"); % load CTD probe measurement (if not available , load "CTD-probe/example_boorsa.mat" to avoid conflict in simulink)
%CTD = load("CTD-probe/example_borsa.mat");
depth_CTD = CTD.Depth;
maxDepth  = 10;


if flag.useCTDprofile % Load and reduce densityProfile array to 50 meters (target depth)
    
    index = interp1(depth_CTD, 1:length(depth_CTD), maxDepth, 'nearest'); % find the nearest index for 50 meter depth 
    depth_CTD = depth_CTD(1:index);
    densityProfile = CTD.Density(1:index);

    rho_water_interval = max(densityProfile) - min(densityProfile);
else % use constant water density
    densityProfile = 999;    % constant water density (flag.useCTDprofile = false)
    fprintf(fid, 'constant water density beeing used: %.1f [kg/m^3]\n', densityProfile);
end
% load dimention parameters in struct param.
parameters
g = 9.81; % [m/s^2]
Cd_top = 1.28;%0.9 % Drag coefficient when floating to surface
Cd_bottom  = 1.28; %Drag coefficient when sinking


% [V_max, V_piston] = calc_preset_volume(flag, param);
% rho_min_vehicle = param.mass/V_max;
if flag.searchPresetVolume
    param.h_bot = 0.050; %Initial guess of outer lid height from vehicle bottom
    if flag.useCTDprofile
        param.h_bot = search_presett_volume(flag, param, min(densityProfile));
    else
        param.h_bot = search_presett_volume(flag, param, densityProfile);
    end
end
fprintf(fid, 'Vehicle presetting of outer lid height from bottom: h_bot = %.1f [mm]\n', param.h_bot*1000);

[V_max, V_piston] = calc_preset_volume(flag, param);
rho_min_vehicle = param.mass/V_max;
rho_max_vehicle = param.mass/(V_max-V_piston);
delta_rho_vehicle = rho_max_vehicle - rho_min_vehicle;
delta_rho_water = max(densityProfile) - min(densityProfile);
h_pis_eq = find_piston_equilibrium(flag, param, min(densityProfile));
fprintf(fid, 'Piston position that result in equilibrium water surface %.2f [mm] \n\t piston_pos < %.2f --> floating \n\t piston_pos > %.2f --> sinking\n', h_pis_eq, h_pis_eq, h_pis_eq);
fprintf(fid, 'h_bot seems to have an deviation of about 1.5mm (sim %.2f == real %.2f) BUT THIS NEDS VERIFICATION, AND MIGHT BE CHANGING WITH NEW BATTERY PACKAGE\n', param.h_bot*1000, (param.h_bot*1000-1.5));
fprintf(fid, 'rho min vehicle: %.2f \t rho max vehicle %.2f \t delta rho vehicle: %.2f \n', rho_min_vehicle, rho_max_vehicle, delta_rho_vehicle);
fprintf(fid, 'rho min water:   %.2f \t rho max water:  %.2f \t delta rho water: %.2f \n',min(densityProfile), max(densityProfile), delta_rho_water);
if (min(densityProfile) ~= densityProfile(1))
   fprintf(fid, 'Warning! minimum water density not at surface. Make sure vehicle are able to float to surface.\n');
end
if rho_max_vehicle < max(densityProfile)
     if flag.useCTDprofile
        index = interp1(densityProfile, 1:length(densityProfile), rho_max_vehicle, 'nearest'); % find the nearest index for rho_max_density
        max_diving_depth = depth_CTD(index);
        fprintf(fid, 'Vehicle are unable to dive deeper then: %.1f meter\n',max_diving_depth);

     end
end
delta_pos_positon_min = (1/51200000); % Linear movement per step


%%%% piston parameters
v_piston_max = 2.44e-3; %%  1e-3; % max linear speed of piston [m/s]
fprintf(fid, '\nmaximum linear speed of piston: %.3f [mm/s]\n', v_piston_max*1000);


%% plot density profile and vehicle density range.
figure(1)
hold on 
if flag.useCTDprofile
    plot(densityProfile, depth_CTD, 'b');
else
   plot(densityProfile*[1,1], [0,maxDepth],'b'); 
end
plot(rho_min_vehicle*[1,1], [0, maxDepth], 'r');
plot(rho_max_vehicle*[1,1], [0, maxDepth], 'y');
hold off
set(gca,'YDir','Reverse');
grid();
title('Density Profile');
xlabel("Density [kg/m^3]");
ylabel("Depth[m]");
legend("water density profile", "min density vehicle", "max density vehicle");

%% PID constants

% optimum PID constant 
Kp = 0.026;
Ki = 0.001;
Kd = 0.1;

alpha = 0.4; % Tuning parameter for EMA filter [0,1]
sampleTime = 0.5;
% sysd = tf([alpha, 0], [1, -(1-alpha)], sampleTime); %Discrete transfer function (EMA filter) 
% sysc = d2c(sysd); % Transfer function for Pressure sensor (EMA filter)
%bode(sysc) % Plot bode diagram of EMA filter
integralTreshold = 10; %threshold in meter when integral is activated.
offsetSensor = 0;


% simulation depth and time
step_depth1 = 5;
step_depth2 = 0;
stepTime = 120;

fprintf(fid, '\nSimulation Parameters:\n');
fprintf(fid, 'Kp = %.4f \nKi = %.4f \nKd = %.4f \nIntegral treshold = %.2f\n', Kp, Ki, Kd, integralTreshold);
fprintf(fid, 'alpha (EMA filter) = %.2f \ntime constant = %.2f [seconds]\n', alpha, sampleTime); 

%% Run Simulation 
tspan = [0 60*2]; % Time span for simulation
%options = simset('MaxStep', 0.5,'MinStep',1e-11, 'AbsTol', 1e-11, 'RelTol', 1e-11);
%set_param('buoyancy2','AlgebraicLoopSolver','LineSearch');
out = sim('buoyancy_simulation.slx',tspan);     %,options);

%% Plot Simulation result 
Depth = out.stepResponse.signals.values(:,2);
step = out.stepResponse.signals.values(:,1);
timeStep = out.stepResponse.time;

pisPos = out.PistonPosition.signals.values;
PIDout = out.PIDout.signals.values;
timePiston = out.PistonPosition.time;


fig5 = figure(5);
subplot(2,1,1);
    % hold on;
    %SEGMENT 2
    hold on;
    plot(timeStep, Depth);
    plot(timeStep, step,':');
    legend('Measured Depth','Target Depth');
    title('Measured Depth');
    xlabel('Time [s]');
    xlim([0 max(timeStep)]);
    ylim([-0.1 (max(Depth)+0.2)]);
    ylabel('Depth [m]');
    set(gca,'YDir','reverse');
    grid();
    hold off;

subplot(2,1,2);
    hold on;
    plot(timePiston, pisPos*1000);
    plot(timePiston, PIDout*1000);
    xlabel('Time [s]');
    xlim([0 max(timeStep)]);
    ylabel('Position [mm]');
    title('Piston position');
    legend("Piston position","PID output");
    grid();
    hold off;

fig6 = figure(6); % PID term contribution
    hold on 
    plot(out.PIDterms.time, out.PIDterms.signals.values(:,1)*1000);
    plot(out.PIDterms.time, out.PIDterms.signals.values(:,2)*1000);
    plot(out.PIDterms.time, out.PIDterms.signals.values(:,3)*1000);
    hold off
    grid()
    xlabel('Time[s]'); ylabel("Piston Position [mm]")
    legend('Kp-term', 'Ki-term', 'Kd-term');
    
    
% Close BV-info.txt
fclose(fid);
%% save simulation results for comparing with measured response
save('../BuoyancyVehiclePlotTestData/simOut.mat','out')


