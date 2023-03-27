%% constants
clear all; close all; clc;
%%


%%%%%%%%%%%
searchPresetVolume = false; % search for preset volume, false use param.h_bot 
useCTDprofile = true;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       ; % boolean variable, (true/false) idicate wheter ctd profile is available.  
%%%%%%%%%%
CTD = load("CTD-probe/example_boorsa.mat"); % load CTD probe measurement
depth_CTD = CTD.Depth;
densityProfile = CTD.Density;

if useCTDprofile
    index = interp1(depth_CTD, 1:length(depth_CTD), 50, 'nearest'); % find the nearest index for 50 meter depth 
    depth_CTD = depth_CTD(1:index);
    densityProfile = densityProfile(1:index);

    rho_water_interval = max(densityProfile) - min(densityProfile)
end
fresh_water_density = 999;    % constant water density (useCTDprofile = false)

% load dimention parameters in struct param.
parameters


g = 9.81; % [m/s^2]
Cd_top = 1.28;%0.9 % Drag coefficient when floating to surface
Cd_bottom  = 1; %Drag coefficient when sinking 
[V_max, V_piston] = calc_preset_volume(param);

rho_min_vehicle = param.mass/V_max;
if searchPresetVolume
    param.h_bot = 0.08; %Initial value of outer lid height from vehicle bottom  
end

while searchPresetVolume %% NEEDS REFACTORY
    if useCTDprofile
       
       [V_max, V_piston] = calc_preset_volume(param);
       rho_min_vehicle = param.mass/V_max;
        if rho_min_vehicle < min(densityProfile)
            param.h_bot = param.h_bot +0.001; % Add one mm (for safety limit)
            searchPresetVolume = false;
            disp("recomended presetting of volume: ");
            disp(param.h_bot);
        end
        param.h_bot = param.h_bot +0.001; % Add one mm
    else
        param.h_bot = param.h_bot +0.001; % Add one mm
        [V_max, V_piston] = calc_preset_volume(param);
        rho_min_vehicle = param.mass/V_max;
        if rho_min_vehicle < fresh_water_density
            param.h_bot = param.h_bot +0.001; % Add one mm (for safety limit)            
            searchPresetVolume = false;
            disp("recomended presetting of volume: ");
            disp(param.h_bot);
        end
        %param.h_bot = param.h_bot +0.001; % Add one mm
    end
        
end

[V_max, V_piston] = calc_preset_volume(param);
rho_min_vehicle = param.mass/V_max;
rho_max_vehicle = param.mass/(V_max-V_piston);
rho_vehicle_interval = rho_max_vehicle - rho_min_vehicle
find_pistion_equilibrium(param, true, 20, rho_max_vehicle, depth_CTD);

delta_pos_positon_min = (1/51200000); % Linear movement per step


%%%% piston parameters
v_piston_max = 2.44e-3; %%  1e-3; % max linear speed of piston [m/s]

% if (rho_vehicle_interval < rho_water_interval) || useCTDprofile
%     disp("Warnign: densiti of water has greater range then vehicle");
% end

%% plot density profile and vehicle density range.
figure(1)
hold on 
if useCTDprofile
    plot(densityProfile, depth_CTD, 'b');
else
   plot(fresh_water_density*[1,1], [1,50]) 
end
plot(rho_min_vehicle*[1,1], [0, 50], 'r');
plot(rho_max_vehicle*[1,1], [0, 50], 'y');
hold off
set(gca,'YDir','Reverse');
grid();
title('Density Profile');
xlabel("Density [kg/m^3]");
ylabel("Depth[m]");
legend("water density profile", "min density vehicle", "max density cehicle");

%% PID constants

% tryout
% Kp = 0.023;
% Ki = 0.005;
% Kd = 0.03;

% Optimum fresh water tank parameters.
Kp = 0.026;
Ki = 0.001;
Kd = 0.05;

alpha = 0.4; % Tuning parameter for EMA filter [0,1]
sampleTime = 0.5;
% sysd = tf([alpha, 0], [1, -(1-alpha)], sampleTime); %Discrete transfer function (EMA filter) 
% sysc = d2c(sysd); % Transfer function for Pressure sensor (EMA filter)
%bode(sysc) % Plot bode diagram of EMA filter
integralTreshold = 2; %threshold in meter when integral is activated.
offsetSensor = 0;


%% run simulation
step_depth1 = 10;
step_depth2 = 1;
stepTime = 180;
%Simulation 
tspan = [0 180*2]; % Time span for simulation
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

fig6 = figure(6) % PID term contribution
    hold on 
    plot(out.PIDterms.time, out.PIDterms.signals.values(:,1)*1000);
    plot(out.PIDterms.time, out.PIDterms.signals.values(:,2)*1000);
    plot(out.PIDterms.time, out.PIDterms.signals.values(:,3)*1000);
    hold off
    grid()
    xlabel('Time[s]'); ylabel("Piston Position [mm]")
    legend('Kp-term', 'Ki-term', 'Kd-term');
%% save simulation results for comparing with measured response
save('../BuoyancyVehiclePlotTestData/simOut.mat','out')


