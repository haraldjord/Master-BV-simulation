function [V_max, V_piston] = calc_preset_volume(param, h_bot)
%% Calculate maximum volume of vehicle, when piston at bottom position. In SI Unit [meter] and [kg] 
V_led = 0.014*(pi*0.0275^2) + 0.013*(pi*0.019^2) + 6*((0.005*pi*0.0035^2 - 0.003*pi*0.002^2));    % volume of LED "house" = V1 + V2 + 6*V_screwHead

% Volume in pressure sensor hole
V_sensor = 0.007*(pi*0.003^2);

%V_pistonBlock is the volume beneath the piston
V_pistonBlock = 0.0115*(pi*0.015^2);

% Volume of outer lid
r_outerLid = param.r_outerLid;
h_outerLid = param.h_outerLid;

A_outerLid = pi*r_outerLid^2;
V_outerLid = A_outerLid * h_outerLid + 0.005*pi*(0.078^2-0.065^2);% volume outer lid + ring with bigger diameter


%Volume of vehicle house, where y is preset value 
if nargin == 1 % get h_bot from parameter 
    h_bot = param.h_bot; % adjustable length from bottom vehicle to outer lid range = [10, 500]mm
end
r_house = param.r_house; % radius of vehicle house

A_vehicle = pi*(r_house^2);
V_vehicle = A_vehicle*h_bot;

%effective volume = V_max = V_vehicle + V_outerLid + V_led - V_pistonBlock
V_max = V_vehicle + V_outerLid + V_led - V_pistonBlock - V_sensor;

%% calculate V_piston
h_pistonMax = param.h_pistonMax;
h_pistonBlock = param.h_pistonBlock;
r_piston = param.r_piston; % TODO measure this diameter!

A_piston = pi*r_piston^2;
V_piston = (h_pistonMax - h_pistonBlock)*A_piston;


end

