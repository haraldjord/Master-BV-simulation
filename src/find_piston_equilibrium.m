function [h_pis_equilibrium] = find_piston_equilibrium(flag, param, waterDensity, h_bot)
%FIND_PISTION_EQUILIBRIUM Summary of this function goes here
%   Finding piston position that yeald equilibrium at given depth and
%   density profile
%   Return piston height in mm.
    
    if nargin ==4
        V_max = calc_preset_volume(flag,param,h_bot);
    elseif nargin == 3
        V_max = calc_preset_volume(flag, param);
    else 
        disp("Warning! function needs 3 or 4 inputs");
    end    
    rho = waterDensity;
    m = param.mass;
    A_pis = param.Areal_piston;

    h_pis_equilibrium = ((V_max - (m/rho))/A_pis);
    h_pis_equilibrium = h_pis_equilibrium *1000; % meter to mm
end

