function [h_pis_equilibrium] = find_piston_equilibrium(param, waterDensity, h_bot)
%FIND_PISTION_EQUILIBRIUM Summary of this function goes here
%   Finding piston position that yeald equilibrium at given depth and
%   density profile
%   Return piston height in mm.
    
    if nargin ==3
        V_max = calc_preset_volume(param,h_bot);
    else
        V_max = calc_preset_volume(param);
    end    
    rho = waterDensity;
    m = param.mass;
    A_pis = param.Areal_piston;

    h_pis_equilibrium = ((V_max - (m/rho))/A_pis);
    h_pis_equilibrium = h_pis_equilibrium *1000; % meter to mmm
end

