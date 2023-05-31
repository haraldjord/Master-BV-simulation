function [h_bot] = search_presett_volume(flag, param, min_water_density)
%SEARCH_PRESETT_VOLUME Summary of this function goes here
%   Search for presetting of outer lid, that result in buoyancy of the
%   first 10 mm piston displacement 
%   input: param - struct with dimentions parameters
%          min_water_density - mimum water density at the specific location
%   output: h_bot - presetting of height between vehicle bottom and outer
%           lid, that result in floating the first 10 mm piston displacament. 

h_bot_temp = param.h_bot;

searching = true;
while searching 

    h_piston_temp = find_piston_equilibrium(flag ,param, min_water_density, h_bot_temp);
    if h_piston_temp < 10 %h_piston_equilibrium less then 10mm 
        h_bot_temp = h_bot_temp + 0.0001; % Add one 0.1mm    
    else
        searching = false;
        
    end
end

h_bot = h_bot_temp;
end

