function [M, S] = mean_std_angle_0to180(in)

% mean_std_angle_0to180 will calculate the mean and std of a set of angles
% (in degrees) between 0 and 180° on polar considerations 

[M360, S360] = mean_std_angle(2*in);
M = 0.5*M360;
S = 0.5*S360; 
