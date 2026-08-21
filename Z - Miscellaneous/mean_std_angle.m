function [M, S] = mean_std_angle(in,dim,sens)
% mean_std_angle will calculate the mean and std of a set of angles (in degrees) based
% on polar considerations.
%
% Usage: [M, S] = mean_std_angle(in,dim,sens)
%
% in is a vector or matrix of angles (in degrees)
% out is the mean of these angles along the dimension dim
%
% If dim is not specified, the first non-singleton dimension is used.
%
% A sensitivity factor is used to determine oppositeness, and is how close
% the mean of the complex representations of the angles can be to zero
% before being called zero.  For nearly all cases, this parameter is fine
% at its default (1e-12), but it can be readjusted as a third parameter if
% necessary:
%
% [out] = meanangle(in,dim,sensitivity)
%
% Written by J.A. Dunne, 10-20-05
%
if nargin<3
    sens = 1e-12;
end
if nargin<2
    ind = min(find(size(in)>1));
    if isempty(ind)
        %This is a scalar
        M = in;S=0;
        return
    end
    dim = ind;
end
in_rad = in * pi/180;
z_in = exp(1i*in_rad);
mid = mean(z_in,dim,'omitnan');
M = wrapTo360(atan2(imag(mid),real(mid))*180/pi);
M(abs(mid)<sens) = nan;
in_rotbyM = wrapTo180(in-wrapTo360(M));
S = std(in_rotbyM,'omitnan');
