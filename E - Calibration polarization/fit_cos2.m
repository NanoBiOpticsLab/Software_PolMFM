function [s, fit] = fit_cos2(I,calib_step_degree,max_angle)
% fit the function A*cos(2x-phi)+B to find intensity and absolute 0°
% input angle parameters are in degree (°)
% polarization 
    x=(0:calib_step_degree:max_angle);
    yu = max(I);
    yl = min(I);
    yr = (yu-yl);                               % Range of ‘y’
    yz = I-yu+(yr/2);
    zx = x(yz .* circshift(yz,[0 1]) <= 0);     % Find zero-crossings
    per = 2*mean(diff(zx));                     % Estimate period
    ym = mean(I);                               % Estimate offset
    fit = @(b,x)  b(1).*(cosd(2.*(x - b(2)))) + b(3);    % Function to fit
    fcn = @(b) sum((fit(b,x) - I).^2);                              % Least-Squares cost function
    s = fminsearch(fcn, [yr;  -1;  ym]);

% 
% xp = linspace(min(x),max(x));
% figure(153)
% hold on
% plot((x),I,'--.',  (xp),fit(s,xp), '-')
% grid

