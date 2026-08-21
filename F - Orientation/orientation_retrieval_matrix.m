function results_orientation = orientation_retrieval_matrix(Kmatrix,I);

% results_orientation = orientation_retrieval(Kmatrix,I);
% Computes the orientation parameters from intensities I in the 6 channels of
% the 4P-MFM setup; 
% INPUTS
% Kmatrix is the polarization propagation matrix of the
% microscope.
% I = cell array of N (N=nb of spots) cells; each cell is a 1x6 array with
% the intensity in each of the 6 channel 
% OUTPUT
% results_orientation contains the orientation parameters for the
% spots given in I, in the following order: [rho, rho (-dz), rho (+dz), delta, delta (-dz), delta (+dz)]
% Computations are performed either combining all the planes of the MFM or
% only with information coming from -dz or +dz.


% reorder Kmatrix to have it in the order
% XX0   YY0   XY0    --> plane -dz
% XX90  YY90  XY90   --> plane -dz
% XX45  YY45  XY45   --> plane 0
% XX135 YY135 XY135  --> plane 0
% XX0   YY0  XY0     --> plane +dz
% XX90  YY90 XY90    --> plane +dz

K_ordered= Kmatrix;

K_ordered_1 = Kmatrix(1:4,:);    %for -dz plane

K_ordered_3 = [Kmatrix(5,:);    %for +dz plane
    Kmatrix(6,:);
    Kmatrix(3,:);
    Kmatrix(4,:)];

%% Compute orientations 

%initialize results 
results_orientation = zeros(size(I,2),6); 

for ii=1:size(I,2)
    
    % find 2nd dipole moment (M vector)
    I_1 = I{ii}(1:4)'; %for -dz plane
    I_3 = [I{ii}(5:6)'; I{ii}(3:4)']; %for +dz plane
    M_1 = pinv(K_ordered_1)*I_1;
    M_3 = pinv(K_ordered_3)*I_3;
    M = pinv(K_ordered)*I{ii}';
    % Compute P factors
    P_uv = 2*M(3)/(M(1)+M(2));
    P_xy   = (M(1)-M(2))/(M(1)+M(2));
    P_xy_1 = (M_1(1)-M_1(2))/(M_1(1)+M_1(2));
    P_xy_3 = (M_3(1)-M_3(2))/(M_3(1)+M_3(2));
    P   = sqrt(P_uv^2+P_xy^2);
    P_1 = sqrt(P_uv^2+P_xy_1^2);
    P_3 = sqrt(P_uv^2+P_xy_3^2);
    
    lambda   = (1-P)/(3-P);
    lambda_1 = (1-P_1)/(3-P_1);
    lambda_3 = (1-P_3)/(3-P_3);
    
    if lambda>=0 && lambda<=1/3
        delta   = 2*acosd((-1+sqrt(9-24*lambda))/2); % ajouter abs dans le acosd ????????
    else
        delta = NaN;
    end
    if lambda_1>=0 && lambda_1<=1/3
        delta_1 = 2*acosd((-1+sqrt(9-24*lambda_1))/2);
    else
        delta_1 = NaN;
    end
    if lambda_3>=0 && lambda_3<=1/3
        delta_3 = 2*acosd((-1+sqrt(9-24*lambda_3))/2);
    else 
        delta_3 = NaN;
    end
    rho = 0.5*wrapTo360(atan2d(P_uv,P_xy));
    rho_1 = 0.5*wrapTo360(atan2d(P_uv,P_xy_1));
    rho_3 = 0.5*wrapTo360(atan2d(P_uv,P_xy_3));
    
    results_orientation(ii,:) = ([delta delta_1 delta_3 rho rho_1 rho_3]); 

end


