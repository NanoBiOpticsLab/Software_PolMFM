% Copyright (C) <2026> <Louise REGNIER and Bassam HAJJ>

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://gnu.org>.
% Use code with caution.
%%-------------------------------------------
clear   
close all 
clc
%% Initialization : Specify file name and path 

root_folder = uigetdir(pwd, 'Select the root folder named "Software_PolMFM"');
if isequal(root_folder,0)
    error('No folder selected. Operation cancelled by user.');
end
% root_folder = '/Users/lyr25/Library/CloudStorage/OneDrive-INSTITUTCURIE/Paper-PolMFM/Soumission Nature Comm/Codes/'; 
addpath(genpath(root_folder));% this adds all the required functions to the Matlab path

% *************** ACQUISITION FILES: raw data .tif file name and folder ***************

path_data = [root_folder '/0 - DATA_Tutorial/DATA_SilicaBeads_SLB_NR_2024-01-05_002_substack/'];
file_data = '002_frames3000to3300_pixelsize=120nm_EM200_30ms.tif';

% *************************************************************************************

% *************** CALIBRATION FILES: .tif file names and folders for the 3 types of calibrations ***************

% i) z-stack with beads calibration
file_calib = 'image_Pos0.ome.tif';
path_calib = [root_folder '/0 - DATA_Tutorial/Calibrations/Fluobeads-200nm_laser-560nm_100nm-step/003/'];

% ii) intensity calibration with lamp (unpolarized)
file_calib_int =    'image_Pos0.ome.tif';
path_calib_int =    [root_folder '/0 - DATA_Tutorial/Calibrations/Calib_Intensity_NotPolar_TL3D/'];

% iii) polarization calibration: path and file name of the first step of the polarizer
path_calib_polar = [root_folder '/0 - DATA_Tutorial/Calibrations/Calib-Polar_5deg-steps_filter617nm/Calib_Polar_2024-01-11_0/images/RAW_DATA/'];
file_calib_polar = 'image_Pos0.ome.tif';

% *************************************************************************************

disp('------DONE selecting images to analyze------');
%% A - Calibration channels : build the calibration matrices 
% Specify the bead and intensity calibration file paths and names
N_channels = 6;

% create the distorsions matrices (FOV calibration) ;
calib_distorsion_channels_superlocver_v2(file_calib,path_calib,N_channels);
% intenstity correction
calib_intensity_correction(file_calib_int,path_calib_int,N_channels);

display ('----DONE Generating/loading calibration files------')
%% B - Reconstruction  
% [file_data,path_data] = uigetfile('*.tif;*.tiff','Folder','E:\Louise\DATA\'); 
% create the new .tiff : save 1 image / channel (total 6 stacks)

[border,file_image_reco,path_image_reco] = reconstruct_channel_images_superlocver_v2(path_data,file_data,fullfile(path_calib,file_calib),fullfile(path_calib_int,file_calib_int),N_channels);

fid = fopen([path_data 'calib_files_4_' file_data '.txt'], 'w');
fprintf(fid, ['Calibration files for stack ' strrep([path_data file_data],'\','\\' ) ' : ']);
fprintf(fid,['\n channels (beads) calibration file = ' strrep([path_calib file_calib],'\','\\' )]);
fprintf(fid,[ '\n Intensity calibration file = ' strrep([path_calib_int file_calib_int],'\','\\' )]);
fclose(fid);
        

% reconstruct MFM images for detection
disp('----DONE reconstruction/selection of channel images------')

%% C - Pre-détection + D - Localization

[microscope_settings, PSF_theo, Detect]=experimental_settings();

% Choose which frames to analyse:
nb_frames_tot=size(dir(fullfile(path_image_reco, '*.tif')),1);
frames = [1 nb_frames_tot]; % frames(1) is the first frame, frame(end) id the last frame to be analysed

if isfile([path_data file_data(1:end-4)  '_results_fit_fr' num2str(frames(1)) 'to' num2str(frames(2)) '.mat'])
    load([path_data file_data(1:end-4)  '_results_fit_fr' num2str(frames(1)) 'to' num2str(frames(2)) '.mat'],'spots_fit','FIT_Result')
    disp('---- Localization file loaded! ----')
else
    [spots_fit_0, spots_detected, FIT_Result] = predetection_fit_3planes(file_image_reco,path_image_reco,microscope_settings,PSF_theo,frames);
    disp('---- 3D-GF Localization DONE ----')
end
% save fit results
file_save = [file_data(1:end-4) '_results_fit_fr' num2str(frames(1)) 'to' num2str(frames(2)) '.mat'];
uiwait(msgbox(['Saving fitted results:' newline newline ...
               '- File: ' file_save newline ...
               '- Location: ' path_data newline ...
               '- Variable saved: spots_fit'], ...
              'Saving', 'modal'));
% convert to table to include headers
spots_fit = array2table(spots_fit_0(:,1:10),'VariableNames',{'frame' 'x [nm]' 'y [nm]' 'z [nm]' 'intensity [u?]' 'bkd[u?]' 'resnorm[u?]' 'sigmax [nm]' 'sigmay [nm]' 'sigmaz [nm]'});
save([path_data file_save],'spots_fit',"-v7.3")

disp('---- Localization file saved! ----')
%% D bis - optional: display localizations
filter_threshold = spots_fit.("sigmax [nm]")>80 &spots_fit.("sigmax [nm]")<200 ...
& spots_fit.("z [nm]")<1000 & spots_fit.("z [nm]")>-500 ...   
& spots_fit.("sigmaz [nm]")<600 & spots_fit.("sigmaz [nm]")>160 ...
;%
figure(), 
scatter3(spots_fit.("x [nm]")(filter_threshold,1),spots_fit.("y [nm]")(filter_threshold,1),spots_fit.("z [nm]")(filter_threshold,1),10,spots_fit.("z [nm]")(filter_threshold,1),'filled')
axis equal, axis image, colormap turbo
xlabel('x [nm]'),ylabel('y [nm]'),zlabel('z [nm]'),
h = colorbar;
ylabel(h,'Z [nm]')
title('Fit result: 3D scatter plot of localizations')
%% E - Calibration polarization 

% create propagation matrix K
Kmatrix = calib_orientation_Kmatrix(file_calib_polar,path_calib_polar);
[up_one, ~, ~] = fileparts(path_calib_polar(1:end-1)); [up_two, ~, ~] = fileparts(up_one); [up_three, ~, ~] = fileparts(up_two);
save([up_three filesep 'propagation_matrix.mat'], 'Kmatrix');
uiwait(msgbox(['Saving propagation matrix:' newline newline ...
               '- File: propagation_matrix.mat' newline ...
               '- Location: ' up_three newline ...
               '- Variable saved: Kmatrix'], ...
              'Saving', 'modal'));

% if correction for dipole emission
dipolar_emission_correction = 0;
if  dipolar_emission_correction
    Kmatrix = K_dipole_generation(Kmatrix,microscope_settings);
end

display ('----DONE generated polarization calibration files------')

%% F - Orientation

% Define box parameters for intensity calculation
    % Loop over all frames
        % Loop over all spots
            % Convert spot fit coordinates from nm to pixels
            % Convert coordinates for each channel using transformation matrices
            % Retrieve the images


% uncomment the following parameters if you want custom ones
box_parameters = struct();
% box_parameters.add_px = 10;
% box_parameters.factor_box = 3;
% box_parameters.box_size = 10;
% box_parameters.min_box_size = 5;

results_all = []; %zeros(size(spots_fit,1),16); 
if  istable(spots_fit)
    spots_fit = table2array(spots_fit);
end
tic

                h_waitbar = waitbar(0, 'Initializing...');
                for fr = frames(1):frames(2)
                    in_frame = find(spots_fit(:,1)==fr);                 
                    if isempty(in_frame)==0
                        % convert localization fit result to each channel with transformation matrix
                        spots_channels = find_spots_channels(file_data,path_data,spots_fit(in_frame,:),file_calib,path_calib,microscope_settings,fr,border);
                        % Compute PSFs intensities in each channel
                        I = find_intensity_spots_v2(file_data,path_data,file_image_reco,path_image_reco,spots_fit(in_frame,:),spots_channels,microscope_settings,fr,box_parameters);
                        % compute orientations
                        results_orientation = orientation_retrieval_matrix(Kmatrix,I);
                        results_all = [results_all ; [spots_fit(in_frame,1:10) ((results_orientation))]];
                    end
                    perc = (fr - frames(1)) / (frames(2) - frames(1));
                    waitbar(perc, h_waitbar, sprintf('Processing frame %d / %d  (%3.1f%%)', fr, frames(2), perc*100));

                end
                waitbar(1, h_waitbar, 'Done!');
                pause(1);
                close(h_waitbar);
                      

toc

dz = microscope_settings.pixel_size.z; % nm
z_threshold = 100; % [nm]
dz_low = dz-z_threshold; dz_high = dz+z_threshold;
centered_z  = (results_all(:,4)>=dz_low & results_all(:,4)<=dz_high); 
z_low = results_all(:,4)<dz_low;
z_high = results_all(:,4)>dz_high; 
results_all = [results_all zeros(size(results_all,1),2)]; % add two columns: col#17 = delta according to z pos; col#18 = rho according to z pos;
results_all(z_low,[17 18]) = results_all(z_low,[12 15]);
results_all(z_high,[17 18]) = results_all(z_high,[13 16]);
results_all(centered_z,[17 18]) = results_all(centered_z,[11 14]);

%offset orientations to match horizontal polarisation to horizontal axis of the camera- microscope dependant
offset = 60; % orientation offset (can be tuned) [degrees]
results_all(:,[14 15 16 18]) = 0.5*wrapTo360((results_all(:,[14 15 16 18])+offset)*2); 
 
% save the results 
if microscope_settings.Photons_convert == 0 
    Results_titles = {'frame' 'x [nm]' 'y [nm]' 'z [nm]' 'intensity [a.u.]' 'bkd[a.u.]' 'resnorm[a.u.]' 'sigmax [nm]' 'sigmay [nm]' 'sigmaz [nm]' 'delta'  'deltaz1' 'deltaz3' 'rho' 'rhoz1' 'rhoz3' 'deltaz' 'rhoz'};
elseif microscope_settings.Photons_convert == 1 
    Results_titles = {'frame' 'x [nm]' 'y [nm]' 'z [nm]' 'intensity [ph]' 'bkd[ph]' 'resnorm[ph]' 'sigmax [nm]' 'sigmay [nm]' 'sigmaz [nm]' 'delta'  'deltaz1' 'deltaz3' 'rho' 'rhoz1' 'rhoz3' 'deltaz' 'rhoz'};
end

%remove NaNs --> now: do not remove NaNs
idx_Nan =isnan(results_all ); 
Results_noNaN = results_all; %results_all(sum(idx_Nan,2)==0,:);
file_name_results= [path_data [file_data(1:end-4) '_locs-and-orientation-results.csv']];
file_name_results_0=file_name_results; ij=0;
while exist(file_name_results)==2
    ij=ij+1;
    file_name_results = [file_name_results_0(1:end-4) '_' num2str(ij) '.csv'];
end
results_table = array2table(Results_noNaN,'VariableNames',Results_titles);
writetable(results_table,file_name_results);

file_name_results= [path_data [file_data(1:end-4) '_locs-and-orientation-results.txt']];
file_name_results_0=file_name_results; ij=0;
while exist(file_name_results)==2
    ij=ij+1;
    file_name_results = [file_name_results_0(1:end-4) '_' num2str(ij) '.txt'];
end
fid=fopen(file_name_results,'w');
[r,s]=size(Results_noNaN);
fprintf(fid,['%g' repmat('\t %g',1,s-1) '\n'],Results_noNaN');
fclose(fid);

disp('Orientation computation: DONE')

uiwait(msgbox(['Orientation computation complete:' newline newline ...
               '- .csv and .txt files saved' newline ...
               '- Location: ' file_name_results(1:end-4)], ...
              'Saving - Done', 'modal'));
%% F-bis (Optional) Scatter plot of the locs and orientations (color-coded)

% filtering bad localisations
filter_locs = results_table.("sigmax [nm]") < 200 &  results_table.("sigmax [nm]") > 80 & ...
     results_table.("sigmaz [nm]")  < 600 &  results_table.("sigmaz [nm]") > 160 & ...
     results_table.("z [nm]")  < 1000 &  results_table.("z [nm]") > -500 & ...
    results_table.delta>20 & results_table.delta < 150 ;

figure(), 
scatter3(results_all(:,2),results_all(:,3),results_all(:,4),5,'filled','MarkerFaceColor',[.8 .8 .8],'MarkerFaceAlpha',0.1)
hold on;
scatter3(results_all(filter_locs,2),results_all(filter_locs,3),results_all(filter_locs,4),10,(results_all(filter_locs,14)),'filled')
axis equal, axis image, 
xlabel('x [nm]'),ylabel('y [nm]'),zlabel('z [nm]'),
colormap hsv  
colorbar, clim([0 180]), 
title({'Results 3D - Color = \rho [°] ','gray = filtered locs'}) 
