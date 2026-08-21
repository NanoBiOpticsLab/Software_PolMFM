function [SPOTS_FIT, SPOTS_DETECTED, FIT_RESULT] = predetection_fit_3planes(file_image_reco,path_image_reco,microscope_settings,PSF_theo,frames)
% predetection_3planes Performs pre-detection of 3 planes MFM single molecules
% Input: file_image & path_image =  location of the first frame to analyze.
% The file should come from a 4P-MFM acquisition and should be already
% reconstructed into the 6 channels the following order (z,p): (-1;0°),(-1;90°),(0;45°),(0;135°),(+1;0°),(+1;90°).
% microscope_settings =  parameters of the microscope 
% PSF_theo =  parameters of the theoretical PSF of the microscope
% Output:
% spots_detected =  predetection results. 
%   col 1 = frame
%   col 2-4 = X-Y-Z (in pixel, round)
%   col 5-10 = xmnin-xmax-ymin-ymax-zmin-zmax for sub_spots mosaic boundaries
% spots_fit = fit gaussian of the pre-detected spots 
%   col 1 = frame
%   col 2-4 = Y-X-Z localization in nm / 
%   col 5-6-7 = amp - bgd - resnorm 
%   col 8-9-10 = sigma x-y-z
%   col 11-12-13 = centroid Y-X-Z
%   col 14-15-16 = localization in nm (not shifted, in the reference of the subs_pots mosaic)
%   col 17 = nb of iterations     
% FIT_Result = structure containing the result of the fit for each spot. 

% Louise Régnier
% Adapted from MFM_piepline
% 20/06/2023

%% create the 3 planes MFM image from the 6 polarized channels

% warning('For now processing is done on 1 frame only #fr')
% Initialize results 
SPOTS_FIT = [];
SPOTS_DETECTED = [];
FIT_RESULT = [];

%% Find the reconstructed files
[~,image_reco_name_1,ext] = fileparts(file_image_reco);
fileList =  dir(fullfile(path_image_reco, '*.tif'));
l=size(fileList,1); 
nb_digits=numel(num2str(l));

%% First step: choose the threshold; test on several frames; Define the parameters for predetection
button = 'Yes';
% options.detect_th  = prctile(image_struct.data_filtered(:),99);
options.detect_th  = 250;

%default values for filtering
defaultValue_2{1} = num2str(4);
defaultValue_2{2} = num2str(4);
defaultValue_2{3} = num2str(.5);
defaultValue_2{4} = num2str(.5);

while strcmp(button,'Yes')
    
    frame_to_test = inputdlg(['Choose frame to test for detection threshold  (between ' num2str(frames(1)) ' and ' num2str(frames(2)) ')'],'Input',[1 35],{'1'});
    fr_test = str2double(frame_to_test{1});
    options.pixel_size = microscope_settings.pixel_size;
    
    %load image
    file_image_fr = fullfile(path_image_reco,[image_reco_name_1(1:end-nb_digits) num2str(fr_test,['%.' num2str(nb_digits) 'd']) ext]);
    tiff_file= file_image_fr;
    h = imfinfo(tiff_file).Height;
    w = imfinfo(tiff_file).Width;
    stack_size=size(imfinfo(tiff_file),1);
    imageMFM=zeros(h,w,stack_size/2); %
    for ii=1:stack_size/2
        imageMFM(:,:,ii) =  0.5*double(imread(tiff_file,(ii-1)*2+1))+0.5*double(imread(tiff_file,(ii-1)*2+2));
    end
    imageMFM_MIP = max(imageMFM,[],3);
    
    image_struct.data = imageMFM;
    
    figure(486786),
    % subplot(141), 
    subplot(9,3,[1 4 7]), imagesc(imadjust(mat2gray(imageMFM(:,:,1)),[0.05 0.7])),axis equal,axis image, title('z = -dz'),
    % subplot(142),  
    subplot(9,3,[10 13 16]), imagesc(imadjust(mat2gray(imageMFM(:,:,2)),[0.05 0.7])),axis equal,axis image, title('z = 0'), 
    % subplot(143),  
    subplot(9,3,[19 22 25]), imagesc(imadjust(mat2gray(imageMFM(:,:,3)),[0.05 0.7])),axis equal,axis image, title('z = +dz'),
    % subplot(144),  
    subplot(9,3,[2 3 5 6 8 9 11 12]), imagesc(imadjust(mat2gray(imageMFM_MIP),[0.05 0.7])),   axis equal, axis image, title('Max Int Proj')
    colormap gray,
    sgtitle(['MFM Image - frame #' num2str(fr_test)])
      
    %filtering param

    dlgTitle = 'Filtering: specify kernel size';
    prompt_2(1) = {'XY Kernel size for Gaussian smoothing to remove Background'};
    prompt_2(2) = {'Z  Kernel size for Gaussian smoothing to remove Background'};
    prompt_2(3) = {'XY Kernel size of Gaussian smoothing for spots enhancement'};
    prompt_2(4) = {'Z  Kernel size of Gaussian smoothing for spots enhancement'};
    
    userValue = inputdlg(prompt_2,dlgTitle,1,defaultValue_2);
    
    if( ~ isempty(userValue))
        kernel_size.bgd_xy = str2double(userValue{1});
        kernel_size.bgd_z  = str2double(userValue{2});
        kernel_size.psf_xy = str2double(userValue{3});
        kernel_size.psf_z   = str2double(userValue{4});
        defaultValue_2{1} = num2str(kernel_size.bgd_xy);
        defaultValue_2{2} = num2str(kernel_size.bgd_z);
        defaultValue_2{3} = num2str(kernel_size.psf_xy);
        defaultValue_2{4} = num2str(kernel_size.psf_z);
    end
    
    flag.output     = 0;
    
    img_filt = img_filter_Gauss_v3(image_struct,kernel_size,flag);
    
    image_struct.data_filtered = uint32(img_filt);
    %filter = kernel;
        
    % %- Show filtered image (maximum projection)
    % figure(4538900),
    % h_plot = imshow(imadjust(mat2gray(max(img_filt,[],3))),[]);
    % title(['Maximum projection of filtered image (Gaussian)- frame #' num2str(fr_test)], 'FontSize',8);
    % colormap(hot)
    
    %predetction param
    regions.xy = 5; %crop
    regions.z  = 0; %crop
    % regions.xy_sep = 5; %for nonmaxsupr; not implemented yet
    % regions.z_sep  = 0;
    options.size_detect = regions;
    
    % define cell_prop
    cell_prop(1).x      = [1 1 w w];
    cell_prop(1).y      = [1 h h 1];
    cell_prop(1).reg_type = 'Rectangle';
    cell_prop(1).reg_pos  = [1 1 w h];
    cell_prop(1).label = 'EntireImage';
    
    
    % define flag_struct
    flag_struct.mode_predetect = 'connectcomp' ; %'nonMaxSupr'
    flag_struct.reg_pos_sep = 0;
    flag_struct.region_smaller = 1;
    flag_struct.output = 1;
    
    options.cell_prop = cell_prop;
    options.box_spots_mosaic.xy = 5;
    options.box_spots_mosaic.z = 2;
    button_threshold = 'No';
    while strcmp(button_threshold,'No')
        [spots_detected, img_mask, CC_GOOD] = spots_predetect_v17(image_struct,options,flag_struct);
        button_threshold = questdlg('Do you validate threshold?');
        if strcmp(button_threshold,'No')
            thr = inputdlg('Detection threshold ?','Input',1,{num2str(options.detect_th)});
            options.detect_th = str2double(cell2mat(thr));
        end
    end
    button = questdlg('Test filtering/thresholding on another frame ? If no: analysis will run over all frames.');
end

%
msgbox({['KERNEL SIZE:'],['bgd (xy) = ' num2str(kernel_size.bgd_xy)],...
    ['bgd (z) = ' num2str(kernel_size.bgd_z)],...
    ['snr (xy) = ' num2str(kernel_size.psf_xy)],...
    ['snr (z) = ' num2str(kernel_size.psf_z)],...
    ['THRESHOLD: ' num2str(options.detect_th)],...
    },'Final parameters','modal')
    
nnf = 0; file_fitparam = [path_image_reco '\fit_parameters_' num2str(nnf) '.txt'];
while isfile(file_fitparam)
    nnf = nnf+1; file_fitparam = [path_image_reco '\fit_parameters_' num2str(nnf) '.txt'];
end
 fid = fopen(file_fitparam, 'w');
 fprintf(fid,[ 'Analyzed frames: ' num2str(frames(1)) ' to ' num2str(frames(2)) ...
     '\n KERNEL SIZE: \n bgd (xy) = ' num2str(kernel_size.bgd_xy) ...
    '\n bgd (z) = ' num2str(kernel_size.bgd_z) ...
    '\n snr (xy) = ' num2str(kernel_size.psf_xy) ...
    '\n snr (z) = ' num2str(kernel_size.psf_z) ...
    '\n THRESHOLD = ' num2str(options.detect_th)...
    '\n \n MICROSCOPE PARAMETERS: \n x-y pixel_size = ' num2str(microscope_settings.pixel_size.xy) ' nm' ...
    '\n z pixel_size = ' num2str(microscope_settings.pixel_size.z) ' nm'...
    '\n excitation filter wavelength = ' num2str(microscope_settings.Ex) ' nm' ...
    '\n emission filter wavelength = ' num2str(microscope_settings.Em) ' nm' ...
    '\n numerical aperture = ' num2str(microscope_settings.NA)   ...
    '\n refractive index = ' num2str(microscope_settings.RI) ...
    '\n camera sensitivity = ' num2str(microscope_settings.Cam_sensitivity) ...
    '\n camera EM gain = ' num2str(microscope_settings.Cam_EM) ...
    '\n Convert photons (1=yes /0=no)' num2str(microscope_settings.Photons_convert)...
    ]);
 fclose(fid);
    


%% Now start loop over all frames
    h_waitbar = waitbar(0,'Initializing ...');
tic
flag_struct.output = 0;
for fr = frames(1):frames(2)
    file_image_fr = fullfile(path_image_reco,[image_reco_name_1(1:end-nb_digits) num2str(fr,['%.' num2str(nb_digits) 'd']) ext]);
    tiff_file=file_image_fr;
    h = imfinfo(tiff_file).Height;
    w = imfinfo(tiff_file).Width;
    stack_size=size(imfinfo(tiff_file),1);
    imageMFM=zeros(h,w,stack_size/2); %
    
    for ii=1:stack_size/2
        imageMFM(:,:,ii) =  0.5*double(imread(tiff_file,(ii-1)*2+1)+imread(tiff_file,(ii-1)*2+2));
    end    
    image_struct.data = imageMFM;
    
    %% Filter the image
    flag.output     = 0;   
    img_filt = img_filter_Gauss_v3(image_struct,kernel_size,flag);
    image_struct.data_filtered = uint32(img_filt);
        
    %% Pre-detect       
    [spots_detected, img_mask, CC_GOOD] = spots_predetect_v17(image_struct,options,flag_struct);
    [sub_spots, sub_spots_filt] = spots_predetect_mosaic_v1(image_struct,img_mask,spots_detected,flag_struct);
    %% Fitting
    
    %=== Fit with 3D Gaussian
    
    % disp('Fitting: STARTED ... ')
    
    %- Some parameters
    
    pixel_size            = microscope_settings.pixel_size;
    flag_struct.parallel  = 0;
    % fit_limits            = handles.fit_limits;
    %
    % warning('set fitlimts')
    % fit_limits.sigma_xy_min = 20;%
    % fit_limits.sigma_z_min = 20;
    % fit_limits.sigma_xy_max = 250;
    % fit_limits.sigma_z_max = 500;
    
    fit_limits.sigma_xy_min = PSF_theo.xy_nm*0.2;%
    fit_limits.sigma_z_min = PSF_theo.z_nm*0.2;
    fit_limits.sigma_xy_max = PSF_theo.xy_nm*2.5;
    fit_limits.sigma_z_max = PSF_theo.z_nm*2;
    
    %== Determine mode of fitting
    mode_fit        = 'sigma_free_xz';
    par_start       = [];
    handles.par_fit = [];
    
    % sigmaxy - sigmaz - centerx - centery - centerz  - amp - bgd
    
    bound.lb = [fit_limits.sigma_xy_min fit_limits.sigma_z_min -inf -inf -inf 0   0  ];
    bound.ub = [fit_limits.sigma_xy_max fit_limits.sigma_z_max inf  inf  inf  inf inf];
    

    %- Call fitting routine
    parameters.pixel_size  = pixel_size;
    parameters.PSF_theo    = PSF_theo;
    parameters.par_start   = par_start;
    parameters.flag_struct = flag_struct;
    parameters.mode_fit    = mode_fit;
    parameters.bound       = bound;
    parameters.box_spots_mosaic.xy = options.box_spots_mosaic.xy;
    parameters.box_spots_mosaic.z  = options.box_spots_mosaic.z;
    parameters.Photons_convert = microscope_settings.Photons_convert;
    parameters.Cam_EM = microscope_settings.Cam_EM;
    parameters.Cam_sensitivity = microscope_settings.Cam_sensitivity; 

    [spots_fit, FIT_Result] = spots_fit_batch_3D_Gauss_v7(spots_detected, sub_spots, parameters);
   
    SPOTS_FIT = [SPOTS_FIT; [fr*ones(size(spots_fit,1),1) spots_fit]];
    SPOTS_DETECTED = [SPOTS_DETECTED; [fr*ones(size(spots_fit,1),1) spots_detected]];
    FIT_RESULT = [FIT_RESULT FIT_Result];

    perc=(fr-frames(1))/frames(2)*100;
    waitbar(perc/100,h_waitbar,sprintf(' processing %3.2f%% ',perc))
    
end

waitbar(100/100,h_waitbar,sprintf(' DONE '))
pause(1);
close(h_waitbar);

toc


