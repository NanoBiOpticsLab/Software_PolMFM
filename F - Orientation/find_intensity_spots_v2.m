function I = find_intensity_spots_v2(file_data,path_data,file_image_reco,path_image_reco,spots_fit,spots_channels,microscope_settings,fr,box_parameters); 

% Compute the intensity factors of each spot in each channel of the frame
% Inputs: 
% file_data,path_data: image (original to analyze)
% file_image_reco,path_image_reco: corresponding reconstructed frame 
% spots_fit: localizations of all the molecules in the frame - localization is given in the reconstructed framework
% spots_channels: corresponding localizations in the original reference (6 locs per spot, for each channel)
% microscope_settings
% fr: frame number to be reconstructed
% Output: I is a structure containing for each spot of the frame the
% intensities in each channel
% Louise Régnier; 04/12/2023

% v2:   - define a minimal size box
%       - reduce the nb of fits 
%       - 2 possibilities for size box computation: 'adaptative' or 'fixed'
%       - the parameters are the following: 
%           *add_px =  window size around the box before intensity
%           computation (center +/- add_px): all molecules closer to the
%           edge of the image are discared. Defines a max box size for
%           adaptative method.
%           *factor_box = only for adaptive box: box = center +/-round(factor_box*PSF_size_px)
%           *min_box_size = only for adaptive box: box cannot be smaller
%           than center +/- min_box_size
%           *box_size = only for fixed box: size box = center +/- box_size
%       - default parameters: add_px = 8; factor_box = 3; box_size = 5; min_box_size = 5;
%% define parameters for box
box_parameters.answer_method_box_intensity = 'Fixed';
answer_method_box_intensity = box_parameters.answer_method_box_intensity; % 'fixed' or 'adaptative'
if isfield(box_parameters,'add_px')
    add_px = box_parameters.add_px;
else
    add_px = 8;
end

if isfield(box_parameters,'factor_box')  % for 'adaptative' box size definition
    factor_box = box_parameters.factor_box;
else 
    factor_box = 3;
end

if isfield(box_parameters,'box_size') % for    'fixed    box size definition
    box_size = box_parameters.box_size;
else 
    box_size = 5;   
end

if isfield(box_parameters,'min_box_size')  % define minimum box parameter
    min_box_size = box_parameters.min_box_size;
else 
    min_box_size = 5;
end

%%

nb_planes = 6;
nb_spots = size(spots_channels{1},1);
I{1} = zeros(1,nb_planes); % intensity factors for each channel for each bead

spots_fit_px = spots_fit(:,2:3)/microscope_settings.pixel_size.xy  + 1; % positions of the spots in reconstructed framework
%%
img = (double(imread(fullfile(path_data,file_data),fr))); % original image frame

%% Find the reconstructed files
[~,image_reco_name_1,ext] = fileparts(file_image_reco);
fileList =  dir(fullfile(path_image_reco, '*.tif'));
l=size(fileList,1); 
nb_digits=numel(num2str(l));


%load image
    file_image_fr = fullfile(path_image_reco,[image_reco_name_1(1:end-nb_digits) num2str(fr,['%.' num2str(nb_digits) 'd']) ext]);
    tiff_file= file_image_fr;
    h = imfinfo(tiff_file).Height;
    w = imfinfo(tiff_file).Width;
    stack_size=size(imfinfo(tiff_file),1);
    I_3d=zeros(h,w,stack_size/2); % build MFM image (3 planes, no polar)
    for ii=1:stack_size/2
        I_3d(:,:,ii) =  0.5*double(imread(tiff_file,(ii-1)*2+1)+imread(tiff_file,(ii-1)*2+2));
    end

 
%% loop on all spots 

parameters.pixel_size = microscope_settings.pixel_size; %nm
parameters.flags.crop =0;
parameters.flags.output = 0;
parameters.par_crop = [];
parameters.par_microscope= microscope_settings; %nm

% figure, imagesc(img), axis equal , colormap gray, axis image
switch answer_method_box_intensity
    case 'Adaptative'
        for s = 1:nb_spots
            if round(spots_fit_px(s,1))-add_px>0 && round(spots_fit_px(s,1))+add_px < size(I_3d,1)+1 && round(spots_fit_px(s,2))-add_px>0 && round(spots_fit_px(s,2))+add_px < size(I_3d,2)+1
                spot =  I_3d(round(spots_fit_px(s,1))-add_px:round(spots_fit_px(s,1))+add_px,round(spots_fit_px(s,2))-add_px:round(spots_fit_px(s,2))+add_px,:);
                for ch = 1:nb_planes

                    % subimage of the bead #s in plane #ch

                    %        im1_bead = im1(round(pos_bead2(2,n_PSF))-add_px:round(pos_bead2(2,n_PSF))+add_px,round(pos_bead2(1,n_PSF))-add_px:round(pos_bead2(1,n_PSF))+add_px);
                    % figure, imagesc(spot)

                    % Compute the 2D GF to obtain the extension of the PSF in each plane (sigma PSF) and save the value
                    if mod(ch,2)==1 % do the fit for each MFM plane (ie every 2 channel)
                        img_PSF2d.data = spot(:,:,round(ch/2));
                        [PSF_fit2d, img_PSF2d] = PSF_2D_Gauss_fit_v8(img_PSF2d,parameters);
                        PSF_size_px = sqrt(2*log(2))*abs(PSF_fit2d.sigma_xy/parameters.pixel_size.xy); %FWHM %
                    end
                    box_x = round(spots_channels{ch}(s,2))-round(factor_box*PSF_size_px):round(spots_channels{ch}(s,2))+round(factor_box*PSF_size_px); %box_x;%
                    box_y = round(spots_channels{ch}(s,1))-round(factor_box*PSF_size_px):round(spots_channels{ch}(s,1))+round(factor_box*PSF_size_px); %box_y;%

                    % if box is too small, minimum box size definition
                    if size(box_y,2)<min_box_size
                        box_y = round(spots_channels{ch}(s,1))-floor(min_box_size/2):round(spots_channels{ch}(s,1))+floor(min_box_size/2);
                    end
                    if size(box_x,2)<min_box_size
                        box_x = round(spots_channels{ch}(s,2))-floor(min_box_size/2):round(spots_channels{ch}(s,2))+floor(min_box_size/2);
                    end
                    % if box too big, maximum box size definition
                    if size(box_y,2)>2*add_px+1
                        box_y = round(spots_channels{ch}(s,1))-add_px:round(spots_channels{ch}(s,1))+add_px;
                    end
                    if size(box_x,2)>2*add_px+1
                        box_x = round(spots_channels{ch}(s,2))-add_px:round(spots_channels{ch}(s,2))+add_px;
                    end
                    if box_x(1)>0 && box_y(1)>0 && box_x(end)<1+size(img,1) && box_y(end)<1+size(img,2)
                        % PSF_ch = img(box_x,box_y)...
                        %     -mean([mean(img(box_x,box_y(1))) mean(img(box_x,box_y(end))) mean(img(box_x(1),box_y(2:end-1)))  mean(img(box_x(end),box_y(2:end-1)))]);
                        PSF_ch = img(box_x,box_y)...
                            -mean([(img(box_x,box_y(1)))' (img(box_x,box_y(end)))' (img(box_x(1),box_y(2:end-1)))  (img(box_x(end),box_y(2:end-1)))]);

                        I{s}(1,ch) = sum(PSF_ch(:)); % compute intensity
                        % hold on, rectangle('Position',[min(box_y)-0.5 min(box_x)-0.5 max(box_y)-min(box_y)+1 max(box_x)-min(box_x)+1],'EdgeColor','g') %need -0.5 to fit the pixel entirely (rather than the corrdinate of the pixel)
                    else
                        % warning('out of window')
                        I{s}(1,ch) = NaN;
                    end
                end

            else
                I{s}=  zeros(1,nb_planes);
            end

        end

    case 'Fixed'
        for s = 1:nb_spots
            if round(spots_fit_px(s,1))-add_px>0 && round(spots_fit_px(s,1))+add_px < size(I_3d,1)+1 && round(spots_fit_px(s,2))-add_px>0 && round(spots_fit_px(s,2))+add_px < size(I_3d,2)+1
                for ch = 1:nb_planes
                    % subimage of the bead #s in plane #ch
                    %        im1_bead = im1(round(pos_bead2(2,n_PSF))-add_px:round(pos_bead2(2,n_PSF))+add_px,round(pos_bead2(1,n_PSF))-add_px:round(pos_bead2(1,n_PSF))+add_px);

                    %fixed box size (2*box_size+1)
                    box_x = round(spots_channels{ch}(s,2)) - box_size :round(spots_channels{ch}(s,2)) + box_size ;
                    box_y = round(spots_channels{ch}(s,1)) - box_size :round(spots_channels{ch}(s,1)) + box_size ;

                    if box_x(1)>0 && box_y(1)>0 && box_x(end)<1+size(img,1) && box_y(end)<1+size(img,2)
                        % PSF_ch = img(box_x,box_y)...
                        %     -mean([mean(img(box_x,box_y(1))) mean(img(box_x,box_y(end))) mean(img(box_x(1),box_y(2:end-1)))  mean(img(box_x(end),box_y(2:end-1)))]);
                        PSF_ch = img(box_x,box_y)...
                            -mean([(img(box_x,box_y(1)))' (img(box_x,box_y(end)))' (img(box_x(1),box_y(2:end-1)))  (img(box_x(end),box_y(2:end-1)))]);

                        I{s}(1,ch) = sum(PSF_ch(:)); % compute intensity
                        % hold on, rectangle('Position',[min(box_y)-0.5 min(box_x)-0.5 max(box_y)-min(box_y)+1 max(box_x)-min(box_x)+1],'EdgeColor','r') %need -0.5 to fit the pixel entirely (rather than the corrdinate of the pixel)
                    else
                        % warning('out of window')
                        I{s}(1,ch) = NaN;
                    end
                end

            else
                I{s}=  zeros(1,nb_planes);
            end

        end

end