function spots_channel = find_spots_channels(file_data,path_data,spots_fit,file_calib,path_calib,microscope_settings,fr,border) 
calibrationFile=[path_calib,[file_calib(1:end-4) '.mat']];
load(calibrationFile,'transformPerSlice','cornerPerSlice');
channelshiftFile=[path_calib,[file_calib(1:end-4) '_channel-shift.mat']];
load(channelshiftFile,'dist_p_xy');

%%
N_spots = size(spots_fit,1);
nb_planes = 6;
clear spots_channel
spots_channel{1} = zeros(N_spots,2);
planes_transfo = [1 4 2 5 3 6];
for k = 1:N_spots
    
    coord = zeros(nb_planes,2);
    for kk=1:nb_planes
        coord_0 = [spots_fit(k,3) spots_fit(k,2)]/microscope_settings.pixel_size.xy  + border -1 + [dist_p_xy(planes_transfo(kk),2) dist_p_xy(planes_transfo(kk),1)];
        coord(kk,:) = coord_0;
%               transformPerSlice{planes_transfo(kk)}.b * coord_0 * transformPerSlice{planes_transfo(kk)}.T + ...
%              -mean(transformPerSlice{planes_transfo(planes_transfo(kk))}.c) ;%...
%             0*[dist_p_xy(planes_transfo(kk),2) dist_p_xy(planes_transfo(kk),1)] +...
%             0*[cornerPerSlice(planes_transfo(kk),2) cornerPerSlice(planes_transfo(kk),1)]+...
%             0*[1 1]; 
        spots_channel{kk}(k,:) = coord(kk,:);
    end
    
end

flag_imageoutput = 0;
if flag_imageoutput==1
    img = (double(imread(fullfile(path_data,file_data),fr)));
    figure, imagesc(img-gaussSmooth(img,fr,'same')), axis equal, axis image
    M = 'rgwmyc';
    hold on, for kk= 1:nb_planes
        plot(spots_channel{kk}(:,1),spots_channel{kk}(:,2),[M(kk) 'o'],'MarkerSize',7)
    end
    imcontrast
    colormap gray
end
% 
% cornerPerSlice(kk,1):cornerPerSlice(kk,1)+size(stack,1)-1,cornerPerSlice(kk,2):cornerPerSlice(kk,2)+size(stack,2)-1
%         stack=zeros(cornerPerSlice(1,3)-cornerPerSlice(1,1)-3,cornerPerSlice(1,4)-cornerPerSlice(1,2)-3,nb_planes,'uint16');

    
%%
% 
% figure, plot(spots_fit(:,2)/microscope_settings.pixel_size.xy + 1, spots_fit(:,1)/microscope_settings.pixel_size.xy + 1,'o')
% axis equal , axis image

% %%
%  im_ch1 = (double(imread([ path 'image_Pos0.ome\image_Pos0.ome_corr_0020.tif'],...
%      1)));
%  
%  im_ch2 = (double(imread([ path 'image_Pos0.ome\image_Pos0.ome_corr_0020.tif'],...
%      2)));
%  
%  im_ch3 = (double(imread([ path 'image_Pos0.ome\image_Pos0.ome_corr_0020.tif'],...
%      3)));
%  
%  im_ch4 = (double(imread([ path 'image_Pos0.ome\image_Pos0.ome_corr_0020.tif'],...
%      4)));
%  
%  im_ch5 = (double(imread([ path 'image_Pos0.ome\image_Pos0.ome_corr_0020.tif'],...
%      5)));
%  
%  im_ch6 = (double(imread([ path 'image_Pos0.ome\image_Pos0.ome_corr_0020.tif'],...
%      6)));
% figure, subplot(231),imagesc(im_ch1-gaussSmooth(im_ch1,15,'same')), axis equal, axis image
% hold on, plot(spots_fit(:,2)/microscope_settings.pixel_size.xy + 1, spots_fit(:,1)/microscope_settings.pixel_size.xy + 1,'rx','MarkerSize',10)
% subplot(232),imagesc(im_ch2-gaussSmooth(im_ch2,15,'same')), axis equal, axis image
% hold on, plot(spots_fit(:,2)/microscope_settings.pixel_size.xy + 1, spots_fit(:,1)/microscope_settings.pixel_size.xy + 1,'rx','MarkerSize',10)
% subplot(233),imagesc(im_ch3-gaussSmooth(im_ch3,15,'same')), axis equal, axis image
% hold on, plot(spots_fit(:,2)/microscope_settings.pixel_size.xy + 1, spots_fit(:,1)/microscope_settings.pixel_size.xy + 1,'rx','MarkerSize',10)
% subplot(234),imagesc(im_ch4-gaussSmooth(im_ch4,15,'same')), axis equal, axis image
% hold on, plot(spots_fit(:,2)/microscope_settings.pixel_size.xy + 1, spots_fit(:,1)/microscope_settings.pixel_size.xy + 1,'rx','MarkerSize',10)
% subplot(235),imagesc(im_ch5-gaussSmooth(im_ch5,15,'same')), axis equal, axis image
% hold on, plot(spots_fit(:,2)/microscope_settings.pixel_size.xy + 1, spots_fit(:,1)/microscope_settings.pixel_size.xy + 1,'rx','MarkerSize',10)
% subplot(236),imagesc(im_ch6-gaussSmooth(im_ch6,15,'same')), axis equal, axis image
% hold on, plot(spots_fit(:,2)/microscope_settings.pixel_size.xy + 1, spots_fit(:,1)/microscope_settings.pixel_size.xy + 1,'rx','MarkerSize',10)

%%
