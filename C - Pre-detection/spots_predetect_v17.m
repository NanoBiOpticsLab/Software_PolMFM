function  [spots_detected, img_mask, CC_GOOD] = spots_predetect_v17(image_struct,options,flag_struct)
% Function returns pixel location of local maximum
%
% spots_detected(:,1:3)   ... xyz position of detected local maximum in image
% spots_detected(:,10:11) ... Intensity in raw image and filtered image
%
% v1 Feb 28, 2010
% - Original implementation
%
% v2 March 1, 2010
% - Restrict search on area of nucleus alone
%
% v3 March 2, 2010
% - Adjust to new definiton of outlines of cell and TS
%
% v4
%- Adjust to new defintion of cell_prop
%
% v5
% - new version which allow thresholding of the score parameter


%- Empty structures in case you don't use connected component
CC_best = {};
CC_GOOD = {};
prop_img_detect = [];

%= Options
size_detect = options.size_detect;
detect_th   = options.detect_th;
cell_prop   = options.cell_prop;
pixel_size  = options.pixel_size;
box_spots_mosaic_xy = options.box_spots_mosaic.xy;
box_spots_mosaic_z = options.box_spots_mosaic.z;


%= Dimension of the image
%image      = image_struct.data;
image_filt   = image_struct.data_filtered;
image        = image_struct.data; 
[dim.Y dim.X dim.Z] = size(image_filt);

%= Set image outside of nucleus to zero
image_filt_mask = image_filt;   
    
    %=== Restrict analysis to a certain sub-region
    
    %- No field defined
    if not(isfield(cell_prop,'mask_cell_3D'))
        
        [X_grid,Y_grid] = meshgrid(1:dim.X,1:dim.Y);
        mask_cell_2D = inpolygon(X_grid,Y_grid,cell_prop.x,cell_prop.y);
        mask_cell_3D = repmat(mask_cell_2D,[1,1,dim.Z]);  


    end
    

    %- Set image outside of sub-region to zero

        %- Entire cell
            image_filt_mask(not(mask_cell_3D)) = 0;
  
    
        
    [dim.Y dim.X dim.Z] = size(image_filt_mask);



%===  Detect local maximum by nonmaximal suppression
% Comment: toolbox might have to be installed again if error message pop up.
% Give 3d coordinates of all identified non-suppressed point locations (n x d)
% Coordinates are integer - no sub-pixel information
% (1,1,1) is pixel in upper left corner on first focal plane


switch flag_struct.mode_predetect

    case 'nonMaxSupr'
%         disp('... non maximum supression ...');
%            
%         if flag_struct.reg_pos_sep == 0
%         
%             rad_detect = round([size_detect.xy size_detect.xy size_detect.z]);
%         else
%             rad_detect = round([size_detect.xy_sep size_detect.xy_sep size_detect.z_sep]);
%             
%         end
%            
%         pos_dum = nonMaxSupr(double(image_filt_mask), rad_detect,detect_th);
%         
%         
%         if size(pos_dum,1) > 1
%         
%             %- Remove spots that are within half the detection radius
%             %  This can occur when small pixels are used
%             pos_sort = sortrows(pos_dum);
% 
%             pos_diff = abs(diff(pos_sort));
% 
%             ind_x_0 = pos_diff(:,2) <= ceil(size_detect.xy);
%             ind_y_0 = pos_diff(:,1) <= ceil(size_detect.xy);
%             ind_z_0 = pos_diff(:,3) <= ceil(size_detect.z);
% 
%             pos_diff(ind_x_0,2) = 0;
%             pos_diff(ind_y_0,1) = 0;
%             pos_diff(ind_z_0,3) = 0;
% 
%             ind_remove = ismember(pos_diff,[0 0 0],'rows');
% 
%             pos_pre_detect = pos_sort;
%             pos_pre_detect(ind_remove,:) = [];
% 
%         else
%             pos_pre_detect = pos_dum;
%         end
% 
%         clearvars pos_dum;

    case 'connectcomp'

        % disp('... connected components ...');
        
        %- Connected components
        par_ccc.conn        = 26;   % Connectivity in 3D
        par_ccc.thresholds  = detect_th;
        [dum, dum, CC]      = multithreshstack_v4(image_filt_mask,par_ccc);
        
        %- Get centroid of each identified region
        CC_best = CC{1};
        S = regionprops(CC_best,'Centroid');
        N_spots = CC_best.NumObjects;

        centroid_linear  = [S.Centroid]';
        centroid_matrix  = round(reshape(centroid_linear,3,N_spots))';
        
        pos_pre_detect(:,1) = centroid_matrix(:,2);
        pos_pre_detect(:,2) = centroid_matrix(:,1);
        pos_pre_detect(:,3) = centroid_matrix(:,3);
        
end

%- Add coordinates if cropped    


%===  Remove spots which are close to edge of image and sort rows
% disp('... remove spots close to the edge ...');

ind_x   = (pos_pre_detect(:,1) > size_detect.xy) & (pos_pre_detect(:,1) <= dim.Y-size_detect.xy);
ind_y   = (pos_pre_detect(:,2) > size_detect.xy) & (pos_pre_detect(:,2) <= dim.X-size_detect.xy);


ind_in_cell = ind_x & ind_y;



%- Get coordinates of good spots
pos_spots_GOOD     = pos_pre_detect(ind_in_cell,:);   
pos_spots_GOOD_lin = sub2ind(size(image_filt), pos_spots_GOOD(:,1),pos_spots_GOOD(:,2),pos_spots_GOOD(:,3));


%- Assign values
spots_detected(:,1:3)  = pos_spots_GOOD;
spots_detected(:,10)   = image(pos_spots_GOOD_lin);
spots_detected(:,11)   = image_filt(pos_spots_GOOD_lin);
    % define boundaries for sub_spots images
spots_detected(:,4) = max(spots_detected(:,1)-box_spots_mosaic_xy,1);
spots_detected(:,5) = min(spots_detected(:,1)+box_spots_mosaic_xy,dim.Y);
spots_detected(:,6) = max(spots_detected(:,2)-box_spots_mosaic_xy,1);
spots_detected(:,7) = min(spots_detected(:,2)+box_spots_mosaic_xy,dim.X);
spots_detected(:,8) = max(spots_detected(:,3)-box_spots_mosaic_z,1);
spots_detected(:,9) = min(spots_detected(:,3)+box_spots_mosaic_z,dim.Z);    
      
%- Get CC only for good spots
if ~isempty(CC_best) && CC_best.NumObjects > 0
        
    %- Correct pixel-lists
    pixel_list_crop = CC_best.PixelIdxList(ind_in_cell);
    pixel_list_full = {};
    
    for i_list = 1:length(pixel_list_crop)
             
       %- Get list 
       list_loop =  pixel_list_crop{i_list};
       
       [y_sub,x_sub,z] = ind2sub(size(image_filt_mask),list_loop);
       
       %- Correct x & y      
       
       y_list = y_sub ;
       x_list = x_sub ;
       
       %- Make new list 
       pixel_list_full{i_list} = sub2ind(size(image_filt), y_list, x_list, z);
           
    end
    
    %- Get CC for best spots
    try 
        CC_GOOD                   = CC_best;
        CC_GOOD.NumObjects        = length(pixel_list_crop);
        CC_GOOD.PixelIdxList      = pixel_list_full;
        CC_GOOD.PixelIdxList_crop = pixel_list_crop;
    catch err
        disp('Error in spots_predetect_v17')
        err
    end
 
end


%=== Plot results of spot detection
%=   Subtract one from each value to center cross in pixel since we don't
%    have sub-pixel pointing accuracy (and the way matlab handles sub-pixel pointing)
%- Masked image for plot
img_mask.max_xy    = max(image_filt_mask,[],3);
img_mask.max_xz    = squeeze(max(image_filt_mask,[],1));
  %flag_struct.output=0;  
    
  if flag_struct.output
      %
      %     figure(4573), imagesc(imadjust(mat2gray(max(image_filt_mask,[],3)))), axis equal, axis image, colormap hot
      %     hold on, plot(pos_pre_detect(:,2),pos_pre_detect(:,1),'gx'), title('MIP Filtered image + predetections')

      %==== Actual plot
      h_fig = figure(486786);


      %     %- All spots
      %     h1 = subplot(2,2,1);
      %     imshow(img_mask.max_xy,[]); hold on
      %     plot(pos_pre_detect(:,2),pos_pre_detect(:,1),'gx','MarkerSize',5)
      %     hold off
      %     title('All detected spots')
      %
      %     h2 = subplot(2,2,2);
      %     imshow(img_mask.max_xz',[],'XData',[0 (dim.X-1)*pixel_size.xy],'YData',[0 (dim.Z-1)*pixel_size.z]); hold on
      %     plot(pos_pre_detect(:,2)*pixel_size.xy-pixel_size.xy,pos_pre_detect(:,3)*pixel_size.z-pixel_size.z,'gx','MarkerSize',5)
      %     hold off
      %     title('All detected spots')


      %- Spots away from edge
      pos_in_img = pos_pre_detect(ind_in_cell,:);

      % h3 = subplot(1,2,1);
      h3 = subplot(9,3,[14 15 17 18 20 21 23 24]);
      imshow(imadjust(mat2gray(img_mask.max_xy)),[]); hold on
      plot(pos_in_img(:,2),pos_in_img(:,1),'ro','MarkerSize',5)
      hold off
      title('XY')

      % h4 = subplot(2,2,2);
      h4 = subplot(9,3,[26 27]);
      imshow(img_mask.max_xz',[],'XData',[0 (dim.X-1)*pixel_size.xy],'YData',[0 (dim.Z-1)*pixel_size.z]); hold on
      plot(pos_in_img(:,2)*pixel_size.xy-pixel_size.xy,pos_in_img(:,3)*pixel_size.z-pixel_size.z,'ro','MarkerSize',5)
      hold off
      title('XZ')

      sgtitle({'MIP Filtered image + predetections (red circle)','Press any key to continue'})

      %     linkaxes([h1,h3], 'xy');
      %     linkaxes([h2,h4], 'xy');

      linkaxes([h3], 'xy');
      linkaxes([h4], 'xy');

      set(h_fig,'Color','w')
      disp('Paused - Press Enter')
      pause()

  end
