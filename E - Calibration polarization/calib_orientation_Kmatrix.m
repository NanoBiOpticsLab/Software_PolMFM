function [K] = calib_orientation_Kmatrix(file_calib_polar,path_calib_polar)

% Experimental acquisition for calibration: white light, nb_fr for each
% wavelength, arbitrary number of wavelength. Write on disk, in folder
% named 'folder_0', 'folder_1', etc...; arbitrary step in degree; go from
% 0° to 180°
% Input: path_calib_polar & file_calib_polar folder and file of the first
% calibration data (polarizer at 0°)
% Output: Propagation matrix K. 
%       K is a 6x3 matrix; 
%       line order is [0°/-dz; 90°/-dz; 45°; 135°; 0°/+dz; 90°/+dz]; 
%       column order is [XX YY XY] 
%
% Louise Régnier; written 04/07/2023

N_channels = 6;

answer = inputdlg({'Polarizer angular step used for calibration (°) ?','Maximum angle used for polarization calibration (°) ?'},...
    'Input',1,{'5','180'});
if isempty(answer)
    error('User cancelled input.');
end
calib_step_degree = str2double(answer{1});
max_angle = str2double(answer{2});
if isnan(calib_step_degree) || isnan(max_angle)
    error('Invalid numeric input for calibration step or maximum angle.');
end

pos_idx = strfind(path_calib_polar,'_0');
max_idx = max_angle/calib_step_degree; 

clear I K Areas Area_bkg
I = zeros(N_channels,max_idx+1);
K = zeros(N_channels,3); % K matrix: 4x3 Matrice + 2x3 matrice for duplicate of 0-90° polar
Areas = zeros(N_channels,4); %Position of all rectangles
Area_bkg = zeros(1,4);

% warning('turn polarizer counterclockwise - comment and replace  lg 41 & 59 & 82 & 92 for clockwise rotation')
ipolar_count=0;
for ipolar = max_idx:-1:0 %  /!\/!\/!\/!\ HERE /!\/!\/!\/!\ counterclockwise
% for ipolar = 0:max_idx % /!\/!\/!\/!\ HERE /!\/!\/!\/!\ clockwise  
    
   fullfile = [path_calib_polar(1:pos_idx-1) '_' num2str(ipolar) path_calib_polar(pos_idx+2:end) file_calib_polar];
   if isfile(fullfile)
   info=imfinfo(fullfile);
   nb_fr = size(info,1);
   
       im = zeros(info(1).Width,info(1).Height,nb_fr);
       
       for fr = 1:nb_fr % get the image for the wavelength selected
           im(:,:,fr) = imread(fullfile,'Index',fr);
       end

       im_avg = mean(im,3);

       if ipolar == max_idx %  /!\/!\/!\/!\ HERE /!\/!\/!\/!\ counterclockwise %if first polarization: select the areas to compute the intensity
           % if ipolar==0 % %  /!\/!\/!\/!\ HERE /!\/!\/!\/!\ clockwise %if first polarization: select the areas to compute the intensity
           figure(4534503); imagesc(im_avg); colormap gray; axis equal; axis image; title({'Polarized calibration (1st image)','Channels and bgd selected'})
           Areas =  [ 43    80    92   139;  47   304    92   139; 213    78    92   139; ...
               214   297    92   139;   379    81    92   139;    379   303    92   139];
           for i_channel = 1:N_channels
               rectangle('Position',Areas(i_channel,:),'EdgeColor','blue','SelectionHighlight','on');
           end
           roi.Position = [0.5 0.5 11 43]; rectangle('Position',roi.Position,'EdgeColor','blue','SelectionHighlight','on');
           Area_bkg = round(roi.Position);
       end

       for ch=1:N_channels
%            I(ch,ipolar+1) = mean(im_avg(Areas(ch,2):Areas(ch,2)+Areas(ch,4)-1,Areas(ch,1):Areas(ch,1)+Areas(ch,3))-1,'all')...
%                - mean(im_avg(Area_bkg(2):Area_bkg(2)+Area_bkg(4)-1,Area_bkg(1):Area_bkg(1)+Area_bkg(3)-1),'all');  %  /!\/!\/!\/!\ HERE /!\/!\/!\/!\ clockwise 
           I(ch,ipolar_count+1) = mean(im_avg(Areas(ch,2):Areas(ch,2)+Areas(ch,4)-1,Areas(ch,1):Areas(ch,1)+Areas(ch,3))-1,'all')...
               - mean(im_avg(Area_bkg(2):Area_bkg(2)+Area_bkg(4)-1,Area_bkg(1):Area_bkg(1)+Area_bkg(3)-1),'all'); %  /!\/!\/!\/!\ HERE /!\/!\/!\/!\ counterclockwise 
       end
       
   
   else
       
        for ch=1:N_channels
%            I(ch,ipolar+1) = 0; %  /!\/!\/!\/!\ HERE /!\/!\/!\/!\ clockwise 
            I(ch,ipolar_count+1) = 0;  %  /!\/!\/!\/!\ HERE /!\/!\/!\/!\ counterclockwise 
       end
       
   end
   ipolar_count=ipolar_count+1;
   
end
%% plot calibration curves &  fit a sinus to find absolute polarization (where is 0° ?)
clear s
clear Phi
      Phi = zeros(1,N_channels);
      close(figure(456153))
       figure(456153);subplot(1,3,[1 2]); 
       hold on; title('Channel polarization calibration curves')
       for ch=1:N_channels
           plot(0:calib_step_degree:max_angle,I(ch,:))
           [s{ch}, fit] = fit_cos2(I(ch,:),calib_step_degree,max_angle);
           xp = linspace(0,max_angle);
           hold on, plot(xp,fit(s{ch},xp),'-.');
           if s{ch}(1)<0 % if coeff1 is negative: phase shift of 90°
               Phi(ch) = rem(wrapTo360(s{ch}(2)+90),180) ;
           else
               Phi(ch) =  rem(wrapTo360(s{ch}(2)),180) ;
           end
           disp(['channel ' num2str(ch) ': phi = ' num2str(Phi(ch))])
       end
       xlabel('Polarization [°]'); ylabel('Intensity')
       legend(['0°- MFM=-dz'],['Fit 0°- MFM=-dz'],['90°- MFM=-dz'],['Fit 90°- MFM=-dz'],['45°- MFM=0'],['Fit 45°- MFM=0'],...
           ['135°- MFM=0'],['Fit 135°- MFM=0'],['0°- MFM=+dz'],['Fit 0°- MFM=+dz'],['90°- MFM=+dz'],['Fit 90°- MFM=+dz'],'Location','eastoutside')
       grid on

  
%% Build K matrix 
% find absolute angle polar 0: 
%     phi0 = mean([-s{1}(2) -s{2}(2) -s{5}(2) -s{6}(2)])
    phi0 = mean([rem(wrapTo360(Phi(1)),180) rem(wrapTo360(Phi(5)),180) ...
        rem(wrapTo360(Phi(2)-90),180) rem(wrapTo360(Phi(6)-90),180) ...
         rem(wrapTo360(Phi(3)+45),180) rem(wrapTo360(Phi(4)+135),180)])
    %polarizer at 0° 
    I01 = fit(s{1},phi0);
    I02 = fit(s{2},phi0);
    I03 = fit(s{3},phi0);
    I04 = fit(s{4},phi0);
    I05 = fit(s{5},phi0);
    I06 = fit(s{6},phi0);
    A2 = I01+I02+I03+I04+I05+I06;
    K(1,1) = I01/A2;
    K(2,1) = I02/A2;
    K(3,1) = I03/A2;
    K(4,1) = I04/A2;
    K(5,1) = I05/A2;
    K(6,1) = I06/A2;
    %polariser at 90°
    I901 = fit(s{1},phi0+90);
    I902 = fit(s{2},phi0+90);
    I903 = fit(s{3},phi0+90);
    I904 = fit(s{4},phi0+90);
    I905 = fit(s{5},phi0+90);
    I906 = fit(s{6},phi0+90);
    A2 = I901+I902+I903+I904+I905+I906;
    K(1,2) = I901/A2;
    K(2,2) = I902/A2;
    K(3,2) = I903/A2;
    K(4,2) = I904/A2;
    K(5,2) = I905/A2;
    K(6,2) = I906/A2;
    %polariser at 45°
    I451 = fit(s{1},phi0+45);
    I452 = fit(s{2},phi0+45);
    I453 = fit(s{3},phi0+45);
    I454 = fit(s{4},phi0+45);
    I455 = fit(s{5},phi0+45);
    I456 = fit(s{6},phi0+45);
    A2 = I451+I452+I453+I454+I455+I456;
    K(1,3) = (2*I451-I01-I901)/A2;
    K(2,3) = (2*I452-I02-I902)/A2;
    K(3,3) = (2*I453-I03-I903)/A2;
    K(4,3) = (2*I454-I04-I904)/A2;
    K(5,3) = (2*I455-I05-I905)/A2;
    K(6,3) = (2*I456-I06-I906)/A2;
figure(456153), subplot(133), imagesc(K), title('2nd order dipole moment "K" matrix'), 
ylabel('Channel'), yticklabels({'0-', '90-', '45', '135', '0+','90+'})
xticks([1 2 3]); xticklabels({'\mu_x \mu_x','\mu_y \mu_y','\mu_x \mu_y',}); axis image, axis equal
colorbar    






