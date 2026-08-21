% code to reconstruct and correct for intensity, of stack of images 
% 11/11/2023: v2 implementation
%   modify the intensity correction : use factors provided by (unpolarization) lamp, computed with calib_intensity_correction.m 
% 01/12/2023: fix the name definition of the reconstructed files

function [border,fullname_fr1,fullpath_reco] = reconstruct_channel_images_superlocver_v2(PathName,FileName,fullfile_calib,fullfile_calib_int,nb_planes)

    % [FileName_transfo, PathName_transfo]=uigetfile('*.mat','Select the transformation matrix',path_calib);
    % calibrationFile=[path_calib,'image_Pos0.ome.mat'];
    calibrationFile = [fullfile_calib(1:end-3) 'mat'];

    % [FileName_intcorr, PathName_intcorr]=uigetfile('*.mat','Select the intensity correction matrix',path_calib_int);
    % correctionfactorFile=[path_calib_int,'image_Pos0.ome_intensity_corr_factors.mat'];
    correctionfactorFile = [fullfile_calib_int(1:end-4) '_intensity_corr_factors.mat'];

    % [FileName, PathName]=uigetfile('*.tif','Select the file to reconstruct',path_data,'MultiSelect','on');


    %% read calibration parameters and correction 

    % load transformation matrix
    load(calibrationFile,'transformPerSlice','cornerPerSlice');
    % load([calibrationFile(1:end-4) '_pos-int-cor.mat'],'correctionfactor')
    
    intensity_corr = 0; % = 0 if from lamp calibration / else if from bead calibration
    if intensity_corr == 0
    % load intensity correction and pick the correction factors corresponding to your emission wavelength
    load(correctionfactorFile,'I_corr_factor')
%     load(correctionfactorFile,'list_wavelength')
%     diff_filter_em = zeros(1,size(list_wavelength,2));
%     for ii = 1:size(list_wavelength,2)
%         em_filter = str2num(list_wavelength{ii}(1:3));
%         diff_filter_em(1,ii) = em_filter-microscope_settings.Em;
%     end
%     [~,idxmin]=min(abs(diff_filter_em));
    correctionfactor = I_corr_factor{1};
    else 
        load(correctionfactorFile,'correctionfactor');
    end
% warning('Intensity correction factors: now from bead calib file - uncomment l.23-32 to correct intensity from lamp acquisition')
    h = waitbar(0,'Initializing ...');
%% read image to reconstruct and initialize reconstrution data
    imFilename=[PathName, FileName];
    info=imfinfo(imFilename);
    numFrames=length(info);
    [~,Name,Ext] = fileparts(FileName);
    idx_ome = strfind(Name,'.ome'); %look for '.ome' in file name
    if idx_ome>0
        Name = Name(1:idx_ome-1); % remove .ome extension
    end
    folder_reco = [Name '_reco'];
    mkdir(char(PathName),folder_reco); % create new folder that will contain reconstructed data
    fullpath_reco=[PathName folder_reco];
    %%
    stack=zeros(cornerPerSlice(1,3)-cornerPerSlice(1,1)-3,cornerPerSlice(1,4)-cornerPerSlice(1,2)-3,nb_planes,'uint16');
    mask=true(size(stack,1),size(stack,2));

    % prompt = {'ch#1','ch#2','ch#3','ch#4','ch#5','ch#6'};
    % dlgtitle = 'Enter camera base level';
    % dims = [1 45];
    % definput = {'150','150','150','150','150','150'};
    % answer = inputdlg(prompt,dlgtitle,dims,definput);
    % averreadnoise=[str2double((answer{1})),str2double((answer{2})),str2double((answer{3})),...
    %     str2double((answer{4})),str2double((answer{5})),str2double((answer{6}))];
    averreadnoise = 150*ones(1,6); % camera background level in each channel; 
    % init name of the reconstructed images for frame #1 
    nb_digits=numel(num2str(numFrames)); 
    fullname_fr1 =  [Name '_' num2str(1,['%.' num2str(nb_digits) 'd']) '.tif'];

    for ff=1:numFrames

        im=imread(imFilename,'Index',ff);
        writeFilename_ff = fullfile(fullpath_reco,[Name '_' num2str(ff,['%.' num2str(nb_digits) 'd']) '.tif']);

        for kk=1:nb_planes
            patch=im(cornerPerSlice(kk,1):cornerPerSlice(kk,1)+size(stack,1)-1,cornerPerSlice(kk,2):cornerPerSlice(kk,2)+size(stack,2)-1);
            [X Y]=meshgrid(0:size(patch,1)-1,0:size(patch,2)-1);%to match peak coordinates from calibrationBEads.m routine
            X=X';
            Y=Y';
            %apply transformation to meshgrid
            aux= ([X(:) Y(:)]-repmat(transformPerSlice{kk}.c(1,:),[size(patch,1)*size(patch,2) 1])) * transformPerSlice{kk}.T'/transformPerSlice{kk}.b;

            XI=reshape(aux(:,1),size(X));
            YI=reshape(aux(:,2),size(Y));
            %interpolate images
            stack(:,:,kk) = uint16(interp2(Y,X,double(patch),YI,XI,'*cubic'));%uint16
            if(ff==1)
                    %find maximum region without zero-border filling
                    %%maskPos=find(XI<cornerPerSlice(1,2) | YI<cornerPerSlice(1,1) | XI>size(patch,1) | YI>size(patch,2))
                    maskPos=find(XI<0 | YI<0 | XI>size(patch,1)-1 | YI>size(patch,2)-1); % change here
                    mask(maskPos)=false;

                end
        end

        perc=ff/numFrames*100;
        waitbar(perc/100,h,sprintf(' processing %3.2f%% ',perc))
        %  waitbar(perc/100,h,sprintf(' processing %3.2f%% of the stack %2d/%2d total',perc,i,numStack))
        if(ff==1)
            %find cropping area to make sure we do not include zero-filling after
            %aligning the images
            border=1;
            % choice = questdlg('Do you want to specify the border value?', ... % One could fix one border value.
            %     'Border value', ...
            %     'Yes','No','No');
            choice = 'No';
            switch choice
                case 'Yes'
                    border=input('Please specify the border value:');
                case 'No'
                    while(sum(mask(border:end-border+1,border:end-border+1)==false)>0)
                        border=border+1;
                    end
            end
        end
        %% Intensity correction
        corr=1./correctionfactor;
        C = zeros(size(stack));
        ij=0;
        for k=[1 4 2 5 3 6]% 1:nb_planes %
            ij = ij+1;
          
            C(:,:,ij)=((stack(:,:,k)-averreadnoise(k))*corr(k));
            if k==1
                imwrite(uint16(C(border:end-border+1,border:end-border+1, 1)),writeFilename_ff, 'Compression', 'none', 'WriteMode', 'overwrite');
                  else
                imwrite(uint16(C(border:end-border+1,border:end-border+1, ij)),writeFilename_ff, 'Compression', 'none', 'WriteMode', 'append');
            end
            %                 end
            %                 figure(10)
            %                 imagesc(C(border:end-border+1,border:end-border+1, k) );
            %MMM(:,:,k)=C(border:end-border+1,border:end-border+1, k) ;
            %                 end


            %% Saving
            %     outF=[outFilename num2str(ff-1,'%.4d') '.tif'];
            %     imwrite(stack(border:end-border+1,border:end-border+1, 1), outF, 'Compression', 'none', 'WriteMode', 'overwrite');
            %     for kk=2:size(stack,3)
            %         imwrite(stack(border:end-border+1,border:end-border+1, kk), outF, 'Compression', 'none', 'WriteMode', 'append');
            %     end
        end


        %%
        %perc=100;
        %             waitbar(perc/100,h,sprintf('Finish. Images in the ''Reconstructed'' subfolder',perc))

    end
    waitbar(100/100,h,sprintf(' DONE '))
    disp(['Cropped ' num2str(border) ' pixels to avoid zero-filling after alignment'])
    border_value_name =[fullpath_reco,'\' ,FileName(1:end-4), '_border_value.mat']
    save(border_value_name,'border')
end