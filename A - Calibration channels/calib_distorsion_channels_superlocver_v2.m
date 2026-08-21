function calib_distorsion_channels_superlocver_v2(file_calib,path_calib,N_channels)
% calib_distorsion_channels  Registrate and calibrate all the channels together;
%   [T]=calib_distorsion_channels(file,path,N_channels)
%
%   Input: file = filename of the image
%          path = path folder
%          N_channels = number of channels to be calibrated in the raw data
%          rk: if the image file given is a stack with several slices, the function will show you all the slices and you have to choose one slice only for the calibration
%   Output: save calib files


tiff_file=fullfile(path_calib,file_calib);
h = imfinfo(tiff_file).Height;
w = imfinfo(tiff_file).Width;
stack_size=size(imfinfo(tiff_file),1);
imageData=zeros(h,w,stack_size);

for ii=1:stack_size
   imageData(:,:,ii) =  double(imread(tiff_file,ii));
end
imageData_MIP = max(imageData,[],3);

figure(531), imagesc(imageData_MIP), colormap gray, axis('image')
% imcontrast;pause();
% 
% Areas = zeros(N_channels,4); %Position of all rectangles
% msgbox(['Select the rectangle n°1 (create the rectangle - then double-click)'],'Hi','modal')
% roi = drawrectangle(); wait(roi);
% Areas(1,:)= round(roi.Position);
% for i_channel = 2:N_channels
%     msgbox(['Select the rectangle n°' num2str(i_channel) ' (drag the rectangle - then double-clik)'],'Hi','modal')
%     wait(roi);
%     Areas(i_channel,:)= round(roi.Position);
%     
% end
% IF MFM well aligned on set up:
hh=ceil(h/3);
ww=ceil(w/2);

Areas = [   1 1 hh ww;
    1+hh 1 hh ww;
    1+2*hh 1 hh ww;
    1 1+ww hh ww;
    1+hh 1+ww hh ww;
    1+2*hh 1+ww hh ww;];



hold on,
for ii=1:N_channels
    rectangle('Position',Areas(ii,:),'EdgeColor','b','LineWidth',3)
end                                 
%% EVERYTHING BELOW IS FROM PIPELINE N PLANES MFM
%% REMOVE BACKGROUND
maxprojection=imageData_MIP;
background = imopen(maxprojection,strel('disk',10));
maxprojectionBG = maxprojection - background;
%% DETECT PEAKS

button = 'No';

while strcmp(button,'No')
    
    
    Center_of_mass = 1; Binarization = abs(1-Center_of_mass); %USER CHOICE
    aa= 10; % window half-size for center detection
    

    figure(531),imagesc(maxprojectionBG), title('Background remove'),
    colormap gray, axis equal, axis image

    hold on,
    for ii=1:N_channels
        rectangle('Position',Areas(ii,:),'EdgeColor','b','LineWidth',3)
    end
    title({'Click on each bead ONCE only.',...
        'Each bead appears in 6 channels but select it in any one channel, not all 6.', ...
        'Contrast can be asjusted'})
    uiwait(msgbox(['Bead selection instructions:' newline newline ...
        '- Click each bead ONCE only' newline ...
        '- Each bead appears in 6 channels - pick any one channel' newline ...
        '- Do NOT select the same bead in multiple channels' newline ...
        ' You can adjust contrast in imconstrast window' newline newline ...
        'Press Enter when done'], ...
       'Instructions', 'modal'));
    imcontrast;
    [xi,yi] = getpts;
    hold on, plot(xi,yi,'ro')
    nb_points_init=size(xi,1);
    coord = zeros(nb_points_init*N_channels,2);
    coord(1:nb_points_init,:) = [round(xi) round(yi)];
    %translate tothe other channels
    r_matching = {'9','35'}; % ADJUST THIS PARAMETER TO MATCH COORDINATES BETWEEN CHANNELS
    % inputdlg({'r1','r2'},'Coordinate matching',1,{'9','35'});%USER CHOICE
    r1=str2double(cell2mat(r_matching(1))); r2=str2double(cell2mat(r_matching(2)));
    for i = 1:nb_points_init
        if((coord(i,1) > Areas(1,1)-1+aa)&&(coord(i,1) < Areas(1,1)+Areas(1,3)-1-aa) && (coord(i,2) < Areas(1,2)+Areas(1,4)-1-aa)&& (coord(i,2) > Areas(1,2)-1+aa)); %first tile
            CH = 1;
            coord((CH-1)*nb_points_init+i,:) = [round(xi(i)) round(yi(i))];
            for ch=[2 3 4 5 6]
                [ch_i, ch_j] = ind2sub([3 2], ch);
                coord((ch-1)*nb_points_init+i,:) =  [coord((CH-1)*nb_points_init+i,1)+(ch_i-1)*(hh-r1) coord((CH-1)*nb_points_init+i,2)+(ch_j-1)*(ww-r2)];
            end
        elseif ((coord(i,1) > Areas(2,1)-1+aa)&&(coord(i,1) < Areas(2,1)+Areas(2,3)-1-aa) && (coord(i,2) < Areas(2,2)+Areas(2,4)-1-aa)&& (coord(i,2) > Areas(2,2)-1+aa)); %second tile
            CH = 2;
            coord((CH-1)*nb_points_init+i,:) = [round(xi(i)) round(yi(i))];
            for ch=[1 3 4 5 6]
                [ch_i, ch_j] = ind2sub([3 2], ch);
                coord((ch-1)*nb_points_init+i,:) =  [coord((CH-1)*nb_points_init+i,1)+(ch_i-2)*(hh-r1) coord((CH-1)*nb_points_init+i,2)+(ch_j-1)*(ww-r2)];
            end
        elseif ((coord(i,1) > Areas(3,1)-1+aa)&&(coord(i,1) < Areas(3,1)+Areas(3,3)-1-aa) && (coord(i,2) < Areas(3,2)+Areas(3,4)-1-aa)&& (coord(i,2) > Areas(3,2)-1+aa)); %third tile
            CH = 3;
            coord((CH-1)*nb_points_init+i,:) = [round(xi(i)) round(yi(i))];
            for ch=[1 2 4 5 6]
                [ch_i, ch_j] = ind2sub([3 2], ch);
                coord((ch-1)*nb_points_init+i,:) =  [coord((CH-1)*nb_points_init+i,1)+(ch_i-3)*(hh-r1) coord((CH-1)*nb_points_init+i,2)+(ch_j-1)*(ww-r2)];
            end
        elseif ((coord(i,1) > Areas(4,1)-1+aa)&&(coord(i,1) < Areas(4,1)+Areas(4,3)-1-aa) && (coord(i,2) < Areas(4,2)+Areas(4,4)-1-aa)&& (coord(i,2) > Areas(4,2)-1+aa));  %fourth tile
            CH = 4;
            coord((CH-1)*nb_points_init+i,:) = [round(xi(i)) round(yi(i))];
            for ch=[1 2 3 5 6]
                [ch_i, ch_j] = ind2sub([3 2], ch);
                coord((ch-1)*nb_points_init+i,:) =  [coord((CH-1)*nb_points_init+i,1)+(ch_i-1)*(hh-r1) coord((CH-1)*nb_points_init+i,2)+(ch_j-2)*(ww-r2)];
            end
        elseif ((coord(i,1) > Areas(5,1)-1+aa)&&(coord(i,1) < Areas(5,1)+Areas(5,3)-1-aa) && (coord(i,2) < Areas(5,2)+Areas(5,4)-1-aa)&& (coord(i,2) > Areas(5,2)-1+aa));  %fifth tile
            CH = 5;
            coord((CH-1)*nb_points_init+i,:) = [round(xi(i)) round(yi(i))];
            for ch=[1 2 3 4 6]
                [ch_i, ch_j] = ind2sub([3 2], ch);
                coord((ch-1)*nb_points_init+i,:) =  [coord((CH-1)*nb_points_init+i,1)+(ch_i-2)*(hh-r1) coord((CH-1)*nb_points_init+i,2)+(ch_j-2)*(ww-r2)];
            end
        elseif ((coord(i,1) > Areas(6,1)-1+aa)&(coord(i,1) < Areas(6,1)+Areas(6,3)-1-aa) & (coord(i,2) < Areas(6,2)+Areas(6,4)-1-aa)& (coord(i,2) > Areas(6,2)-1+aa));  %sixt tile
            CH = 6;
            coord((CH-1)*nb_points_init+i,:) = [round(xi(i)) round(yi(i))];
            for ch=[1 2 3 4 5]
                [ch_i, ch_j] = ind2sub([3 2], ch);
                coord((ch-1)*nb_points_init+i,:) =  [coord((CH-1)*nb_points_init+i,1)+(ch_i-3)*(hh-r1) coord((CH-1)*nb_points_init+i,2)+(ch_j-2)*(ww-r2)];
            end
        else
            coord(i,:) = [0 0];
        end
    end
    % remove zeros (comes from points too close to the edge)
    coordWithZeros = (coord(:,1)~=0);
    all_coord = coord(coordWithZeros,:);
    nb_points = size(all_coord,1)/N_channels;
    hold on, plot(all_coord(:,1),all_coord(:,2),'go')
    out=zeros(nb_points*N_channels,2);
    for i=1:nb_points*N_channels
        wd_width = 7;
        window = maxprojection(all_coord(i,2)-wd_width:round(all_coord(i,2))+wd_width,round(all_coord(i,1))-wd_width:round(all_coord(i,1))+wd_width);
        [~,ind]=max(window(:));
        [out(i,1), out(i,2)]=ind2sub([wd_width*2+1 wd_width*2+1],ind);
    end
    centroids=[out(:,2)+all_coord(:,1)-(wd_width+1),out(:,1)+all_coord(:,2)-(wd_width+1)];
    
    %%%%%%%%%%%%%%%%%%
    
    [~,order] = sort(centroids,1);
    centroids=centroids(order(:,1),:);
    figure(3345453),
    imagesc(maxprojection),colormap gray
    hold on
    scatter(centroids(:,1),centroids(:,2), 'ro')
    axis equal, axis image; title({'You can Adjust the Contrast','Press Enter for next step'})
    imcontrast; pause()
    
    %question
    button = questdlg('Do you validate the bead selection ? Choose no if you want to select beads again');
end
%% COORDINATES MATCHING


%% group and split the localizations between the 9 tiles and remove the one near the edge
ws= 3;  %window half-size for gaussian fitting

index1 =find((centroids(:,1) > Areas(1,1)-1+ws)&(centroids(:,1) < Areas(1,1)+Areas(1,3)-1-ws) & (centroids(:,2) < Areas(1,2)+Areas(1,4)-1-ws)& (centroids(:,2) > Areas(1,2)-1+ws)); %first tile
index2 =find((centroids(:,1) > Areas(2,1)-1+ws)&(centroids(:,1) < Areas(2,1)+Areas(2,3)-1-ws) & (centroids(:,2) < Areas(2,2)+Areas(2,4)-1-ws)& (centroids(:,2) > Areas(2,2)-1+ws)); %second tile
index3 =find((centroids(:,1) > Areas(3,1)-1+ws)&(centroids(:,1) < Areas(3,1)+Areas(3,3)-1-ws) & (centroids(:,2) < Areas(3,2)+Areas(3,4)-1-ws)& (centroids(:,2) > Areas(3,2)-1+ws)); %third tile
index4 =find((centroids(:,1) > Areas(4,1)-1+ws)&(centroids(:,1) < Areas(4,1)+Areas(4,3)-1-ws) & (centroids(:,2) < Areas(4,2)+Areas(4,4)-1-ws)& (centroids(:,2) > Areas(4,2)-1+ws));  %fourth tile
index5 =find((centroids(:,1) > Areas(5,1)-1+ws)&(centroids(:,1) < Areas(5,1)+Areas(5,3)-1-ws) & (centroids(:,2) < Areas(5,2)+Areas(5,4)-1-ws)& (centroids(:,2) > Areas(5,2)-1+ws));  %fifth tile
index6 =find((centroids(:,1) > Areas(6,1)-1+ws)&(centroids(:,1) < Areas(6,1)+Areas(6,3)-1-ws) & (centroids(:,2) < Areas(6,2)+Areas(6,4)-1-ws)& (centroids(:,2) > Areas(6,2)-1+ws));  %sixt tile

    u = cell(0);
    u{1} = index1;
    u{2} = index2;
    u{3} = index3;
    u{4} = index4;
    u{5} = index5;
    u{6} = index6;
    
%             V=get(u,'Value')
%             idxactive=find([V{:}]~=0);
    
            for iii=1:N_channels
%                 HH=get(u(idxactive(iii)),'userdata')
%                 centroids_1bis0{iii}= centroids(HH,:);
                centroids_1bis0{iii}= centroids(u{iii},:);
                [a_sorted,order] = sort(centroids_1bis0{iii},1);
                centroids_1bis0{iii}=centroids_1bis0{iii}(order(:,1),:);
                %centroids_1bis0{4}= centroids(index8,:);
           
            end
            if ~isdeployed
    assignin('base','centroids_1bis0',centroids_1bis0);
    assignin('base','maxprojection',maxprojection);
            end
            %% pick the planes of interrest
            
            idxref = 5; % FIXED for now
            translation=[
                Areas(idxref,1)-Areas(1,1)-r1	Areas(idxref,2)-Areas(1,2)-r2;
                Areas(idxref,1)-Areas(2,1)     Areas(idxref,2)-Areas(2,2)-r2;
                Areas(idxref,1)-Areas(3,1)+r1  Areas(idxref,2)-Areas(3,2)-r2;
                Areas(idxref,1)-Areas(4,1)-r1  Areas(idxref,2)-Areas(4,2);
                Areas(idxref,1)-Areas(5,1)     Areas(idxref,2)-Areas(5,2);
                Areas(idxref,1)-Areas(6,1)+r1  Areas(idxref,2)-Areas(6,2);];
            
            
            if ~isdeployed
                assignin('base','translation',translation)
         end
         
         
         % aa=str2num(get(handles.edit_window_size,'String'));
         for ii=1:N_channels
             idx = u{ii};
             centroids_1{ii}=centroids(idx,:);
             [a_sorted,order] = sort(centroids_1{ii},1);
             centroids_1{ii}=centroids_1{ii}(order(:,1),:);
            centroids_1b{ii}=centroids_1{ii}+translation(ii,:);
             
         end
         
      
         hold on
         Markers={'c<','rs','go','m>','bv','r^','g*','yh','rd','y.','mx','yp'};


for ik=1:N_channels
    

hold on

 plot(centroids_1 {ik}(:,1),centroids_1 {ik}(:,2),Markers{ik},'MarkerSize',16)
end

% figure
% for ik=1:N_channels
% hold on 
% plot(centroids_1b {ik}(:,1),centroids_1b {ik}(:,2),Markers{ik})
% end
% set(gca,'Ydir','reverse'),
% axis equal
assignin('base','centroids_1',centroids_1)
              
% hold on
% scatter(centroids_1{idxref}(:,1),centroids_1{idxref}(:,2),70,'rs','filled')

        
     
     assignin('base','centroids_1b',centroids_1b)
for i=1:N_channels
c = xcorr2(centroids_1b{idxref},centroids_1b{i});
%[ssr,snd] = max(abs(c(:)));
% figure(1)
% plot(c(:))
% title('Cross-Correlation')
% hold on
% plot(snd,ssr,'or')
% 
% text(snd*1.05,ssr,'Maximum')
% hold off
[max_cc,imax] = max(abs(c(:)));
[ypeak,xpeak] = ind2sub(size(c),imax(1));
corr_offset = [(ypeak-size(centroids_1b{idxref},1)) (xpeak-size(centroids_1b{idxref},2))];

centroids_1b{i}=centroids_1b{i}-corr_offset;
% figure(453153)
% hold on 
% plot(centroids_1b {i}(:,1),centroids_1b {i}(:,2),Markers{i})
% title('xcorr applied')


end

Nref=size(centroids_1{idxref},1);
XY1=squeeze (centroids_1b{idxref}(:,:))  ;
clear XY0 XY2
XYo=[];
XY2=[];
for ii=1:N_channels
%                 if ii~=1 % WHYYYYYY??????????
    XYo=[XYo;squeeze( centroids_1{ii}(:,:))];
    Npts(ii)=size(centroids_1b{ii},1);
    XY2=[ XY2;squeeze(centroids_1b{ii}(:,:))];
%                 end
    
end

        
%% clustering the coordinates to find the centers of a groupped positions
 %XY3=[XY1;XY2];%% all the coordinates of all the active planes in an arry 
 XY3=[XY2];%% all the coordinates of all the active planes in an arry 

 [idx,C] = kmeans(XY3,Nref,'Start',XY1);% Kmean clustering; C containes the center of the clusters

% figure
% h=gscatter(XY3(:,1),XY3(:,2),idx)
% hold on
% plot(C(:,1),C(:,2),'ko')
% hold on 
% plot(XY1(:,1),XY1(:,2),'>r')
% hold on
% set(gca,'Ydir','reverse')
%% translating XY2 coordinates to the center


%% remove the one that don't have a match in the other planes
assignin('base','XY1',XY1);
assignin('base','XY2',XY2);
assignin('base','C',C);
clear idx

idx = [];
% button = 'No';
% while strcmp(button,'No')
    % r can be adjusted to match the coordinates between channels
    % r_answer = inputdlg('distance for matching coordinates ?','Input',1,{'12'});
    r= 12; %str2double(cell2mat(r_answer));
    Centro={};
    idx= rangesearch(XY1,XY2,r);
    m=find(cellfun(@(v)isempty(v),idx)==1);
    for ij=1:length(m)
        idx{m(ij)}=0;
    end
    idx1 = zeros(sum(Npts),1);
    for ik=1:size(idx,1)
        if length(idx{ik})>1
            idx1(ik)=0;
        else
            idx1(ik)=idx{ik};
        end
    end
    %     idx=[idx{:}]; %convert from cell to matrix

    k=1;
    for i=1:length(centroids_1 {1}(:,1))
        index=(find(idx1==i));

        if length(index) == N_channels 
            for j=1:N_channels
                Centro{j}(k,:)=XYo(index(j),:);
            end
            k=k+1;
        end

    end
    Centro=Centro(~cellfun('isempty',Centro));

    figure, imagesc((maxprojection));hold on
    for u=1:size(Centro,2)
        if ~isempty(Centro{u})
            scatter(Centro{u}(:,1),Centro{u}(:,2),70,'o')
        else
            errordlg('no matching coordinates found within the specified range');

        end
        set(gca,'Ydir','reverse')
    end
    axis image, axis equal
    title('Matching Coordinates - one color per channel')

    % %question
    % button = questdlg('Do you validate bead channels assignation?');

% end
%% fit with a gaussian

for k=1:N_channels
    for i=1:size(Centro{k},1)
        start_location=round(Centro{k}(i,:));
        coordinates{k}(i,:)=gaussian_fitting(maxprojection,start_location,ws);
    end
end
 
 
%% CALIBRATION
filename= fullfile(path_calib,file_calib);
assignin('base','maxprojection',maxprojection)
z_step_answer = inputdlg('Z scan step size (nm)?','Input',1,{'100'});
z_step= str2double(cell2mat(z_step_answer));


global coordinates_finals
c_f=coordinates_finals;
coord_final={};
% readjust coordinates to be able to work with
for i=1:size(Centro{1},1)
    for j=1:N_channels
        coord_final{i}(j,:)= coordinates {j} (i,:);%[cellfun(@(v)v(i,1),coordinates),cellfun(@(v)v(i,2),coordinates)]; 
    end
end
%assignin('base','c_f',c_f)
%% Z_calibration and Transformation Matrix toguether in one function
param.Z_Calibration = 1; % USER DEFINE: do the calibration
param.Generate_transformation_matrix = 1;% USER DEFINE: generate transfo matrix
param.Areas = Areas;

Calibration6planes_4PMFM_v1(filename,maxprojection,coord_final,N_channels,z_step,ws,param);

end