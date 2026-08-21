function Calibration6planes_4PMFM_v1(filename,maxprojection,coord_final,N_channels,z_step,ws,param)

%function Generate_transformation_matrix (filename,maxprojection,coord,nb)
% From CalibrationNplanes.m written byBassam Hajj
% 31/12/2013

% function to generate transformation matrix between different panels
% it finds the transformation to a mean value between all panels
%( it should be changed to the central positions only)

% inputs:
%
%     - filename
%     - corrd : fitted coordinates
%     - z-step: in nm
%     - ws : window half size
%     - nb : number of planes
% outputs : generates a transformation matrix

% Modified By Yasmina Fedala
% 18/07/2019
% merge Z_calibration function and

% Adapted to 4P_MFM 6 planes by Louise Régnier 
% 30/05/2023
%% Channel numbering correspondance

PolMFM_ch_nb = [1 3 5 2 4 6]; % in this function, the channels are numbered from left to right, top and then bottom (in the order  0°z1 ; 45°z2 ; 0°z3 ; 90°z1 ;  135°z2 ; 90°z32)
% however in the PolMFM we number them in the order: 0°z1 ; 90°z1 ; 45°z2 ; 135°z2 ; 0°z3 ; 90°z32
% The array PolMFM_ch_nb makes the correspondance between these two
% numbering when we display the results we keep the PolMFM numbering and
% the user is not lost
%% first step Z_calibration
info = imfinfo(filename);
num_images = numel(info); % size of the z-stack
% 
% imageMIP = imread(filename, 1);
% for ii=2:num_images
%     imageMIP = imageMIP+imread(filename, ii);
% end
% imageMIP = imageMIP/num_images;

% Calculate the mean value for the nine slices
x=0:1:num_images-1;
x=x*z_step;


avmaxint=zeros(N_channels,1);
avZpos=zeros(N_channels/2,1);
%%
for k=1:length(coord_final)  %loop over all the beads %input points  containes the (x,y) coordinates for each spot in each plane
    Zinfo{k}=zeros(num_images,N_channels);
    for j=1:num_images  % loop over the whole z stack
        B = imread(filename, j);
        for i=1:N_channels % loop over all the channels 
            % define a patch of size ws around each (x,y)coordinate
            %       in each plane(i=1:nb) at a
            %       image(j=1:num_images)
            patch2= B(floor(coord_final{k}(i,2)-ws):floor(coord_final{k}(i,2)+ws),floor(coord_final{k}(i,1)-ws):floor(coord_final{k}(i,1)+ws)); % window around the bead/slice/channel
            Zinfo{k}(j,i)=mean2(patch2);%% mean intensity of each patch
        end
    end
    
    % Gaussian fitting
    
    for i=1:N_channels %loop over all the channels
        [g h]=max(Zinfo{k}(:,i));
        c=min(Zinfo{k}(:,i));
        %define the Gaussian 1D function
        myfun = @(A) A(1)*exp(-(((x(:)-A(2)).^2/(2*A(3)^2)))) + A(4) -Zinfo{k}(:,i); % gaussian 1D function: remove the mean intensity of each plane 
        % Initial guess parameters
        A0=[g;x(h);300;c];
        options = optimset('Algorithm',{'levenberg-marquardt',0.005},'Display', 'off');
        fitparams2{k}(i,:) = (lsqnonlin(myfun,A0,[],[],options))';
        gaussianvalues=fitparams2{k}(i,1)*exp(-(((x(:)-fitparams2{k}(i,2)).^2/(2*fitparams2{k}(i,3)^2)))) + fitparams2{k}(i,4);
        figure(4556+k);
        subplot(3,6,i)
        plot(x,Zinfo{k}(:,i),x,gaussianvalues)
        grid on, grid minor
%         ylim([0 4500])
        title(['Channel #' num2str(PolMFM_ch_nb(i))])
        figure(45323),subplot(3,6,i),hold on
        plot(x,Zinfo{k}(:,i),x,gaussianvalues), xlabel('Z (nm)'), ylabel('Intensity')
        grid on, grid minor
%         ylim([0 4500])
        title(['Channel #' num2str(PolMFM_ch_nb(i))])
        
    end
    
%     f=fitparams2{k};
%     [f,index] = sortrows(f,2,'ascend');%% to be sure that the planes are correctly ordered
%     fitparams2{k}= f;
    f = fitparams2{k};
    save([filename(1:end-4) 'point' num2str(k) '.mat'], 'f')%amplitude - z position - width - offset
    
    Z{k}= [fitparams2{k}(1,2)+fitparams2{k}(4,2); fitparams2{k}(2,2)+fitparams2{k}(5,2); fitparams2{k}(3,2)+fitparams2{k}(6,2)]/2;
   
    figure(4556+k);
    subplot(3,6,[7 15])
    bar(Z{k})
%     ylim([0 2600])
    xlabel('Z-plane')
    ylabel('Z [nm]')
    grid on, grid minor
    
    avmaxint=avmaxint+fitparams2{k}(:,1);
    avZpos=avZpos+Z{k};
    
    X=1:N_channels/2;
    [a,b]=fit(X',Z{k},'poly1');
    a1(k)=a.p1;
    Y=a.p1*X+a.p2;
        subplot(3,6,[10 18])
    plot(X,Z{k},'--rs',X,Y,'b')
%     ylim([1400 2600])
    xlabel('Z-plane'); ylabel('Z [nm]'), xlim([0.9 3.1])
    legend(['Z_p_o_s'],['linear fit'],'Location','northwest')
    grid on, grid minor
    title(['Bead #' num2str(k)])
%     coord_final{k}=coord_final{k}(index,:); % Reorder the coordinates in the sense of Zposition
    coord_final{k}(:,3)=repmat(Z{k},2,1);
    sgtitle(['Bead #' num2str(k)])
    
    figure(45323),subplot(3,6,[10 18]), hold on
    plot(X,Z{k},'--s',X,Y,'-')
    sgtitle('All beads')
    legend()

end
  assignin('base','coordReorder',coord_final)

%%
for k=1:length(coord_final)
    display(['The average separation for point number ' num2str(k) ' is ' num2str(a1(k)) ' nm'])
end
delta_z_mean = mean(a1(:));
display (['the total average separation is ' num2str(delta_z_mean) ' nm'])
%% saving average max intensity and average z position into .mat file
avmaxint=avmaxint/length(coord_final);
avZpos=avZpos/length(coord_final);
imax=find(avmaxint==max(avmaxint));
correctionfactor=avmaxint/avmaxint(imax);
assignin('base','avZpos',avZpos)
assignin('base','avmaxint',avmaxint)
assignin('base','correctionfactor',correctionfactor)

if (param.Z_Calibration==1)
save([filename(1:end-4)  '_pos-int-cor.mat'], 'avZpos', 'avmaxint','correctionfactor')
%% saving average separation
fid = fopen([filename(1:end-4) '_Sep_estimation.txt'], 'w');
for k=1:length(coord_final)
    fprintf(fid, ' The average z seperation for point number %d is %5.3f nm\n',k, a1(k));
end
fprintf(fid, ' The mean z seperation is %5.3f nm\n',delta_z_mean);
fclose(fid);
end

global int_corr_pathname
int_corr_pathname=[filename(1:end-4)  '_pos-int-cor.mat'];



%% calculate transformation Matrix
% filename=Filename;
% nb=nb_planes;


numberofbeads=length(coord_final)*N_channels;

peaks=[];%zeros(numberofbeads,2);
% peaks(:,1)= coord(:,2);%4
% peaks(:,2)= coord(:,1);%2
for k=1:length(coord_final)
    %peaks=[peaks; coord{k}(:,:)];%4
    peaks=[peaks; coord_final{k}(:,2),coord_final{k}(:,1),coord_final{k}(:,3) ];
    %peaks(:,2)=[peaks(:,2), coord{k}(:,1)];%4
end

% figure
% scatter(peaks(:,2),peaks(:,1)),
%%
%decompose peaks in the 6 different channels
% yy = [1 min([param.Areas(2,1) param.Areas(5,1)]) min([param.Areas(3,1) param.Areas(6,1)])  size(maxprojection,1)];
% xx = [1 min([param.Areas(4,2) param.Areas(5,2) param.Areas(6,2)]) size(maxprojection,2)];


%% rearrange
numBeadsPerSlice=numberofbeads/N_channels;
peaksPerSlice=zeros(numBeadsPerSlice,2,N_channels);%contains the (x,y) coordinate of each peak with respect to corner reference (xMin,yMin)

cornerPerSlice=zeros(N_channels,4);%contains the bottom left corner (xMin,yMin) for each slice
    
  
    P.XYdata=[];
    P.Zdata=[];
    P.corners=[];
    P.MeanZ=[];
    %nb_planes=nb
    
    
    
    for s=1:N_channels
%         ii=comb(s,1);
%         jj=comb(s,2);
%         pos=find(peaks(:,1)>=xx(ii) & peaks(:,1)<xx(ii+1) & peaks(:,2)>=yy(jj) & peaks(:,2)<yy(jj+1));
        pos = find(peaks(:,1)>=param.Areas(s,2) & peaks(:,1)<param.Areas(s,2)+param.Areas(s,4) & peaks(:,2)>=param.Areas(s,1) & peaks(:,2)<param.Areas(s,1)+param.Areas(s,3));
        nn= min(length(pos),numBeadsPerSlice);
        assignin('base','numBeadsPerSlice',numBeadsPerSlice)
        %numBeadsPerSlice
        %length(pos)
        pos=pos(1:nn,:);
        if(length(pos)~=numBeadsPerSlice)
            pos(end+1:numBeadsPerSlice,:)=nan;
            error(['Slice ' num2str(s) ' does not contain the right number of beads']);
        end
%         peaksPerSlice(:,:,s)=peaks(pos,1:2)-repmat([xx(ii) yy(jj)],[numBeadsPerSlice 1]);
        peaksPerSlice(:,:,s)=peaks(pos,1:2)-repmat([param.Areas(s,2) param.Areas(s,1)],[numBeadsPerSlice 1]);
        [ss pp]=sortrows(peaksPerSlice(:,:,s));
        peaksPerSlice(:,:,s)=peaksPerSlice(pp,:,s);
        cornerPerSlice(s,:)=[param.Areas(s,2) param.Areas(s,1) param.Areas(s,2)+param.Areas(s,4)-1 param.Areas(s,1)+param.Areas(s,3)-1];
%         [ss pp]=sortrows(peaksPerSlice(:,:,s));
%         peaksPerSlice(:,:,s)=peaksPerSlice(pp,:,s);
%         cornerPerSlice(s,:)=[param.Areas(s,2) param.Areas(s,1) param.Areas(s,2)+param.Areas(s,4)-1 param.Areas(s,1)+param.Areas(s,3)-1];
        P.XYdata=[P.XYdata; {peaksPerSlice(:,:,s)}];
        P.Zdata=[P.Zdata;[peaks(pos,3)]'];
        P.corners=[P.corners;{cornerPerSlice(s,:)}];
        %P.MeanZ=[P.MeanZ, mean(P.Zdata{s})];
        P.MeanZ(s)=mean([peaks(pos,3)]);
        
    end

%align all the slices using teh Generalized Procrustes Algorithm (fancy
%term for an iterative method that finds translation, rotation and scaling)
%http://en.wikipedia.org/wiki/Procrustes_distance
tolD=1e-3;

%find mean location of points
Pref=mean(peaksPerSlice,3);

dTotal=1e6;
dTotalOld=dTotal*10;

iter=0;
while(abs(dTotalOld-dTotal)>tolD && iter<1000)
    dTotalOld=dTotal;
    %calculate transformation for each slice
    dTotal=0;
    PrefAux=0;
    for kk=1:N_channels
        [d, Z_transformed] = procrustes(Pref,peaksPerSlice(:,:,kk),'Reflection',false);
        dTotal=dTotal+sum(sqrt(sum((Pref-Z_transformed).^2,2)));
        PrefAux=PrefAux+Z_transformed;
    end
    Pref=PrefAux/N_channels;
    dTotal=dTotal/(N_channels*numBeadsPerSlice);
    %disp(['Iter=' num2str(iter) '.Average distance is ' num2str(dTotal) ])
    iter=iter+1;
end
disp(['Number of iteration is ' num2str(iter) ])
disp(['Average distance is ' num2str(dTotal) ])

%calculate transformation for each frame
transformPerSlice=cell(N_channels,1);
for kk=1:N_channels
    [d, Z, tr] = procrustes(Pref,peaksPerSlice(:,:,kk),'Reflection',false);
    transformPerSlice{kk}=tr;
end

%%
%debugging: check that everything aligns perfectly
dist_p_xy = zeros(N_channels,2);
figure;
zz=lines(N_channels);
plot(Pref(:,2),Pref(:,1),'ko');
for kk=1:N_channels
    pp=transformPerSlice{kk}.b * peaksPerSlice(:,:,kk) * transformPerSlice{kk}.T + transformPerSlice{kk}.c;
    hold on;plot(pp(:,2),pp(:,1),'+','Color',zz(kk,:));
    dist_p(:,kk)=sqrt((Pref(:,2)-pp(:,2)).^2+(Pref(:,1)-pp(:,1)).^2); %distances to the reference point
    dist_p_xy(kk,2)= mean((peaksPerSlice(:,2,kk)+param.Areas(kk,1))-Pref(:,2));
    dist_p_xy(kk,1)= mean((peaksPerSlice(:,1,kk)+param.Areas(kk,2))-Pref(:,1));
    list_p(:,:,kk)=pp; %list of tranformed points
end
hold off;
title 'All crosses should be superimpose one on top of the other around the black circle'
std_dist=std(dist_p,0,2);
mean_std=mean(std_dist);
% display('standard deviation of the beads allignement is ' )
% display(num2str(std_dist));
% display ('and there mean is ');
% display(num2str(mean_std))
center_position=mean(list_p,3);
std_position=std(list_p,0,3);
mean(std_position,1);

%--------------------

%%
%saves transformations
%disp(['Saving transformations at ' [pathIm basenameIm '.mat']])
%save([pathIm basenameIm '.mat'],'transformPerSlice','cornerPerSlice','peaksPerSlice');
if  (param.Generate_transformation_matrix==1)
uiwait(msgbox(['Saving transformation matrices between channels: ' filename(1:end-4) '.mat ------------'...
    'Saving channels tranlsation vectors from original image: ' filename(1:end-4) '_channel-shift.mat'], 'Saving', 'modal'));
% disp(['Saving transformations at ' [filename(1:end-4) '.mat']]);
save([filename(1:end-4) '.mat'],'transformPerSlice','cornerPerSlice','peaksPerSlice');
% disp(['Saving channel shift from original image  at ' [filename(1:end-4) '_channel-shift.mat']]);
save([filename(1:end-4) '_channel-shift.mat'],'dist_p_xy');


%% saving allignement error standard deviation results
save([filename(1:end-4) '_position.mat'],'list_p','dist_p','std_dist','mean_std');
end
%%
global trans_matrix_pathname
trans_matrix_pathname=[filename(1:end-4) '.mat'];





end