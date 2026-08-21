clear all
close all

% written by Bassam Hajj - 07/11/2023
% x , y and rho are the positions and angles in degrees
%%
rho=0:0.1:2*pi;
x=3*cos(rho); 
y=3*sin(rho);
rho=rho*180/pi;%to degrees

figure,plot(x,y,'o')
axis equal

%% Choose file
[file,path] = uigetfile({'*.csv'},'Select the localization file','F:\Louise\DATA\DATA2023');
filename=fullfile(path,file);
%% Load data   x , y and rho are the positions and angles in degrees

D=importdata(filename);
preamp=1; EM=1;
% preamp=15.4; EM=40;
px_size_x = 120; px_size_z = 376; %nm

N_photons = D.data(:,5).*D.data(:,8).*D.data(:,9).*D.data(:,10)/(px_size_x*px_size_x*px_size_z)*(2*pi)^(3/2)*(preamp/EM);% filtering
filter_on = 1;
if filter_on==1
filter_threshold = D.data(:,8)<70 | D.data(:,8)>180 ...
    |  D.data(:,10)>600 | D.data(:,10)<150 ...
    | D.data(:,5)>1e10 | D.data(:,5)<0 ...
     | D.data(:,11)>180 | D.data(:,11)<0 ...
     | D.data(:,2)>1e9 | D.data(:,2)<0 ...
     | D.data(:,3)>1e9 | D.data(:,3)<0 ...
     | N_photons>10000 | N_photons<500 ...
     ;
D.data(filter_threshold,:) = [];
display(['filer ' num2str(sum(filter_threshold)) ' points'])
end

index = isnan(D.data(:,14)); %find NaN data in the orientation and remove the corresponding rows
display(['removing ' num2str(sum(index)) ' points with NAN orientation values'])
D.data(index,:)=[];
x=D.data(:,2);
y=D.data(:,3);
z=D.data(:,4);
%%
figure, subplot(121),histogram(N_photons((~filter_threshold)),'BinWidth',250,'Normalization','probability'), title('N_p_h - thresholded')
subplot(122), histogram(N_photons,0:500:50000),title('N_p_h - all')
% xlim([0 1e5])

%% visual parameters
l_segment=200; %length of the segment in nm
nb=1000; %length(x)  %number of pts to show 
if nb>size(D.data,1)
    nb=size(D.data,1);
end
rho_0=+60; % angle offset in degrees to correctly define the angle compared to the visual representation of the data 
size_points=15;
rho=D.data(:,14);
rho=rho+rho_0;

dx=cosd(rho)*l_segment/2;
dy=sind(rho)*l_segment/2;

u1=x-dx;%lower x limit of the segment
u2=x+dx; %higher x limit of the segment

v1=y-dy;%lower y limit of the segment
v2=y+dy; %higher y limit of the segment

% colormap(jet(256));
figure();
%ax=figure
ax1 = axes;
% %define the color code
cmap = hsv(255); 
ax1.Colormap = cmap;
%Orientation_normalized=rho/pi; %normalize the orientation between 0 and pi
Orientation_normalized = (rem(rho,180))/180;%normalize the orientation between 0 and pi
vecColorIdx = round(Orientation_normalized * (size(cmap,1)-1)) + 1;

display(['ploting'])


for i=1:nb
    plot([u1(i),u2(i)],[v1(i),v2(i)],'Color', cmap(vecColorIdx(i),:))
    hold on 
end
axis equal

% view(2)
ax2 = axes;
scatter(x(1:nb),y(1:nb),size_points,z(1:nb),'filled')
axis equal

linkaxes([ax1,ax2])
ax2.Visible = 'off';


hold off

% colorbar (ax1)
% colorbar (ax2)
%%Give each one its own colormap
colormap(ax1,'hsv')
colormap(ax2,'cool')
%%Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.17 .11 .685 .815]);
cb1 = colorbar(ax1,'Position',[.05 .11 .0675 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
display(['Done'])
