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
clear all
%% load results 
if exist('results_all','var') == 0
    %choose file
    [file_results,path_results] = uigetfile({'*.txt;*.csv'},'Select result file to display');
table_results=readtable(fullfile(path_results,file_results),'PreserveVariableNames',true);
results_all = table2array(table_results);
answer_method_orientation = questdlg('Method for orientation retrieval in the selected file?','please answer','Propagation matrix', 'Polarization factors', 'Cancel', 'Cancel');
end


%% Figure: scatter plot + color code for rho orientation
% figure(),
% scatter(results_all(:,2),results_all(:,3),[],(results_all(:,14)),'filled')
% axis equal, axis image
% xlabel('x [nm]'),ylabel('y [nm]'),caxis([0 180])
% title(['Results - Color = \rho - method = ' answer_method_orientation])

figure(), 
scatter3(results_all(:,2),results_all(:,3),results_all(:,4),5,(results_all(:,1)),'filled')
axis equal, axis image, 
xlabel('x [nm]'),ylabel('y [nm]'),zlabel('z [nm]'),
% caxis([0 180])
colormap jet
title(['Results 3D - Color = \rho - method = ' answer_method_orientation])

%% Figure: histogram of orientations 
size_bin = 5;
figure(),sgtitle(['Histogram of in plane angles \rho - method = ' answer_method_orientation])
subplot(2,2,[1 3])
histogram(results_all(:,14),0:size_bin:180)
xlabel('\rho [°]'),ylabel('# of occurences'), grid on, grid minor
title('\rho_a_l_l')
subplot(2,2,2)
histogram(results_all(:,15),0:size_bin:180)
xlabel('\rho [°]'),ylabel('# of occurences'), grid on, grid minor
title('\rho_-_1')
subplot(2,2,4)
histogram(results_all(:,16),0:size_bin:180), grid on, grid minor
xlabel('\rho [°]'),ylabel('# of occurences')
title('\rho_+_1')


figure(),sgtitle({['Histogram of wobbling angles \delta - method = ' answer_method_orientation],[num2str(size(results_all,1)) ' detected spots']})
subplot(2,2,[1 3])
histogram(results_all(:,11),0:size_bin:180)
xlabel('\delta [°]'),ylabel('# of occurences'), grid on, grid minor
title('\delta_a_v_g = mean(\delta_+_1 , \delta_-_1)')
subplot(2,2,2)
histogram(results_all(:,12),0:size_bin:180)
xlabel('\delta [°]'),ylabel('# of occurences'), grid on, grid minor
title('\delta_-_1')
subplot(2,2,4)
histogram(results_all(:,13),0:size_bin:180), grid on, grid minor
xlabel('\delta [°]'),ylabel('# of occurences')
title('\delta_+_1')

%% Choose ROI & plot histograms

figure(), 
sgtitle(['Orientation analysis in ROI - method = ' answer_method_orientation])
subplot(221)
scatter(results_all(:,2),results_all(:,3),5,(results_all(:,14)),'filled')
axis equal, axis image, 
xlabel('x [nm]'),ylabel('y [nm]'),caxis([0 180]),colormap jet
title(['Results 2D - Color = \rho - method = ' answer_method_orientation]), colorbar

h=impoly(gca)
position=wait(h);
position=[position; squeeze(position(1,1)), squeeze(position (1,2))];
in=inpolygon(results_all(:,2),results_all(:,3),position(:,1),position(:,2));
results_ROI=results_all(in,:);

hold on 
scatter(results_ROI(:,2),results_ROI(:,3),20,(results_ROI(:,14)),'filled'), colormap jet
plot(position(:,1),position(:,2))
hold off
clear h position in

subplot(222)
histogram(results_ROI(:,14),0:size_bin:180)
xlabel('\rho [°]'),ylabel('# of occurences'), grid on, grid minor
title('\rho_a_v_g = mean(\rho_+_1 , \rho_-_1)')
subplot(2,2,3)
histogram(results_ROI(:,15),0:size_bin:180)
xlabel('\rho [°]'),ylabel('# of occurences'), grid on, grid minor
title('\rho_-_1')
subplot(2,2,4)
histogram(results_ROI(:,16),0:size_bin:180), grid on, grid minor
xlabel('\rho [°]'),ylabel('# of occurences')
title('\rho_+_1')


%%  Choose rho & delta from z-loc 
thr_z0 = 200; % at z>thr_z0 or z<-thr_z0: choose orientation result from single plane 
results_all_b = zeros(size(results_all,1),size(results_all,2)+2);
results_all_b(:,1:size(results_all,2)) = results_all; 
for ii = 1:size(results_all,1)
    if results_all(ii,4)>thr_z0
        results_all_b(ii,18)=results_all(ii,16);
        results_all_b(ii,17)=results_all(ii,13);
    elseif results_all(ii,4)<=-thr_z0
        results_all_b(ii,18)=results_all(ii,15);
        results_all_b(ii,17)=results_all(ii,12);
    else
        results_all_b(ii,17)=results_all(ii,14);
        results_all_b(ii,18)=results_all(ii,11);
    end
end
%% save z loc discrimination
VarNames = table_results.Properties.VariableNames;
VarNames{end+1} = 'delta_z'; VarNames{end+1} = 'rho_z';

[filename, pathname]=uiputfile({'*.txt' },'Save merged data',path_results);
file_name_full=fullfile(pathname,[filename(1:end-4),'.txt']);
fid=fopen(file_name_full,'w');
display('Saving data .txt')
[r,s]=size(results_all_b);
fprintf(fid,['%g' repmat(['\t %g'],1,s-1) '\n'],results_all_b');
fclose(fid)
display('done')
file_name_full=fullfile(pathname,[filename(1:end-4),'.csv']);
display('Saving data .csv')
results_all_b_table = array2table(results_all_b,'VariableNames',VarNames);
writetable(results_all_b_table,file_name_full);
display('done')


%% plot histograms 

figure(),sgtitle(['Histogram of in plane angles \rho - choosen with z - method = ' answer_method_orientation])
subplot(2,2,[1 3])
histogram(results_all_b(:,18),0:2:180)
xlabel('\rho [°]'),ylabel('# of occurences'), grid on, grid minor
title('\rho_a_v_g = (\rho_+_1 or \rho_-_1)')
subplot(2,2,2)
histogram(results_all_b(:,15),0:2:180)
xlabel('\rho [°]'),ylabel('# of occurences'), grid on, grid minor
title('\rho_-_1')
subplot(2,2,4)
histogram(results_all_b(:,16),0:2:180), grid on, grid minor
xlabel('\rho [°]'),ylabel('# of occurences')
title('\rho_+_1')


figure(),sgtitle({['Histogram of wobbling angles \delta - choosen with z - method = ' answer_method_orientation],[num2str(size(results_all_b,1)) ' detected spots']})
subplot(2,2,[1 3])
histogram(results_all_b(:,17),0:5:180)
xlabel('\delta [°]'),ylabel('# of occurences'), grid on, grid minor
title('\delta_a_v_g = (\delta_+_1 or \delta_-_1)')
subplot(2,2,2)
histogram(results_all_b(:,12),0:5:180)
xlabel('\delta [°]'),ylabel('# of occurences'), grid on, grid minor
title('\delta_-_1')
subplot(2,2,4)
histogram(results_all_b(:,13),0:5:180), grid on, grid minor
xlabel('\delta [°]'),ylabel('# of occurences')
title('\delta_+_1')
