clear all
close all
%%

calib_0 = 'F:\Louise\DATA\DATA2024\2024-01-11_SilicaBeads_Lipids_NR\Calibration\Calib-Polar_5deg-steps_filter617nm\propagation_matrix_RevertPolarRot_corr.mat';
date_0 = '2024-01-11';
% load(calib_0);
% clear figure(53406) 
% figure(53406)
% for ii = 1:3
%     for jj= 1:6
%         hold on, scatter(jj+(ii-1)*6,K{1,1}(jj,ii),15,'filled','r')
%     end
% end

calib_1 = 'F:\Louise\DATA\DATA2024\2024-01-16_CalibrationPolar_Whitelight-vs-TL\TL3D_CalibPolar_5degstep_filter617nm\propagation_matrix_RevertPolarRot_corr.mat';
date_1 = '2024-01-16';
% load(calib_1);
% figure(53406)
% for ii = 1:3
%     for jj= 1:6
%         hold on, scatter(jj+(ii-1)*6,K{1,1}(jj,ii),15,'filled','g')
%     end
% end

calib_2 = 'F:\Louise\DATA\DATA2024\2024-01-18_Polarized-Beads_PolMFM\Calibration-Polar_TL3D_channels617nm-676nm\propagation_matrix_RevertPolarRot_corr.mat';
date_2 = '2024-01-18';
% load(calib_2);
% figure(453054), subplot(243),imagesc(abs(K{1,1})), colormap jet, clim([-0 0.35]), title(date_2)
% figure(53406)
% for ii = 1:3
%     for jj= 1:6
%         hold on, scatter(jj+(ii-1)*6,K{1,1}(jj,ii),15,'filled','b')
%     end
% end

calib_3 = 'F:\Louise\DATA\DATA2024\2024-01-29_SilicaBeads_Lipids_NR\Calibration\calibPolar_5degsteps_TL3d-617nmfilter\propagation_matrix_RevertPolarRot.mat';
date_3 = '2024-01-29';
% load(calib_3);
% figure(53406)
% for ii = 1:3
%     for jj= 1:6
%         hold on, scatter(jj+(ii-1)*6,K{1,1}(jj,ii),15,'filled','c')
%     end
% end

calib_4 = 'F:\Louise\DATA\DATA2024\2024-02-09_RPE1_DNA_SytoxOrange\Calibration_Polar_5degsteps_TL3D_617nm\propagation_matrix_RevertPolarRot.mat';
date_4 = '2024-02-09';
% load(calib_4);
% figure(53406)
% for ii = 1:3
%     for jj= 1:6
%         hold on, scatter(jj+(ii-1)*6,K{1,1}(jj,ii),15,'filled','m')
%     end
% end

calib_5 = 'F:\Louise\DATA\DATA2024\2024-02-19_CricketTestis_SytoxOrange_PolMFM\Calibration\Polarization_5degsteps_617nmfilter\propagation_matrix_revertPolarRot.mat';
date_5 = '2024-02-19';
% load(calib_5);
% figure(453054), subplot(246),imagesc(abs(K{1,1})), colormap jet, clim([0 0.35]), title(date_5)
% figure(53406)
% for ii = 1:3
%     for jj= 1:6
%         hold on, scatter(jj+(ii-1)*6,K{1,1}(jj,ii),15,'filled','y')
%     end
% end

calib_6 = 'F:\Louise\DATA\DATA2024\2024-02-27_CricketTestes_SytoxOrange_PolMFM\Results\Calibration\CalibPolar_filter617nm_5degsteps_TL3D\propagation_matrix_RevertPolarRot.mat';
date_6 = '2024-02-27';
% load(calib_6);
% figure(53406)
% for ii = 1:3
%     for jj= 1:6
%         hold on, scatter(jj+(ii-1)*6,K{1,1}(jj,ii),15,'filled','k')
%     end
% end

list_calib = {calib_0,calib_1,calib_2,calib_3,calib_4,calib_5,calib_6 }; 
list_date = {date_0,date_1,date_2,date_3,date_4,date_5,date_6};
K_all = zeros(6,3,size(list_calib,2));
for c = 1:size(list_calib,2)
    load(list_calib{c});
    figure(453054), subplot(2,4,c), 
    imagesc(abs(K{1,1})), colormap jet, clim([-0 0.35]), title(list_date{c})
    K_all(:,:,c) = K{1,1}; 
end
%%
figure()
for ii = 1:3
    for jj= 1:6
        hold on, errorbar((jj+(ii-1)*6),mean(squeeze(K_all(jj,ii,:))),std(squeeze(K_all(jj,ii,:))))
    end
end
%%
c = colormap("cool"); nc = size(c,1); 
clr_map = c(1:round(nc/size(list_calib,2)):nc,:);
K_temp = permute(K_all,[2 1 3]);
K_all_array = reshape(K_temp,[18 7]);
figure, 
boxchart(K_all_array' - mean(K_all_array',1))
xlabel('calibration matrix index')
ylabel('Coefficient - mean(coeff)')
figure, s = scatter(1:18,K_all_array-mean(K_all_array,2),40);
for c = 1:size(list_calib,2)
    s(c).MarkerEdgeColor = clr_map(c,:);
end
hold on, errorbar(1:18,zeros(1,18),std(K_all_array,[],2),'.','Color',[0.5 0.5 0.5])
xlabel('calibration matrix index')
ylabel('Coefficient - mean(Coeff)')
legend(list_date)
