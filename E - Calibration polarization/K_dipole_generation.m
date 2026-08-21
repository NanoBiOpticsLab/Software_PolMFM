function [K_matrix] = K_dipole_generation(K_input,microscope_settings)

%% Input: experimental K
% /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\  Definition of channels 0-90 / 45-135 are
% not the same !???
K_exp = K_input; %K_exp(1,:)=K_input(2,:);K_exp(2,:)=K_input(1,:);K_exp(5,:)=K_input(6,:);K_exp(6,:)=K_input(5,:);
figure, subplot(221),imagesc(K_exp),clim([-1 1]),colormap jet,title('K_e_x_p'), axis equal, axis image

%% Ktheo = K theoretical matrix
n_2 = 1.33;                  % refractive index before objective: water
n_1 = microscope_settings.RI; %515; %1.33;                % refractive index after objective
n = n_1/n_2;
NA = microscope_settings.NA;                % numerical aperture
theta_1max = asin(NA/n_1); % opening angle of objective
theta_1c = asin(n_2/n_1); % critical angle
factor = 1 ;% aperture reduction: here value is 1 because no aperture reduction (2d orientation) 
theta2 = theta_1c; %theta_1max;  %theta_1c when not accounting for SAF
theta1 = asin(factor*sin(theta2)); % reduced channel
z_dipole = 0; %height of dipole over sample plane
wave_len = microscope_settings.Em*1e-3;       %emission wavelength in [um]

% General parameters
cos2 = @(p1,t1) sqrt(1-(n_1/n_2.*sin(t1)).^2);
sin_2phi1 = @(p1,t1) sin(2*p1);
rho = @(p1,t1) sin(t1);
t_p = @(p1,t1) 2*n_2*cos2(p1,t1)./(n_2.*cos(t1)+n_1.*cos2(p1,t1));             % Fresnel transmission coefficient for p-polar
t_s = @(p1,t1) 2*n_2*cos2(p1,t1)./(n_2.*cos2(p1,t1)+n_1.*cos(t1));             % Fresnel transmission coefficient for s-polar
psy_depth = @(p1,t1) 2*pi*z_dipole*n_2/wave_len*cos2(p1,t1);

% Electric fields to integrate
E_x_mux = @(p1,t1) n* ((cos(t1)./cos2(p1,t1)).*t_s(p1,t1).*sin(p1).^2 + t_p(p1,t1).*cos(p1).^2.*sqrt(1-rho(p1,t1).^2)).* exp(1i*psy_depth(p1,t1));
E_x_muy = @(p1,t1) -n/2*sin_2phi1(p1,t1).* ((cos(t1)./cos2(p1,t1)).*t_s(p1,t1) - t_p(p1,t1).*sqrt(1-rho(p1,t1).^2)).* exp(1i*psy_depth(p1,t1));
E_x_muz = @(p1,t1) -n^2*(cos(t1)./cos2(p1,t1)).*t_p(p1,t1).*rho(p1,t1).*cos(p1).* exp(1i*psy_depth(p1,t1));

E_y_mux = @(p1,t1) -n/2*sin_2phi1(p1,t1) .* ( (cos(t1)./cos2(p1,t1)).*t_s(p1,t1) - t_p(p1,t1).*sqrt(1-rho(p1,t1).^2) ).* exp(1i*psy_depth(p1,t1));
E_y_muy = @(p1,t1) n* ((cos(t1)./cos2(p1,t1)).*t_s(p1,t1).*cos(p1).^2 + t_p(p1,t1).*sin(p1).^2.*sqrt(1-rho(p1,t1).^2)).* exp(1i*psy_depth(p1,t1));
E_y_muz = @(p1,t1) -n^2*(cos(t1)./cos2(p1,t1)).*t_p(p1,t1).*rho(p1,t1).*sin(p1).* exp(1i*psy_depth(p1,t1));

E_45_mux = @(p1,t1) (E_x_mux(p1,t1)+E_y_mux(p1,t1))/sqrt(2);
E_45_muy = @(p1,t1) (E_x_muy(p1,t1)+E_y_muy(p1,t1))/sqrt(2);
E_45_muz = @(p1,t1) (E_x_muz(p1,t1)+E_y_muz(p1,t1))/sqrt(2);

E_135_mux = @(p1,t1) (-E_x_mux(p1,t1)+E_y_mux(p1,t1))/sqrt(2);
E_135_muy = @(p1,t1) (-E_x_muy(p1,t1)+E_y_muy(p1,t1))/sqrt(2);
E_135_muz = @(p1,t1) (-E_x_muz(p1,t1)+E_y_muz(p1,t1))/sqrt(2);

% Factors of the matrix
XXx = @(p1,t1) abs(E_x_mux(p1,t1)).^2.*sin(t1);
XXy = @(p1,t1) abs(E_y_mux(p1,t1)).^2.*sin(t1);

YYx = @(p1,t1) abs(E_x_muy(p1,t1)).^2.*sin(t1);
YYy = @(p1,t1) abs(E_y_muy(p1,t1)).^2.*sin(t1);

ZZx = @(p1,t1) abs(E_x_muz(p1,t1)).^2.*sin(t1);
ZZy = @(p1,t1) abs(E_y_muz(p1,t1)).^2.*sin(t1);

XX45 = @(p1,t1) abs(E_45_mux(p1,t1)).^2.*sin(t1);
XX135 = @(p1,t1) abs(E_135_mux(p1,t1)).^2.*sin(t1);

YY45 = @(p1,t1) abs(E_45_muy(p1,t1)).^2.*sin(t1);
YY135 = @(p1,t1) abs(E_135_muy(p1,t1)).^2.*sin(t1);

ZZ45 = @(p1,t1) abs(E_45_muz(p1,t1)).^2.*sin(t1);
ZZ135 = @(p1,t1) abs(E_135_muz(p1,t1)).^2.*sin(t1);

XYx = @(p1,t1) 2*real(conj(E_x_mux(p1,t1)) .* E_x_muy(p1,t1)).*sin(t1);
XYy = @(p1,t1) 2*real(conj(E_y_mux(p1,t1)) .* E_y_muy(p1,t1)).*sin(t1);

XY45 = @(p1,t1) 2*real(conj(E_45_mux(p1,t1)) .* E_45_muy(p1,t1)).*sin(t1);
XY135 = @(p1,t1) 2*real(conj(E_135_mux(p1,t1)) .* E_135_muy(p1,t1)).*sin(t1);

XZx = @(p1,t1) 2*real(conj(E_x_mux(p1,t1)) .* E_x_muz(p1,t1)).*sin(t1);
XZy = @(p1,t1) 2*real(conj(E_y_mux(p1,t1)) .* E_y_muz(p1,t1)).*sin(t1);

XZ45 = @(p1,t1) 2*real(conj(E_45_mux(p1,t1)) .* E_45_muz(p1,t1)).*sin(t1);
XZ135 = @(p1,t1) 2*real(conj(E_135_mux(p1,t1)) .* E_135_muz(p1,t1)).*sin(t1);

YZx = @(p1,t1) 2*real(conj(E_x_muy(p1,t1)) .* E_x_muz(p1,t1)).*sin(t1);
YZy = @(p1,t1) 2*real(conj(E_y_muy(p1,t1)) .* E_y_muz(p1,t1)).*sin(t1);

YZ45 = @(p1,t1) 2*real(conj(E_45_muy(p1,t1)) .* E_45_muz(p1,t1)).*sin(t1);
YZ135 = @(p1,t1) 2*real(conj(E_135_muy(p1,t1)) .* E_135_muz(p1,t1)).*sin(t1);

% intensity integrals high NA channels
XX0 = integral2(XXx,0,2*pi,0,theta1,'method','iterated');% Ex(mux)*Ex(mux)
YY0 = integral2(YYx,0,2*pi,0,theta1,'method','iterated');% Ex(muy)*Ex(muy)
ZZ0 = integral2(ZZx,0,2*pi,0,theta1,'method','iterated');% Ex(muz)*Ex(muz)
XY0 = integral2(XYx,0,2*pi,0,theta1,'method','iterated');% Ex(mux)*Ex(muy)
XZ0 = integral2(XZx,0,2*pi,0,theta1,'method','iterated');
YZ0 = integral2(YZx,0,2*pi,0,theta1,'method','iterated');

XX90 = integral2(XXy,0,2*pi,0,theta1,'method','iterated');
YY90 = integral2(YYy,0,2*pi,0,theta1,'method','iterated');
ZZ90 = integral2(ZZy,0,2*pi,0,theta1,'method','iterated');
XY90 = integral2(XYy,0,2*pi,0,theta1,'method','iterated');
XZ90 = integral2(XZy,0,2*pi,0,theta1,'method','iterated');
YZ90 = integral2(YZy,0,2*pi,0,theta1,'method','iterated');

% intensity integrals low NA channels
XX45 = integral2(XX45,0,2*pi,0,theta2,'method','iterated');
YY45 = integral2(YY45,0,2*pi,0,theta2,'method','iterated');
ZZ45 = integral2(ZZ45,0,2*pi,0,theta2,'method','iterated');
XY45 = integral2(XY45,0,2*pi,0,theta2,'method','iterated');
XZ45 = integral2(XZ45,0,2*pi,0,theta2,'method','iterated');
YZ45 = integral2(YZ45,0,2*pi,0,theta2,'method','iterated');

XX135 = integral2(XX135,0,2*pi,0,theta2,'method','iterated');
YY135 = integral2(YY135,0,2*pi,0,theta2,'method','iterated');
ZZ135 = integral2(ZZ135,0,2*pi,0,theta2,'method','iterated');
XY135 = integral2(XY135,0,2*pi,0,theta2,'method','iterated');
XZ135 = integral2(XZ135,0,2*pi,0,theta2,'method','iterated');
YZ135 = integral2(YZ135,0,2*pi,0,theta2,'method','iterated');

sumXX = 2*(XX0 + XX90) + XX45 + XX135; % this is needed to normalize the Z column in the experimental K's /!\ channels 0&90° are repeated twice


Ktheo = [XX0 YY0 XY0; XX90 YY90 XY90; XX45 YY45 XY45; XX135 YY135 XY135; XX0 YY0 XY0; XX90 YY90 XY90;]./sumXX;

subplot(222),imagesc(Ktheo), caxis([-1 1]); colormap jet, title('K_t_h_e_o'), axis equal, axis image

%% STEP 1 : intensity from a rotating dipole, from which we'll deduced the dipole-K matrix Kmu 
alphaStep=5;
alpha =  (0:alphaStep:180-alphaStep)*pi/180;
Nalpha = length(alpha); 
mux = cos(alpha); muy = sin(alpha); muz = 0;
% intensity from a dipole projected along 0° and 90° directions
Imu0 = XX0*mux.^2 + YY0*muy.^2 + ZZ0*muz.^2 + XY0*muy.*mux + XZ0*mux.*muz + YZ0*muy.*muz ;
Imu90 = XX90*mux.^2 + YY90*muy.^2 + ZZ90*muz.^2 + XY90*muy.*mux + XZ90*mux.*muz + YZ90*muy.*muz ;
Imu = Imu90 + Imu0;
Imu0 = Imu0./sum(Imu);
Imu90 = Imu90./sum(Imu);
% reduced channel along 45° and 135° directions
Imu45red = XX45*mux.^2 + YY45*muy.^2 + ZZ45*muz.^2 + XY45*muy.*mux + XZ45*mux.*muz + YZ45*muy.*muz ;
Imu135red = XX135*mux.^2 + YY135*muy.^2 + ZZ135*muz.^2 + XY135*muy.*mux + XZ135*mux.*muz + YZ135*muy.*muz ;
Imu45red = Imu45red./sum(Imu);
Imu135red = Imu135red./sum(Imu);

%matrix Kmu extracted from the mu rotation : just a verification that this leads to Ktheo
Itot = 2*(Imu0 + Imu90) + Imu45red + Imu135red; % channels 0&90° are repeated twice in PolMFM
%1st column XX : polarizer at 0° - idx 1 of alpha
XXmu0 = Imu0(1)/Itot(1); XXmu90 = Imu90(1)/Itot(1); XXmu45 = Imu45red(1)/Itot(1); XXmu135 = Imu135red(1)/Itot(1) ;
%2nd column : polarizer at 90° idx 19
YYmu0 = Imu0(19)/Itot(19); YYmu90 = Imu90(19)/Itot(19); YYmu45 = Imu45red(19)/Itot(19); YYmu135 = Imu135red(19)/Itot(19) ;
%4th column XY : polarizer at 45° idx 10
XYmu0 = 2*Imu0(10)/Itot(10)-XXmu0-YYmu0; XYmu90 = 2*Imu90(10)/Itot(10)-XXmu90-YYmu90;
XYmu45 = 2*Imu45red(10)/Itot(10)-XXmu45-YYmu45; XYmu135 = 2*Imu135red(10)/Itot(10)-XXmu135-YYmu135;
% not useful  polarizer at 135° idx 28 
% XYmu0_2 = 2*Imu0(28)/Itot(28)-XXmu0-YYmu0; XYmu90_2 = 2*Imu90(28)/Itot(28)-XXmu90-YYmu90;
% XYmu45_2 = 2*Imu45red(28)/Itot(28)-XXmu45-YYmu45; XYmu135_2 = 2*Imu135red(28)/Itot(28)-XXmu135-YYmu135;
% XYmu0=mean(XYmu0,XYmu0_2); XYmu90=mean(XYmu90,XYmu90_2);XYmu45=mean(XYmu45,XYmu45_2);XYmu135=mean(XYmu135,XYmu135_2);
% 3rd columne :  taken from theoretical column
ZZmu0 = ZZ0./sumXX ; ZZmu90 = ZZ90./sumXX ; ZZmu45 = ZZ45./sumXX ; ZZmu135 = ZZ135./sumXX ; 

Kmu = [XXmu0 YYmu0 XYmu0; XXmu90 YYmu90 XYmu90; XXmu45 YYmu45 XYmu45; XXmu135 YYmu135 XYmu135;XXmu0 YYmu0 XYmu0; XXmu90 YYmu90 XYmu90]
subplot(224),imagesc(Kmu),caxis([-1 1]),colormap jet,title('K_d_i_p_o_l_e'), axis equal, axis image

%% STEP 2 : intensity from rotating polarizer with a depolarizing sample (nanobeads, lamp)
% intensity from a depolarized emitter followed by a polarizer in the BFP
Ipol0 = (XX0 + YY0 + ZZ0)*cos(alpha).^2 ;
Ipol90 = (XX90 + YY90 + ZZ90)*sin(alpha).^2 ;
Ipol = Ipol0 + Ipol90;
Ipol0 = Ipol0./sum(Ipol);
Ipol90 = Ipol90./sum(Ipol);
% reduced channel, mimics the 45/135 deg channels
Ipol45red = (XX45 + YY45 + ZZ45)*cos(alpha+3*pi/4).^2 ;
Ipol135red = (XX135 + YY135 + ZZ135)*sin(alpha+3*pi/4).^2 ;
% Ipolred = Ipol45red + Ipol135red;
Ipol45red = Ipol45red./sum(Ipol);
Ipol135red = Ipol135red./sum(Ipol);

% measured matrix K from the polarizer rotation
Itotpol = 2*(Ipol0 + Ipol90) + Ipol45red + Ipol135red;  % channels 0&90° are repeated twice in PolMFM
%polarizer at 0° idx 1
XXpol0 = Ipol0(1)/Itotpol(1); XXpol90 = Ipol90(1)/Itotpol(1); XXpol45 = Ipol45red(1)/Itotpol(1); XXpol135 = Ipol135red(1)/Itotpol(1) ;
%polarizer at 90° idx 19
YYpol0 = Ipol0(19)/Itotpol(19); YYpol90 = Ipol90(19)/Itotpol(19); YYpol45 = Ipol45red(19)/Itotpol(19); YYpol135 = Ipol135red(19)/Itotpol(19) ;
%polarizer at 45° idx 10, 135° idx 28 
XYpol0 = 2*Ipol0(10)/Itotpol(10)-XXpol0-YYpol0; XYpol90 = 2*Ipol90(10)/Itotpol(10)-XXpol90-YYpol90;
XYpol45 = 2*Ipol45red(10)/Itotpol(10)-XXpol45-YYpol45; XYpol135 = 2*Ipol135red(10)/Itotpol(10)-XXpol135-YYpol135;
% theoretical column
ZZpol0 = ZZ0./sumXX ; ZZpol90 = ZZ90./sumXX ; ZZpol45 = ZZ45./sumXX ; ZZpol135 = ZZ135./sumXX ; 
%correction factor between z and x,y channels
% Kzz=Kzz/sum(Kzz); %normalized to 1 for simulated values
%  Kzz_factor=sum(K_calc(:,3))/sum(K_calc(:,1));% normalized Kzz to sum(Kzz)=1
%  Kzz=Kzz*Kzz_factor; %normalized to proportion z component
% Kzz/sum(Kzz)*sum(K_calc(:,3))/sum(K_calc(:,1))

Kpol = [XXpol0 YYpol0 XYpol0; XXpol90 YYpol90 XYpol90; XXpol45 YYpol45 XYpol45; XXpol135 YYpol135 XYpol135;XXpol0 YYpol0 XYpol0; XXpol90 YYpol90 XYpol90]
% sumXXpol = XXpol0 + XXpol90 + XXpol45 + XXpol135
% sumYYpol = YYpol0 + YYpol90 + YYpol45 + YYpol135
% sumZZpol = ZZpol0 + ZZpol90 + ZZpol45 + ZZpol135
% sumYYpol = XYpol0 + XYpol90 + XYpol45 + XYpol135
% fpol = sumZZpol/sumXXpol

 subplot(223),imagesc(Kpol),caxis([-1 1]),colormap jet,title('K_p_o_l'), axis equal, axis image

FigInt=figure('name', 'rotation dipole vs polarizer');
      plot(alpha*180/pi,Ipol0(:)', 'r', 'linewidth', 2); hold on ;  
      plot(alpha*180/pi,Ipol90(:)', 'b', 'linewidth', 2); hold on ;  
      plot(alpha*180/pi,Ipol45red(:)', 'g', 'linewidth', 2); hold on ;  
      plot(alpha*180/pi,Ipol135red(:)', 'c', 'linewidth', 2); hold on ;  
      plot(alpha*180/pi,Imu0(:)','--r', 'linewidth', 2); hold on ; 
      plot(alpha*180/pi,Imu90(:)','--b', 'linewidth', 2); hold on ;  
      plot(alpha*180/pi,Imu45red(:)','--g', 'linewidth', 2); hold on ; 
      plot(alpha*180/pi,Imu135red(:)','--c', 'linewidth', 2); hold on ;  
legend('Ipol0', 'Ipol90','Ipol45', 'Ipol135','Imu0','Imu90','Imu45','Imu135');
xlabel('\alpha','fontsize',10); ylabel('Int','fontsize',10);

%% STEP 3 : correcting factors to apply to K for mimicking a dipole rotation
% fit parameters of the responses, for both dipole rotation and polarizer rotation
    C2=cos(2*alpha);
    S2=sin(2*alpha);
    Int1 = Imu0; Int2 = Imu90; Int3 = Imu45red; Int4 = Imu135red;  Int5 = Imu0; Int6 = Imu90; 
    A0_1=mean(Int1,2);  A0_2=mean(Int2,2);  A0_3=mean(Int3,2);  A0_4=mean(Int4,2);  A0_5=mean(Int5,2);  A0_6=mean(Int6,2);
    A2_1=2*mean(Int1.*C2,2); A2_2=2*mean(Int2.*C2,2);A2_3=2*mean(Int3.*C2,2);A2_4=2*mean(Int4.*C2,2);A2_5=2*mean(Int5.*C2,2);A2_6=2*mean(Int6.*C2,2);
    B2_1=2*mean(Int1.*S2,2); B2_2=2*mean(Int2.*S2,2);B2_3=2*mean(Int3.*S2,2);B2_4=2*mean(Int4.*S2,2);B2_5=2*mean(Int5.*S2,2);B2_6=2*mean(Int6.*S2,2);
    A2n_1=A2_1/A0_1;A2n_2=A2_2/A0_2;A2n_3=A2_3/A0_3;A2n_4=A2_4/A0_4;A2n_5=A2_5/A0_5;A2n_6=A2_6/A0_6;
    B2n_1=B2_1/A0_1;B2n_2=B2_2/A0_2;B2n_3=B2_3/A0_3;B2n_4=B2_4/A0_4;B2n_5=B2_5/A0_5;B2n_6=B2_6/A0_6;

    Int1pol = Ipol0; Int2pol = Ipol90; Int3pol = Ipol45red; Int4pol = Ipol135red; Int5pol = Ipol0; Int6pol = Ipol90;
    A0_1pol=mean(Int1pol,2);  A0_2pol=mean(Int2pol,2);  A0_3pol=mean(Int3pol,2);  A0_4pol=mean(Int4pol,2);  A0_5pol=mean(Int5pol,2);  A0_6pol=mean(Int6pol,2); 
    A2_1pol=2*mean(Int1pol.*C2,2); A2_2pol=2*mean(Int2pol.*C2,2);A2_3pol=2*mean(Int3pol.*C2,2);A2_4pol=2*mean(Int4pol.*C2,2);A2_5pol=2*mean(Int5pol.*C2,2);A2_6pol=2*mean(Int6pol.*C2,2);
    B2_1pol=2*mean(Int1pol.*S2,2); B2_2pol=2*mean(Int2pol.*S2,2);B2_3pol=2*mean(Int3pol.*S2,2);B2_4pol=2*mean(Int4pol.*S2,2);B2_5pol=2*mean(Int5pol.*S2,2);B2_6pol=2*mean(Int6pol.*S2,2);
    A2n_1pol=A2_1pol/A0_1pol;A2n_2pol=A2_2pol/A0_2pol;A2n_3pol=A2_3pol/A0_3pol;A2n_4pol=A2_4pol/A0_4pol;A2n_5pol=A2_5pol/A0_5pol;A2n_6pol=A2_6pol/A0_6pol;
    B2n_1pol=B2_1pol/A0_1pol;B2n_2pol=B2_2pol/A0_2pol;B2n_3pol=B2_3pol/A0_3pol;B2n_4pol=B2_4pol/A0_4pol;B2n_5pol=B2_5pol/A0_5pol;B2n_6pol=B2_6pol/A0_6pol;
 
% transformation of Ipol into Imu : factors to apply to Ipol to make it look like Imu
%     corrA0_1 = A0_1/A0_1pol;     corrA0_2 = A0_2/A0_2pol;     corrA0_3 = A0_3/A0_3pol;    corrA0_4 = A0_4/A0_4pol;
%     corrA2_1 = A2_1/A2_1pol;     corrA2_2 = A2_2/A2_2pol;     corrA2_3 = A2_3/A2_3pol;    corrA2_4 = A2_4/A2_4pol;
%     corrB2_1 = B2_1/B2_1pol;     corrB2_2 = B2_2/B2_2pol;     corrB2_3 = B2_3/B2_3pol;    corrB2_4 = B2_4/B2_4pol;
% we noticed that there can be a problem of 0 division: therefore we add some epsilon term
    corrA0_1 = (A0_1+0.00001)/(A0_1pol+0.00001);     corrA0_2 = (A0_2+0.00001)/(A0_2pol+0.00001);     corrA0_3 = (A0_3+0.00001)/(A0_3pol+0.00001);    corrA0_4 = (A0_4+0.00001)/(A0_4pol+0.00001);     corrA0_5 = (A0_5+0.00001)/(A0_5pol+0.00001);    corrA0_6 = (A0_6+0.00001)/(A0_6pol+0.00001);
    corrA2_1 = (A2_1+0.00001)/(A2_1pol+0.00001);     corrA2_2 = (A2_2+0.00001)/(A2_2pol+0.00001);     corrA2_3 = (A2_3+0.00001)/(A2_3pol+0.00001);    corrA2_4 = (A2_4+0.00001)/(A2_4pol+0.00001);     corrA2_5 = (A2_5+0.00001)/(A2_5pol+0.00001);    corrA2_6 = (A2_6+0.00001)/(A2_6pol+0.00001);
    corrB2_1 = (B2_1+0.00001)/(B2_1pol+0.00001);     corrB2_2 = (B2_2+0.00001)/(B2_2pol+0.00001);     corrB2_3 = (B2_3+0.00001)/(B2_3pol+0.00001);    corrB2_4 = (B2_4+0.00001)/(B2_4pol+0.00001);     corrB2_5 = (B2_5+0.00001)/(B2_5pol+0.00001);    corrB2_6 = (B2_6+0.00001)/(B2_6pol+0.00001);
   
    %% STEP 4 : application of the correction to K_exp
alphaStep=5;
alpha =  (0:alphaStep:180-alphaStep)*pi/180;
Nalpha = length(alpha); 
Ipolexp0_1 = K_exp(1,1)*cos(alpha).^2 + K_exp(1,2)*sin(alpha).^2  + K_exp(1,3)*sin(alpha).*cos(alpha) ; %MFM plane 1 (z=dz)
Ipolexp90_1 = K_exp(2,1)*cos(alpha).^2 + K_exp(2,2)*sin(alpha).^2  + K_exp(2,3)*sin(alpha).*cos(alpha) ;%MFM plane 1 (z=dz)
Ipolexp45 = K_exp(3,1)*cos(alpha).^2 + K_exp(3,2)*sin(alpha).^2  + K_exp(3,3)*sin(alpha).*cos(alpha) ;%MFM plane 2 (z=0)
Ipolexp135 = K_exp(4,1)*cos(alpha).^2 + K_exp(4,2)*sin(alpha).^2 + K_exp(4,3)*sin(alpha).*cos(alpha) ;%MFM plane 2 (z=0)
Ipolexp0_3 = K_exp(5,1)*cos(alpha).^2 + K_exp(5,2)*sin(alpha).^2 + K_exp(5,3)*sin(alpha).*cos(alpha) ;%MFM plane 3 (z=+dz)
Ipolexp90_3 = K_exp(6,1)*cos(alpha).^2 + K_exp(6,2)*sin(alpha).^2  + K_exp(6,3)*sin(alpha).*cos(alpha) ;%MFM plane 3 (z=+dz)

IpolexpTOT = sum(Ipolexp0_1)+sum(Ipolexp90_1);
Ipolexp0_1 = Ipolexp0_1./IpolexpTOT;
Ipolexp90_1 = Ipolexp90_1./IpolexpTOT;
Ipolexp45 = Ipolexp45./IpolexpTOT;
Ipolexp135 = Ipolexp135./IpolexpTOT;
Ipolexp0_3 = Ipolexp0_3./IpolexpTOT;
Ipolexp90_3 = Ipolexp90_3./IpolexpTOT;

    Int1polexp = Ipolexp0_1; Int2polexp = Ipolexp90_1; Int3polexp = Ipolexp45; Int4polexp = Ipolexp135; Int5polexp = Ipolexp0_3; Int6polexp = Ipolexp90_3;
    A0_1polexp=mean(Int1polexp,2);  A0_2polexp=mean(Int2polexp,2);  A0_3polexp=mean(Int3polexp,2);  A0_4polexp=mean(Int4polexp,2);   A0_5polexp=mean(Int5polexp,2);  A0_6polexp=mean(Int6polexp,2);
    A2_1polexp=2*mean(Int1polexp.*C2,2); A2_2polexp=2*mean(Int2polexp.*C2,2);A2_3polexp=2*mean(Int3polexp.*C2,2);A2_4polexp=2*mean(Int4polexp.*C2,2);A2_5polexp=2*mean(Int5polexp.*C2,2);A2_6polexp=2*mean(Int6polexp.*C2,2);
    B2_1polexp=2*mean(Int1polexp.*S2,2); B2_2polexp=2*mean(Int2polexp.*S2,2);B2_3polexp=2*mean(Int3polexp.*S2,2);B2_4polexp=2*mean(Int4polexp.*S2,2);B2_5polexp=2*mean(Int5polexp.*S2,2);B2_6polexp=2*mean(Int6polexp.*S2,2);
    A2n_1polexp=A2_1polexp/A0_1polexp;A2n_2polexp=A2_2polexp/A0_2polexp;A2n_3polexp=A2_3polexp/A0_3polexp;A2n_4polexp=A2_4polexp/A0_4polexp;A2n_5polexp=A2_5polexp/A0_5polexp;A2n_6polexp=A2_6polexp/A0_6polexp;
    B2n_1polexp=B2_1polexp/A0_1polexp;B2n_2polexp=B2_2polexp/A0_2polexp;B2n_3polexp=B2_3polexp/A0_3polexp;B2n_4polexp=B2_4polexp/A0_4polexp;B2n_5polexp=B2_5polexp/A0_5polexp;B2n_6polexp=B2_6polexp/A0_6polexp;

    % application of the factors to Ipol
    A0_1exp = A0_1polexp*corrA0_1;     A0_2exp = A0_2polexp*corrA0_2;     A0_3exp = A0_3polexp*corrA0_3;     A0_4exp = A0_4polexp*corrA0_4;     A0_5exp = A0_5polexp*corrA0_5;     A0_6exp = A0_6polexp*corrA0_6; 
    A2_1exp = A2_1polexp*corrA2_1;     A2_2exp = A2_2polexp*corrA2_2;     A2_3exp = A2_3polexp*corrA2_3;     A2_4exp = A2_4polexp*corrA2_4;     A2_5exp = A2_5polexp*corrA2_5;     A2_6exp = A2_6polexp*corrA2_6; 
    B2_1exp = B2_1polexp*corrB2_1;     B2_2exp = B2_2polexp*corrB2_2;     B2_3exp = B2_3polexp*corrB2_3;     B2_4exp = B2_4polexp*corrB2_4;     B2_5exp = B2_5polexp*corrB2_5;     B2_6exp = B2_6polexp*corrB2_6; 
    % normalization of factors
    A2n_1exp=A2_1exp/A0_1exp;A2n_2exp=A2_2exp/A0_2exp;A2n_3exp=A2_3exp/A0_3exp;A2n_4exp=A2_4exp/A0_4exp;A2n_5exp=A2_5exp/A0_5exp;A2n_6exp=A2_6exp/A0_6exp;
    B2n_1exp=B2_1exp/A0_1exp;B2n_2exp=B2_2exp/A0_2exp;B2n_3exp=B2_3exp/A0_3exp;B2n_4exp=B2_4exp/A0_4exp;B2n_5exp=B2_5exp/A0_5exp;B2n_6exp=B2_6exp/A0_6exp; 

    % corrected intensity 
    IntFit_1=A0_1exp+A2_1exp*cos(2*alpha)+B2_1exp*sin(2*alpha);
    IntFit_2=A0_2exp+A2_2exp*cos(2*alpha)+B2_2exp*sin(2*alpha);
    IntFit_3=A0_3exp+A2_3exp*cos(2*alpha)+B2_3exp*sin(2*alpha);
    IntFit_4=A0_4exp+A2_4exp*cos(2*alpha)+B2_4exp*sin(2*alpha);
    IntFit_5=A0_5exp+A2_5exp*cos(2*alpha)+B2_5exp*sin(2*alpha);
    IntFit_6=A0_6exp+A2_6exp*cos(2*alpha)+B2_6exp*sin(2*alpha);

Figexpcorr =figure('name', 'exp corrected vs uncorrected');
       plot(alpha*180/pi,Ipolexp0_1(:)','r', 'linewidth', 2); hold on ; 
       plot(alpha*180/pi,Ipolexp90_1(:)','b', 'linewidth', 2); hold on ;  
       plot(alpha*180/pi,Ipolexp45(:)','g', 'linewidth', 2); hold on ; 
       plot(alpha*180/pi,Ipolexp135(:)','c', 'linewidth', 2); hold on ;  
       plot(alpha*180/pi,Ipolexp0_3(:)','m', 'linewidth', 1); hold on ; 
       plot(alpha*180/pi,Ipolexp90_3(:)','k', 'linewidth', 1); hold on ;  
       plot(alpha*180/pi,IntFit_1(:)', '--r', 'linewidth', 2); hold on ;  
       plot(alpha*180/pi,IntFit_2(:)', '--b', 'linewidth', 2); hold on ;  
       plot(alpha*180/pi,IntFit_3(:)', '--g', 'linewidth', 2); hold on ;  
       plot(alpha*180/pi,IntFit_4(:)', '--c', 'linewidth', 2); hold on ; 
       plot(alpha*180/pi,IntFit_5(:)', '--m', 'linewidth', 1); hold on ;  
       plot(alpha*180/pi,IntFit_6(:)', '--k', 'linewidth', 1); hold on ; 
legend('Ipolexp0_1', 'Ipolexp90_1','Ipolexp45','Ipolexp135','Ipolexp0_3', 'Ipolexp90_3','Icorrected0_1', 'Icorrected90_1','Icorrected45', 'Icorrected135','Icorrected0_3', 'Icorrected90_3');
xlabel('\alpha','fontsize',10); ylabel('Int','fontsize',10);


%% STEP 5 : deduced K_exp matrix from transformed intensity response IntFit
% measured matrix K from the polarizer rotation
Intfittot = IntFit_1 + IntFit_2 + IntFit_3 + IntFit_4 + IntFit_5 + IntFit_6 ;
%polarizer at 0° idx 1
XXpol0 = IntFit_1(1)/Intfittot(1); XXpol90 = IntFit_2(1)/Intfittot(1); XXpol45 = IntFit_3(1)/Intfittot(1); XXpol135 = IntFit_4(1)/Intfittot(1) ; XXpol0_3 = IntFit_5(1)/Intfittot(1); XXpol90_3 = IntFit_6(1)/Intfittot(1);
%polarizer at 90° idx 19
YYpol0 = IntFit_1(19)/Intfittot(19); YYpol90 = IntFit_2(19)/Intfittot(19); YYpol45 = IntFit_3(19)/Intfittot(19); YYpol135 = IntFit_4(19)/Intfittot(19) ;YYpol0_3 = IntFit_5(19)/Intfittot(19); YYpol90_3 = IntFit_6(19)/Intfittot(19); 
%polarizer at 45° idx 10, 135° idx 28 
XYpol0 = 2*IntFit_1(10)/Intfittot(10)-XXpol0-YYpol0; XYpol90 = 2*IntFit_2(10)/Intfittot(10)-XXpol90-YYpol90;
XYpol45 = 2*IntFit_3(10)/Intfittot(10)-XXpol45-YYpol45; XYpol135 = 2*IntFit_4(10)/Intfittot(10)-XXpol135-YYpol135;
XYpol0_3 = 2*IntFit_5(10)/Intfittot(10)-XXpol0_3-YYpol0_3; XYpol90_3 = 2*IntFit_6(10)/Intfittot(10)-XXpol90_3-YYpol90_3;
% theoretical column
% ZZpol0 = ZZ0./sumXX ; ZZpol90 = ZZ90./sumXX ; ZZpol45 = ZZ45./sumXX ; ZZpol135 = ZZ135./sumXX ; 

Kpolfit = [XXpol0 YYpol0 XYpol0; XXpol90 YYpol90 XYpol90; XXpol45 YYpol45 XYpol45; XXpol135 YYpol135 XYpol135;XXpol0_3 YYpol0_3 XYpol0_3; XXpol90_3 YYpol90_3 XYpol90_3]

figure, imagesc(Kpolfit),caxis([-1 1]),colormap jet,title('K_c_o_r_r'), axis equal, axis image

%% /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ put again K matrix in the order of the setup

K_matrix = Kpolfit;% K_matrix(1,:)=Kpolfit(2,:);K_matrix(2,:)=Kpolfit(1,:);K_matrix(5,:)=Kpolfit(6,:);K_matrix(6,:)=Kpolfit(5,:);

end