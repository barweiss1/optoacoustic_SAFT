% --------------- Comments ---------------
% This test adds noise to the sinogram and looks at SAFT performance under
% this condition

% This test file runs 2D scans on the image and tests SAFT performace on 2D
% scans. 
% The SNR and FWHM are tested for all SAFT methods

% --------------- Code ---------------

%% clear Workspace
clear all
close all
clc

%% add paths to required folders

addpath('SIR_Functions');
addpath('Create_Image_Patterns');
addpath('SAFT_Methods');


%% parameters definition

% Parameter definition. All units are SI. The focus of the transducer is
% always at the center of the image.

% These are the parameters from the paper - "Improved optoacoustic microscopy through
% three-dimensional spatial impulse response synthetic aperture focusing
% technique", Jake Turner et al to recreate Fig. 2
D = 24.6e-3; % sensor diameter - NA = 45 degrees
F = 12.3e-3; % focal length 
c = 1480 ; % speed of sound 
time_res = 5 ; % time to depth sample factor (Nt = Nz*time_res)
space_res = 25e-6; % resolution in z-axis for smooth time derivative of the SIR
image_width =8.4e-6*c; % width that gives 2.8us of depth in time
image_height =4.48e-6*c; % height that gives 2.8us of depth in time

Nz = 750; % Number of pixels in the horizontal direction. 
Ny = 400; % Number of pixels in the vertical diraction
detector_num = Ny; % number of detector positions


line_index=Ny/2-1; % index for cross section plot

% create the imaged object pattern 

% % create points image
% I=zeros(Ny,Nz);
% width=1;
% height=1;
% x_center=Nz/2-floor(width/2);
% y_center = Ny/2;
% space=floor(Nz/15);
% rect_num=13;
% I=create_rect_image(I,rect_num,space,x_center,y_center,width,height);
% points_pos=x_center-space*((rect_num-1)/2):space:x_center+space*((rect_num-1)/2);

% create angled line image
I=zeros(Ny,Nz);
width=1;
height=200;
line_num=5;
x_center=Nz/2-floor(width/2);
y_center=Ny/2;
angle=20;
space=floor(Nz/5);
I = create_line_image(I,line_num,space,x_center,y_center,height,angle);
points_pos=x_center-space*((line_num-1)/2):space:x_center+space*((line_num-1)/2);


% I = create_circle_image(I,rect_num,space,x_center,y_center,radius);
% radius = 15;




%% plot original image

figure(1);
imagesc(I);
colormap gray
colorbar
impixelinfo;
title('original image');
xlabel('z');
ylabel('y');

% create time vectors
[t,y,z] = create_tyz_vectors(F,c,time_res,image_width,image_height,Nz,Ny);
t_plot=t-F/c;
dt = t(2)-t(1);

%% Sinogram - Long Run Time

% create sinogram 
[s,detector_pos]=sinogram(I,D,F,c,time_res,space_res,image_width,image_height,detector_num,Nz,Ny);

%% calculate SIR (SAFT weighting map)

% the SIR should be calculated at the same dimensions of the sinogram scan 
sir = Sph_SIR_map_wrapper(D, F, c, t, detector_pos(detector_num/2+1:end), c*t_plot ,space_res);

%% Plot SIR map

figure(24);
imagesc(t_plot,detector_pos,sir);
title('SIR for SAFT-SIR');
xlabel('t[sec]');
ylabel('y[m]');

%% Add Gaussian noise 

sigma_per=15; % std of the noise in [%] from maximum of sinogram
max_s=max(max(s));
s_noise=s+(randn(size(s))*max_s*(sigma_per/100)); % adding gaussian noise 


%% plot sinogram and noisy sinogram 

dz = z(2) - z(1);
point_distance = space*dz;

figure(2);
subplot(2,1,1);
imagesc(t_plot(1:end),detector_pos,s_noise);
% ylim([-0.6e-3 0.6e-3]);
% xlim([-0.8e-6 0.8e-6]);
colormap bone
colorbar
xlabel('t[sec]','FontSize',16);
ylabel('detector position[m]','FontSize',16);
mytitle = ['noisy sinogram \sigma=', num2str(sigma_per), '% of maximum'];
title(mytitle,'Interpreter','tex','FontSize',20);
subplot(2,1,2);
plot(t_plot(1:end-1),s_noise(line_index,:));
xlabel('t[sec]','FontSize',16);
ylabel('amplitude','FontSize',16);
% xlim([-0.8e-6 0.8e-6]);
title('Center Cross Section','FontSize',20);


figure(3);
subplot(2,1,1);
imagesc(t_plot(1:end),detector_pos,s);
% ylim([-0.6e-3 0.6e-3]);
% xlim([-0.8e-6 0.8e-6]);
colormap bone
colorbar
xlabel('t[sec]','FontSize',16);
ylabel('detector position[m]','FontSize',16);
title('Clean Sinogram','FontSize',20);
subplot(2,1,2);
plot(t_plot(1:end-1),s(line_index,:));
xlabel('t[sec]','FontSize',16);
ylabel('amplitude','FontSize',16);
% xlim([-0.8e-6 0.8e-6]);
title('Center Cross Section','FontSize',20);

% SNR calculation
s_SNR=sinogram_SNR(s_noise(line_index,:),points_pos,space,time_res);
% this section is for debugging the singoram_SNR function
% sigma=max_s*(sigma_per/100);
% [pks,locs]=findpeaks(s(Nz/2-1,:));
% real_SNR=mag2db(pks/sigma);
% % plot SNR
% figure(18);
% subplot(2,1,1);
% plot(points_pos*dz,s_SNR);
% hold on
% plot(points_pos*dz,real_SNR);
% xlabel('position[m]');
% ylabel('SNR[dB]');
% title('SNR as a function of depth');
% legend('calculated SNR','actual SNR');
% hold off
% subplot(2,1,2);
% plot(t(1:end-1)*c,s(Nz/2-1,:));
% xlabel('position[m]');
% ylabel('Amplitude');
% title('Clean Sinogram Central Cross Section');

% Add Raw Singoram SNR to SNR plot
figure(20);
plot(z(1)+ points_pos*dz,s_SNR,'DisplayName','Singoram','LineWidth',1.5);
xlabel('position[m]','FontSize',16);
ylabel('SNR[dB]','FontSize',16);
title('SNR as a function of depth','FontSize',20);
legend

% Add FWHM to Plot
s_FWHM=sinogram_FWHM(s(line_index,:),points_pos,space,time_res,dz);
figure(21);
plot(z(1)+ points_pos*dz,s_FWHM,'DisplayName','Sinogram','LineWidth',1.5);
xlabel('position[m]','FontSize',16);
ylabel('FWHM[m]','FontSize',16);
title('FWHM as a function of depth','FontSize',20);
legend

%% SAFT

saft_n=SAFT(s_noise,t,detector_pos,c,F,D,'rect');
saft=SAFT(s,t,detector_pos,c,F,D,'rect');

figure(4);
subplot(2,1,1);
imagesc(t_plot(1:end-1),detector_pos,saft);
% ylim([-0.6e-3 0.6e-3]);
% xlim([-0.8e-6 0.8e-6]);
colormap bone
colorbar
xlabel('t[sec]','FontSize',16);
ylabel('detector position[m]','FontSize',16);
title('SAFT Clean','FontSize',20);
subplot(2,1,2);
plot(t_plot(1:end-1),saft(line_index,:));
xlabel('t[sec]','FontSize',16);
ylabel('amplitude','FontSize',16);
% xlim([-0.8e-6 0.8e-6]);
title('Center Cross Section','FontSize',16);

figure(5);
subplot(2,1,1);
imagesc(t_plot(1:end-1),detector_pos,saft_n);
% ylim([-0.6e-3 0.6e-3]);
% xlim([-0.8e-6 0.8e-6]);
colormap bone
colorbar
xlabel('t[sec]','FontSize',16);
ylabel('detector position[m]','FontSize',16);
title('SAFT Noisy','FontSize',20);
subplot(2,1,2);
plot(t_plot(1:end-1),saft_n(line_index,:));
xlabel('t[sec]','FontSize',16);
ylabel('amplitude','FontSize',16);
% xlim([-0.8e-6 0.8e-6]);
title('Center Cross Section','FontSize',20);

% Add Raw Singoram SNR to SNR plot
saft_SNR=sinogram_SNR(saft_n(line_index,:),points_pos,space,time_res);
figure(20);
hold on
plot(z(1)+ points_pos*dz,saft_SNR,'DisplayName','SAFT','LineWidth',1.5);
drawnow
hold off

% Add FWHM to Plot
saft_FWHM=sinogram_FWHM(saft(line_index,:),points_pos,space,time_res,dz);
figure(21);
hold on
plot(z(1)+ points_pos*dz,saft_FWHM,'DisplayName','SAFT','LineWidth',1.5);
drawnow
hold off


%% SAFT-SIR 

saft_sir_n=SAFT_SIR(s_noise,t,detector_pos,c,F,D,'rect');
saft_sir=SAFT_SIR(s,t,detector_pos,c,F,D,'rect');


figure(8);
subplot(2,1,1);
imagesc(t_plot(1:end-1),detector_pos,saft_sir_n);
% ylim([-0.6e-3 0.6e-3]);
% xlim([-0.8e-6 0.8e-6]);
colormap bone
colorbar
xlabel('t[sec]','FontSize',16);
ylabel('detector position[m]','FontSize',16);
title('SAFT-SIR Noisy','FontSize',20);
subplot(2,1,2);
plot(t_plot(1:end-1),saft_sir_n(line_index,:));
xlabel('t[sec]','FontSize',16);
ylabel('amplitude','FontSize',16);
% xlim([-0.8e-6 0.8e-6]);
title('Center Cross Section','FontSize',20);

figure(9);
subplot(2,1,1);
imagesc(t_plot(1:end-1),detector_pos,saft_sir);
% ylim([-0.6e-3 0.6e-3]);
% xlim([-0.8e-6 0.8e-6]);
colormap bone
colorbar
xlabel('t[sec]','FontSize',16);
ylabel('detector position[m]','FontSize',16);
title('SAFT-SIR Clean','FontSize',20);
subplot(2,1,2);
plot(t_plot(1:end-1),saft_sir(line_index,:));
xlabel('t[sec]','FontSize',16);
ylabel('amplitude','FontSize',16);
% xlim([-0.8e-6 0.8e-6]);
title('Center Cross Section','FontSize',20);

% Add Raw Singoram SNR to SNR plot
saft_sir_SNR=sinogram_SNR(saft_sir_n(line_index,:),points_pos,space,time_res);
figure(20);
hold on
plot(z(1)+ points_pos*dz,saft_sir_SNR,'DisplayName','SAFT-SIR','LineWidth',1.5);
drawnow
hold off

% Add FWHM to Plot
saft_sir_FWHM=sinogram_FWHM(saft_sir(line_index,:),points_pos,space,time_res,dz);
figure(21);
hold on
plot(z(1)+ points_pos*dz,saft_sir_FWHM,'DisplayName','SAFT-SIR','LineWidth',1.5);
drawnow
hold off

%% SAFT Spherical-SIR

saft_sir_sph_n=SAFT_SIR_sph(s_noise,t,detector_pos,c,F,D,sir);
saft_sir_sph=SAFT_SIR_sph(s,t,detector_pos,c,F,D,sir);


figure(10);
subplot(2,1,1);
imagesc(t_plot(1:end-1),detector_pos,saft_sir_sph_n);
% ylim([-0.6e-3 0.6e-3]);
% xlim([-0.8e-6 0.8e-6]);
colormap bone
colorbar
xlabel('t[sec]','FontSize',16);
ylabel('detector position[m]','FontSize',16);
title('SAFT Spherical-SIR Noisy','FontSize',20);
subplot(2,1,2);
plot(t_plot(1:end-1),saft_sir_sph_n(line_index,:));
xlabel('t[sec]','FontSize',16);
ylabel('amplitude','FontSize',16);
% xlim([-0.8e-6 0.8e-6]);
title('Center Cross Section','FontSize',20);

figure(11);
subplot(2,1,1);
imagesc(t_plot(1:end-1),detector_pos,saft_sir_sph);
% ylim([-0.6e-3 0.6e-3]);
% xlim([-0.8e-6 0.8e-6]);
colormap bone
colorbar
xlabel('t[sec]','FontSize',16);
ylabel('detector position[m]','FontSize',16);
title('SAFT Spherical-SIR Clean','FontSize',20);
subplot(2,1,2);
plot(t_plot(1:end-1),saft_sir_sph(line_index,:));
xlabel('t[sec]','FontSize',16);
ylabel('amplitude','FontSize',16);
% xlim([-0.8e-6 0.8e-6]);
title('Center Cross Section','FontSize',20);

% Add Raw Singoram SNR to SNR plot
saft_sir_sph_SNR=sinogram_SNR(saft_sir_sph_n(line_index,:),points_pos,space,time_res);
figure(20);
hold on
plot(z(1)+ points_pos*dz,saft_sir_sph_SNR,'DisplayName','SAFT Spherical-SIR','LineWidth',1.5);
drawnow
hold off

% Add FWHM to Plot
saft_sir_sph_FWHM=sinogram_FWHM(saft_sir_sph(line_index,:),points_pos,space,time_res,dz);
figure(21);
hold on
plot(z(1)+ points_pos*dz,saft_sir_sph_FWHM,'DisplayName','SAFT Spherical-SIR','LineWidth',1.5);
drawnow
hold off


%% Pointwise SAFT - a different implementation of SAFT that is easier to debug

saft_n=pointwise_SAFT(s_noise,t,detector_pos,c,F,D,'rect');
saft=pointwise_SAFT(s,t,detector_pos,c,F,D,'rect');

figure(6);
subplot(2,1,1);
imagesc(t_plot(1:end-1),detector_pos,saft);
ylim([-0.6e-3 0.6e-3]);
xlim([-0.8e-6 0.8e-6]);
colormap bone
colorbar
xlabel('t[sec]');
ylabel('detector position[m]');
title('SAFT Clean');
subplot(2,1,2);
plot(t_plot(1:end-1),saft(line_index,:));
xlabel('t[sec]');
ylabel('amplitude');
xlim([-0.8e-6 0.8e-6]);
title('Center Cross Section');

figure(7);
subplot(2,1,1);
imagesc(t_plot(1:end-1),detector_pos,saft_n);
ylim([-0.6e-3 0.6e-3]);
xlim([-0.8e-6 0.8e-6]);
colormap bone
colorbar
xlabel('t[sec]');
ylabel('detector position[m]');
title('SAFT Noisy');
subplot(2,1,2);
plot(t_plot(1:end-1),saft_n(line_index,:));
xlabel('t[sec]');
ylabel('amplitude');
xlim([-0.8e-6 0.8e-6]);
title('Center Cross Section');



%% SAFT-CF 

[saft_cf_n,CFn]=SAFT_CF(s_noise,t,detector_pos,c,F,D,'rect');
[saft_cf,CF]=SAFT_CF(s,t,detector_pos,c,F,D,'rect');


figure(6);
subplot(2,1,1);
imagesc(t_plot(1:end-1),detector_pos,saft_cf_n);
ylim([-0.6e-3 0.6e-3]);
xlim([-0.8e-6 0.8e-6]);
colormap bone
colorbar
xlabel('t[sec]');
ylabel('detector position[m]');
title('SAFT-CF Noisy');
subplot(2,1,2);
plot(t_plot(1:end-1),saft_cf_n(line_index,:));
xlabel('t[sec]');
ylabel('amplitude');
xlim([-0.8e-6 0.8e-6]);
title('Center Cross Section');

figure(7);
subplot(2,1,1);
imagesc(t_plot(1:end-1),detector_pos,saft_cf);
ylim([-0.6e-3 0.6e-3]);
xlim([-0.8e-6 0.8e-6]);
colormap bone
colorbar
xlabel('t[sec]');
ylabel('detector position[m]');
title('SAFT-CF Clean');
subplot(2,1,2);
plot(t_plot(1:end-1),saft_cf(line_index,:));
xlabel('t[sec]');
ylabel('amplitude');
xlim([-0.8e-6 0.8e-6]);
title('Center Cross Section');

%% SAFT-SIR-CF 

[saft_sir_cf_n,SIR_CFn]=SAFT_SIR_CF(s_noise,t,detector_pos,c,F,D,'rect');
[saft_sir_cf,SIR_CF]=SAFT_SIR_CF(s,t,detector_pos,c,F,D,'rect');


figure(10);
subplot(2,1,1);
imagesc(t_plot(1:end-1),detector_pos,saft_sir_cf_n);
ylim([-0.6e-3 0.6e-3]);
xlim([-0.8e-6 0.8e-6]);
colormap bone
colorbar
xlabel('t[sec]');
ylabel('detector position[m]');
title('SAFT-SIR-CF Noisy');
subplot(2,1,2);
plot(t_plot(1:end-1),saft_sir_cf_n(Nz/2-1,:));
xlabel('t[sec]');
ylabel('amplitude');
xlim([-0.8e-6 0.8e-6]);
title('Center Cross Section');

figure(11);
subplot(2,1,1);
imagesc(t_plot(1:end-1),detector_pos,saft_sir_cf);
ylim([-0.6e-3 0.6e-3]);
xlim([-0.8e-6 0.8e-6]);
colormap bone
colorbar
xlabel('t[sec]');
ylabel('detector position[m]');
title('SAFT-SIR-CF Clean');
subplot(2,1,2);
plot(t_plot(1:end-1),saft_sir_cf(Nz/2-1,:));
xlabel('t[sec]');
ylabel('amplitude');
xlim([-0.8e-6 0.8e-6]);
title('Center Cross Section');

%% SAFT delay check


saft_delay=SAFT_delay(s,t,detector_pos,c,F,D);

figure(12);
subplot(2,1,1);
imagesc(t_plot(1:end),y,s);
colormap bone
colorbar
xlabel('t[sec]');
ylabel('y[m]');
title('sinogram');
subplot(2,1,2);
imagesc(t_plot(1:end),y,squeeze(saft_delay(line_index,:,:)));
colormap bone
colorbar
xlabel('t[sec]');
ylabel('y[m]');
title('SAFT delay sinogram');

