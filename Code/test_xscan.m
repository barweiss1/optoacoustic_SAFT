% --------------- Comments ---------------
% This test is the same as Test_Noise but here we do scans in the x-axis
% instead of y-axis there, meaning the ortogonal lateral direction 

% --------------- Code ---------------

%% clear Workspace
clear all
close all
clc

%% parameters definition

% Parameter definition. All units are SI. The focus of the transducer is
% always at the center of the image.

% These are the parameters from the paper - "Improved optoacoustic microscopy through
% three-dimensional spatial impulse response synthetic aperture focusing
% technique", Jake Turner et al to recreate Fig. 2
D = 24.6e-3; % sensor diameter - NA = 45 degrees
F = 12.3e-3; 
c = 1480 ; 
time_res = 5 ;
space_res = 25e-6; % resolution in z-axis for smooth time derivative of the SIR
image_width = 6.6e-6*c; % width that gives 2.8us of depth in time
image_height = 4.4e-6*c; % height that gives 2.8us of depth in time
x_width = 4.4e-6*c; % width of scan in x direction

Nz = 450; % Number of pixels in the horizontal direction. 
Ny = 300; % Number of pixels in the vertical diraction
Nx = 300; % Number of pixels in x-scan direction
detector_num = Nx; % number of detector positions
x = linspace(-x_width/2,x_width/2,Nx);

line_index=Nx/2;

% create rectangle image
I=zeros(Ny,Nz);
width=1;
height=60;
x_center=Nz/2-floor(width/2);
% y_center=Ny/2-floor(height/2);
y_center = Ny/2;
space=Nz/6;
rect_num=7;
line_num=5;
radius = 6;
angle=20;


%% plot original image

% I=create_rect_image(I,rect_num,space,x_center,y_center,width,height);
I=create_circle_image(I,rect_num,space,x_center,y_center,radius);
% I = create_line_image(I,line_num,space,x_center,y_center,height,angle);

points_pos=x_center-space*((rect_num-1)/2):space:x_center+space*((rect_num-1)/2);
% points_pos=x_center-space*((line_num-1)/2):space:x_center+space*((line_num-1)/2);

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

%% Sinogram - Long Run Time

% create sinogram 
[s,detector_pos]=x_scan(I,D,F,c,time_res,space_res,image_width,image_height,x_width,detector_num,Nz,Ny,Nx);

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
xlabel('t[sec]');
ylabel('detector position[m]');
mytitle = ['noisy sinogram \sigma=', num2str(sigma_per), '% of maximum'];
title(mytitle,'Interpreter','tex');
subplot(2,1,2);
plot(t_plot(1:end-1),s_noise(line_index,:));
xlabel('t[sec]');
ylabel('amplitude');
% xlim([-0.8e-6 0.8e-6]);
title('Center Cross Section');


figure(3);
subplot(2,1,1);
imagesc(t_plot(1:end),detector_pos,s);
% ylim([-0.6e-3 0.6e-3]);
% xlim([-0.8e-6 0.8e-6]);
colormap bone
colorbar
xlabel('t[sec]');
ylabel('detector position[m]');
title('Clean Sinogram');
subplot(2,1,2);
plot(t_plot(1:end-1),s(line_index,:));
xlabel('t[sec]');
ylabel('amplitude');
% xlim([-0.8e-6 0.8e-6]);
title('Center Cross Section');

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
plot(z(1)+ points_pos*dz,s_SNR,'DisplayName','Singoram');
xlabel('position[m]');
ylabel('SNR[dB]');
title('SNR as a function of depth');
legend

% Add FWHM to Plot
s_FWHM=sinogram_FWHM(s(line_index,:),points_pos,space,time_res,dz);
figure(21);
plot(z(1)+ points_pos*dz,s_FWHM,'DisplayName','Sinogram');
xlabel('position[m]');
ylabel('FWHM[m]');
title('FWHM as a function of depth');
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
xlabel('t[sec]');
ylabel('detector position[m]');
title('SAFT Clean');
subplot(2,1,2);
plot(t_plot(1:end-1),saft(line_index,:));
xlabel('t[sec]');
ylabel('amplitude');
% xlim([-0.8e-6 0.8e-6]);
title('Center Cross Section');

figure(5);
subplot(2,1,1);
imagesc(t_plot(1:end-1),detector_pos,saft_n);
% ylim([-0.6e-3 0.6e-3]);
% xlim([-0.8e-6 0.8e-6]);
colormap bone
colorbar
xlabel('t[sec]');
ylabel('detector position[m]');
title('SAFT Noisy');
subplot(2,1,2);
plot(t_plot(1:end-1),saft_n(line_index,:));
xlabel('t[sec]');
ylabel('amplitude');
% xlim([-0.8e-6 0.8e-6]);
title('Center Cross Section');

% Add Raw Singoram SNR to SNR plot
saft_SNR=sinogram_SNR(saft_n(line_index,:),points_pos,space,time_res);
figure(20);
hold on
plot(z(1)+ points_pos*dz,saft_SNR,'DisplayName','SAFT');
drawnow
hold off

% Add FWHM to Plot
saft_FWHM=sinogram_FWHM(saft(line_index,:),points_pos,space,time_res,dz);
figure(21);
hold on
plot(z(1)+ points_pos*dz,saft_FWHM,'DisplayName','SAFT');
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
xlabel('t[sec]');
ylabel('detector position[m]');
title('SAFT-SIR Noisy');
subplot(2,1,2);
plot(t_plot(1:end-1),saft_sir_n(line_index,:));
xlabel('t[sec]');
ylabel('amplitude');
% xlim([-0.8e-6 0.8e-6]);
title('Center Cross Section');

figure(9);
subplot(2,1,1);
imagesc(t_plot(1:end-1),detector_pos,saft_sir);
% ylim([-0.6e-3 0.6e-3]);
% xlim([-0.8e-6 0.8e-6]);
colormap bone
colorbar
xlabel('t[sec]');
ylabel('detector position[m]');
title('SAFT-SIR Clean');
subplot(2,1,2);
plot(t_plot(1:end-1),saft_sir(line_index,:));
xlabel('t[sec]');
ylabel('amplitude');
% xlim([-0.8e-6 0.8e-6]);
title('Center Cross Section');

% Add Raw Singoram SNR to SNR plot
saft_sir_SNR=sinogram_SNR(saft_sir_n(line_index,:),points_pos,space,time_res);
figure(20);
hold on
plot(z(1)+ points_pos*dz,saft_sir_SNR,'DisplayName','SAFT-SIR');
drawnow
hold off

% Add FWHM to Plot
saft_sir_FWHM=sinogram_FWHM(saft_sir(line_index,:),points_pos,space,time_res,dz);
figure(21);
hold on
plot(z(1)+ points_pos*dz,saft_sir_FWHM,'DisplayName','SAFT-SIR');
drawnow
hold off

