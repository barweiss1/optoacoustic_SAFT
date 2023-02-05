% --------------- Comments ---------------
% This test file runs 3D scans on the image and tests SAFT performace on 3D
% scans.


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


D = 24.6e-3; % transducer diameter - NA = 45 degrees
F = 12.3e-3; % Focal length 
c = 1480 ; % speed of sound 
time_res = 3 ; % time to depth sample factor (Nt = Nz*time_res)
space_res = 25e-6; % resolution in z-axis for smooth time derivative of the SIR
image_width = 5.5e-6*c; % width that gives 5.5us of depth in time
image_height = 2.75e-6*c; % height that gives 2.7us of depth in time
x_width = 2.75e-6*c; % width of scan in x direction

Nz = 300; % Number of pixels in the depth direction. 
Ny = 150; % Number of pixels in the vertical direction
Nx = 150; % Number of pixels in x-scan direction
detector_num = Nx; % number of detector positions
x = linspace(-x_width/2,x_width/2,Nx);
dx = x(2)-x(1);

line_index=Nx/2; % index for cross section plot

% create the imaged object pattern 

% create angled line image
I=zeros(Ny,Nz);
width=1;
height=140;
x_center=Nz/2-floor(width/2);
% y_center=Ny/2-floor(height/2);
y_center = Ny/2;
space=floor(Nz/4);
line_num=3;
angle=30;
I = create_line_image(I,line_num,space,x_center,y_center,height,angle);
points_pos=x_center-space*((line_num-1)/2):space:x_center+space*((line_num-1)/2);


% % create horizontal line image
% I=zeros(Ny,Nz);
% width=4;
% height=2;
% x_center=Nz/2-floor(width/2);
% % y_center=Ny/2-floor(height/2);
% y_center = floor(Ny/2);
% space=floor(Nz/8);
% rect_num=5;
% I=create_rect_image(I,rect_num,space,x_center,y_center,width,height);
% points_pos=x_center-space*((rect_num-1)/2):space:x_center+space*((rect_num-1)/2);

% create circle image 
% I=create_circle_image(I,rect_num,space,x_center,y_center,radius);


%% plot original image


% create time vectors
[t,y,z] = create_tyz_vectors(F,c,time_res,image_width,image_height,Nz,Ny);
t_plot=t-F/c;

figure(1);
imagesc(z,[-1*flip(y) y],I);
colormap gray
colorbar
impixelinfo;
title('original image','FontSize',20);
xlabel('z[m]','FontSize',16);
ylabel('y[m]','FontSize',16);


%% Sinogram - Long Run Time

% create sinogram 
tic
[s,detector_pos_x,detector_pos_y]=scan_3D(I,D,F,c,time_res,space_res,image_width,image_height,x_width,detector_num,Nz,Ny,Nx);
% get run time
toc

%% Add Gaussian noise 

sigma_per=15; % std of the noise in [%] from maximum of sinogram
max_s=max(max(max(s)));
s_noise=s+(randn(size(s))*max_s*(sigma_per/100)); % adding gaussian noise 


%% plot sinogram and noisy sinogram 

dz = z(2) - z(1);
dx = x(2) - x(1);
point_distance = space*dz;

% MAP - maximum amplitude projection
s_map = squeeze(max(abs(s),[],1));
s_n_map = squeeze(max(abs(s_noise),[],1));

figure(2);
subplot(2,1,1);
imagesc(t_plot(1:end),detector_pos_y,s_n_map);
% ylim([-0.6e-3 0.6e-3]);
% xlim([-0.8e-6 0.8e-6]);
colormap bone
colorbar
xlabel('t[sec]');
ylabel('detector position[m]');
mytitle = ['noisy sinogram \sigma=', num2str(sigma_per), '% of maximum'];
title(mytitle,'Interpreter','tex');
subplot(2,1,2);
plot(t_plot(1:end-1),s_n_map(line_index,:));
xlabel('t[sec]');
ylabel('amplitude');
% xlim([-0.8e-6 0.8e-6]);
title('Center Cross Section');


figure(3);
subplot(2,1,1);
imagesc(t_plot(1:end),detector_pos_y,s_map);
% ylim([-0.6e-3 0.6e-3]);
% xlim([-0.8e-6 0.8e-6]);
colormap bone
colorbar
xlabel('t[sec]','FontSize',16);
ylabel('detector position[m]','FontSize',16);
title('Raw Signal','FontSize',20);
subplot(2,1,2);
plot(t_plot(1:end-1),s_map(line_index,:));
xlabel('t[sec]');
ylabel('amplitude');
% xlim([-0.8e-6 0.8e-6]);
title('Center Cross Section');


%% plot raw data cross sections 

figure(17);
subplot(3,1,1);
imagesc(c*t_plot(1:end-1),detector_pos_y,abs(squeeze(s(75,:,:))));
colormap bone
title_text = ['Raw Data Cross section in X direction at the focus plane 0 Degrees'];
colorbar
ax = gca;
ax.FontSize = 16;
xlabel('z[m]','FontSize',16);
ylabel('detector position[m]','FontSize',16);
title(title_text,'FontSize',20);
subplot(3,1,2);
imagesc(t_plot(1:end-1),detector_pos_y,abs(squeeze(s(85,:,:))));
colormap bone
title_text = ['Raw Data Cross section in X direction at ' , num2str(dx*10*1e3), '[mm] distance from the focus plane'];
colorbar
xlabel('t[sec]');
ylabel('detector position[m]');
title(title_text);
subplot(3,1,3);
imagesc(t_plot(1:end-1),detector_pos_y,abs(squeeze(s(95,:,:))));
colormap bone
title_text = ['Raw Data Cross section in X direction at ' , num2str(dx*20*1e3), '[mm] distance from the focus plane'];
colorbar
xlabel('t[sec]');
ylabel('detector position[m]');
title(title_text);


%% SAFT-SIR-Y

saft_y_n=SAFT_SIR_3D_y(s_noise,t,detector_pos_y,c,F,D,'rect');
saft_y=SAFT_SIR_3D_y(s,t,detector_pos_y,c,F,D,'rect');

% MAP - maximum amplitude projection
saft_y_map = squeeze(max(abs(saft_y),[],1));
saft_y_n_map = squeeze(max(abs(saft_y_n),[],1));

figure(6);
subplot(2,1,1);
imagesc(t_plot(1:end-1),detector_pos_y,saft_y_n_map);
% ylim([-0.6e-3 0.6e-3]);
% xlim([-0.8e-6 0.8e-6]);
colormap bone
colorbar
xlabel('t[sec]');
ylabel('detector position[m]');
title('SAFT-SIR-Y Noisy');
subplot(2,1,2);
plot(t_plot(1:end-1),saft_y_n_map(line_index,:));
xlabel('t[sec]');
ylabel('amplitude');
% xlim([-0.8e-6 0.8e-6]);
title('Center Cross Section');

figure(7);
subplot(2,1,1);
imagesc(t_plot(1:end-1),detector_pos_y,saft_y_map);
% ylim([-0.6e-3 0.6e-3]);
% xlim([-0.8e-6 0.8e-6]);
colormap bone
colorbar
xlabel('t[sec]');
ylabel('detector position[m]');
title('SAFT-SIR-Y Clean');
subplot(2,1,2);
plot(t_plot(1:end-1),saft_y_map(line_index,:));
xlabel('t[sec]');
ylabel('amplitude');
% xlim([-0.8e-6 0.8e-6]);
title('Center Cross Section');




%% SAFT-SIR-X

saft_x_n=SAFT_SIR_3D_x(s_noise,t,detector_pos_x,c,F,D,'rect');
saft_x=SAFT_SIR_3D_x(s,t,detector_pos_x,c,F,D,'rect');

% MAP - maximum amplitude projection
saft_x_map = squeeze(max(abs(saft_x),[],1));
saft_x_n_map = squeeze(max(abs(saft_x_n),[],1));


figure(8);
subplot(2,1,1);
imagesc(c*t_plot(1:end-1),detector_pos_y,saft_x_n_map);
% ylim([-0.6e-3 0.6e-3]);
% xlim([-0.8e-6 0.8e-6]);
colormap bone
colorbar
ax = gca;
ax.FontSize = 16;
xlabel('z[m]','FontSize',16);
ylabel('detector position[m]','FontSize',16);
title('SAFT-SIR-X Noisy 30 Degrees','FontSize',20);
subplot(2,1,2);
plot(t_plot(1:end-1),saft_x_n_map(line_index,:));
xlabel('t[sec]');
ylabel('amplitude');
% xlim([-0.8e-6 0.8e-6]);
title('Center Cross Section');


figure(9);
subplot(2,1,1);
imagesc(t_plot(1:end-1),detector_pos_y,saft_x_map);
% ylim([-0.6e-3 0.6e-3]);
% xlim([-0.8e-6 0.8e-6]);
colormap bone
colorbar
xlabel('t[sec]');
ylabel('detector position[m]');
title('SAFT-SIR-X Clean');
subplot(2,1,2);
plot(t_plot(1:end-1),saft_x_map(line_index,:));
xlabel('t[sec]');
ylabel('amplitude');
% xlim([-0.8e-6 0.8e-6]);
title('Center Cross Section');



%% plot cross sections 

figure(15);
subplot(3,1,1);
imagesc(c*t_plot(1:end-1),detector_pos_y,abs(squeeze(saft_x_n(75,:,:))));
colormap bone
title_text = ['SAFT-SIR-X Cross section in X direction at the focus plane 30 Degrees'];
colorbar
ax = gca;
ax.FontSize = 16;
xlabel('z[m]','FontSize',16);
ylabel('detector position[m]','FontSize',16);
title(title_text,'FontSize',20);
subplot(3,1,2);
imagesc(t_plot(1:end-1),detector_pos_y,abs(squeeze(saft_x_n(85,:,:))));
colormap bone
title_text = ['SAFT-SIR-X Cross section in X direction at ' , num2str(dx*10*1e3), '[mm] distance from the focus plane'];
colorbar
xlabel('t[sec]');
ylabel('detector position[m]');
title(title_text);
subplot(3,1,3);
imagesc(t_plot(1:end-1),detector_pos_y,abs(squeeze(saft_x_n(95,:,:))));
colormap bone
title_text = ['SAFT-SIR-X Cross section in X direction at ' , num2str(dx*20*1e3), '[mm] distance from the focus plane'];
colorbar
xlabel('t[sec]');
ylabel('detector position[m]');
title(title_text);


%% Displacement Test  

% calculate ideal line location
dy = y(2) - y(1);
y_len = height*dy*cos(angle*pi/180); % height is actualy the line length
ideal_y = [ -y_len/2 y_len/2];
ideal_x = [0 0];


% divide time to 5 segments
T = length(t_plot);
extra_pixels = T - Nz*time_res; % the time vector represnts a longer depth than the z vector so we extract the extra pixels to find the lines' coordinates in the time vectir 
points_pos=x_center-space*((line_num-1)/2):space:x_center+space*((line_num-1)/2);
line_pos_t = floor(points_pos*time_res + extra_pixels/2); % convert points_pos to time vector positions 
%line_pos_t = [220 320 420 520 620];
pix_range = floor(-space*time_res/2:space*time_res/2);

% plot MAPs for each line and its relevant time segment
figure(18);
for i = 1:line_num
    
    % calculate MAP for each segment in the depth direction 
    seg_map = squeeze(max(abs(saft_x_n(:,:,line_pos_t(i)+pix_range)),[],3));
    
    % plot
    subplot(2,2,i); % change to fit number of lines in the image
    imagesc(detector_pos_y,detector_pos_x,seg_map);
    hold on
    line = plot(ideal_y,ideal_x,'--r','LineWidth',3.0);
    line.Color(4) = 0.4; % set line transperancy
    hold off
    colormap bone
    title_text = ['SAFT-SIR-X MAP in Z direction for the ', iptnum2ordinal(i) ,' line'];
    colorbar
    ax = gca;
    ax.FontSize = 12;
    xlabel('y[m]');
    ylabel('x[m]');
    title(title_text);
    
end

%% plot depth direction MAP for each line

% divide time to 5 segments
T = length(t_plot);
extra_pixels = T - Nz*time_res; % the time vector represnts a longer depth than the z vector so we extract the extra pixels to find the lines' coordinates in the time vectir 
points_pos=x_center-space*((line_num-1)/2):space:x_center+space*((line_num-1)/2);
line_pos_t = floor(points_pos*time_res + extra_pixels/2); % convert points_pos to time vector positions 
%line_pos_t = [220 320 420 520 620];
pix_range = floor(-space*time_res/2:space*time_res/2);

% plot MAPs for each line and its relevant time segment
figure(16);
for i = 1:line_num
    
    % calculate MAP for each segment in the depth direction 
    seg_map = squeeze(max(abs(saft_x_n(:,:,line_pos_t(i)+pix_range)),[],3));
    
    % plot
    subplot(2,2,i);
    imagesc(detector_pos_y,detector_pos_x,seg_map);
    colormap bone
    title_text = ['SAFT-SIR-X MAP in Z direction for the ', iptnum2ordinal(i) ,' line'];
    colorbar
    ax = gca;
    ax.FontSize = 12;
    xlabel('y[m]');
    ylabel('x[m]');
    title(title_text);
    
end


%% SAFT-SIR-cross

saft_cross_n=SAFT_SIR_3D_cross(s_noise,t,detector_pos_x,detector_pos_y,c,F,D,'rect');
saft_cross=SAFT_SIR_3D_cross(s,t,detector_pos_x,detector_pos_y,c,F,D,'rect');

% MAP - maximum amplitude projection
saft_cross_map = squeeze(max(abs(saft_cross),[],1));
saft_cross_n_map = squeeze(max(abs(saft_cross_n),[],1));

figure(10);
subplot(2,1,1);
imagesc(t_plot(1:end-1),detector_pos_y,saft_cross_n_map);
% ylim([-0.6e-3 0.6e-3]);
% xlim([-0.8e-6 0.8e-6]);
colormap bone
colorbar
xlabel('t[sec]');
ylabel('detector position[m]');
title('SAFT-SIR-Cross Noisy');
subplot(2,1,2);
plot(t_plot(1:end-1),saft_cross_n_map(line_index,:));
xlabel('t[sec]');
ylabel('amplitude');
% xlim([-0.8e-6 0.8e-6]);
title('Center Cross Section');

figure(11);
subplot(2,1,1);
imagesc(t_plot(1:end-1),detector_pos_y,saft_cross_map);
% ylim([-0.6e-3 0.6e-3]);
% xlim([-0.8e-6 0.8e-6]);
colormap bone
colorbar
xlabel('t[sec]');
ylabel('detector position[m]');
title('SAFT-SIR-Cross Clean');
subplot(2,1,2);
plot(t_plot(1:end-1),saft_cross_map(line_index,:));
xlabel('t[sec]');
ylabel('amplitude');
% xlim([-0.8e-6 0.8e-6]);
title('Center Cross Section');

% Add Raw Singoram SNR to SNR plot
saft_cross_SNR=sinogram_SNR(saft_cross_n_map(line_index,:),points_pos,space,time_res);
figure(20);
hold on
plot(z(1)+ points_pos*dz,saft_cross_SNR,'DisplayName','SAFT-SIR-Cross');
drawnow
hold off

% Add FWHM to Plot
saft_cross_FWHM=sinogram_FWHM(saft_cross_map(line_index,:),points_pos,space,time_res,dz);
figure(21);
hold on
plot(z(1)+ points_pos*dz,saft_cross_FWHM,'DisplayName','SAFT-SIR-Cross');
drawnow
hold off


