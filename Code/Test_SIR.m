% --------------- Comments ---------------
% A test for the SIR itself -  test different distances, focuses and
% diameters

% --------------- Code ---------------


clear all
close all
clc

%% calculte SIR - 10mm image

% Parameter definition. All units are SI. The focus of the transducer is
% always at the center of the image.
D = 13e-3; % sensor diameter
%F = 25.4e-3; % sensor focal length
F = 6e-3;
c = 1480 ; 
time_res = 5 ;
image_width = 10e-3;

Nz=250;
detector_num = 250; % number of detector positions

% create rectangle image
I=zeros(Nz);
x_center=Nz/2;
y_center=Nz/2;
space=Nz/25;
rect_num=9;
width=1;
height=1;
I=create_rect_image(I,rect_num,space,x_center,y_center,width,height);



Nz = 250; % Number of pixels in the axial direction.

[t_10mm,y_10mm,z_10mm] =create_tyz_vectors(F,c,time_res,image_width,Nz);
dt_10mm = t_10mm(2) - t_10mm(1);

% Calculation of the SIR.
h_10mm = Sph_F_exact_wrapper(D, F, c, t_10mm, y_10mm, z_10mm ) ;
dh_10mm = time_deriv_im(h_10mm,dt_10mm);


%% Find amplification scheme

SIR_center=squeeze(dh_10mm(:,:,125));
peaks=max(abs(SIR_center));

figure;
plot(z_10mm,peaks);
xlabel('z[m]');
ylabel('Peak Value');
title('Peak Values per distace');


%% calculate sinogram 

[s_10mm,detector_pos_10mm]=sinogram2(I,D,F,c,time_res,image_width,detector_num,Nz);

%% calculte SIR - 2mm image

% Parameter definition. All units are SI. The focus of the transducer is
% always at the center of the image.
D = 13e-3; % sensor diameter
%F = 25.4e-3; % sensor focal length
F = 25.4e-3;
c = 1480 ; 
time_res = 5 ;
image_width = 2e-3;

Nz = 250; % Number of pixels in the axial direction.

[t_2mm,y_2mm,z_2mm] =create_tyz_vectors(F,c,time_res,image_width,Nz);
dt_2mm = t_2mm(2) - t_2mm(1);

% Calculation of the SIR.
h_2mm = Sph_F_exact_wrapper(D, F, c, t_2mm, y_2mm, z_2mm ) ;
dh_2mm = time_deriv_im(h_2mm,dt_2mm);

%% calculate sinogram 

[s_2mm,detector_pos_2mm]=sinogram2(I,D,F,c,time_res,image_width,detector_num,Nz);


%% plot SIRs - image witdhs

figure(1);
subplot(3,1,1);
plot(t_10mm(1:end-1),squeeze(dh_10mm(:,125,125)));
title("SIR for 10mm at z=" + num2str(z_10mm(125)) + " y=" + num2str(y_10mm(1)));
xlabel("t");
subplot(3,1,2);
plot(t_10mm(1:end-1),squeeze(dh_10mm(:,125,100)));
title("SIR for 10mm at z=" + num2str(z_10mm(100)) + " y=" + num2str(y_10mm(1)));
xlabel("t");
subplot(3,1,3);
plot(t_10mm(1:end-1),squeeze(dh_10mm(:,125,150)));
title("SIR for 10mm at z=" + num2str(z_10mm(150)) + " y=" + num2str(y_10mm(1)));
xlabel("t");

figure(2);
subplot(3,1,1);
plot(t_2mm(1:end),squeeze(h_2mm(:,126,125)));
title("SIR for 2mm at z=" + num2str(z_2mm(125)) + " y=" + num2str(y_2mm(1)));
xlabel("t");
subplot(3,1,2);
plot(t_2mm(1:end),squeeze(h_2mm(:,126,100)));
title("SIR for 2mm at z=" + num2str(z_2mm(100)) + " y=" + num2str(y_2mm(1)));
xlabel("t");
subplot(3,1,3);
plot(t_2mm(1:end),squeeze(h_2mm(:,126,150)));
title("SIR for 2mm at z=" + num2str(z_2mm(150)) + " y=" + num2str(y_2mm(1)));
xlabel("t");

%% plot sinograms 

figure(3);
subplot(2,1,1);
imagesc(t_10mm(1:end-1),y_10mm,s_10mm);
colorbar
title('sinogram for 10mm image width');
xlabel('t[sec]');
ylabel('y[m]');

subplot(2,1,2);
imagesc(t_2mm(1:end-1),y_2mm,s_2mm);
colorbar
title('sinogram for 2mm image width');
xlabel('t[sec]');
ylabel('y[m]');

figure(4);
imagesc(s_10mm-s_2mm);
colorbar
title('difference between the sinograms');
xlabel('t[sec]');
ylabel('y[m]');

%% calculte SIR - F 12mm 

% Parameter definition. All units are SI. The focus of the transducer is
% always at the center of the image.
D = 13e-3; % sensor diameter
%F = 25.4e-3; % sensor focal length
F = 12e-3;
c = 1480 ; 
time_res = 5 ;
image_width = 4e-3;

detector_num = 250; % number of detector positions

% create rectangle image
I=zeros(Nz);
x_center=Nz/2;
y_center=Nz/2;
space=Nz/25;
rect_num=9;
width=1;
height=1;
I=create_rect_image(I,rect_num,space,x_center,y_center,width,height);



Nz = 250; % Number of pixels in the axial direction.

[t_F12mm,y_F12mm,z_F12mm] =create_tyz_vectors(F,c,time_res,image_width,Nz);
dt_F12mm = t_F12mm(2) - t_F12mm(1);

% Calculation of the SIR.
h_F12mm = Sph_F_exact_wrapper(D, F, c, t_F12mm, y_F12mm, z_F12mm ) ;
dh_F12mm = time_deriv_im(h_F12mm,dt_F12mm);


%% Find amplification scheme

SIR_center=squeeze(dh_F12mm(:,:,125));
peaks=max(abs(SIR_center));

figure;
plot(z_F12mm,peaks);
xlabel('z[m]');
ylabel('Peak Value');
title('Peak Values per distace');

%% calculte SIR - F 6mm 

% Parameter definition. All units are SI. The focus of the transducer is
% always at the center of the image.
D = 13e-3; % sensor diameter
%F = 25.4e-3; % sensor focal length
F = 6e-3;
c = 1480 ; 
time_res = 5 ;
image_width = 4e-3;

detector_num = 250; % number of detector positions

% create rectangle image
I=zeros(Nz);
x_center=Nz/2;
y_center=Nz/2;
space=Nz/25;
rect_num=9;
width=1;
height=1;
I=create_rect_image(I,rect_num,space,x_center,y_center,width,height);


Nz = 250; % Number of pixels in the axial direction.

[t_F6mm,y_F6mm,z_F6mm] =create_tyz_vectors(F,c,time_res,image_width,Nz);
dt_F6mm = t_F6mm(2) - t_F6mm(1);

% Calculation of the SIR.
h_F6mm = Sph_F_exact_wrapper(D, F, c, t_F6mm, y_F6mm, z_F6mm ) ;
dh_F6mm = time_deriv_im(h_F6mm,dt_F6mm);

%% Find amplification scheme

SIR_center=squeeze(dh_F12mm(:,:,125));
peaks=max(abs(SIR_center));

figure;
plot(z_F12mm,peaks);
xlabel('z[m]');
ylabel('Peak Value');
title('Peak Values per distace');


%% plot SIRs - Focal Lengths

figure(5);
subplot(3,1,1);
plot(t_F12mm(1:end-1),squeeze(dh_F12mm(:,125,125)));
title("SIR for F12mm at z=" + num2str(z_F12mm(125)) + " y=" + num2str(y_F12mm(1)));
xlabel("t");
subplot(3,1,2);
plot(t_F12mm(1:end-1),squeeze(dh_F12mm(:,125,100)));
title("SIR for F12mm at z=" + num2str(z_F12mm(100)) + " y=" + num2str(y_F12mm(1)));
xlabel("t");
subplot(3,1,3);
plot(t_F12mm(1:end-1),squeeze(dh_F12mm(:,125,150)));
title("SIR for F12mm at z=" + num2str(z_F12mm(150)) + " y=" + num2str(y_F12mm(1)));
xlabel("t");

figure(6);
subplot(3,1,1);
plot(t_F6mm(1:end),squeeze(h_F6mm(:,125,125)));
title("SIR for F6mm at z=" + num2str(z_F6mm(125)) + " y=" + num2str(y_F6mm(1)));
xlabel("t");
subplot(3,1,2);
plot(t_F6mm(1:end),squeeze(h_F6mm(:,125,100)));
title("SIR for F6mm at z=" + num2str(z_F6mm(100)) + " y=" + num2str(y_F6mm(1)));
xlabel("t");
subplot(3,1,3);
plot(t_F6mm(1:end),squeeze(h_F6mm(:,125,150)));
title("SIR for F6mm at z=" + num2str(z_F6mm(150)) + " y=" + num2str(y_F6mm(1)));
xlabel("t");

