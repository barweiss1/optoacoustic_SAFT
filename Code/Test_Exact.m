% --------------- Comments ---------------
% Test of the SIR itself 

% --------------- Code ---------------


clear all
close all
clc

%% calculte SIR

% Parameter definition. All units are SI. The focus of the transducer is
% always at the center of the image.
D = 13e-3; % sensor diameter
%F = 25.4e-3; % sensor focal length
F = 25.4e-3;
c = 1480 ; 
time_res = 5 ;
image_width = 10e-3;

Nz = 120; % Number of pixels in the axial direction.

z = linspace(-image_width/2, image_width/2, Nz) ;

Ny = floor( Nz / 2 ); % Number of pixels in the perpendicular direction.
% The SIR is calculated for the positive half plane in front of the
% transducer and then mirrored ( in the wrapper below ).

y = z( Ny + 1 : end ) ;

%Definition of the time step and time vector
dx = z(2) - z(1) ; 
dt = dx / c ;
dt = dt / time_res ;
fs = 1/dt;

t1 = F -  ( sqrt( 2 ) / 2 ) * image_width ;
tend = F + ( sqrt( 2 ) / 2 ) * image_width ;

t = t1/c : dt : tend/c ;

% Calculation of the SIR.
h = Sph_F_exact_wrapper(D, F, c, t, y, z ) ;
dh = time_deriv_im(h,dt);


%% plot SIR for several points 

figure;
subplot(3,1,1);
plot(t(1:end-1),squeeze(dh(:,71,60)));
title("SIR time derivative at z=" + num2str(z(60)) + " y=" + num2str(y(1)));
xlabel("t");
subplot(3,1,2);
plot(t(1:end-1),squeeze(dh(:,61,40)));
title("SIR time derivative at z=" + num2str(z(40)) + " y=" + num2str(y(20)));
xlabel("t");
subplot(3,1,3);
plot(t(1:end-1),squeeze(dh(:,61,80)));
title("SIR time derivative at z=" + num2str(z(80)) + " y=" + num2str(y(20)));
xlabel("t");


figure;
subplot(3,1,1);
plot(t(1:end),squeeze(h(:,61,60)));
title("SIR at z=" + num2str(z(60)) + " y=" + num2str(y(1)));
xlabel("t");
subplot(3,1,2);
plot(t(1:end),squeeze(h(:,61,40)));
title("SIR at z=" + num2str(z(40)) + " y=" + num2str(y(1)));
xlabel("t");
subplot(3,1,3);
plot(t(1:end),squeeze(h(:,61,80)));
title("SIR at z=" + num2str(z(80)) + " y=" + num2str(y(1)));
xlabel("t");
