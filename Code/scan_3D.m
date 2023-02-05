function [s,detector_pos_x,detector_pos_y] = scan_3D(I,D,F,c,time_res,space_res,image_width,image_height,x_width,detector_num,Nz,Ny,Nx)
% Function that creates a sinogram from a 2D desired source image
% and the acoustic lence parameters 
% we create a 3D scan by assuming that other than the 2D image inputed the
% rest of the space has no emitting sources.
% we use a 2D image to save computations by using 2D sums instead of 3D.
% 
% we use the y_scan function to create planes of scans in y direction. and
% each time move the SIR in the x direction and pass the moved SIR to
% y_scan to simulate movement of transducer in that direction. that way we
% get a 3D scan.
% 
%
% INPUTS (all units are in SI)
% I(y,z) - the 2D image used for the scan - the image is parallel to one of
% the main axis of the transducer
% D = transducer diameter
% F = transducer focal length
% c = speed of sound
% time_res = number of time steps per dx 
% space_res = resolution in z-axis for smooth time derivative of the SIR
% image_width = the width of the image
% image_height = the height of the image
% x_width = width of scan in x direction 
% detector_num = the number of the places to put the detector at
% ,for full image scan 
% Nz = number of pixels in z direction
% Ny = number of pixels in y direction
% Nx = number of pixels in x direction scan 
% 
% 
% 
% 
% OUTPUT
%
% s(x,y,t) = the sinogram of the image given the lence parameters
% detector_pos = vector of all the detector positions 

% define scan and image direction vectors
[t,y,z] = create_tyz_vectors(F,c,time_res,image_width,image_height,Nz,Ny);
dt = t(2)-t(1);
detector_num_up = floor((detector_num+1)/2);
detector_pos = linspace(0,x_width/2,detector_num_up);
detector_pos_y = linspace(-image_height/2,image_height/2,Ny); 

s=zeros(detector_num_up,Ny,length(t)-1); %initialize sinogram

% for each detector - move the detector, calculate the SIR cross-section that meats the image and add it
% to the sinogram at the relevant x-value
for i=1:detector_num_up
  
    
    % calculate r coordinates on the transducer plane for the SIR
    % calculation
    r = sqrt(detector_pos(i).^2 + y.^2);
    
    % calculate the SIR for the needed area for the current scan
    h = Sph_F_exact_wrapper(D, F, c, t, r, z );
    dh = time_deriv_im(h,dt,c,space_res); % taking the time derivative to get the pressure wave
    dh = -1*dh; % the formula for pressure has a - sign
    
    % add y_scan line to the sinogram 
    [s_y,~] = y_scan(I,F,c,time_res,image_width,image_height,Nz,Ny,dh);
    s(i,:,:) = s_y;
    
end

% SIR is symmetric so we can duplicate the it symmetrically to save
% computations 

detector_pos_x = [-1*flip(detector_pos(2:end)) detector_pos];
s = [flip(s,1) ; s(2:end,:,:)];


end