function [s,detector_pos] = x_scan(I,D,F,c,time_res,space_res,image_width,image_height,x_width,detector_num,Nz,Ny,Nx)
% Function that creates a sinogram from an image of a desired source image
% and the acoustic lence parameters in the x direction
% we scan here a 2D slice of the image where the transducer moves with
% respect to the image
% 
%
% INPUTS (all units are in SI)
% I(y,z) - image we want to recieve
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
% s(x,t) = the sinogram of the image given the lence parameters
% detector_pos = vector of all the detector positions 

% define scan and image direction vectors
[t,y,z] = create_tyz_vectors(F,c,time_res,image_width,image_height,Nz,Ny);
dt = t(2)-t(1);
detector_num_up = floor((detector_num+1)/2);
detector_pos = linspace(0,x_width/2,detector_num_up);

% check if I is a sparse matrix to improve compute time
non_zero_per = nnz(I)/numel(I);
if non_zero_per <= 0.05 % if more than 5% are non zero than the matrix is not sparse
    I = sparse(I);
end


s=zeros(detector_num_up,length(t)-1); %initialize sinogram

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

    % sum up signals from every location to get detected signal
    for tt=1:length(t)-1
        p_xyt=I.*squeeze(dh(tt,:,:)); % pressure at each point at time index tt
        s(i,tt)=sum(sum(p_xyt));
    end
    
end

% SIR is symmetric so we can duplicate the it symmetrically to save
% computations 

detector_pos = [-1*flip(detector_pos(2:end)) detector_pos];
s = [flip(s,1) ; s(2:end,:)];


end