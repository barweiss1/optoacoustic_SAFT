function [s,detector_pos] =y_scan(I,F,c,time_res,image_width,image_height,Nz,Ny,sir)
% Function that creates a sinogram from an image of a desired source image
% and the acoustic lence parameters
% this function is used to create the 3D scan
% 
%
% INPUTS (all units are in SI)
% I(y,z) - image used for the scan 
% F = transducer focal length
% c = speed of sound
% time_res = number of time steps per dx 
% image_width = the width of the image
% image_height = the height of the image, for full image scan 
% Nz = number of pixels in z direction
% Ny = number of pixels in y direction
% sir = 2D sir used for the scan (acording to the x location of the
% transducer)
% 
% 
% OUTPUT
%
% s(y,t) = the sinogram of the image given the lence parameters
% detector_pos = vector of all the detector positions 


[Y,~] = size(I);
[t,y,z] = create_tyz_vectors(F,c,time_res,image_width,image_height,Nz,Ny);

% check if I is a sparse matrix to improve compute time
is_I_sparse = true;
non_zero_per = nnz(I)/numel(I);
if non_zero_per > 0.05 % if more than 5% are non zero than the matrix is not sparse
    is_I_sparse = false;
end


% define detector vectors
detector_num = Ny;
delta_D = floor(Ny/detector_num); % number of pixel to move the detector each time
pixel_tran = (Y/2 : -delta_D : -Y/2)+1; % vector of translations of the image for the scan
detector_pos = linspace(-image_height/2,image_height/2,detector_num); 


s=zeros(detector_num,length(t)-1); %initialize sinogram

% for each detector - move the image, calculate recieved signal and add it
% to the sinogram
for i=1:detector_num
   
    I_ref=imtranslate(I,[0, pixel_tran(i)]); % image moved to imitate the scan 
    % if I is sparse than I_ref is sparse, so convert to sparse
    % representation to improve compute time
    if is_I_sparse
       I_ref=sparse(I_ref); 
    end
    % sum up signals from every location to get detected signal
    for tt=1:length(t)-1
         p_xyt=I_ref.*squeeze(sir(tt,:,:)); % pressure at each point at time index tt
        s(i,tt)=sum(sum(p_xyt));
    end
    
end



end