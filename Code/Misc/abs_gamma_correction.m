function I_gamma = abs_gamma_correction(I,gamma)

% Function that does gamma correction for an image with negative values
% the idea is to do gamma correction in the absolute value instead of the
% normalized image and then after the correction return the sign to each
% pixel
%
%
% INPUTS
% I = original image 
% gamma = gamma value
%
% OUTPUT
%
% I_gamma = image after the gamma correction 

I_abs=abs(I); %create absolute value image for the gamma correction
I_sign=sign(I); % save the sign of the pixel to reconstruct 

I_norm=im_double_norm(I_abs); % convert to double [0 1] representation

I_abs_gamma=I_norm.^gamma; % gamma correction

I_gamma=I_abs_gamma.*I_sign; % correct sign of each pixel

end