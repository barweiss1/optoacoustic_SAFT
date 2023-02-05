function I_norm = im_double_norm(I)
% Function normalizes a general image to double [0 1] range
%
% INPUTS
% I = original image 
% 
%
% OUTPUT
%
% I_norm = image normalized to [0 1] range

min_val=min(min(I));
max_val=max(max(I));

I_norm = (I - min_val)/(max_val - min_val);

end