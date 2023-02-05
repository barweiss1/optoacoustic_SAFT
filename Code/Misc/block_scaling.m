function scaled_im = block_scaling(I,block_width,block_height)

% Function that rescales blocks in an image for easy presentation purposes
%
% Inputs:
% I = image to rescale
% block_width = width of the scale block 
% block_height = height of scale block
%
% Outputs:
% scaled_im = rescaled image block by block 

[X,Y] = size(I);

scaled_im=zeros(size(I));

max_val=max(max(I));

zero_im=(abs(I)<0.01*max_val);

for j=1:block_height:X-block_height
    for i=1:block_width:Y-block_width

        scaled_im(j:j+block_height,i:i+block_width)=im_double_norm(I(j:j+block_height,i:i+block_width));
        
    end
end

%scaled_im(zero_im==1)=0.5;


end