function dh= time_deriv_im(h,dt,c,space_res)
% Function that takes the discrete time derivative of a sequence of images
%
% INPUTS
% h(t,y,z) = the image sequence that the time derivative is calculated on
% dt = time difference 
% c = speed of sound
% space_res = resolution in z-axis for smoothing of the derivative
%
%
% OUTPUT
%
% dh(t,y,z) = the discrete time derivative of the sequence

% create gaussian derivative filter - this is to avoid artifacts caused by
% high frequency amplification of 2 sample derivative

sigma = space_res/c;
t = 0:dt:2*sigma;
t = [ -t(end:-1:2) t];
gaus = (1/sqrt(2*pi*sigma^2))*exp(-t.^2/(2*sigma^2)); % define gaussian
deriv = [-1 1]; % define derivative filter
gaus_deriv=conv(gaus,deriv,'valid'); % create gaussian derivative filter


[T,Y,Z]= size(h);
dh=zeros(T,Y,Z);
for y=1:Y
    for z=1:Z
        
        % take the smooth time derivative for a point (z,y)
        dh(:,y,z) = conv(h(:,y,z),gaus_deriv,'same');
        
    end
    
end

