function dh = Sph_SIR_map_wrapper(D, F, c, t, y, z ,space_res)

% Function to calculate the SIR of a spherically focused round transducer,
% within a given ROI - here the SIR represents the reciever sensitivity at
% each position in the ROI
%
% INPUT (all units in SI)
% D = transducer diameter
% F = transducer focal length
% c = speed of sound
% t = time vector in which the SIR is calculated. It should be bigger than
% (see following) the ROI
% y = points perpendicular to the axis of the transducer where the SIR 
% is to be computed. The vectors yz define a half plane (ROI) with only positive
% values of y.
% z = points along the axis of the transducer where the SIR is to be
% computed. The focus of the transducer is located at z == 0 (i.e., for z
% vectors that do not include z == 0, the focus is not included in the
% calculation).
%
% OUTPUT
%
% h(t,y,z) = SIR of a spherically focused detector.

Nz = length( z ) ;
Ny = length( y ) ;
odd = rem( Nz,2 );

dh = zeros( Ny , Nz );

% define gaussian derivative filter
dt = t(2) - t(1);
sigma = space_res/c;
t_gaus = 0:dt:2*sigma;
t_gaus = [ -t_gaus(end:-1:2) t_gaus];
gaus = (1/sqrt(2*pi*sigma^2))*exp(-t_gaus.^2/(2*sigma^2)); % define gaussian
deriv = [-1 1]; % define derivative filter
gaus_deriv=conv(gaus,deriv,'valid'); % create gaussian derivative filter

for ii = 1:Nz

    for jj = 1:Ny
    
        h0 = Sph_F_exact(D, F, c, t, y( jj ), z( ii ) ) ;
        dh0 = conv(h0,gaus_deriv,'same');
        dh( jj, ii ) = max(abs(dh0));
          
    end

end 

% If the number of pixels in the ROI is uneven, the amplitude at the focus
% is interpolated from the amplitudes of the SIRs along the axis.


dh = [flip(dh,1) ; dh];
