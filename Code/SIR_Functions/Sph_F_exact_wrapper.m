function h = Sph_F_exact_wrapper(D, F, c, t, y, z )

% Function to calculate the SIR of a spherically focused round transducer,
% within a given ROI.
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

h = zeros( length( t ) , Ny , Nz );

for ii = 1:Nz

    for jj = 1:Ny
    
        h0 = Sph_F_exact(D, F, c, t, y( jj ), z( ii ) ) ;
        h( :, jj, ii ) = h0;
          
    end

end 

% If the number of pixels in the ROI is uneven, the amplitude at the focus
% is interpolated from the amplitudes of the SIRs along the axis.
if odd
   
   h = intp_h( h , z );

end

h = mirror_h( h );
