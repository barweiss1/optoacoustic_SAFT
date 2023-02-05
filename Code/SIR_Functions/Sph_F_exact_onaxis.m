function h = Sph_F_exact_onaxis(D, F, c, t, z ) 

% Function to calculate the SIR of a spherically focused round transducer
% at a given point z along its axis. Cf. "Transient fields of concave
% annular arrays", Arditi, M. et al., Ultrasonic Imaging 3, 37-61 (1981)
% for details on the calculation.
%
% INPUT (all units in SI)
%
% D = transducer diameter
% F = transducer focal length
% c = speed of sound
% t = time vector in which the SIR is calculated. It should be bigger than
% (see following) the ROI
% z = points along the axis of the transducer where the SIR is to be
% computed. The focus of the transducer is located at z == 0.
%
% OUTPUT
%
% h(t,z) = SIR of a spherically focused detector.

d = F.*(1 - sqrt(1 - (D./F).^2) );

r1 = sqrt( D.^2 + (F - d + z ).^2 );
r2 = r1;
r = abs(z);

ct = c.*t ;

if z<=0
    r0 = F - r;
elseif z>0
    r0 = F + r;
end

if ( z == 0 ) 
    
    h = zeros(size(ct)) ;
    h( min( abs(ct - F) ) == ( abs( ct - F ) ) ) = d ;
    
else
    M = ( r0 + r1 )/2;

    Az = r1 - r0;

    up = ct - M;

    % argum = abs( up ./Az );

    rect = zeros(size(t));

    rect( abs( up./Az )  > 1/2 ) = 0;
%     rect( abs( up./Az )  == 1/2 ) = .5;
    rect( abs( up./Az )  < 1/2 ) = 1;
    
    if nnz(rect) == 0
        rect(min(abs(up./Az)) == abs(up./Az)) = 1;
    end

    h = c .* F ./ ( abs( z ) ) .* rect ; 
end