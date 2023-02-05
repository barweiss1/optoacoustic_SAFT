function h = Sph_F_exact(D, F, c, t, y, z ) 

% Function to calculate the SIR of a spherically focused round transducer
% at a given point. Cf. "Transient fields of concave
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
% y = point perpendicular to the axis of the transducer where the SIR 
% is to be computed.
% z = points along the axis of the transducer where the SIR is to be
% computed. The focus of the transducer is located at z == 0.
%
% OUTPUT
%
% h(t,y,z) = SIR of a spherically focused detector.


D = D/2;
r = sqrt ( y.^2 + z.^2 );

d = F.*(1 - sqrt(1 - (D./F).^2) );

theta_tmp = myasin( abs(z), y );

ct = c.*t ;

theta_tmp(isnan(theta_tmp)) = 0;
thetad = asin( D./F );


if (abs(theta_tmp) < thetad)
    R = 1 ; 
else
    R = 2 ; 
end

theta = myasin(z, y ) ;

if z <= 0
    r0 = F - r;
elseif z > 0
    r0 = F + r;
end

r1 = sqrt( ( D - (y) ).^2 + (F - d + (z) ).^2 );
r2 = sqrt( ( D + (y) ).^2 + (F - d + (z) ).^2 );

if r2 < r1
    tmp = r1;
    r1 = r2;
    r2 = tmp;
end

% Define eta(t)

    eta1 = ( 1 - d./F ) ./ (sin( theta )+eps);
    eta2 = ( F.^2 + r.^2 - ct.^2 );
    eta3 = eta2 ./ ( 2.*r.*F ) ;
   
    eta2 = eta3 ./ ( tan( theta ) + eps ) ;
    eta = eta1 + eta2 ;
    eta = F .* eta ;
    eta ( imag( eta ) ~=0 ) = 0;
   
% Define sigma(t)
    
    sigma = eta3.^2;
    sigma = sqrt( 1 - sigma ) ;
    sigma = F .* sigma ;
    sigma( imag( sigma ) ~= 0 ) = 0 ;
    
% SIR

    h1 = c.*F./r ;
    h1(isnan(h1)) = 0;
    h2 = 1/pi .* acos( eta ./ (sigma + eps) ) ;
    h2( imag(h2) ~= 0 ) = 0 ;
    h2 = h1 .* h2 ;
    h2(isnan(h2)) = 0;
    
    H1 = zeros( size( t ) ) ;
    H2 = H1 ;
    
    if ( z == 0 ) && ( y == 0 ) 
        
    H1( min( abs(ct - F) ) == ( abs( ct - F ) ) ) = d ;
    
    elseif ( z ~= 0 ) && ( y == 0 )
        
        H1 = Sph_F_exact_onaxis(D, F, c, t, z ) ;
       
    else

        if R == 1
             if z <= 0

                H1 = h1.*(( r0 < ct) .* ( ct < r1 )) ;
                H2 = h2.*(( r1 < ct) .* ( ct < r2 )) ;

             elseif z > 0

                H1 = h1.*(( r2 < ct) .* ( ct < r0 )) ;
                H2 = h2.*(( r1 < ct) .* ( ct < r2 )) ;

             end

        elseif R == 2

                H2 = h2.*(( ct < r2 ) .* ( r1 < ct )) ;
              
        end
    end

h = ( H1 + H2 ) ;

    

