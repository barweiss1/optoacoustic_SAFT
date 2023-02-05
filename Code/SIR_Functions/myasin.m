function as = myasin( x, y ) 
% Auxiliary function to make sure that the arcsine of the cartesian
% coordinates (x,y) lies within the interval [0º , 360º ).

r = sqrt( x.^2 + y.^2 );
as = asin( y./r );
if (x < 0) && (y > 0)
    as = pi - as ;
elseif (x < 0) && (y < 0 )
    as = pi + abs(as) ;
elseif (x > 0) && (y < 0 );
    as = 2*pi - abs(as) ;
end