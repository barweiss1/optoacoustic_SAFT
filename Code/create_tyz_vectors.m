function [t,y,z] = create_tyz_vectors(F,c,time_res,image_width,image_height,Nz,Ny)

% Function that creates the t,y,z vectors based on the given parameters
%
% INPUTS
% Parameter definition. All units are SI. The focus of the transducer is
% F = transducer focal length
% c = speed of sound
% time_res = number of time steps per dx 
% image_width = the width of the image
% image_height = the height of the image
% Nz = number of pixels in z direction
% Ny = number of pixels in y direction
%
% 
% OUTPUT
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


z = linspace(-image_width/2, image_width/2, Nz) ;

%Ny = floor( Nz / 2 ); % Number of pixels in the perpendicular direction.
% The SIR is calculated for the positive half plane in front of the
% transducer and then mirrored ( in the wrapper below ).

%y = z( Ny + 1 : end ) ;
y = linspace(-image_height/2, image_height/2, Ny);
y = y(floor(Ny/2) +1 :end);

%Definition of the time step and time vector
dx = z(2) - z(1) ; 
dt = dx / c ;
dt = dt / time_res ;

t1 = F -  ( sqrt( 2 ) / 2 ) * image_width ;
tend = F + ( sqrt( 2 ) / 2 ) * image_width ;

% t1 = F -  image_width/2 ;
% tend = F + image_width/2 ;

t = t1/c : dt : tend/c ;


end