function h = intp_h( h, z )

% Auxiliary function to interpolate the amplitude of the SIR at focus.

Nz = size( h, 3 );

% Find the axis of the transducer.
field = squeeze( max( h, [], 1 ) );
ax = field( : , 1 );
z0 = ceil( Nz / 2 );

% Divide the SIR amplitudes between the ones before and after the focus.
ax1 = ax( 1 : (z0 - 1) );
z1 = z( 1 : (z0 - 1) );

ax2 = ax( (z0 + 2) : Nz );
z2 = z( (z0 + 2) : Nz );

% Find the interpolant for the amplitudes for z < F and z > F.
p1 = fit( z1(:) , ax1(:) , 'cubicspline' );
p2 = fit( z2(:) , ax2(:) , 'cubicspline' );

% Extrapolate the amplitude at the focus from both sides.
z1_0 = p1( z( z0 ) );
z2_0 = p2( z( z0 ) );

% Average the two.
z0_0 = (z1_0 + z2_0) / 2 ;

% Substitute amplitude value.
h0 = h( :, z0, 1 ) ;
h0 = h0/max( h0 ) .* z0_0 ;
h( :, z0, 1 ) = h0 ; 



