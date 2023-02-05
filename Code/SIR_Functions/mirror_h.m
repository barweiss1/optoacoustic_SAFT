function h2 = mirror_h( h )

% Auxiliary function to mirror the SIR to the y<0 half-plane.

% Extract dimensions of the SIR
len = size( h, 1 ) ;
Ny = size( h, 2 );
Nz = size( h, 3 );
odd = rem(Nz,2);

% Define the dimensions of the new SIR
if odd
    h2 = zeros( len, 2*Ny - 1, Nz );
else
    h2 = zeros( len, 2*Ny, Nz );
end

% Mirror the SIR, tim-sample by time-sample.
for ii = 1:len
    
    tmp = squeeze( h(ii, : , : ) ) ;
    
    if odd
        
        tmp2 = tmp( :, 2:Ny ) ;
        tmp = [ fliplr(tmp2) tmp(:,1) tmp2 ] ;
        
    else
        
        tmp = [ flipud(tmp); tmp ] ;
        
    end
    
    h2(ii, : , : ) = tmp ;
    
end

