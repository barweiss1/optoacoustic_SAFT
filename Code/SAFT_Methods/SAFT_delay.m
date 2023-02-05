function saft_delay = SAFT_delay(s,t,y,c,F,D)

% Function implements just the delay part of the saft algorithm 
%
% INPUTS (all units are in SI)
% s(y,t) = the original sinogram
% t = time vector of the sinogram 
% y = represents the recievers axis in the sinogram
% c = speed of sound
% F = transducer focal length
% D =  Diameter of the transducer
% 
% 
% 
% OUTPUT
% saft = the sinogram after the SAFT delay only correction for each
% detector line
% 

[Y,T] = size(s);

saft_delay=zeros(Y,Y,T);

% calculate distances for the delay vector
dt=t(2)-t(1);
dy=y(2)-y(1);
z=t*c;
r=abs(z-F); % tune dimensions to fit


%implement on each t at a time the SAFT delay and sum
for k =1:Y % for each detector position we create the delay image
    
    saft_delay(k,:,:)=s; % initialize to no delay
    for i=1:T


        % calculte N (the number of lines to sum up) based on the the distance
        % from the focus 
        D_eff=D*abs(F-z(i))/F; % this is the width of the aperture in point z
        N=ceil(D_eff/dy); % discretize to get number of detector lines 

        %calculate delay vector for given time
        y_dist=y(1:N)-y(1); % vector of distances in y axis (x in the paper) from the VD
        r_prime=(r(i).^2+y_dist.^2).^(1/2); % sythesized point distance from VD
        delta_t=sign(z(i)-F)*((r(i)-r_prime)/c); % delays calculated as defined in the paper
        delta_t_index=round(delta_t/dt); % discretize delays and convert to samples 
      
        for j=2:N

            % check that the delay is legal
            if(i - delta_t_index(j) < 1 || i - delta_t_index(j) > T )
                continue;
            end
            % if in bounds delay appropriately 
            if( j+k < Y )
                saft_delay(k,j+k,i)=s(j+k,i-delta_t_index(j));
            end
            if( k-j > 1 )
                saft_delay(k,k-j,i)=s(k-j,i-delta_t_index(j));
            end


        end

    end
end

end