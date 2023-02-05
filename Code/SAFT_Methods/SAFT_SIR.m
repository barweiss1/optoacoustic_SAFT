function saft = SAFT_SIR(s,t,y,c,F,D,window_name)

% Function implements SAFT with a use of a normalized window for the
% summation, window type can also we chosen
%
%
% The number of elements in the summation changes with the depth based on
% the transducer's natural aperture
%
% INPUTS (all units are in SI)
% s(y,t) = the original sinogram
% t = time vector of the sinogram 
% y = represents the recievers axis in the sinogram
% c = speed of sound
% F = transducer focal length
% D =  Diameter of the transducer
% window_name =  the type of desired window , supported windows are:
% 'Hamming', 'Hann', 'rect', 'Bartlett'
% 
% 
% OUTPUT
% saft = the sinogram after the SAFT correction
% 
saft=zeros(size(s));

[Y,T] = size(s);

% calculate distances for the delay vector
dt=t(2)-t(1);
dy=y(2)-y(1);
z=t*c;
r=abs(z-F); % tune dimensions to fit
N_min =6;

%implement on each t at a time the SAFT delay and sum
for i=1:T
    
    
    % calculte N (the number of lines to sum up) based on the the distance
    % from the focus 
    D_eff=D*abs(F-z(i))/F; % this is the width of the aperture in point z
    N=round(D_eff/dy)+N_min; % discretize to get number of detector lines 
    % check if N is valid
    if (N > Y/2)
      N = floor(Y/2); 
    end
    M=2*N+1; % total number of points in the summation - N above, N below and the 1 in the center
    win=get_window(M,window_name); % get window type for the filtering 
    
    % add central summation line with delay
    mid_index = N+1; % index for the middle of the window
    saft(:,i)=win(mid_index)*s(:,i);
    
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
        % SAFT delay and sum
        % add contribuitions from rows above
        saft(j:end,i)=saft(j:end,i)+win(mid_index+j)*s(1:end-j+1,i-delta_t_index(j));
        % add contribuitions from rows below
        saft(1:end-j+1,i)=saft(1:end-j+1,i)+win(mid_index-j)*s(j:end,i-delta_t_index(j));
        
        
    end
    
end


end