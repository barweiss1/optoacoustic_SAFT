function saft = pointwise_SAFT(s,t,y,c,F,D,window_name)

% Function implements basic SAFT on a given sinogram according to the introduction of
% "Delay-multiply-and-sum-based synthetic aperture focusing in
% photoacoustic microscopy" ,Jongin Park,a Seungwan Jeon,a Jing Meng,b Liang
% Song,c Jin S. Lee,a and Chulhong Kima
% This implementation is easier to debug since we can see the summed signal
% and mask
%
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
r=abs(z-F); 

% calculte N (the number of lines to sum up)  
N=floor(0.5*D/dy); % discretize to get number of detector lines 
if (N > Y/5)
    N = floor(Y/5); 
end
M=2*N+1; % total number of points in the summation - N above, N below and the 1 in the center
win=get_window(M,window_name); % get window type for the filtering 


%implement on each t at a time the SAFT delay and sum
for i=1:T
    
    
    %calculate delay vector for given time
    y_dist=y(1:N)-y(1); % vector of distances in y axis (x in the paper) from the VD
%    y_dist=dy*s(0:N-1);
    r_prime=sqrt(r(i).^2+y_dist.^2); % sythesized point distance from VD
    delta_t=sign(z(i)-F)*((r(i)-r_prime)/c); % delays calculated as defined in the paper
    delta_t_index=round(delta_t/dt); % discretize delays and convert to samples 
    dt_index_sym = [delta_t_index(end:-1:1) 0 delta_t_index]; % symmetric version of delta_t_index
    
    for yy=1:Y
        
        
        %intialize delayed vector
        delay_sig = zeros(M,1);
        delay_mask = zeros(size(s));
        for jj = -N:N 
            
            % check bounds of delay and signal
            if(i - dt_index_sym(jj+N+1) < 1 || i - dt_index_sym(jj+N+1) > T || yy+jj > Y || yy+jj < 1)
                 continue;
            end
            % save the required point 
            delay_sig(jj+N+1)=s(yy+jj,i-dt_index_sym(jj+N+1));
            delay_mask(yy+jj,i-dt_index_sym(jj+N+1))=1;
        end
        
        % create mask for debugging summation - add breaking point for
        % debug and run below line
        % delay_mask_s = imfuse(im_double_norm(s),delay_mask,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
        % sum all delayed points
        saft(yy,i)=sum(win.*delay_sig);
    end
    
    
end



end