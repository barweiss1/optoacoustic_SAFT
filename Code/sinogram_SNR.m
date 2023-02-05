function SNR = sinogram_SNR(scan_line,points_pos,space,time_res)

% Function recieves a single scan line of a sinogram and returns the SNR
% for each point
%
%
% INPUTS 
% scan line = time signal representing the scan line
% points_pos = vector representing the positions of the points 
% space = space (in discretized time points) between the peaks 
% time_res = time oversampling factor relatively to space
% 
% 
% OUTPUT
% SNR = SNR vector at each peak point 
% 


% convert to time resolution with time_res factor
space_t=time_res*space;
points_pos_t=time_res*points_pos;

% 
point_num = length(points_pos);
%SNR = zeros(point_num,1);
sigmas = zeros(point_num,1);
peaks = zeros(point_num,1);
len=length(scan_line);
segment_len=floor(space_t/4);


for i=1:point_num
    
   % check valid bounds for signal segment
   start_i = points_pos_t(i) - segment_len;
   if(points_pos_t(i)-segment_len < 1)
     start_i=1;
   end
   
   end_i = points_pos_t(i)+segment_len;
   if(points_pos_t(i)+segment_len > len)
     end_i=len;
   end
   
   % seperate the signal around the peak point
   sig=scan_line(start_i:end_i);
   % estimate the clean signal with a basic moving average filter
   sig_clean=movmean(sig,15);
   peaks(i)=max(abs(sig_clean)); % calculate peak value
   sigmas(i)=std(sig-sig_clean); % (sig - sig_clean) estimates the noise to compute its magnitude
    
end

SNR=mag2db(peaks./sigmas); % 20*log(snr)

end