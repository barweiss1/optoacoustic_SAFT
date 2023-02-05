function FWHM = sinogram_FWHM(scan_line,points_pos,space,time_res,dz)

% Function recieves a single scan line of a sinogram and returns the FWHM
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
% FWHM = FWHM vector at each peak point 
% 


% convert to time resolution with time_res factor
space_t=time_res*space;
points_pos_t=time_res*points_pos;

% 
point_num = length(points_pos);
FWHM = zeros(point_num,1);
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
   % find the maximum in absolute value and calculate the half maximum
   peak_abs = max(abs(sig));
   peak_sign = 2*( (max(sig) > abs(min(sig))) - 0.5); % get the sign of the maximum for the calculation
   left_index = find(sig*peak_sign >= peak_abs/2,1 ,'first');
   right_index = find(sig*peak_sign >= peak_abs/2,1 ,'last');
   FWHM(i)=(right_index-left_index)*dz/time_res;
   
end


end