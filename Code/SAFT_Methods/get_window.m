function win = get_window(len,window_name)

% this function returns a window function of a specified type and length
%
% Inputs
% window_name =  the type of window needed, supported windows are:
% 'Hamming', 'Hann', 'rect', 'Bartlett'
% default window is 'rect'
% len = the number of samples in the window
%
%
% Outputs
% win = a window according to the specified parameters

if (~exist('window_name','var'))
   win = rectwin(len); 
else
    if (strcmp(window_name,'rect'))
       win = rectwin(len);

    elseif (strcmp(window_name,'Hann'))
        win = hann(len);

    elseif (strcmp(window_name,'Hamming'))
        win = hamming(len);

    elseif (strcmp(window_name,'Bartlett'))
        win = bartlett(len);

    else 
        win = rectwin(len);
    end

end

% normalize
%win=sqrt(len/sum(win.^2))*win;
%win=sqrt(len)*(win/sum(win));
win=len*(win/sum(win));


end