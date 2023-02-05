function new_image = add_angled_line(orig_image, x, y, len,angle)

% Function that adds a rectangle to an image with specified properties
%
% INPUTS
% orig_image = the image that the rectangle is added to in double
% representation
% x,y = line x,y coordinates (top left corner)
% len = length of the line in pixels
% angle = angle of the line (in degrees) with respect to y
%
% OUTPUT
%
% new_image = the original image with a rectangle added at the specified
% position

new_image=orig_image;
% check line bounds 
[im_height,im_width]=size(new_image);

angle_rad = (angle/180)*pi; % convert to radians for the tan function
height = floor(len * cos(angle_rad));

y_start = y - round(height/2);
y_end = y + round(height/2);
if(im_height+1 < y_end)
    y_end=im_height;
end
if(1 > y_start)
    y_start=1;
end
y_vec = y_start:y_end;

x_vec = round((y_vec - y)*tan(angle_rad) + x);


% add the line
for i = 1:length(y_vec)
    new_image(y_vec(i),x_vec(i))=1;
end

end
