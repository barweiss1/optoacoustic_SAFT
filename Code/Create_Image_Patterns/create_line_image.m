function new_image = create_line_image(orig_image, line_num, space, x_center, y_center, len,angle)

% Function that adds angled lines to an image with specified properties
%
% INPUTS
% orig_image = the image that the lines are added to in double
% representation
% line_num = total rectangle number added to the image
% space = space from one line to the next
% x_center,y_center = center line x,y coordinates (top left corner)
% len = length of the line in pixels
% angle = angle of the line (in degrees) with respect to y
%
% OUTPUT
%
% new_image = the original image with lines added according to the specified pattern

new_image=orig_image;

x = x_center - floor(line_num/2)*(1+space);
y = y_center;

for i=0:line_num-1
    new_image=add_angled_line(new_image,x+i*(1+space),y,len,angle);
end

end