function new_image = create_circle_image(orig_image, circle_num, space, x_center, y_center, radius)

% Function that adds a circles to an image with specified properties
%
% INPUTS
% orig_image = the image that the rectangles are added to in double
% representation
% rect_num = total rectangle number added to the image
% space = space from one circle to the next
% x_center,y_center = center x,y coordinates
% radius = radius of the circles
%
% OUTPUT
%
% new_image = the original image with circles added according to the specified pattern

new_image=orig_image;

x = x_center - floor(circle_num/2)*(radius+space);
y = y_center;

for i=0:circle_num-1
    new_image=add_circle(new_image,x+i*(radius+space),y,radius);
end

end