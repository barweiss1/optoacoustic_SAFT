function new_image = add_circle(orig_image, x, y, radius)

% Function that adds a rectangle to an image with specified properties
%
% INPUTS
% orig_image = the image that the rectangle is added to in double
% representation
% x,y = rectangle x,y coordinates (top left corner)
% width,height = x and y rectange size respectively
%
% OUTPUT
%
% new_image = the original image with a rectangle added at the specified
% position


[im_height,im_width]=size(orig_image);

x_vec = 1:im_width;
y_vec = 1:im_height;

dist_im = (x_vec - x).^2 + (y_vec - y).^2';
circle_im = dist_im < radius^2;



% add the rectangle
%new_image = orig_image + sqrt(dist_im).*circle_im; % circle with gradient 
new_image = orig_image | circle_im;

end
