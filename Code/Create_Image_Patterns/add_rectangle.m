function new_image = add_rectangle(orig_image, x, y, width, height)

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

new_image=orig_image;
% check rectangle bounds 
[im_height,im_width]=size(new_image);
if(im_width < x+width-1)
    x_vec=x:im_width;
else
    x_vec=x:x+width-1;
end

if(im_height < y+height-1)
    y_vec=y:im_height;
else
    y_vec=y:y+height-1;
end

% add the rectangle
new_image(y_vec,x_vec)=1;

end
