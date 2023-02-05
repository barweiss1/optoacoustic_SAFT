function new_image = create_rect_image(orig_image, rect_num, space, x_center, y_center, width, height)

% Function that adds a rectangle to an image with specified properties
%
% INPUTS
% orig_image = the image that the rectangles are added to in double
% representation
% rect_num = total rectangle number added to the image
% space = space from one rectangle to the next
% x,y = first (highest) rectangle x,y coordinates (top left corner)
% width,height = x and y rectange size respectively
%
% OUTPUT
%
% new_image = the original image with rectangles added according to the specified pattern

new_image=orig_image;

x = x_center - floor(rect_num/2)*(width+space);
y = y_center;

for i=0:rect_num-1
    new_image=add_rectangle(new_image,x+i*(width+space),y,width,height);
end

end