function [newX,newY] = rotatePath(posx,posy,tAngle)

% multiply by a rotation matrix
newX = posx * cos(tAngle) - posy * sin(tAngle); %rotate posx
newY = posx * sin(tAngle) + posy * cos(tAngle); % rotate posy
disp('')