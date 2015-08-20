function [ p ] = sample3d(xyz)
%SAMPLE Summary of this function goes here
%   Detailed explanation goes here

x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);

x_max = max(x);
y_max = max(y);
z_max = max(z);

p = zeros(x_max, y_max, z_max);
    
for i=1:length(x)
    xi = x(i,1);
    yi = y(i,1);
    zi = z(i,1);
    p(xi, yi, zi) = p(xi, yi, zi) + 1.0;        
end

p = p / double(length(x));

end

