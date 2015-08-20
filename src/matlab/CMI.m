function [ r ] = CMI( pxyz )
%CMI Conditional Mutual Information MI(X;Y|Z)
%   Input is a 3-dimensional array p(x,y,z)
%   Calculation is:
%   MI(X;Y|Z) = \sum p(x,y,z) log (p(x,y|z) / (p(x|z) * p(y|z)))
%   p(x,y|z)  = p(x,y,z) / p(z)
%   p(x|z)    = sum_y p(x,y,z) / p(z)
%   p(y|z)    = sum_x p(x,y,z) / p(z)
%   p(z)      = sum_{x,y} p(x,y,z)

r = 0;

eps = 0.00000001;

pz  = sum(sum(pxyz, 2), 1);
pxz = sum(pxyz, 2);
pyz = sum(pxyz, 1);

for x = 1:size(pxyz,1)
    for y = 1:size(pxyz,2)
        for z = 1:size(pxyz,3)
            if abs(pxyz(x,y,z)) > eps && abs(pz(1,1,z)) > eps && abs(pxz(x,1,z)) > eps && abs(pyz(1,y,z)) > eps
                r = r + pxyz(x,y,z) * ( log(pxyz(x,y,z) / pz(1,1,z)) - log(pxz(x,1,z) / pz(1,1,z) * pyz(1,y,z) / pz(1,1,z)));
            end
        end
    end
end

end