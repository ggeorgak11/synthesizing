function n = fitPlaneLSQ(uvw)
% Fit a plane to the given subset in a least squares sense
% w = n(1)*u + n(2)*v + n(3)

A = uvw;
b = uvw(:,3);
A(:,3) = 1;

n = A \ b;

end