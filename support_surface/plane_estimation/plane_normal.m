function [n, d] = plane_normal(s)
% Extract a normal for a given plane
n = -(s.n) / norm(s.n);
d = 1 / norm(s.n);
