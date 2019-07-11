function [overlap_ratios, overlap_areas, edge_lengths] = ComputeCircumcircleOverlap (SN, CC, RCC)
% ComputeCircumcircleOverlap : Computes the level of overlap between the
% circumcircles of adjacent triangles.
% Inputs: 
%       SN - ntri x 3 indicating triangle adjacency using NaNs for
%       neighbors that don't exist
%
%       CC - ntri x 2 array of circumcircle centers
%
%       RCC - ntri x 1 vector of circumcircle radii

overlap_areas  = zeros(size(SN));
overlap_ratios = zeros(size(SN));
edge_lengths   = zeros(size(SN));

% For each neighbor
for i=1:3
    index = SN(:,i);
    
    % Assign the NaN indices some real value these values won't matter
    % anyway
    index(isnan(index)) = 1;
    
    % Get the radii of the adjacent circumcircle
    r2 = RCC(index);
    
    % r <= R
    r = min(RCC, r2);
    R = max(RCC, r2);
    
    % Compute distance between the circumcircle centers
    delta = CC - CC(index,:);
    
    d = sqrt(sum(delta.^2,2));
    
    % Compute angles using Cosine rule
    c1 = ((R.^2 + d.^2) - r.^2) ./ (2 * R .* d);
    c2 = ((r.^2 + d.^2) - R.^2) ./ (2 * r .* d);
    
    % Apply clipping to c1 and c2 to eliminate crazy values caused by
    % numerical errors - including situations where d == 0
    % Note that c1 should always be positive but c2 may be -ve
    c1 = max(-1, min(c1, 1));
    c2 = max(-1, min(c2, 1));
    
    % Compute theta1 and theta2
    theta1 = acos(c1);
    theta2 = acos(c2);
    
    % Compute sines from the cosines
    s1 = sqrt(1 - c1.^2);
    s2 = sqrt(1 - c2.^2);
    
    % The following expression should work even in cases where c2 < 0 =>
    % theta2 > (pi/2). In this case the triangle area is added not
    % subtracted.
    overlap_areas(:,i) = (theta1 - (s1.*c1)).*(R.^2) + (theta2 - (s2.*c2)).*(r.^2);
    
    % Overlap ratio - ratio of overlap area to area of the smaller circumcircle.
    overlap_ratios(:,i) = overlap_areas(:,i) ./ (pi * r.^2);
    
    % Edge length in pixels 2 * R * sin(theta1)
    edge_lengths(:,i) = 2*(s1.*R);
    
    % Find the cases where the two circumcenters are on the same side of
    % the edge. This is signaled by c2 < 0. We set the edge_length -ve to
    % signal this case.
    % t = (c2 < 0);
    % edge_lengths(t, i) = -edge_lengths(t, i);

    % Handle the case where d=0. This can happen because of integral edgel
    % coordinates and it indicates that the two circumcircles are
    % coincident so the overlap ratio is 1.
    
    t = (d == 0);
    
    overlap_areas(t,i)  = pi * (R(t)).^2;
    overlap_ratios(t,i) = 1.0;
    edge_lengths(t,i)   = 2*R(t);
end