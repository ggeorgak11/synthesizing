function d2 = compute_squared_edge_distances (SN, A)
% compute_squared_edge_distances : Computes the similarity between neighboring
% triangles based on their Atttribute values stored in an ntri x ndim matrix with
% one row for each triangle. Note that the result is the square of the
% Euclidean norm between the attribute vectors of neighboring nodes in the
% Delaunay graph. We leave it as a square to avoid the sqrt which may be
% unnecessary.

missing_edges = isnan(SN);

% Fill in missing edges with 1 to avoid indexing problems
SN(missing_edges) = 1;

d2 = zeros(size(SN));

ndims = size(A,2);

for i = 1:ndims
    % Pull out the channel
    channel = A(:,i);
    
    % Nifty trick to compute the difference between a pixel and its
    % neighbors in each channel
    delta = bsxfun(@minus, channel(SN), channel);
    
    d2 = d2 + delta.^2;
end

% missing edges get zero weights
d2(missing_edges) = 0;