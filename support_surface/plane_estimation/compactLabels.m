function out = compactLabels (labels)
% compactLabels : takes an input label image and reorders the labels to
% avoid gaps in the sequence. These gaps may occur due to rasterization
% where some triangle labels don't show up due to triangle fighting.

% Use unique to compute the unique index for each pixel
[~, ~, labels2] = unique(labels(:));

out  = reshape(labels2, size(labels));