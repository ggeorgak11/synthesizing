function out = colorImageSegments(Irgb, labels, edge_color)
% colorImageSegments : Takes a labeled image and assigns each segment the
% mean color in that segment. The edge_color is used to label edge pixels.
% it should be a 3x1 vector of floats. The output will be uint8

Irgb = double(Irgb)/255;

red    = Irgb(:,:,1);
green  = Irgb(:,:,2);
blue   = Irgb(:,:,3);

count = accumarray (labels(:), 1);

red_sums   = accumarray (labels(:), red(:));
green_sums = accumarray (labels(:), green(:));
blue_sums  = accumarray (labels(:), blue(:));

lut = [red_sums./count, green_sums./count, blue_sums./count];

lut = [edge_color; lut];

% Find all of the places that are near label changes
edges = (labels ~= imerode(labels, ones(3))) | (labels ~= imdilate(labels, ones(3)));

labels = labels+1;
labels(edges) = 1;

% This works but label2rgb may complain about lut entries duplicating the
% zerocolor of [1 1 1]
% out = label2rgb(labels, lut);

% Note that if the first argument to ind2rgb is not floating point it adds
% a 1 which screws things up.

out = ind2rgb(double(labels), lut);

end

