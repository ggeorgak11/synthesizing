
% Georgios Georgakis 2016

function m = sample_depth(depth, pos)

rad=10;
pos_x = pos(2);
pos_y = pos(1);
c = 0;
hypo_depth = 0;
% sample around the position for the median depth
for w = pos_x-rad : pos_x+rad
    for h = pos_y-rad : pos_y+rad
        if w<1 || h<1 || w>size(depth,1) || h>size(depth,2) || depth(w,h)==0
            continue;
        else
            c = c + 1;
            hypo_depth(c) = depth(w,h);
        end
    end
end

%if ~exist('hypo_depth','var'), keyboard, end

m = median(hypo_depth);
