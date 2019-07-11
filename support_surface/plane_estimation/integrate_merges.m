function lut = integrate_merges(n, old_labels, new_labels)
% integrate_merges : This function integrates a series of label
% aggregations into a final lookup table

lut = 1:n;

% The next array stores the pointer to the next index in the list 0 marks
% the end of the list
next = zeros(n,1);

% The tail array stores the indices of the tail of the list. For labels
% that have been merged the tail is zero.
tail = 1:n;

if (length(old_labels) ~= length(new_labels))
    error ('old_labels and new_labels not the same size');
end

%% Merge labels

nmerges = length(old_labels);

for i = 1:nmerges
    old_label = old_labels(i);
    new_label = new_labels(i);
    
    % Merge the lists
    next(tail(new_label)) = old_label;
    tail(new_label) = tail(old_label);
    tail(old_label) = 0;
end


%% Fill in lut

for i = 1:n
    if (tail(i))
        idx = i;
        
        while(idx)
            lut(idx) = i;
            idx = next(idx);
        end
    end
end


%% You can normalize entries in the lut using unique
% [~, ~, normalized_lut] = unique(lut);

% Old version - short but inefficient
% for i = 1:nmerges
%     lut(lut == old_labels(i)) = new_labels(i);
% end
