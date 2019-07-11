
function show_seg(im, depth, seg)

load('classMapping40.mat'); % load class mapping of nyu % loads className

subplot(1,3,1);
imagesc(im); 
axis image;

subplot(1,3,2);
imagesc(depth);
axis image;

subplot(1,3,3);
imagesc(seg);
axis image;

lbls = unique(seg);
lbls = lbls + 1;

for j=1:size(lbls,1)
    %[lbls(j) className{lbls(j)}]
    fprintf('%d %s\n', lbls(j)-1, className{lbls(j)});
end
disp(' ');

