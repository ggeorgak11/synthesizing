

function seg = segment_object(img, mask)

seg = zeros(size(img,1), size(img,2), size(img,3));
seg1 = zeros(size(img,1), size(img,2));
seg2 = zeros(size(img,1), size(img,2));
seg3 = zeros(size(img,1), size(img,2));

r=img(:,:,1); g=img(:,:,2); b=img(:,:,3);
ind = find(mask>0);

seg1(ind) = r(ind);
seg2(ind) = g(ind);
seg3(ind) = b(ind);

seg(:,:,1) = seg1;
seg(:,:,2) = seg2;
seg(:,:,3) = seg3;
%figure; imagesc(uint8(seg));

