function [ normSpxl ] = fitPlaneAffine( dxc, dyc, dzc )
%
%

%         pxls = find(lim == k); 
%         % get the depth coordinates of these pxls
%         dxc = dx(pxls);
%         dyc = dy(pxls);
%         dzc = dz(pxls);
        samp_thresh = 500; % 1000

        sp_valid = find(dzc ~= 0);
        sample = 1;
        if length(sp_valid) > 1000
            sample = round(length(sp_valid)/samp_thresh); % reduce this to get more sampled points
        end
        
        indval = sp_valid(1:sample:end);
        ws = 1./dzc(indval);
        us = dxc(indval).*ws;
        vs = dyc(indval).*ws;
        A2 = [us vs ones(length(indval),1) ws];

%            A2 = [dxc(indval) dyc(indval) dzc(indval) ones(length(indval),1)];

        [u,s,v] = svd(A2);
        normSpxl = zeros(1,4);
        normSpxl(1:4) = v(:,end)';
        conf = cond(A2); % s(4,4);
        
        nrm = norm(normSpxl(1:3));
        normSpxl(1:3) = normSpxl(1:3)/nrm;
        normSpxl(4) = normSpxl(4)/nrm;


end

