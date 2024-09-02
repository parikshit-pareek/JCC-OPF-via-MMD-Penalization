%% Construction of 'h' matrix : h is the matrix having sigma^2 values for kernel calculation


% Parikshit Pareek, NTU Sg


% h = [hgg_ga [Hga ; Haa]];
% This code exploits the symmetry of the 'h' matrix to save computational time.

function [h] = h_making(const,Ak_all)
hgg_ga = ones(2*const.ngen,const.ngen); % [hgg ; hga]

hga = zeros(const.ngen,const.ngen);

for m = 1:const.ngen
    for l = 1:m
     hga(m,l) = meddistance(Ak_all(:,:,m)'*Ak_all(:,:,const.ngen+l))^2;
    end
end

for m = 1:const.ngen
    for l = 1:m
     haa(m,l) = meddistance(Ak_all(:,:,const.ngen+m)'*Ak_all(:,:,const.ngen+l))^2;
    end
end
Hga= hga+hga' - diag(diag(hga));
Haa= haa+haa' - diag(diag(haa));

h = [hgg_ga [Hga ; Haa]];
h(h == 0) = 1;
end