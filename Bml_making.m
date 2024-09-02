%% Construction of Bml matrix : Matrix for quadratic term of MMD
% Quad Term : g'*Bml*g;


% Parikshit Pareek, NTU Sg

% This code exploits the symmetry of the 'Bml' matrix to save computational time.
% Bml = [B_gg B_ga; B_ag B_aa];

% Symmetric Submatrics: B_gg and B_aa
% B_ga : Non-symmetric (zero in some systems (118,39) not zero in some systems (14,30). Why ?)
% B_ag : Non-symmetric

function [Bml] = Bml_making(h,const,Ak_all,rd_set_coeff)
B_ga = zero(const.ngen,const.ngen); % [hgg ; hga]
B_ag = zero(const.ngen,const.ngen); % [hgg ; hga]

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

Bml = [B_gg B_ga; B_ag B_aa];
end