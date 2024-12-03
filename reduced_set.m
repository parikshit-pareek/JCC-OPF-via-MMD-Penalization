% Code for obtaining reduced set method randomly selecting vectors 

function [RS] = reduced_set(Des_DCOPF,target_size)
[big_sample_size_row,~, ~]=size(Des_DCOPF.xs_pen);
idex=randi(big_sample_size_row, round(big_sample_size_row/target_size) ,target_size)';
alp=1*ones(big_sample_size_row,1)./big_sample_size_row;
op = struct();
kerx= myProcessOptions(op, 'mmd_kernel', KGaussian(meddistance(Des_DCOPF.xs_pen)^2)); 
kxx = kerx.eval(Des_DCOPF.xs_pen',Des_DCOPF.xs_pen');
xx = Des_DCOPF.xs_pen ;
parfor j = 1:round(big_sample_size_row/target_size)
        rd_set = xx(idex(:,j),:,:); % Reduced samples for only bus 8 and 15
        op = struct();
        ker= myProcessOptions(op, 'mmd_kernel', KGaussian(meddistance(rd_set)^2)); 
        kz=ker.eval(rd_set', rd_set');
        kzx=ker.eval(rd_set', Des_DCOPF.xs_pen');
        rd_set_coeff(:,j)=(pinv(kz)*kzx*alp); % Weight of the reduced set 
        Gap(j) = alp'*kxx*alp + rd_set_coeff(:,j)'*kz*rd_set_coeff(:,j) - 2* alp'*kzx'*rd_set_coeff(:,j);
end

RS.Gap = Gap;
RS.rd_set_coeff = rd_set_coeff;
RS.idex = idex;
% Optimal RS
[RS.min_gap,j_min]= min(Gap);
RS.idex_opt = idex(:,j_min);
RS.rd_set_coeff_opt = rd_set_coeff(:,j_min);

end

