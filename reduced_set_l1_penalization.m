%% Code for selecting redueced set using l_1 Penalization 



function [z,b_opt,rd_set,K,sol] = reduced_set_l1_penalization(solver_name,xs,lambda,K)

yalmip clear

nt_r = size(xs,1);
if isempty(K)
    op = struct();
    ker= myProcessOptions(op, 'mmd_kernel', KGaussian(meddistance(xs)^2)); 
    K=ker.eval(xs', xs');
end
alp=1*ones(nt_r,1)/nt_r;

K_alp = K*repmat(alp,[1,nt_r]);

b_p = sdpvar(nt_r,1);
b_n = sdpvar(nt_r,1);
cj = ones(nt_r,1);

obj_1 = [b_p-b_n]'*K*[b_p-b_n];
obj_2 = lambda*cj'*(b_p+b_n) - 2 * sum(K_alp)*[b_p-b_n];
cons = [b_p >= 0; b_n>=0];


ops = sdpsettings('solver',solver_name,'verbose',0);
sol = optimize(cons,obj_1+obj_2,ops);


b = round(value(b_p-b_n),3);

z = b~=0;

rd_set = xs(z,:); 
op = struct();
ker= myProcessOptions(op, 'mmd_kernel', KGaussian(meddistance(rd_set)^2)); 
kzx = ker.eval(rd_set', xs');
kz=ker.eval(rd_set', rd_set');
b_opt = (pinv(kz)*kzx*alp)./nt_r; % Weight of the reduced set 

b_opt = b_opt/sum(b_opt);








%% --------------------------------------------------------------------------- Local COdes
end