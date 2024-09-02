

clear
clc
load('test_data_weibull.mat')

parfor i = 1:10


nt_r_test = 10000;

[~,x_test] = xs_generate_DG(data,rbus,nt_r_test,length(rbus),rated_cap,dist_type);
[Resa{i}] = sample_approximation(data,Des_DCOPF,x_test,solver_name,rbus,rated_cap,dist_type);

end

function [SA_DC_CCOPF] = sample_approximation(data,Des_DCOPF,x_test,solver_name,rbus,rated_cap,dist_type)

% solver_name = 'mosek';

tic
[SA_DC_CCOPF] = desired_dist_DCOPF_sample(data,Des_DCOPF,solver_name,rbus,rated_cap,dist_type);

SA_DC_CCOPF.time = toc;

%% Testing via Samples 
nt_r_x = 10000;
g_x = SA_DC_CCOPF.g;
for j= 1:nt_r_x
[A1_test(:,:,j),Ao_test(:,j)] = dcopf_a1ao_dhatr(x_test(j,:)',data,rbus);
 F_x(:,j)   = A1_test(:,:,j)*g_x + Ao_test(:,j);
end
 
F_x > 10^(-3);
SA_DC_CCOPF.vio_individual = sum(ans,2);
SA_DC_CCOPF.eps_individual = SA_DC_CCOPF.vio_individual/nt_r_x;
[SA_DC_CCOPF.max_eps_individual,SA_DC_CCOPF.max_c] = max(SA_DC_CCOPF.eps_individual);


SA_DC_CCOPF.vio_joint = sum(ans);
SA_DC_CCOPF.vio_idx = SA_DC_CCOPF.vio_joint >=1;
SA_DC_CCOPF.eps_joint = (sum(SA_DC_CCOPF.vio_idx)/nt_r_x);




function[SA_DCOPF] = desired_dist_DCOPF_sample(data,Des_DCOPF,solver_name,rbus,rated_cap,dist_type)
% data = case30;
% base =100;
const = ex_extract_ccDCOPF(data); %% Extract information from the MPC structureconst.Hr = const.H(:,rbus);
const.Hr = const.H(:,rbus);

SA_DCOPF.target_size = round((2/0.05)*(log(1/10^(-4))+2*const.ngen),0); % Number of samples used for SAA
[~,xs] = xs_generate_DG(data,rbus,SA_DCOPF.target_size,length(rbus),rated_cap,dist_type);



total_var = Des_DCOPF.total_var;
%% Formulate the cc-DCOPF Problem via SAA
% Detailed Formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% min_{g,s,eta} c(g)    
% s.t. sum(g) = sum(d)
%      f_l <= Hg*g - Hd*(d_hat) <= f_u
%      g_l <= g <= g_u
%      sum(eta) = 1
%      Prob( f_l <= Hg*[g+sum(d_err)*eta] - Hd*(d_hat+d_err) <= f_u,
%           and g_l <= g+sum(d_err)*eta <= g_u) >= 1 - epsilon
%  Here s: Load shed variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Formulate the Problem Using YALMIP and CCC
g = sdpvar(const.ngen*2,1,'full'); % generation and participation factor [pg;eta]
% eta = sdpvar(const.ngen,1,'full'); % affine corrective control policy


objective = (g(1:const.ngen))'*const.Q*(g(1:const.ngen)) + const.c_g'*(g(1:const.ngen))+sum(const.c_o)...
      + g(const.ngen+1:end)'*(const.Q*total_var*100)*g(const.ngen+1:end);

constr_det = [sum(g(1:const.ngen)) == sum(const.d_hat);
    const.f_l <= const.Hg*(g(1:const.ngen)) - const.Hd*(const.d_hat) <= const.f_u;
    const.g_l <= g(1:const.ngen) <= const.g_u;
%      zeros(const.nload,1) <= s <= const.d_hat;
    sum(g(const.ngen+1:end)) == 1; 
    zeros(const.ngen,1) <= g(const.ngen+1:end) <= ones(const.ngen,1)];

% Sample based constraints for total N samples of xi

const_sample = [const.g_l <= g(1:const.ngen) + sum(xs(1,:))*g(const.ngen+1:end) <= const.g_u; % Generator Limits
               const.f_l <= const.Hg*(g(1:const.ngen)+...
               sum(xs(1,:))*g(const.ngen+1:end)) - const.Hd*(const.d_hat)- const.Hr*xs(1,:)' <= const.f_u]; 
for s = 2:length(xs(:,1))
 const_sample = [ const_sample; 
               const.g_l <= g(1:const.ngen)+ sum(xs(s,:))*g(const.ngen+1:end) <= const.g_u; % Generator Limits
               const.f_l <= const.Hg*(g(1:const.ngen)+...
               sum(xs(s,:))*g(const.ngen+1:end)) - const.Hd*(const.d_hat)- const.Hr*xs(1,:)' <= const.f_u]; % Power Flow limits
end

ops = sdpsettings('solver',solver_name,'verbose',0);

sol = optimize([constr_det;const_sample],objective,ops);

SA_DCOPF.sol = sol;
SA_DCOPF.obj=value(objective);
SA_DCOPF.pg = value(g(1:const.ngen));
SA_DCOPF.eta = value(g(const.ngen+1:end));
SA_DCOPF.g = value(g);




end

function [A_1,A_o]= dcopf_a1ao_dhatr(xi,data,rbus)
const = ex_extract_ccDCOPF(data); %% Extract information from the MPC structure
const.Hr = const.H(:,rbus);
I = eye(const.ngen,const.ngen);
Mg_up = [I sum(xi)*I];
A_1 = [Mg_up;-Mg_up;const.Hg*Mg_up;-const.Hg*Mg_up]; % Mg_down = -Mg_up;
A_o = [-const.g_u;...%-const.g_u-I*(const.d_hatg+xi(const.genbuses));
        const.g_l;...
       -const.Hd*(const.d_hat)-const.Hr*xi-const.f_u;...
        const.Hd*(const.d_hat)+const.Hr*xi+const.f_l];

end
end
