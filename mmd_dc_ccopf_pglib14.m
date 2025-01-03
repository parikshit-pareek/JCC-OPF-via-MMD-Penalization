%% This is the final code to be used for the MMD Penalization based JCC-OPF. 

% Works for Pglib 14 and 57 bus test cases and comparion
% Significant Modification: 1) JCC-OPF as vectorized constraints to allow tractability
%                           2) Decision variable sagragation of MMD to allow the convergence
%                           3) DCOPF now has both generator bus load and load bus node is uncertain



% Parikshit 
% 3 Sept 2021


% Study : 
%         A) Effect of weights
%         b) Effect of number of samples


%-----------------------------------------------------------------------------------------------------------
clear
clc

format short g


data = ext2int(case24_ieee_rts);
% data.gen(11,9) = 170;
solver_name = 'mosek';
dist_type = 'normal';
frac_pen = 0.01;
w_j = 1; 
w_mmd_pg = 0.0002;
w_mmd_alpha = 50;
target_size = 100; % Number of samples used for mmd construction 


[Det_DCOPF] = deterministic_DCOPF(data,solver_name);
Res.Det_DCOPF = Det_DCOPF;
% DC_det = rundcopf(data);

nt_r = 100; % Total number of samples nt_r >>>>> N



[Des_DCOPF] = desired_dist_DCOPF_local(data,target_size,nt_r,frac_pen,dist_type,solver_name);
Res.Des_DCOPF = Des_DCOPF;
if Des_DCOPF.sol.problem  ~=0
    return
end
tic

[z,rd_set_coeff,rd_set,sol] = reduced_set_l1_penalization_local(solver_name,Des_DCOPF.xs_pen,8,Des_DCOPF.K);
target_size = length(rd_set(:,1));

%%  Formulating the CC-DCOPF constraints as MMD  
const = ex_extract_ccDCOPF(data); %% Extract information from the MPC structure
% ----------------------------------------------------------------------------------------------------------------------
% Calulating Value of F_star = A1 * u_star + Ao; at different DCOPF samples
% F_star = zeros(2*const.ngen+2*const.nline,target_size);
xi_all = Des_DCOPF.xs_pen;
g_dc = Des_DCOPF.g;
for j =1:target_size
    [A1_dc,Ao_dc] = dcopf_a1ao_dhatg(xi_all(j,:)',data);
    F_star(:,j) = A1_dc * g_dc  + Ao_dc;  % All values must be less than zero to have a correct feasible DCOPF solution
end
% ----------------------------------------------------------------------------------------------------------------------
% Construction of constraint function in matrix form as F: A_1*u + A_o
% A1_rd_set = zeros(2*const.ngen+2*const.nline,2*const.ngen,target_size);
% Ao_rd_set = zeros(2*const.ngen+2*const.nline,target_size);
for j= 1:target_size
[A1_rd_set(:,:,j),Ao_rd_set(:,j)] = dcopf_a1ao_dhatg(rd_set(j,:)',data);
end
t = size(A1_rd_set,1);
for k = 1:2*const.ngen
    Ak_all(:,:,k) = reshape(A1_rd_set(:,k,:),[t,target_size]);
end
%% ----------------------------------------------------------------------------------------------------------------------
% Construction of MMD objective: <mu_d - mu(g)>
% <mu_d - mu(g)> = < mu_d , mu_d > -2 <mu_d , mu(g)  > + < mu(g) , mu(g) >
% <mu_d - mu(g)> = A + B + C
% ----------------------------------------------------------------------------------------------------------------------
% < mu_d , mu_d > : Inner product of DCOPF constraint set (F_star) in RKHS 
% mmd_A = < mu_d , mu_d >
op = struct();

kerF= myProcessOptions(op, 'mmd_kernel', KGaussian(meddistance(F_star)^2)); 
K_FsFs=kerF.eval(F_star, F_star);
mmd_A = Des_DCOPF.desired_beta' * K_FsFs * Des_DCOPF.desired_beta ;  
% ----------------------------------------------------------------------------------------------------------------------

% 2 <mu_d , mu(g)>: Inner product of DCOPF constraint set (F_star) and parametrized set mu 
%  < mu_d , mu(g) > = < mu_d , mu_Ao > + < mu_d , mu_A1) > = mmd_B_l(g) + mmd_B_c
desired_beta = Des_DCOPF.desired_beta;
hd = ones(1,2*const.ngen);
parfor k =const.ngen+1:2*const.ngen
 Ak = reshape(A1_rd_set(:,k,:),[size(A1_rd_set,1),target_size]);
hd(k) =meddistance(F_star'*Ak)^2;
end
hd(hd==0)=1;

%% The Ak cooresponding to Pg valus is same becoz uncertainty donot affect the unit commitment. Thus, mmd_B_1 = 0 for Pg values
 mmd_B_1 = zeros(2*const.ngen,1); 
 s = size(A1_rd_set,1);
parfor k =const.ngen+1:2*const.ngen % Only for alpha values need to be calculated
Ak = reshape(A1_rd_set(:,k,:),[s,target_size]);
kerFAk = myProcessOptions(op, 'mmd_kernel', KGaussian(hd(k))); 
K_FsAk = kerFAk.eval(F_star, Ak);
mmd_B_1(k,1) = desired_beta' * K_FsAk *rd_set_coeff; % mmd_B_l = mmd_B_1'*g
end

% mmd_B_2
kerFAo = myProcessOptions(op, 'mmd_kernel', KGaussian(meddistance(F_star'*Ao_rd_set)^2)); % this is important to calculate correctly
K_FsAo = kerFAo.eval(F_star, Ao_rd_set);
mmd_B_c = Des_DCOPF.desired_beta' * K_FsAo * rd_set_coeff;

%%  < mu(g) , mu(g) > = Quadra  term + Linear term + Constant term = mmd_C_q + mmd_C_l + mmd_C_c;
% ---------------------------------------------------- Quadratic  Term -----------------------------------------------
tau_K_AkAl_tau = zeros(2*const.ngen,2*const.ngen);
h = h_making(const,Ak_all);
  

for k =1:const.ngen
  for l = 1:const.ngen
kerAkl = myProcessOptions(op, 'mmd_kernel', KGaussian(h(k,l))); % this is important to calculate correctly
% B_ml = kerAkl.eval(Ak_all(:,:,k), Ak_all(:,:,l));
tau_K_AkAl_tau (k,l) = rd_set_coeff' * kerAkl.eval(Ak_all(:,:,k), Ak_all(:,:,l)) * rd_set_coeff;
  end
end

for k =1:2*const.ngen
  for l = const.ngen+1:2*const.ngen
kerAkl = myProcessOptions(op, 'mmd_kernel', KGaussian(h(k,l))); % this is important to calculate correctly
% B_ml = kerAkl.eval(Ak_all(:,:,k), Ak_all(:,:,l));
tau_K_AkAl_tau (k,l) = rd_set_coeff' * kerAkl.eval(Ak_all(:,:,k), Ak_all(:,:,l)) * rd_set_coeff;
  end
end
% mmd_C_q = g'*tau_K_AkAl_tau*g;

tau_K_AoAk_tau = zeros(2*const.ngen,1);
% % ----------------------------------------------------  Linear term
parfor k =const.ngen+1:2*const.ngen % Needed only for alpha as Ak for Pg is same for all samples, K_AoAk is zero for Pg
kerAoAk = myProcessOptions(op, 'mmd_kernel', KGaussian(meddistance(Ao_rd_set'*Ak_all(:,:,k))^2)); % this is important to calculate correctly
tau_K_AoAk_tau(k,1) = rd_set_coeff'*kerAoAk.eval(Ak_all(:,:,k), Ao_rd_set)*rd_set_coeff;
end
% mmd_C_l = g'*tau_K_AoAk_tau;

% Constant term
kerAoAo = myProcessOptions(op, 'mmd_kernel', KGaussian(meddistance(Ao_rd_set)^2)); % this is important to calculate correctly
tau_K_AoAo_tau = rd_set_coeff'*kerAoAo.eval(Ao_rd_set, Ao_rd_set)*rd_set_coeff;
mmd_C_c = tau_K_AoAo_tau;


%----------------------
g = sdpvar(const.ngen*2,1,'full'); % generation and par ipation factor [pg;eta]

J_g  = g(1:const.ngen)'*const.Q*g(1:const.ngen) + const.c_g'*g(1:const.ngen)+sum(const.c_o)... 
       + g(const.ngen+1:end)'*(const.Q*Des_DCOPF.total_var*100)*g(const.ngen+1:end); % Cost function 

% MMD Term Related to Pg and constant
mmd_pg = g(1:const.ngen)'*(tau_K_AkAl_tau(1:const.ngen,1:const.ngen))*g(1:const.ngen) + ...
         g(1:const.ngen)'*(tau_K_AoAk_tau(1:const.ngen) - 2*mmd_B_1(1:const.ngen))+ mmd_C_c - 2*(mmd_B_c) + mmd_A;
 
% MMD Term Related to alpha       
mmd_alpha = g(const.ngen+1:end)'*(tau_K_AkAl_tau(const.ngen+1:end,const.ngen+1:end))*g(const.ngen+1:end) + ...
            g(const.ngen+1:end)'*(tau_K_AoAk_tau(const.ngen+1:end) - 2*mmd_B_1(const.ngen+1:end));

mmd = w_mmd_pg*mmd_pg + w_mmd_alpha*mmd_alpha;% only for reference  calculation;

objective =  w_j*J_g+w_mmd_pg*mmd_pg + w_mmd_alpha*mmd_alpha;


sync = find(const.g_u==0);
% The base point constraints will remain as it is. Basically, this means we want to safisfy these constraints at the mean prediction point 
constr_det = [sum(g(1:const.ngen)) == sum(const.d_hat);
    const.f_l <= const.Hg*(g(1:const.ngen)) - const.Hd*(const.d_hat) <= const.f_u;
    const.g_l <= g(1:const.ngen) <= const.g_u;
    g(const.ngen+sync) == 0;  % removing alphas for synchronous condensors 
%     mmd_alpha >=0;
    sum(g(const.ngen+1:end)) == 1; zeros(const.ngen,1) <= g(const.ngen+1:end) <= ones(const.ngen,1)];


ops = sdpsettings('solver','ipopt','verbose',0);
Res.sol = optimize(constr_det,objective,ops);

g_x = value(g);
Res.pg  = g_x(1:end/2);
Res.eta = g_x(const.ngen+1:end);
Res.gencost = value(J_g);
Res.time = toc+Des_DCOPF.sol.solvertime;
Flow = value(const.Hg*(g(1:const.ngen)) - const.Hd*(const.d_hat));
%% Sample based testing of Chance constraint DCOPF 
nt_r_x =10000;
x_test_for_compare_118_57;
 Des_DCOPF.xs_pen = x_test;
  g_xb = g_x;
  g_xb(1:const.ngen)=g_xb(1:const.ngen);
parfor j= 1:nt_r_x
[A1_test(:,:,j),Ao_test(:,j)] = dcopf_a1ao_dhatg(x_test(j,:)',data);
                   F_x(:,j)   =  A1_test(:,:,j)*g_xb + Ao_test(:,j);           
end
% Line_limit = repmat(const.f_u,[1,nt_r_x]);
% flow_vp = abs(Flow)-Line_limit;
F_x > 10^(-3);
Res.vio_individual = sum(ans,2);
Res.eps_individual = Res.vio_individual/nt_r_x;
[Res.max_eps_individual,~] = max(Res.eps_individual);

Res.vio_joint = sum(ans);
Res.vio_idx = Res.vio_joint >=1;
Res.eps_joint = (sum(Res.vio_idx)/nt_r_x);
Res.costs = [ value(w_mmd_pg*mmd_pg)  value(w_mmd_alpha*mmd_alpha) w_j*value(J_g)  ];%((value(J_g)-DC_det.f)/DC_det.f)*100


Res.w_j = w_j; 
Res.w_mmd_pg = w_mmd_pg;
Res.w_mmd_alpha = w_mmd_alpha;

rbus = const.loadbuses;
Res
 %
% Sample approximation code
% tic
% [Res.SA_DC_CCOPF] = sample_approximation(data,Des_DCOPF,solver_name,rbus);
%  Res.SA_DC_CCOPF.time = toc;
% % eps_obs_joint = SA_DC_CCOPF.eps_joint
% % SA_DC_CCOPF.obj
% 
% Res.SA_DC_CCOPF


function [A_1,A_o]= dcopf_a1ao_dhatg(xi,data)
const = ex_extract_ccDCOPF(data); %% Extract information from the MPC structure
I = eye(const.ngen,const.ngen);
Mg_up = [I sum(xi)*I];
A_1 = [Mg_up;-Mg_up;const.Hg*Mg_up;-const.Hg*Mg_up]; % Mg_down = -Mg_up;
A_o = [-const.g_u;...%-const.g_u-I*(const.d_hatg+xi(const.genbuses));
        const.g_l;...
       -const.Hd*(const.d_hat+xi)-const.f_u;...
        const.Hd*(const.d_hat+xi)+const.f_l];

end





function[Des_DCOPF] = desired_dist_DCOPF_local(data,target_size,nt_r,frac_pen,dist_type,solver_name)
% data = case30;
% base =100;
const = ex_extract_ccDCOPF(data); %% Extract information from the MPC structure
const.g_u = const.g_u*0.9;

% target_size = 100; % Number of samples used for SAA

rated_cap = sum(data.bus(const.loadbuses,3)*frac_pen)/length(const.loadbuses); % Total Rated capacity is fraction of total load at load buses
% dist_type = 'beta';
[~,xs_pen] = xs_generate_DG(data,const.loadbuses,nt_r,length(const.loadbuses),rated_cap,dist_type);

% tau = 0.01;
% % % % Constructing the two dataset given by the Ashley M Hou et.al 2020
% for ii = 1:nt_r
% A = randraw('uniform', [-1, 1], [length(const.loadbuses),length(const.loadbuses)]);
% Abar = A*A';
% for k = 1:length(const.loadbuses)
%     for j = 1:length(const.loadbuses)
%         sig(k,j,ii) = tau*(Abar(k,j)/sqrt(Abar(k,k)*Abar(j,j)))*sqrt(const.d_hat(k)*const.d_hat(j));
%     end
% end
% xi(:,ii) = mvnrnd(zeros(length(const.loadbuses),1),sig(:,:,ii))';
% end
% xs_pen = xi';

% xi = load('data_pglib_118_tau0a05.mat');
% xs_pen = xi';

[~,Des_DCOPF.desired_beta,Des_DCOPF.xi,Des_DCOPF.K,Des_DCOPF.sol_rd] = reduced_set_l1_penalization_local(solver_name,xs_pen,8,[]);


Des_DCOPF.xs_pen = xs_pen;
total_var = sum( var(xs_pen));
Des_DCOPF.total_var = total_var;

xi_all = Des_DCOPF.xi;
target_size = length(xi_all(:,1));
% iz = find(round(xs_pen(1,:),6) == 0);
% xs_nz = xs_pen;
% xs_nz(:,iz) = [];
% Des_DCOPF.xs_nonzero = xs_nz;
% Des_DCOPF.xs_pen = xs_pen;
% total_var = sum( var(xs_pen));
% Des_DCOPF.total_var = total_var;
% 
% [big_sample_size_row,~, ~]=size(Des_DCOPF.xs_pen);
% indices=randi(big_sample_size_row, 1 ,target_size);
% Des_DCOPF.xi_all=Des_DCOPF.xs_pen(indices,:,:);
% Des_DCOPF.xi = xs_pen(indices,:,:);
% op = struct();
% xi = Des_DCOPF.xi;
% xi_all = Des_DCOPF.xi_all;
% % for i=1:no_o
% %reduced_sample(:,:,x)=sample(indices(x,:),:,x);
% %reduced_sample1(:,:,x)=sample1(indices(x,:),:,x);
% ker= myProcessOptions(op, 'mmd_kernel', KGaussian(meddistance(xi)^2)); 
% kz=ker.eval(xi', xi');
% kzx=ker.eval(xi', xs_pen');
% alp=1*ones(big_sample_size_row,1);
% Des_DCOPF.desired_beta=(pinv(kz)*kzx*alp)./big_sample_size_row; % Weight of the reduced set 
% end


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
    const.f_l <= const.Hg*(g(1:const.ngen))- const.Hd*(const.d_hat) <= const.f_u;
    const.g_l <= g(1:const.ngen) <= const.g_u;
%      zeros(const.nload,1) <= s <= const.d_hat;
    sum(g(const.ngen+1:end)) == 1; 
    zeros(const.ngen,1) <= g(const.ngen+1:end) <= ones(const.ngen,1)];

% Sample based constraints for total N samples of xi

const_sample = [const.g_l <= g(1:const.ngen) + sum(xi_all(1,:))*g(const.ngen+1:end) <= const.g_u; % Generator Limits
               const.f_l <= const.Hg*(g(1:const.ngen)+...
               sum(xi_all(1,:))*g(const.ngen+1:end)) - const.Hd*((const.d_hat+xi_all(1,:)')) <= const.f_u]; 
for j = 2:target_size
 const_sample = [ const_sample; 
               const.g_l <= g(1:const.ngen)+ sum(xi_all(j,:))*g(const.ngen+1:end) <= const.g_u; % Generator Limits
               const.f_l <= const.Hg*(g(1:const.ngen)+...
               sum(xi_all(j,:))*g(const.ngen+1:end)) - const.Hd*((const.d_hat+xi_all(j,:)')) <= const.f_u]; % Power Flow limits
end

ops = sdpsettings('solver',solver_name,'verbose',0);

sol = optimize([constr_det;const_sample],objective,ops);

Des_DCOPF.sol = sol;
Des_DCOPF.obj=value(objective);
Des_DCOPF.pg = value(g(1:const.ngen));
Des_DCOPF.eta = value(g(const.ngen+1:end));
Des_DCOPF.g = value(g);

%


end


%% Large Sample Based DCOPF and Testing of Violation probability 

function [SA_DC_CCOPF] = sample_approximation(data,Des_DCOPF,solver_name,rbus)

% solver_name = 'mosek';


[SA_DC_CCOPF] = desired_dist_DCOPF_sample(data,Des_DCOPF,solver_name);



%% Testing via Samples 
nt_r_x = 1000;
x_test = Des_DCOPF.xs_pen; 
g_x = SA_DC_CCOPF.g;
parfor j= 1:nt_r_x
[A1_test(:,:,j),Ao_test(:,j)] = dcopf_a1ao(x_test(j,:)',data,rbus);
 F_x(:,j)   = A1_test(:,:,j)*g_x + Ao_test(:,j);
end
 
F_x > 10^(-3);
SA_DC_CCOPF.vio_individual = sum(ans,2);
SA_DC_CCOPF.eps_individual = SA_DC_CCOPF.vio_individual/nt_r_x;
[SA_DC_CCOPF.max_eps_individual,SA_DC_CCOPF.max_c] = max(SA_DC_CCOPF.eps_individual);


SA_DC_CCOPF.vio_joint = sum(ans);
SA_DC_CCOPF.vio_idx = SA_DC_CCOPF.vio_joint >=1;
SA_DC_CCOPF.eps_joint = (sum(SA_DC_CCOPF.vio_idx)/nt_r_x);




function[SA_DCOPF] = desired_dist_DCOPF_sample(data,Des_DCOPF,solver_name)
% data = case30;
% base =100;
const = ex_extract_ccDCOPF(data); %% Extract information from the MPC structure
SA_DCOPF.target_size = round((2/0.05)*(log(1/10^(-4))+2*const.ngen),0); % Number of samples used for SAA

[big_sample_size_row,~, ~]=size(Des_DCOPF.xs_pen);
indices=randi(big_sample_size_row, 1 ,SA_DCOPF.target_size);
xs_pen = Des_DCOPF.xs_pen(indices,:,:);

total_var = Des_DCOPF.total_var;
target_size = SA_DCOPF.target_size;
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

const_sample = [const.g_l <= g(1:const.ngen) + sum(xs_pen(1,:))*g(const.ngen+1:end) <= const.g_u; % Generator Limits
               const.f_l <= const.Hg*((g(1:const.ngen))+ sum(xs_pen(1,:))*g(const.ngen+1:end)) - const.Hd*((const.d_hat+xs_pen(1,:)')) <= const.f_u]; 
for jj = 2:target_size
 const_sample = [ const_sample; 
               const.g_l <= g(1:const.ngen)+ sum(xs_pen(jj,:))*g(const.ngen+1:end) <= const.g_u; % Generator Limits
               const.f_l <= const.Hg*((g(1:const.ngen))+ sum(xs_pen(jj,:))*g(const.ngen+1:end)) - const.Hd*((const.d_hat+xs_pen(jj,:)')) <= const.f_u]; % Power Flow limits
end

ops = sdpsettings('solver',solver_name,'verbose',0);

sol = optimize([constr_det;const_sample],objective,ops);

SA_DCOPF.sol = sol;
SA_DCOPF.obj=value(objective);
SA_DCOPF.pg = value(g(1:const.ngen));
SA_DCOPF.eta = value(g(const.ngen+1:end));
SA_DCOPF.g = value(g);




end


end


%% ------------------------------------------------------------------------------------------------------------
function [z,b_opt,rd_set,K,sol] = reduced_set_l1_penalization_local(solver_name,xs,lambda,K)

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


b = round(value(b_p-b_n));

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