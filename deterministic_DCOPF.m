function[Det_DCOPF] = deterministic_DCOPF(data,solver_name)
% data = case30;
% base =100;
const = ex_extract_ccDCOPF(data); %% Extract information from the MPC structure

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
g = sdpvar(const.ngen,1,'full'); % generation and participation factor [pg;eta]
% eta = sdpvar(const.ngen,1,'full'); % affine corrective control policy


objective = (g(1:const.ngen))'*const.Q*(g(1:const.ngen)) + const.c_g'*(g(1:const.ngen))+sum(const.c_o);

constr_det = [sum(g(1:const.ngen)) == sum(const.d_hat);
    const.f_l <= const.Hg*(g(1:const.ngen))- const.Hd*(const.d_hat) <= const.f_u;
    const.g_l <= g(1:const.ngen) <= const.g_u];

% Sample based constraints for total N samples of xi

ops = sdpsettings('solver',solver_name,'verbose',0);

sol = optimize([constr_det],objective,ops)

Det_DCOPF.sol = sol;
Det_DCOPF.obj=value(objective);
Det_DCOPF.pg = value(g(1:const.ngen));
Det_DCOPF.eta = value(g(const.ngen+1:end));
Det_DCOPF.g = value(g);


end
