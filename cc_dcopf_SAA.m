clear; close all; 
clc;

disp('begin to solve the problem');
% opt_yalmip = sdpsettings('usex0',1);

%% Overall settings
using_optimizer = 0; 
nMC = 1;
ops.beta = 10^(-3);
% distribution = 'gaussian';
distribution = 'beta';
casename = 'ex_case24_ieee_rts';



datapath = ['C:\Users\pare0001\OneDrive - Nanyang Technological University\Codes\PSCC_2022_CC-OPF_in_RKHS\ConvertChanceConstraint-ccc-master\extras\',casename,'\'];
resultpath = ['~/Documents/gdrive/Results-cc-DCOPF/results/',casename,'/',distribution,'/'];


yalmip('clear');




mpc = loadcase(casename);

%% Extract information from the MPC structure
const = ex_extract_ccDCOPF(mpc);

epsilons_all = 0.05;
for ieps = 1:length(epsilons_all)

% common settings
ops.epsilon = epsilons_all(ieps); 
ops.verbose = 1;


ops.method = 'sample average approximation';
% ops.method = 'robust counterpart';

switch ops.method

    case 'sample average approximation'  
        ops.type = 'sampling and discarding';
        PMAX = 9;
        ops.M = sum(mpc.gen(:,PMAX))*5;
    case 'robust counterpart'
        assert( using_optimizer == 0 ); % not compatiable with YALMIP optimizer 
        ops.type = 'ball';
        disp(ops.type);
    otherwise
        error('unknown method');
end
disp([ops.method,': ',ops.type]);


%% Formulate the cc-DCOPF Problem
% Detailed Formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% min_{g,s,eta} c(g) + c(s)     
% s.t. sum(g) = sum(d) - sum(s)
%      f_l <= Hg*g - Hd*(d_hat-s) <= f_u
%      g_l <= g <= g_u
%      sum(eta) = 1
%      Prob( f_l <= Hg*[g+sum(d_err)*eta] - Hd*(d_hat+d_err-s) <= f_u,
%           and g_l <= g+sum(d_err)*eta <= g_u) >= 1 - epsilon
%  Here s: Load shed variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Formulate the Problem Using YALMIP and CCC
g = sdpvar(const.ngen,1,'full'); % generation
s = sdpvar(const.nload,1,'full'); % load shed
eta = sdpvar(const.ngen,1,'full'); % affine corrective control policy
d_err = sdpvar(const.nload,1,'full'); % forecast errors of loads

objective = g'*const.Q*g + const.c_g'*g + const.c_s*sum(s) + const.c_g'*eta;

constr_det = [sum(g) == sum(const.d_hat) - sum(s);
    const.f_l <= const.Hg*g - const.Hd*(const.d_hat-s) <= const.f_u;
    const.g_l <= g <= const.g_u;
     zeros(const.nload,1) <= s <= const.d_hat;
%     s == 0;
    sum(eta) == 1;
    zeros(const.ngen,1) <= eta <= ones(const.ngen,1)];
constr_inner = [const.g_l <= g + sum(d_err)*eta <= const.g_u;
    const.f_l <= const.Hg*(g+ sum(d_err)*eta) - const.Hd*(const.d_hat+d_err-s) <= const.f_u];
% d_err_data = mvnrnd(zeros(nload,1), 0.05*diag(d_hat),N)';


N_trains = 50;

for iN = 1:length(N_trains)
N_train = N_trains(iN);
N_test = 10^4;
 
%% Solve the Problem
iMC = 1;
        
        traindata = load([casename,...
            '-traindata-N=',num2str(N_train),'-iMC=',num2str(iMC),'.mat'] );
        d_err_data = traindata.d_err;
        constr_prob = prob(constr_inner,d_err, ops.epsilon, d_err_data,ops);
        
        yalmip_ops = sdpsettings('solver','mosek','verbose',0); % sdpt3
        sol = optimize([constr_det; constr_prob],objective,yalmip_ops);    
        disp(sol.info);
       
        if sol.problem ~= 0
            disp('Unsolved Problem');
            continue;
        else
            g = value(g)
            eta = value(eta)
        end
        
        testdata = load([casename,...
            '-testdata-N=',num2str(N_test),'.mat'] );
        
        eps_ofs = estimate_violation_probability(constr_inner, d_err, testdata.d_err, ops);
        disp(['out of sample violation probability is: ', num2str(eps_ofs)]);
        disp(ops.epsilon);
        if strcmp(ops.method, 'sample average approximation')
            [~,violated] = estimate_violation_probability(constr_inner, d_err, traindata.d_err, ops);
            results(iMC).violated = violated;
            disp([num2str(length(violated.indices)),' scenarios being violated']);
        end
        % save results
%         results(iMC).epsilon = ops.epsilon;
%         results(iMC).N_train = N_train;
%         results(iMC).N_test = N_test;
%         results(iMC).eps_ofs = eps_ofs;
%         results(iMC).obj = value(objective);
%         results(iMC).g = value(g);
%         results(iMC).s = value(s);
%         results(iMC).eta = value(eta);
%         results(iMC).solvetime = sol.solvertime;
%         results(iMC).ops = ops;
%       eps_test(iMC,iN) = eps_ofs;
    end
    if strcmp(ops.method, 'scenario approach')
        save([resultpath,casename,'-',ops.method,'-results-N=',num2str(N_train),'.mat'],'results');
    else
        save([resultpath,casename,'-',ops.method,'-',ops.type,...
            '-results-N=',num2str(N_train),'-epsilon=',num2str(ops.epsilon),'.mat'],'results');
    end    
end





