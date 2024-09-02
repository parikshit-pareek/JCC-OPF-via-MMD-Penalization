function [A_1,A_o]= dcopf_a1ao(xi,data,rbus)
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