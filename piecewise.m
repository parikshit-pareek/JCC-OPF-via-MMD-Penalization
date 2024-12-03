%  load('IEEE_30_Bus_data.mat')
% data=IEEE_30_Bus_data;
% ng=size(data.gen(:,1));
% ng=ng(1);
%% Piecewise linearization of the cost function
% CfS:= the cost multipliers for piecewise linearize DCOPF model 
% S := slope of lines
% FPg:= Cost at different values of Pg's 
% FPg:= [F(Pg_min) F(Pg1) F(Pg2) F(Pg_max)]
% First column of FPg represents the F(Pg_min)
function [Pgf,FPg,CfS]=piecewise(data,ng,nbus,Pg,interval,h)
Cf1=data.gencost(:,6);
Cf2=data.gencost(:,5);
Cf0=data.gencost(:,7);
ramp_up=data.gen(:,11);
ramp_down=data.gen(:,12);
if nargin <=3
    Pmin=data.gen(:,10);
    Pmax=data.gen(:,9);
elseif h==1
    Pmin=data.gen(:,10);
    Pmax=data.gen(:,9);
else 
    Pmin=max(data.gen(:,10),Pg(:,h-1)-ramp_down*interval);
    Pmax=min(data.gen(:,9),Pg(:,h-1)+ramp_up*interval);
end
% Cost at any Pgen

Pgf=zeros(ng,4);
% First column of FPg represents the F(Pm_min)
FPg=zeros(ng,4);
S=zeros(ng,3);
for i=1:ng
Pgf(i,:)=linspace(Pmin(i),Pmax(i),4);
    for j=1:4
        FPg(i,j)=Cf0(i)+Cf1(i).*Pgf(i,j)+Cf2(i).*(Pgf(i,j).^2);
    end
end
%  Calculation of slope of lines
for i=1:ng
    for j=1:3
        if (Pgf(i,j+1)-Pgf(i,j))==0
            S(i,j)=0;
        else
        S(i,j)=(FPg(i,j+1)-FPg(i,j))/(Pgf(i,j+1)-Pgf(i,j));
        end
    end
end
% CfS: the cost multipliers for piecewise linearize DCOPF model 
k=1;
CfS=zeros(1,3*ng+nbus);
for i=1:ng
     CfS(1,k:k+2)=S(i,:);
 k=k+3;
end
end

