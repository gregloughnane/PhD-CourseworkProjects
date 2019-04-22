
% This file was used to compare results of Table 6.4.2 In Haftka.

clear all;close all;clc
format compact
syms A1 A2 A3 A4 A5 A6 A7 A8 A9 A10

tic
[K,E,M,rho,rhow,L,c,s] = TenBarTruss();
x0 = 5*ones(1,10);
epsilon = 1e-6;
for i =1:100
    
%%% Invert K Matrix @ Appropriate Design Point
if i == 1
    K_reduced = subs(K,{A1,A2,A3,A4,A5,A6,A7,A8,A9,A10},x0);
    Kinv = K_reduced\eye(length(K_reduced));
    [stress,dstdA,w,dwdA] = Sensitivity(K,Kinv,M,E,L,c,s,rho,rhow,x0);
    M_reduced = subs(M,{A1,A2,A3,A4,A5,A6,A7,A8,A9,A10},x0);
else
    K_reduced = subs(K,{A1,A2,A3,A4,A5,A6,A7,A8,A9,A10},x);
    Kinv = K_reduced\eye(length(K_reduced));
    [stress,dstdA,w,dwdA] = Sensitivity(K,Kinv,M,E,L,c,s,rho,rhow,x);
    M_reduced = subs(M,{A1,A2,A3,A4,A5,A6,A7,A8,A9,A10},x);
end

%%% Move Limits
if i <= 11
    MLl = 0.75 - 0.05*i;
    MLu = 0.75 - 0.05*i;
elseif i>12 && i<=15
    MLl = 0.1; MLu = 0.1;
elseif i>15 && i<=20
    MLl = 0.05; MLu = 0.05;
else
    MLl = 0.02; MLu = 0.02;
end
if i == 1;
    LB = x0 - MLl*x0; UB = x0 + MLu*x0;
else
    LB = x - MLl*x; UB = x + MLu*x;
    x0 = x; clear x;
end

%%% Enforce Side Bounds Based on Move Limits
if min(LB) <= 0.1
   for j = 1:length(LB)
       if LB(j) <= 0.1
            LB(j) = 0.1;
       end
   end
end
if max(UB) >= 10
   for j = 1:length(UB)
       if UB(j) >= 10
           UB(j) = 10;
       end
   end
end

%%% Min Weight s.t. Stress Constraints
A = []; B = [];Aeq=[]; Beq=[];
options = optimset('Display','off','TolFun',1e-6,'Algorithm','SQP','PlotFcns',@optimplotfval);
[x,fval,exitflag,output]=fmincon(@(x) weightapprox(x,x0,w,dwdA),x0,A,B,Aeq,Beq,LB,UB,@ (x) StressApprox(x,x0,stress,dstdA),options);
f(i,1) = fval;
%%% Convergence Check
checkweight(i,1) = fval;
    if i <= 2
        pch(1:2,1) = 1;
    else
        pch(i,1) = abs((checkweight(i,1) - checkweight(i-1,1))/checkweight(i-1,1)*100);
    end
    if i >=3 
        if pch(i,1) < epsilon
            if pch(i-1,1) < epsilon
                if pch(i-2,1) < epsilon
                    break
                end
            end
        end
    end
end
toc
fprintf('\n');fprintf('Optimum Areas');
fprintf('\n');fprintf('%3.2f ',x);fprintf('\n');
fprintf('\n');fprintf('Optimum Weight');
fprintf('\n');fprintf('%3.2f ',fval);fprintf('\n');
f