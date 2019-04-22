
% This file was used to compare results of Table 6.7.4 In Haftka, Case A.

clear all;close all;clc
format compact
syms A1 A2 A3 A4 A5 A6 A7 A8 A9 A10

tic
[K,E,M,rhow,L,c,s] = TenBarTruss();
x0 = 20*ones(1,10);
epsilon = 1e-3;
for i =1:30
%%% Invert K Matrix @ Appropriate Design Point
if i == 1
    K_reduced = subs(K,{A1,A2,A3,A4,A5,A6,A7,A8,A9,A10},x0);
    Kinv = K_reduced\eye(length(K_reduced));
    [stress,dstdA,D,dQdA] = Sensitivity(K,Kinv,E,L,c,s);
else
    K_reduced = subs(K,{A1,A2,A3,A4,A5,A6,A7,A8,A9,A10},x);
    Kinv = K_reduced\eye(length(K_reduced));
    [stress,dstdA,D,dQdA] = Sensitivity(K,Kinv,E,L,c,s);
end

%%% Move Limits
if i <= 8
    MLl = 0.9 - 0.1*i;
    MLu = 0.9 - 0.1*i;
elseif i>8 && i<=12
    MLl = 0.1; MLu = 0.1;
elseif i>12 && i<=20
    MLl = 0.05; MLu = 0.05;
else
    MLl = 0.01; MLu = 0.01;
end
if i == 1;
    LB = x0 - MLl*x0; UB = x0 + MLu*x0;
else
    LB = x - MLl*x; UB = x + MLu*x;
    x0 = x;
end

%%% Enforce Side Bounds Based on Move Limits
if min(LB) <= 0.1
   for j = 1:length(LB)
       if LB(j) <= 0.1
            LB(j) = 0.1;
       end
   end
end
if max(UB) >= 40
   for j = 1:length(UB)
       if UB(j) >= 40
           UB(j) = 40;
       end
   end
end

%%% Min Weight s.t. Stress Constraints
if i == 1
    [coeff,dcon] = DispApprox(x0,D,dQdA);
else
    [coeff,dcon] = DispApprox(x,D,dQdA);
    clear x;
end
A = []; B = [];
Aeq = []; Beq = [];
options = optimset('Display','off','TolFun',1e-6,'Algorithm','SQP','PlotFcns',@optimplotfval);
[x,fval,exitflag,output]=fmincon(@(x) weight(x,L,rhow),x0,A,B,Aeq,Beq,LB,UB,@ (x) StressApprox(x,x0,stress,dstdA),options);

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