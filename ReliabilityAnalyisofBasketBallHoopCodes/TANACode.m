% TANA
% Solve Nonlinearity Index r
options = optimset('Display','iter','TolFun',1e-10); 
x0 = 1;
f = @(r)(LSO - (LS(i) + (1/r) *((xt0^(1-r)*dLSdE0*(mut^r - xt0^r)) + ...
                              (xF0^(1-r)*dLSdF0*(muF^r - xF0^r)) + ...
                              (xEb0^(1-r)*dLSdEb0*(muEb^r - xEb0^r)) + ...
                              (xvb0^(1-r)*dLSdvb0*(muvb^r - xvb0^r)) + ...
                              (xEs0^(1-r)*dLSdEs0*(muEs^r - xEs0^r)) + ...
                              (xvs0^(1-r)*dLSdvs0*(muvs^r - xvs0^r)) ) ) )
[r,fval] = fsolve(f,x0,options)
% Plug r into TANA expansion equation to determine g(x)
syms x1 x2 x3 x4 x5 x6 gTANA
limitstate = sym(LSO,'f')
xT = sym(xT0,'f')
xf = sym(xF0,'f')
xeb = sym(xEb0,'f')
xvb = sym(xvb0,'f')
xes = sym(xEs0,'f')
xvs = sym(xvs0,'f')
dgde = sym(dDdE0,'f')
dgdf = sym(dDdF0,'f')
dgda = sym(dDdA0,'f')
dgdl = sym(dDdL0,'f')
dgdi = sym(dDdI0,'f')
gTANA = symfun(limitstate + (1/r)*((xe^(1-r)*dgde*(x1^r - xe^r)) + ...
                                   (xf^(1-r)*dgdf*(x2^r - xf^r)) + ...
                                   (xa^(1-r)*dgda*(x3^r - xa^r)) + ...
                                   (xl^(1-r)*dgdl*(x4^r - xl^r)) + ...
                                   (xi^(1-r)*dgdi*(x5^r - xi^r)) ),...
                                   [x1 x2 x3 x4 x5])
% Convert g(x) to g(u)
syms u1 u2 u3 u4 u5 gTANAu
x1n = muE + u1 * sigmaE
x2n = muF + u2 * sigmaF
x3n = muA + u3 * sigmaA
x4n = muL + u4 * sigmaL
x5n = muI + u5 * sigmaI
gTANAu = vpa(subs(gTANA,{x1,x2,x3,x4,x5},[x1n,x2n,x3n,x4n,x5n]),4)
fhandle = @gTANAu;
% Solve beta minimization on g(u)
x0 = [0 0 0 0 0]; LB = 2*[-1;-1;-1;-1;-1]; UB = 2*[1 ;1 ;1 ;1 ;1];
A = [ ]; B = [ ]; 
Aeq = []; Beq = []; options = optimset('Display','iter')
[x,b1,exitflag,output] = fmincon(@myfun,x0,A,B,Aeq,Beq,LB,UB,fhandle,options)                                
e1 = abs(( b1 - b0 ) / b0)
xE1 = subs((muE + u1 * sigmaE),{u1,u2,u3,u4,u5},x)
xF1 = subs((muF + u2 * sigmaF),{u1,u2,u3,u4,u5},x)
xA1 = subs((muA + u3 * sigmaA),{u1,u2,u3,u4,u5},x)
xL1 = subs((muL + u4 * sigmaL),{u1,u2,u3,u4,u5},x)
xI1 = subs((muI + u5 * sigmaI),{u1,u2,u3,u4,u5},x)