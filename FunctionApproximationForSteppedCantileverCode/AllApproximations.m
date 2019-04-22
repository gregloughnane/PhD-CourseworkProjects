% Greg Loughnane

%% Initialize
clear all; close all; clc ; format compact
syms x1 x2 x3

%% Tip-Displacement Constraint, Initial Guess, and Corresponding Function Value 
g = -(61/(x1.^3) + 37/(x2.^3) + 19/(x3.^3) - 1);
x1p = [6.5 6.5 6.5];                  % Previous Point
x2p = [9 9 9];                  % Current Point
gx1 = subs(g,{x1,x2,x3},x1p);   % Function Value @ Previous Point
gx2 = subs(g,{x1,x2,x3},x2p);   % Function Value @ Previous Point
xi = [x1 x2 x3];
%% Gradient Information
gradgxi = [diff(g,x1)
           diff(g,x2)
           diff(g,x3)];
dgdxi_x1 = subs(gradgxi,{x1,x2,x3},x1p);
dgdxi_x2 = subs(gradgxi,{x1,x2,x3},x2p);

grad2gxi = [diff(gradgxi(1),x1)   0 0;
            0  diff(gradgxi(2),x2)  0;
            0 0   diff(gradgxi(3),x3)];
d2gdxi2_x1 = subs(diag(grad2gxi),{x1,x2,x3},x1p);

%% Linear, Reciprocal, Conservative, Quadratic
gL = expand(gx1 + sum((xi - x1p).*dgdxi_x1'));
gR = expand(gx1 + sum((xi - x1p).*(x1p./xi).*dgdxi_x1'));
for i = 1:length(x1p)
    if x1p(i)*dgdxi_x1(i) <= 0
        Gi = 1;
    else
        Gi = x1p(i)/x1;
    end
    G(:,i) = Gi;
end
gC = gx1 + sum((xi - x1p).*G.*dgdxi_x1');
gQ = gx1 + sum((xi - x1p).*dgdxi_x1') + 1/2*sum((xi - x1p).^2.*d2gdxi2_x1');

%% TANA
xi1 = x1p; xi2 = x2p;
% Solve Nonlinearity Index r
options = optimset('Display','off','TolFun',1e-10); x0 = 1; % options = optimset('Display','iter');
f = @(r)(gx1 - (gx2 + (1/r)*sum(xi2.*(1-r).*(xi1.^r - xi2.^r).*dgdxi_x2')));
[r,fval] = fsolve(f,x0,options);
% Plug r into TANA expansion equation
xi = [x1 x2 x3];
gTANA = (gx2 + (1/r)*sum(xi2.*(1-r).*(xi.^r - xi2.^r).*dgdxi_x2'));

%% TANA1
% Solve Nonlinearity Indices pi
p0 = [-1;-1;-1]; % Initial Guess
options = optimset('Display','off','TolFun',1e-10);  % options = optimset('Display','iter');
h = @(p)[dgdxi_x2(1) - (dgdxi_x1(1)*(xi2(1)/xi1(1)))^(p(1)-1);...
         dgdxi_x2(2) - (dgdxi_x1(2)*(xi2(2)/xi1(2)))^(p(2)-1);...
         dgdxi_x2(3) - (dgdxi_x1(3)*(xi2(3)/xi1(3)))^(p(3)-1)];
[p,hval] = fsolve(h,p0,options);
% Compute Epsilon Residual
epsilon = gx2 - (gx1 +(xi1(1)^(1-p(1))/p(1)*(xi2(1).^p(1) - xi1(1)^p(1))*dgdxi_x1(1) +...
                       xi1(2)^(1-p(2))/p(2)*(xi2(2).^p(2) - xi1(2)^p(2))*dgdxi_x1(2) +...
                       xi1(3)^(1-p(3))/p(3)*(xi2(3).^p(3) - xi1(3)^p(3))*dgdxi_x1(3)));
% Plug p values into TANA1 expansion equation and correct with epsilon
xi = [x1 x2 x3];
gTANA1 = (gx1 +(xi1(1)^(1-p(1))/p(1)*(xi(1).^p(1) - xi1(1)^p(1))*dgdxi_x1(1) +...
                xi1(2)^(1-p(2))/p(2)*(xi(2).^p(2) - xi1(2)^p(2))*dgdxi_x1(2) +...
                xi1(3)^(1-p(3))/p(3)*(xi(3).^p(3) - xi1(3)^p(3))*dgdxi_x1(3))) +...
                epsilon;

%% TANA2
% Solve Nonlinearity Indices pi
z0 = [1;1;1;0]; % Initial Guess ~ Note that Epsilon is the 4th term
options = optimset('Display','off','TolFun',1e-10);  % options = optimset('Display','iter');
f = @(z)[gx1 - (gx2 +  (xi2(1)^(1-z(1))/z(1)*(xi1(1).^z(1) - xi2(1)^z(1))*dgdxi_x2(1) +...
                        xi2(2)^(1-z(2))/z(2)*(xi1(2).^z(2) - xi2(2)^z(2))*dgdxi_x2(2) +...
                        xi2(3)^(1-z(3))/z(3)*(xi1(3).^z(3) - xi2(3)^z(3))*dgdxi_x2(3) +...
                1/2*z(4)*((xi1(1)-xi2(1))^2 + (xi1(2)-xi2(2))^2)+ (xi1(3)-xi2(3))^2));
                dgdxi_x1(1) - (dgdxi_x2(1)*(xi1(1)/xi2(1)))^(z(1)-1);...
                dgdxi_x1(2) - (dgdxi_x2(2)*(xi1(2)/xi2(2)))^(z(2)-1);...
                dgdxi_x1(3) - (dgdxi_x2(3)*(xi1(3)/xi2(3)))^(z(3)-1)]; 
[z,fval] = fsolve(f,z0,options);
% Plug p and epsilon values into TANA2 expansion equation
xi = [x1 x2 x3];
gTANA2 = (gx2 +  (xi2(1)^(1-z(1))/z(1)*(xi(1).^z(1) - xi2(1)^z(1))*dgdxi_x2(1) +...
                  xi2(2)^(1-z(2))/z(2)*(xi(2).^z(2) - xi2(2)^z(2))*dgdxi_x2(2) +...
                  xi2(3)^(1-z(3))/z(3)*(xi(3).^z(3) - xi2(3)^z(3))*dgdxi_x2(3) +...
         1/2*z(4)*((xi(1)-xi2(1))^2 + (xi(2)-xi2(2))^2)+ (xi(3)-xi2(3))^2)); 

%% Least Squares Regression
% Define factor levels for comparison and code variables (FOR REGRESSION)
low = x1p; high = x2p;
k = length(xi); % # D.V.s

X = csvread('2_Level.csv');
X(:,2:4) = X; X(:,1) = ones(length(X),1);
for i = 1:length(X)
    for j = 1:k
    if X(i,j+1) == -1
        X(i,j+1) = low(1);
    elseif X(i,j+1) == 1
        X(i,j+1) = high(1);
    elseif X(i,j+1) == -0.33
        X(i,j+1) = (low(1) + 2/3);
    elseif X(i,j+1) == 0.33
        X(i,j+1) = (high(1) - 2/3);
    end
    Y(i,:) = subs(g,{x1,x2,x3},X(i,2:4));
    end
end
B = (real((X'*X)\eye(4)))*(X'*Y); % Solve for beta coefficients
B0 = B(1); Bi = B(2:4); xji = [x1 x2 x3]; 
Yj = B0 + sum(Bi.*xji'); % Fit Linear Model
gRSM = vpa(real(Yj),4);

X4 = csvread('4_Level.csv');
X4(:,2:4) = X4; X4(:,1) = ones(length(X4),1);
for i = 1:length(X4)
    for j = 1:k
    if X4(i,j+1) == -1
        X4(i,j+1) = low(1);
    elseif X4(i,j+1) == 1
        X4(i,j+1) = high(1);
    elseif X4(i,j+1) == -0.33
        X4(i,j+1) = (low(1) + 2/3);
    elseif X4(i,j+1) == 0.33
        X4(i,j+1) = (high(1) - 2/3);
    end
    Y4(i,:) = subs(g,{x1,x2,x3},X4(i,2:4));
    end
end
B4 = (real((X4'*X4)\eye(4)))*(X4'*Y4); % Solve for beta coefficients
B04 = B4(1); Bi4 = B4(2:4); xji = [x1 x2 x3]; 
Yj4 = B04 + sum(Bi4.*xji'); % Fit Linear Model
gRSM4 = vpa(real(Yj4),4);

% Run Monte Carlo Simulation and Compute r^2 value
% points = input('How many points? : ')
points = 20;
point = x1p; 
for i = 1:points;
    add_vector(i,:) = 1/points*(x2p-x1p).*rand(1,3);
    point(i+1,:) = x1p + (i)*add_vector(i,:);
    Yi(i,:) = double(subs(g,{x1,x2,x3},point(i,:)));
    Yihat(i,:) = subs(gRSM,{x1,x2,x3},point(i,:));
    Yihat4(i,:) = subs(gRSM4,{x1,x2,x3},point(i,:));
end
Ybar = mean(Yihat);Ybar4 = mean(Yihat4);
rsquared = 1 - sum((Yi - Yihat).^2)./sum((Yi - ones(size(Yi),1)*Ybar).^2);
rsquared4 = 1 - sum((Yi - Yihat4).^2)./sum((Yi - ones(size(Yi),1)*Ybar4).^2);

%% Percent Error For Each Approximation
n = 20; % Number of Xi points to compute function approximations at
t = -1; % Start t at -1
for i=1:n+1
    xnew(i,:) = (t + 1/2).*x2p + (1/2 - t).*x1p;
    gx(i,:) = double(subs(g,{x1,x2,x3},xnew(i,:)));     dg(i,:) = (gx(i,:)-gx(i,:))/gx(i,:);
    gLx(i,:) = double(subs(gL,{x1,x2,x3},xnew(i,:)));   dgL(i,:) = (gLx(i,:)-gx(i,:))/gx(i,:)/100;
    gRx(i,:) = double(subs(gR,{x1,x2,x3},xnew(i,:)));   dgR(i,:) = (gRx(i,:)-gx(i,:))/gx(i,:)/100;
    gCx(i,:) = double(subs(gC,{x1,x2,x3},xnew(i,:)));   dgC(i,:) = (gCx(i,:)-gx(i,:))/gx(i,:)/100;
    gQx(i,:) = double(subs(gQ,{x1,x2,x3},xnew(i,:)));   dgQ(i,:) = (gQx(i,:)-gx(i,:))/gx(i,:)/100;
    gTx(i,:) = subs(gTANA,{x1,x2,x3},xnew(i,:));        dgT(i,:) = (gTx(i,:)-gx(i,:))/gx(i,:)/100;
    gT1x(i,:) = subs(gTANA1,{x1,x2,x3},xnew(i,:));      dgT1(i,:) = (gT1x(i,:)-gx(i,:))/gx(i,:)/100;
    gT2x(i,:) = subs(gTANA2,{x1,x2,x3},xnew(i,:));      dgT2(i,:) = (gT2x(i,:)-gx(i,:))/gx(i,:)/100;
    gRSMx(i,:) = subs(gRSM,{x1,x2,x3},xnew(i,:));       dgRSM(i,:) = (gRSMx(i,:)-gx(i,:))/gx(i,:)/100;
    gRSM4x(i,:) = subs(gRSM,{x1,x2,x3},xnew(i,:));      dgRSM4(i,:) = (gRSM4x(i,:)-gx(i,:))/gx(i,:)/100;
    step(i,:) = t; t = t + 2/n;
end

%% Plot Error For All Methods Against One Another
disp('--------------------------------------------------------------------------------------------')
fprintf('   x1   x2    x3     g(x)   gL(x)  gR(x)  gC(x)  gQ(x)  gT(x) gT1(x) gT2(x) gRSM2(x) gRSM4(x)\n');
disp('--------------------------------------------------------------------------------------------')
for i = 1:length(xnew)
    fprintf('%6.3f %6.3f %6.3f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f \n', ...
  xnew(i,1),xnew(i,2),xnew(i,3),gx(i,1),gLx(i,1),gRx(i,1),gCx(i,1),gQx(i,1),gTx(i,1),gT1x(i,1),gT2x(i,1),gRSMx(i,1),gRSM4x(i,1));
end
figure;plot(step,dg,'o-',step,dgL,step,dgR,step,dgC,step,dgQ,step,dgT,step,dgT1,step,dgT2,'d-',step,dgRSM,'+-',step,dgRSM4,'x-')
legend('g','gL','gR','gC','gQ','gT','gT1','gT2','gRSM2','gRSM4');xlabel('t');ylabel('Percent Error')
title('Percent Error About Xi = (t + 1/2)X2 + (1/2 - t)X1 Points');

%% Plot Function Value vs. Function Approximations for x1=x2=x3
figure;plot(xnew(:,1),gx,'o-',xnew(:,1),gLx,xnew(:,1),gRx,xnew(:,1),gCx,xnew(:,1),gQx,xnew(:,1), ...
    gTx,xnew(:,1),gT1x,xnew(:,1),gT2x,'d-',xnew(:,1),gRSMx,'+-',xnew(:,1),gRSM4x,'x-')
legend('g','gL','gR','gC','gQ','gT','gT1','gT2','gRSM2','gRSM4');
xlabel('Design Point x=x1=x2=x3');ylabel('Function/Approximation Value');title('Comparison of Function Approximations');
