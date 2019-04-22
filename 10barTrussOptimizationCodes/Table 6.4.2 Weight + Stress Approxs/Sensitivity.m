function [stress,dstdA,w,dwdA] = Sensitivity(K,Kinv,M,E,L,c,s,rho,rhow,x)
%% Sensitivity Analysis ~ Adjoint Method
% Element Stresses (psi)
syms A1 A2 A3 A4 A5 A6 A7 A8 A9 A10
%% Displacement @ Joints (in)
p = 100000;
P = [0 0 0 -p 0 0 0 -p 0 0 0 0]'; P_reduced = P(1:8,:);
D = Kinv*P_reduced; 
Q = D; D(9:12) = 0; % Note B.C.s
%% Total Weight
w = rhow*x*L';
%% Element Stresses
u1 = D(1); v1 = D(2);   u2 = D(3); v2 = D(4);   u3 = D(5); v3 = D(6);
u4 = D(7); v4 = D(8);   u5 = D(9); v5 = D(10);  u6 = D(11); v6 = D(12);
st1 = E/L(1)*(c(1)*(u3-u5) + s(1)*(v3-v5)); st2 = E/L(2)*(c(2)*(u1-u3) + s(2)*(v1-v3));
st3 = E/L(3)*(c(3)*(u4-u6) + s(3)*(v4-v6)); st4 = E/L(4)*(c(4)*(u2-u4) + s(4)*(v2-v4));
st5 = E/L(5)*(c(5)*(u3-u4) + s(5)*(v3-v4)); st6 = E/L(6)*(c(6)*(u1-u2) + s(6)*(v1-v2));
st7 = E/L(7)*(c(7)*(u4-u5) + s(7)*(v4-v5)); st8 = E/L(8)*(c(8)*(u3-u6) + s(8)*(v3-v6));
st9 = E/L(9)*(c(9)*(u2-u3) + s(9)*(v2-v3)); st10 = E/L(10)*(c(10)*(u1-u4) + s(10)*(v1-v4));
stress = [st1 st2 st3 st4 st5 st6 st7 st8 st9 st10]';
clear u1 u2 u3 u4 u5 u6 v1 v2 v3 v4 v5 v6
syms u1 u2 u3 u4 u5 u6 v1 v2 v3 v4 v5 v6
st1 = E/L(1)*(c(1)*(u3-u5) + s(1)*(v3-v5)); st2 = E/L(2)*(c(2)*(u1-u3) + s(2)*(v1-v3));
st3 = E/L(3)*(c(3)*(u4-u6) + s(3)*(v4-v6)); st4 = E/L(4)*(c(4)*(u2-u4) + s(4)*(v2-v4));
st5 = E/L(5)*(c(5)*(u3-u4) + s(5)*(v3-v4)); st6 = E/L(6)*(c(6)*(u1-u2) + s(6)*(v1-v2));
st7 = E/L(7)*(c(7)*(u4-u5) + s(7)*(v4-v5)); st8 = E/L(8)*(c(8)*(u3-u6) + s(8)*(v3-v6));
st9 = E/L(9)*(c(9)*(u2-u3) + s(9)*(v2-v3)); st10 = E/L(10)*(c(10)*(u1-u4) + s(10)*(v1-v4));
% Derivatives of Stiffness w.r.t. Areas
dkdA1 = diff(K,A1);dkdA2 = diff(K,A2);dkdA3 = diff(K,A3);dkdA4 = diff(K,A4);
dkdA5 = diff(K,A5);dkdA6 = diff(K,A6);dkdA7 = diff(K,A7);dkdA8 = diff(K,A8);
dkdA9 = diff(K,A9);dkdA10 = diff(K,A10);
% Derivatives of Force Vector w.r.t. Areas
dFdA1(8,1)=0;dFdA2(8,1)=0;dFdA3(8,1)=0;dFdA4(8,1)=0;dFdA5(8,1)=0;
dFdA6(8,1)=0;dFdA7(8,1)=0;dFdA8(8,1)=0;dFdA9(8,1)=0;dFdA10(8,1)=0;
% Partial Derivatives of Displacements w.r.t. Areas
dQdA1 = Kinv*(dFdA1 - dkdA1*Q);dQdA2 = Kinv*(dFdA2 - dkdA2*Q);
dQdA3 = Kinv*(dFdA3 - dkdA3*Q);dQdA4 = Kinv*(dFdA4 - dkdA4*Q);
dQdA5 = Kinv*(dFdA5 - dkdA5*Q);dQdA6 = Kinv*(dFdA6 - dkdA6*Q);
dQdA7 = Kinv*(dFdA7 - dkdA7*Q);dQdA8 = Kinv*(dFdA8 - dkdA8*Q);
dQdA9 = Kinv*(dFdA9 - dkdA9*Q);dQdA10 = Kinv*(dFdA10 - dkdA10*Q);
dQdA = [dQdA1 dQdA2 dQdA3 dQdA4 dQdA5 dQdA6 dQdA7 dQdA8 dQdA9 dQdA1];
displacement_sensitivities = vpa(dQdA,4);
% Partial Derivatives of Element Stresses w.r.t. Displacements
z1 = [diff(st1,u1); diff(st1,v1);diff(st1,u2);diff(st1,v2);diff(st1,u3);diff(st1,v3);diff(st1,u4);diff(st1,v4)];
z2 = [diff(st2,u1); diff(st2,v1);diff(st2,u2);diff(st2,v2);diff(st2,u3);diff(st2,v3);diff(st2,u4);diff(st2,v4)];
z3 = [diff(st3,u1); diff(st3,v1);diff(st3,u2);diff(st3,v2);diff(st3,u3);diff(st3,v3);diff(st3,u4);diff(st3,v4)];
z4 = [diff(st4,u1); diff(st4,v1);diff(st4,u2);diff(st4,v2);diff(st4,u3);diff(st4,v3);diff(st4,u4);diff(st4,v4)];
z5 = [diff(st5,u1); diff(st5,v1);diff(st5,u2);diff(st5,v2);diff(st5,u3);diff(st5,v3);diff(st5,u4);diff(st5,v4)];
z6 = [diff(st6,u1); diff(st6,v1);diff(st6,u2);diff(st6,v2);diff(st6,u3);diff(st6,v3);diff(st6,u4);diff(st6,v4)];
z7 = [diff(st7,u1); diff(st7,v1);diff(st7,u2);diff(st7,v2);diff(st7,u3);diff(st7,v3);diff(st7,u4);diff(st7,v4)];
z8 = [diff(st8,u1); diff(st8,v1);diff(st8,u2);diff(st8,v2);diff(st8,u3);diff(st8,v3);diff(st8,u4);diff(st8,v4)];
z9 = [diff(st9,u1); diff(st9,v1);diff(st9,u2);diff(st9,v2);diff(st9,u3);diff(st9,v3);diff(st9,u4);diff(st9,v4)];
z10 = [diff(st10,u1); diff(st10,v1);diff(st10,u2);diff(st10,v2);diff(st10,u3);diff(st10,v3);diff(st10,u4);diff(st10,v4)];
% Adjoint Vectors
lamda1 = Kinv*z1; lamda2 = Kinv*z2; lamda3 = Kinv*z3; lamda4 = Kinv*z4; lamda5 = Kinv*z5;
lamda6 = Kinv*z6; lamda7 = Kinv*z7; lamda8 = Kinv*z8; lamda9 = Kinv*z9; lamda10 = Kinv*z10;
% Derivatives of Element Stresses w.r.t. Areas8
lamda = lamda1;
dstdA1 = lamda'*(dFdA1 - dkdA1*Q); dstdA2 = lamda'*(dFdA2 - dkdA2*Q);
dstdA3 = lamda'*(dFdA3 - dkdA3*Q); dstdA4 = lamda'*(dFdA4 - dkdA4*Q);
dstdA5 = lamda'*(dFdA5 - dkdA5*Q); dstdA6 = lamda'*(dFdA6 - dkdA6*Q);
dstdA7 = lamda'*(dFdA7 - dkdA7*Q); dstdA8 = lamda'*(dFdA8 - dkdA8*Q);
dstdA9 = lamda'*(dFdA9 - dkdA9*Q); dstdA10 = lamda'*(dFdA10 - dkdA10*Q);
dst1dA = vpa([dstdA1;dstdA2;dstdA3;dstdA4;dstdA5;dstdA6;dstdA7;dstdA8;dstdA9;dstdA10],4);
lamda = lamda2;
dstdA1 = lamda'*(dFdA1 - dkdA1*Q); dstdA2 = lamda'*(dFdA2 - dkdA2*Q);
dstdA3 = lamda'*(dFdA3 - dkdA3*Q); dstdA4 = lamda'*(dFdA4 - dkdA4*Q);
dstdA5 = lamda'*(dFdA5 - dkdA5*Q); dstdA6 = lamda'*(dFdA6 - dkdA6*Q);
dstdA7 = lamda'*(dFdA7 - dkdA7*Q); dstdA8 = lamda'*(dFdA8 - dkdA8*Q);
dstdA9 = lamda'*(dFdA9 - dkdA9*Q); dstdA10 = lamda'*(dFdA10 - dkdA10*Q);
dst2dA = vpa([dstdA1;dstdA2;dstdA3;dstdA4;dstdA5;dstdA6;dstdA7;dstdA8;dstdA9;dstdA10],4);
lamda = lamda3;
dstdA1 = lamda'*(dFdA1 - dkdA1*Q); dstdA2 = lamda'*(dFdA2 - dkdA2*Q);
dstdA3 = lamda'*(dFdA3 - dkdA3*Q); dstdA4 = lamda'*(dFdA4 - dkdA4*Q);
dstdA5 = lamda'*(dFdA5 - dkdA5*Q); dstdA6 = lamda'*(dFdA6 - dkdA6*Q);
dstdA7 = lamda'*(dFdA7 - dkdA7*Q); dstdA8 = lamda'*(dFdA8 - dkdA8*Q);
dstdA9 = lamda'*(dFdA9 - dkdA9*Q); dstdA10 = lamda'*(dFdA10 - dkdA10*Q);
dst3dA = vpa([dstdA1;dstdA2;dstdA3;dstdA4;dstdA5;dstdA6;dstdA7;dstdA8;dstdA9;dstdA10],4);
lamda = lamda4;
dstdA1 = lamda'*(dFdA1 - dkdA1*Q); dstdA2 = lamda'*(dFdA2 - dkdA2*Q);
dstdA3 = lamda'*(dFdA3 - dkdA3*Q); dstdA4 = lamda'*(dFdA4 - dkdA4*Q);
dstdA5 = lamda'*(dFdA5 - dkdA5*Q); dstdA6 = lamda'*(dFdA6 - dkdA6*Q);
dstdA7 = lamda'*(dFdA7 - dkdA7*Q); dstdA8 = lamda'*(dFdA8 - dkdA8*Q);
dstdA9 = lamda'*(dFdA9 - dkdA9*Q); dstdA10 = lamda'*(dFdA10 - dkdA10*Q);
dst4dA = vpa([dstdA1;dstdA2;dstdA3;dstdA4;dstdA5;dstdA6;dstdA7;dstdA8;dstdA9;dstdA10],4);
lamda = lamda5;
dstdA1 = lamda'*(dFdA1 - dkdA1*Q); dstdA2 = lamda'*(dFdA2 - dkdA2*Q);
dstdA3 = lamda'*(dFdA3 - dkdA3*Q); dstdA4 = lamda'*(dFdA4 - dkdA4*Q);
dstdA5 = lamda'*(dFdA5 - dkdA5*Q); dstdA6 = lamda'*(dFdA6 - dkdA6*Q);
dstdA7 = lamda'*(dFdA7 - dkdA7*Q); dstdA8 = lamda'*(dFdA8 - dkdA8*Q);
dstdA9 = lamda'*(dFdA9 - dkdA9*Q); dstdA10 = lamda'*(dFdA10 - dkdA10*Q);
dst5dA = vpa([dstdA1;dstdA2;dstdA3;dstdA4;dstdA5;dstdA6;dstdA7;dstdA8;dstdA9;dstdA10],4);
lamda = lamda6;
dstdA1 = lamda'*(dFdA1 - dkdA1*Q); dstdA2 = lamda'*(dFdA2 - dkdA2*Q);
dstdA3 = lamda'*(dFdA3 - dkdA3*Q); dstdA4 = lamda'*(dFdA4 - dkdA4*Q);
dstdA5 = lamda'*(dFdA5 - dkdA5*Q); dstdA6 = lamda'*(dFdA6 - dkdA6*Q);
dstdA7 = lamda'*(dFdA7 - dkdA7*Q); dstdA8 = lamda'*(dFdA8 - dkdA8*Q);
dstdA9 = lamda'*(dFdA9 - dkdA9*Q); dstdA10 = lamda'*(dFdA10 - dkdA10*Q);
dst6dA = vpa([dstdA1;dstdA2;dstdA3;dstdA4;dstdA5;dstdA6;dstdA7;dstdA8;dstdA9;dstdA10],4);
lamda = lamda7;
dstdA1 = lamda'*(dFdA1 - dkdA1*Q); dstdA2 = lamda'*(dFdA2 - dkdA2*Q);
dstdA3 = lamda'*(dFdA3 - dkdA3*Q); dstdA4 = lamda'*(dFdA4 - dkdA4*Q);
dstdA5 = lamda'*(dFdA5 - dkdA5*Q); dstdA6 = lamda'*(dFdA6 - dkdA6*Q);
dstdA7 = lamda'*(dFdA7 - dkdA7*Q); dstdA8 = lamda'*(dFdA8 - dkdA8*Q);
dstdA9 = lamda'*(dFdA9 - dkdA9*Q); dstdA10 = lamda'*(dFdA10 - dkdA10*Q);
dst7dA = vpa([dstdA1;dstdA2;dstdA3;dstdA4;dstdA5;dstdA6;dstdA7;dstdA8;dstdA9;dstdA10],4);
lamda = lamda8;
dstdA1 = lamda'*(dFdA1 - dkdA1*Q); dstdA2 = lamda'*(dFdA2 - dkdA2*Q);
dstdA3 = lamda'*(dFdA3 - dkdA3*Q); dstdA4 = lamda'*(dFdA4 - dkdA4*Q);
dstdA5 = lamda'*(dFdA5 - dkdA5*Q); dstdA6 = lamda'*(dFdA6 - dkdA6*Q);
dstdA7 = lamda'*(dFdA7 - dkdA7*Q); dstdA8 = lamda'*(dFdA8 - dkdA8*Q);
dstdA9 = lamda'*(dFdA9 - dkdA9*Q); dstdA10 = lamda'*(dFdA10 - dkdA10*Q);
dst8dA = vpa([dstdA1;dstdA2;dstdA3;dstdA4;dstdA5;dstdA6;dstdA7;dstdA8;dstdA9;dstdA10],4);
lamda = lamda9;
dstdA1 = lamda'*(dFdA1 - dkdA1*Q); dstdA2 = lamda'*(dFdA2 - dkdA2*Q);
dstdA3 = lamda'*(dFdA3 - dkdA3*Q); dstdA4 = lamda'*(dFdA4 - dkdA4*Q);
dstdA5 = lamda'*(dFdA5 - dkdA5*Q); dstdA6 = lamda'*(dFdA6 - dkdA6*Q);
dstdA7 = lamda'*(dFdA7 - dkdA7*Q); dstdA8 = lamda'*(dFdA8 - dkdA8*Q);
dstdA9 = lamda'*(dFdA9 - dkdA9*Q); dstdA10 = lamda'*(dFdA10 - dkdA10*Q);
dst9dA = vpa([dstdA1;dstdA2;dstdA3;dstdA4;dstdA5;dstdA6;dstdA7;dstdA8;dstdA9;dstdA10],4);
lamda = lamda10;
dstdA1 = lamda'*(dFdA1 - dkdA1*Q); dstdA2 = lamda'*(dFdA2 - dkdA2*Q);
dstdA3 = lamda'*(dFdA3 - dkdA3*Q); dstdA4 = lamda'*(dFdA4 - dkdA4*Q);
dstdA5 = lamda'*(dFdA5 - dkdA5*Q); dstdA6 = lamda'*(dFdA6 - dkdA6*Q);
dstdA7 = lamda'*(dFdA7 - dkdA7*Q); dstdA8 = lamda'*(dFdA8 - dkdA8*Q);
dstdA9 = lamda'*(dFdA9 - dkdA9*Q); dstdA10 = lamda'*(dFdA10 - dkdA10*Q);
dst10dA = vpa([dstdA1;dstdA2;dstdA3;dstdA4;dstdA5;dstdA6;dstdA7;dstdA8;dstdA9;dstdA10],4);
dstdA = double([dst1dA';  dst2dA'; dst3dA'; dst4dA'; dst5dA'; dst6dA'; dst7dA'; dst8dA'; dst9dA'; dst10dA']);
stress_sensitivities = vpa(dstdA,4);
% Derivatives of Mass w.r.t. Areas
dmdA1 = diff(M,A1);dmdA2 = diff(M,A2);dmdA3 = diff(M,A3);dmdA4 = diff(M,A4);dmdA5 = diff(M,A5);
dmdA6 = diff(M,A6);dmdA7 = diff(M,A7);dmdA8 = diff(M,A8);dmdA9 = diff(M,A9);dmdA10 = diff(M,A10);
dwtdA1=0;dwtdA2=0;dwtdA3=0;dwtdA4=0;dwtdA5=0;dwtdA6=0;dwtdA7=0;dwtdA8=0;dwtdA9=0;dwtdA10=0;
for i = 1:length(M)
    for j = 1:length(M)
        dwdA1 = dmdA1(i,j); dwtdA1 = dwtdA1 + dwdA1;
        dwdA2 = dmdA2(i,j); dwtdA2 = dwtdA2 + dwdA2;
        dwdA3 = dmdA3(i,j); dwtdA3 = dwtdA3 + dwdA3;
        dwdA4 = dmdA4(i,j); dwtdA4 = dwtdA4 + dwdA4;
        dwdA5 = dmdA5(i,j); dwtdA5 = dwtdA5 + dwdA5;
        dwdA6 = dmdA6(i,j); dwtdA6 = dwtdA6 + dwdA6;
        dwdA7 = dmdA7(i,j); dwtdA7 = dwtdA7 + dwdA7;
        dwdA8 = dmdA8(i,j); dwtdA8 = dwtdA8 + dwdA8;
        dwdA9 = dmdA9(i,j); dwtdA9 = dwtdA9 + dwdA9;
        dwdA10 = dmdA10(i,j); dwtdA10 = dwtdA10 + dwdA10;
    end
end
dwdA = double([dwtdA1 dwtdA2 dwtdA3 dwtdA4 dwtdA5 dwtdA6 dwtdA7 dwtdA8 dwtdA9 dwtdA10]/rho*rhow);
weight_sensitivities = vpa(dwdA',4);
end
