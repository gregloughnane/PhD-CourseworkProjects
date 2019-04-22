function [K,E,M,rhow,L,c,s] = TenBarTruss()
% This program computes mass & stiffness matrices.
%
% SYNTAX: [[K,E,M,rhow,L,c,s] = TenBarTruss(x)] = TenBarTruss(x)
%   The input is x, corresponding to a vector of the initial areas, A1-A10.
%   The outputs are:
%       K -         Stiffness Matrix
%       M -         Mass Matrix
%       E -         Elastic Modulus (psi)
%       rhow -      Weight Density (lbm)

%% Define Variables
syms A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 K k real
A = [A1 A2 A3 A4 A5 A6 A7 A8 A9 A10];
E = 1e7;           %psi
rho = 0.1 / (32.172 * 12);  % lbm to lbf conversion --> %lb/in^3
rhow = 0.1;        %rho for weight computations
p = 100000;        %lbs
Length = 360;

%% Stiffness & Mass Matrix Assemblies
% K and M for each truss element (local coordinates)
[k1,m1,c1,s1,L1] =  TrussElement1(0,Length,Length,Length);
[k2,m2,c2,s2,L2] =  TrussElement2(Length,Length,2*Length,Length);
[k3,m3,c3,s3,L3] =  TrussElement3(0,0,Length,0);
[k4,m4,c4,s4,L4] =  TrussElement4(Length,0,2*Length,0);
[k5,m5,c5,s5,L5] =  TrussElement5(Length,0,Length,Length);
[k6,m6,c6,s6,L6] =  TrussElement6(2*Length,0,2*Length,Length); 
[k7,m7,c7,s7,L7] =  TrussElement7(0,Length,Length,0);
[k8,m8,c8,s8,L8] =  TrussElement8(0,0,Length,Length);
[k9,m9,c9,s9,L9] =  TrussElement9(Length,Length,2*Length,0);
[k10,m10,c10,s10,L10] = TrussElement10(Length,0,2*Length,Length);
L = [L1 L2 L3 L4 L5 L6 L7 L8 L9 L10];
c = [c1 c2 c3 c4 c5 c6 c7 c8 c9 c10];
s = [s1 s2 s3 s4 s5 s6 s7 s8 s9 s10];
k1=subs(k1); m1=subs(m1); c1=subs(c1); s1=subs(s1); k2=subs(k2); m2=subs(m2); c2=subs(c2); s2=subs(s2);
k3=subs(k3); m3=subs(m3); c3=subs(c3); s3=subs(s3); k4=subs(k4); m4=subs(m4); c4=subs(c4); s4=subs(s4);
k5=subs(k5); m5=subs(m5); c5=subs(c5); s5=subs(s5); k6=subs(k6); m6=subs(m6); c6=subs(c6); s6=subs(s6);
k7=subs(k7); m7=subs(m7); c7=subs(c7); s7=subs(s7); k8=subs(k8); m8=subs(m8); c8=subs(c8); s8=subs(s8);
k9=subs(k9); m9=subs(m9); c9=subs(c9); s9=subs(s9); k10=subs(k10); m10=subs(m10); c10=subs(c10); s10=subs(s10);
% Assign node numbering to each element
nn1 = [5 3]; nn2 = [3 1]; nn3 = [6 4]; nn4 = [4 2]; nn5 = [4 3]; 
nn6 = [2 1]; nn7 = [5 4]; nn8 = [6 3]; nn9 = [3 2]; nn10 = [4 1];
% Assign global DOF (position in stiffness matrix) to each elemental node
indices = [ 2*nn1(1)-1, 2*nn1(1), 2*nn1(2)-1, 2*nn1(2);
            2*nn2(1)-1, 2*nn2(1), 2*nn2(2)-1, 2*nn2(2);
            2*nn3(1)-1, 2*nn3(1), 2*nn3(2)-1, 2*nn3(2);
            2*nn4(1)-1, 2*nn4(1), 2*nn4(2)-1, 2*nn4(2);
            2*nn5(1)-1, 2*nn5(1), 2*nn5(2)-1, 2*nn5(2);
            2*nn6(1)-1, 2*nn6(1), 2*nn6(2)-1, 2*nn6(2);
            2*nn7(1)-1, 2*nn7(1), 2*nn7(2)-1, 2*nn7(2);
            2*nn8(1)-1, 2*nn8(1), 2*nn8(2)-1, 2*nn8(2);
            2*nn9(1)-1, 2*nn9(1), 2*nn9(2)-1, 2*nn9(2);
            2*nn10(1)-1, 2*nn10(1), 2*nn10(2)-1, 2*nn10(2)];
% K for the structure (global)
K = sym(zeros(12,12));
k = [k1 k2 k3 k4 k5 k6 k7 k8 k9 k10]; % Define elemental stiffnesses in a single matrix for ease of assembly
for i = 1:10
    index = indices(i,:);
    start = 4*i-3;
    stop = 4*i;
    K(index,index) = K(index,index) + k(:,start:stop);
end
K_reduced = K(1:8,1:8);K = K_reduced;    
% M for the structure (global)
M = sym(zeros(12,12));
m = [m1 m2 m3 m4 m5 m6 m7 m8 m9 m10]; % Define elemental stiffnesses in a single matrix for ease of assembly
for i = 1:10
    index = indices(i,:);
    start = 4*i-3;
    stop = 4*i;
    M(index,index) = M(index,index) + m(:,start:stop);
end
M_reduced = M(1:8,1:8);M = M_reduced;             
end

