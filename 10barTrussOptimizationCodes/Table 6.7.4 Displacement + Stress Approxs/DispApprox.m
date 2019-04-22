function [dcon] = DispApprox(x,x0,D,dQdA)
clear A1 A2 A3 A4 A5 A6 A7 A8 A9 A10
syms A1 A2 A3 A4 A5 A6 A7 A8 A9 A10
xi = [A1 A2 A3 A4 A5 A6 A7 A8 A9 A10];
%% Linear Constraint Approximations For Displacements
% subtract = [ D(2,1) + sum((x - x0).*dQdA(2,:)'); %v1
%              D(4,1) + sum((x - x0).*dQdA(4,:)'); %v2
%              D(6,1) + sum((x - x0).*dQdA(6,:)'); %v3
%              D(8,1) + sum((x - x0).*dQdA(8,:)')]; %v4
subtract = [ D(2,1) + sum((xi - x).*dQdA(2,:)); %v1
             D(4,1) + sum((xi - x).*dQdA(4,:)); %v2
             D(6,1) + sum((xi - x).*dQdA(6,:)); %v3
             D(8,1) + sum((xi - x).*dQdA(8,:))]; %v4
dcon = subs(subtract,{A1,A2,A3,A4,A5,A6,A7,A8,A9,A10},[0 0 0 0 0 0 0 0 0 0]);
end


