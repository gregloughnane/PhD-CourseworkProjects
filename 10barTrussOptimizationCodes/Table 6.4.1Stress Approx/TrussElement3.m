function[k_element,m_element,c,s,L] = TrussElement3(x1,y1,x2,y2)
% Determines elemental stiffness matrices for truss elements in local
% coordinates 

syms A3 E I theta rho beta

% Length of bar
L = sqrt((x2-x1)^2 + (y2-y1)^2);

% Cosine of bar angle
c = (x2-x1)/L;

% Sine of bar angle
s = (y2-y1)/L;

% Elemental stiffness
k_element = A3*E/L* [c^2  c*s  -c^2  -c*s
                    c*s  s^2  -c*s  -s^2
                    -c^2  -c*s  c^2  c*s
                    -c*s  -s^2  c*s  s^2];                
% Elemental mass
m_element = rho*A3*L/6* [2 0 1 0
                        0 2 0 1
                        1 0 2 0
                        0 1 0 2];
                
end