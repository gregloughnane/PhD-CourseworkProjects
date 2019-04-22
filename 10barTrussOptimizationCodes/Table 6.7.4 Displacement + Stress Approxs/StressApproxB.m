function [con,ceq] = StressApproxB(x,x0,stress,dstdA,D,dQdA)
con(1) = double(stress(1,1) + sum((x - x0).*dstdA(1,:)));
con(2) = double(stress(2,1) + sum((x - x0).*dstdA(2,:)));
con(3) = double(stress(3,1) + sum((x - x0).*dstdA(3,:)));
con(4) = double(stress(4,1) + sum((x - x0).*dstdA(4,:)));
con(5) = double(stress(5,1) + sum((x - x0).*dstdA(5,:)));
con(6) = double(stress(6,1) + sum((x - x0).*dstdA(6,:)));
con(7) = double(stress(7,1) + sum((x - x0).*dstdA(7,:)));
con(8) = double(stress(8,1) + sum((x - x0).*dstdA(8,:)));
con(9) = double(stress(9,1) + sum((x - x0).*dstdA(9,:)));
con(10) = double(stress(10,1) + sum((x - x0).*dstdA(10,:)));
con(11) = double(D(2,1) + sum((x - x0).*dQdA(2,:)));
con(12) = double(D(4,1) + sum((x - x0).*dQdA(4,:)));
con(13) = double(D(6,1) + sum((x - x0).*dQdA(6,:)));
con(14) = double(D(8,1) + sum((x - x0).*dQdA(8,:)));
con = [ (con(1) - 25e3)/25e3
        (con(2) - 25e3)/25e3
        (con(3) - 25e3)/25e3
        (con(4) - 25e3)/25e3
        (con(5) - 25e3)/25e3
        (con(6) - 25e3)/25e3
        (con(7) - 25e3)/25e3
        (con(8) - 25e3)/25e3
        (con(9) - 25e3)/25e3
        (con(10) - 25e3)/25e3
        (-con(1) - 25e3)/25e3
        (-con(2) - 25e3)/25e3
        (-con(3) - 25e3)/25e3
        (-con(4) - 25e3)/25e3
        (-con(5) - 25e3)/25e3
        (-con(6) - 25e3)/25e3
        (-con(7) - 25e3)/25e3
        (-con(8) - 25e3)/25e3
        (-con(9) - 25e3)/25e3
        (-con(10) - 25e3)/25e3
        (con(11) - 2)/2
        (con(12) - 2)/2
        (con(11) - 2)/2
        (con(12) - 2)/2
        (-con(11) - 2)/2
        (-con(12) - 2)/2
        (-con(11) - 2)/2
        (-con(12) - 2)/2];
ceq=[];
end

