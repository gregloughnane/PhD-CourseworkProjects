function [con,ceq] = StressApproxA(x,x0,stress,dstdA,D,dQdA)
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
ceq(1) = double(D(2,1) + sum((x - x0).*dQdA(2,:)));
ceq(2) = double(D(6,1) + sum((x - x0).*dQdA(6,:)));
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
        (-con(10) - 25e3)/25e3];
ceq=[   (ceq(1) + 2)
        (ceq(2) + 1)];
end

