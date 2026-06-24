function dxdt = odefuneuler(t,x,c1,c2,c3)
    dxdt = zeros(7,1);
    dxdt(1) = c1*x(2)*x(3);
    dxdt(2) = c2*x(3)*x(1);
    dxdt(3) = c3*x(1)*x(2);
    q = zeros(4,1);
    q(1) = x(4);
    q(2) = x(5);
    q(3) = x(6);
    q(4) = x(7);
    q = q./norm(q);
    qmul = @(q1,q2)([q1(1),-q1(2),-q1(3),-q1(4); ...
                 q1(2), q1(1),-q1(4),q1(3); ...
                 q1(3), q1(4), q1(1),-q1(2); ...
                 q1(4),-q1(3), q1(2), q1(1)]*q2);
    dqdt = qmul(q,[0;x(1)/2;x(2)/2;x(3)/2]);
    dxdt(4) = dqdt(1);
    dxdt(5) = dqdt(2);
    dxdt(6) = dqdt(3);
    dxdt(7) = dqdt(4);

end