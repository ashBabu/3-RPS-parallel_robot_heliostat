function F=for_mainprog_3RPS_test(x,A,O1G)
global  rp 

n1=x(1);n3=x(2);o3=x(3);    xc=O1G(1);yc=O1G(2);

n2=-yc/rp;
o1=n2;
o2=n1-(2*xc/rp);

F=[ n1^2    + n2^2    + n3^2    - 1;
    o1^2    + o2^2    + o3^2    - 1;
    A(1)*n1 + A(2)*n2 + A(3)*n3;
%     n1*o1   + n2*o2   + n3*o3;
%     A(1)*o1 + A(2)*o2 + A(3)*o3
    ];



