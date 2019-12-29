function F=Refl_ray(r,OS,N)


r1=r(1);r2=r(2);r3=r(3);


F=[ r1^2 + r2^2 + r3^2 - 1;
    cross(OS,N)- cross(N,r);
    dot(OS,N)- dot(N,r);
    ];



