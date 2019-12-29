function P0p=top_platform(beta,rp)
 P0p=zeros(4,1,3);
  P0p(1,1,1)=rp;         P0p(2,1,1)=0;      P0p(4,1,1)=1;
  P0p(1,1,2)=rp*cosd(beta);         P0p(2,1,2)=rp*sind(beta);   P0p(4,1,2)=1;
  P0p(1,1,3)=rp*cosd(2*beta);         P0p(2,1,3)=rp*sind(2*beta);     P0p(4,1,3)=1;
 