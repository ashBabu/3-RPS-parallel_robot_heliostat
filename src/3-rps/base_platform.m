function B0b=base_platform(alpha,rb)
B0b=zeros(4,1,3);
  B0b(1,1,1)=rb;         B0b(2,1,1)=0;      B0b(4,1,1)=1;
  B0b(1,1,2)=rb*cosd(alpha);         B0b(2,1,2)=rb*sind(alpha);   B0b(4,1,2)=1;
  B0b(1,1,3)=rb*cosd(2*alpha);         B0b(2,1,3)=rb*sind(2*alpha);     B0b(4,1,3)=1;
 