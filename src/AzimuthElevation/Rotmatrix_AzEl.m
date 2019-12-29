function [rotmat,rotmatZ,eta,Zeta]=Rotmatrix_AzEl(alpha,A,ang)
% a=input('enter the x-cordinate of the receiver\n');
% b=input('enter the y-cordinate of the receiver\n');
% c=input('enter the z-cordinate of the receiver\n');
% a=-60;b=-60;c=120;
%  a=0;b=-6;c=12;
rotmat=zeros(3,3,length(alpha));
[xs, ys, zs]=sunposition(alpha',A');
global  OP O1M OO1 Rot_gam
OM=Rot_gam*O1M+OO1; % wrt gcs
for i=1:length(A)
    OS =[xs(1,i);ys(1,i);zs(1,i)]; % wrt gcs
    MS=OS; % wrt gcs ; unit sunray
    MP=(OP-OM)/norm((OP-OM),2); % wrt gcs ; unit reflected ray
    n = MS + MP;
    normal(:,:,i)=n/norm(n,2); % wrt gcs ;  unit normal
    xn(i)=normal(1,1,i);   yn(i)=normal(2,1,i);       zn(i)=normal(3,1,i); qq(i)=sqrt(xn(i)^2+yn(i)^2);
    eta(i)=atan2d(yn(i),xn(i));                             % azimuth
    Zeta(i)=asind(zn(i)/sqrt(xn(i)^2+yn(i)^2+zn(i)^2));   % elevation
    R1 = [cosd(eta(i)) -sind(eta(i)) 0;
          sind(eta(i))  cosd(eta(i)) 0;
        0            0         1];
    R2 = [cosd(90-Zeta(i)) 0 sind(90-Zeta(i));
                0          1      0  ;
        -sind(90-Zeta(i))  0  cosd(90-Zeta(i))];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    axis = [xn(i);yn(i);zn(i)];
    Skew_axis_matrix = [   0      -axis(3)    axis(2);
                         axis(3)     0        -axis(1);
                        -axis(2)   axis(1)       0];
    Rz = expm((ang*pi/180)*Skew_axis_matrix);
    Ak = R1*R2;
    rotmatZ(:,:,i) = Rz*Ak;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    rotmat(:,:,i) = R1*R2;
    %     rotmat(:,:,i)=[cosd(eta(i))*sind(Zeta(i)) -sind(eta(i)) cosd(eta(i))*cosd(Zeta(i));sind(eta(i))*sind(Zeta(i)) cosd(eta(i)) sind(eta(i))*cosd(Zeta(i));-cosd(Zeta(i)) 0 sind(Zeta(i))];
    %     rotmat(:,:,i)=[normal2(1,1,i)*normal2(3,1,i) -normal2(2,1,i) normal2(1,1,i) ;normal2(2,1,i)*normal2(3,1,i) normal2(1,1,i)  normal2(2,1,i) ;-(normal2(1,1,i)^2+normal2(2,1,i)^2) 0 normal2(3,1,i)];
    %     rotmat(:,:,i)=[xn(i)*zn(i)/qq(i) -yn(i)/qq(i) xn(i) ;yn(i)*zn(i)/qq(i) xn(i)/qq(i)  yn(i);-(xn(i)^2+yn(i)^2)/qq(i) 0 zn(i)];
end

