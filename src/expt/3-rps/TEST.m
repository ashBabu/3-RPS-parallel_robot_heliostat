close all
clear all
clc
home
global OP t zc Rot rp ard

% The base of the Az-El heliostat is the origin of the global co-ordinate
% system.
rb=0.400;
rp=0.250;
xl = 6.385;   yl= -1.21; zl = 2.56;  % LASER co-ordinates one fixed on the window frame
zc=1.64;
tx = 10*.3048; ty =0.0; tz=0; t=[tx;ty;tz];
xl_3rps = xl-tx; yl_3rps = yl-ty; zl_3rps = zl-zc;

% [xs ,ys, zs]=sunposition(elevation, azimuth);  % xs,ys,zs are co-ordinates of the sun

OS1=[xl_3rps;yl_3rps;zl_3rps];
unit_sun = OS1/norm(OS1,2);

O1R1=[rb;0;0];O1R2=[-.5*rb;sqrt(3)*rb/2;0];O1R3=[-.5*rb;-sqrt(3)*rb/2;0];
O1Rs=[O1R1 O1R2 O1R3];

GS1=[rp;0;0];   GS2=[-.5*rp;sqrt(3)*rp/2;0];   GS3=[-.5*rp;-sqrt(3)*rp/2;0];
GSs=[GS1 GS2 GS3];

gamma=0;   % the angle x_b(base coordinate system) makes with local-east
Rot=[ cosd(gamma) -sind(gamma) 0;
    sind(gamma)  cosd(gamma) 0;
    0            0       1];
x0 = [0.8;0.61;0.5;-0.01;0.015];
options = optimoptions(@fsolve,'Algorithm','levenberg-marquardt');
% options = optimoptions('fsolve','Display','iter','TolX',1e-9);
% options = optimoptions('fsolve','Jacobian','off');

%%

% Coordintes of the Point to be tracked.
xp = -3.545; yp = 2.74;  zp = 3.10;
OP = [xp; yp; zp];
Act_req =[];
l10 = 1.645; l20 = 1.653;  l30 = 1.650;
% O1P = OP- t;
for i=1:1%length(OS_d)
    %     fhandle=@(x)for_mainprog_3RPS(x,OS_d(:,i));
    fhandle=@(x)for_mainprog_3RPS(x,unit_sun);
    
    [xval(:,i),fval(:,i)] = fsolve(fhandle,x0,options);
    x0=xval(:,i);
    n1=xval(1,i);n3=xval(2,i);o3=xval(3,i);xc=xval(4,i);yc=xval(5,i);
    B0t(:,i)=[xc;yc;zc];
    O1G=[xc;yc;zc];   % w.r.t base coordinate system (bcs)
    O1P_bcs=Rot'*(OP-t);      % O1P_bcs = O1P described in bcs
    g=O1P_bcs-O1G;
    GP=g/norm(g,2);
    %     GS=Rot'*OS_d(:,i);         % Sun vector in bcs
    GS=Rot'*unit_sun;         % Sun vector in bcs
    normal=(GS+GP)/norm((GS+GP),2);
    a1=normal(1,1);   a2=normal(2,1);       a3=normal(3,1);      % wrt global coordinate system
    n2=-yc/rp;
    o1=n2;
    o2=n1-(2*xc/rp);
    Trans_matr = [n1 o1 a1 xc;                               % transformation which takes mirror to base coordinate system
                n2 o2 a2 yc;
                n3 o3 a3 zc;
                0   0  0  1] ;
    B0p1 = Trans_matr*[GS1;1];   % w.r.t base coordinate system
    l1_dir = B0p1-[O1R1;1];
    l1(i,:)=norm(l1_dir,2);
    
    B0p2 = Trans_matr*[GS2;1];
    l2_dir = B0p2-[O1R2;1];
    l2(i,:)=norm(l2_dir,2);
    
    B0p3 = Trans_matr*[GS3;1];
    l3_dir = B0p3-[O1R3;1];
    l3(i,:)=norm(l3_dir,2);
    
    %     Act_req = cat(1, Act_req, [tme(i), 1000*(l1(i)-l10) 1000*(l2(i)-l20) 1000*(l3(i)-l30)]);
    Act_req = cat(1, Act_req, [1000*(l1(i)-l10) 1000*(l2(i)-l20) 1000*(l3(i)-l30)]);
end
   Act_tme = [Act_req(:,1)/2.4 Act_req(:,2)/2.4 Act_req(:,3)/2.4 ];
L = [l1 l2 l3]

act = Act_tme;
[rw, cl] = size(act);
[sort_act, Ind_act] = sort(abs(act),2);



delete(instrfind({'PORT'},{'COM7'}))
ard = arduino('COM7');
% Pins 4,5,6,7,8,9 are PWM pins. Motor 1 uses 4 and 5, Motor 2 uses 6 and 7
% and Motor 3 uses 8 and 9.
% Pins 31,32,33,34,35 are digital pins. Motor 1 uses 31 and NOT gate, Motor
% 2 uses 32 and 33 and Motor 3 uses 34 and 35
% Pins 2,3,18,19,20,21 are used for encoder input signals. Motor 1 uses 2
% and 3, Motor 2 uses 18 and 19 and Motor 3 uses 20 and 21
ard.pinMode(4,'OUTPUT');
ard.pinMode(5,'OUTPUT');
ard.pinMode(6, 'OUTPUT');
ard.pinMode(7,'OUTPUT');
ard.pinMode(8,'OUTPUT');
ard.pinMode(9, 'OUTPUT');

ard.pinMode(30,'OUTPUT');
ard.pinMode(31,'OUTPUT');
ard.pinMode(32,'OUTPUT');
ard.pinMode(33, 'OUTPUT');
ard.pinMode(34,'OUTPUT');
ard.pinMode(35,'OUTPUT');


for i=1:rw
   tic
   while(toc <= sort_act(i,1))
      RunMotor1(act(i,1)) 
      RunMotor2(act(i,2)) 
      RunMotor3(act(i,3)) 
   end
   MOTORSTOP
   tic
   while (toc <= (sort_act(i,2)-sort_act(i,1)))
        if Ind_act(i,2) == 1
            RunMotor1(act(i,1))
        elseif Ind_act(i,2) == 2
            RunMotor2(act(i,2))
        elseif Ind_act(i,2) == 3
            RunMotor3(act(i,3))
        end
        
        if Ind_act(i,3) == 1
            RunMotor1(act(i,1))
        elseif Ind_act(i,3) == 2
            RunMotor2(act(i,2))
        elseif Ind_act(i,3) == 3
            RunMotor3(act(i,3))
        end
   end
   MOTORSTOP
   tic
   while (toc <= (sort_act(i,3)-sort_act(i,2)))
        if Ind_act(i,3) == 1
            RunMotor1(act(i,1))
        elseif Ind_act(i,3) == 2
            RunMotor2(act(i,2))
        elseif Ind_act(i,3) == 3
            RunMotor3(act(i,3))
        end
   end
   MOTORSTOP
   pause
end
