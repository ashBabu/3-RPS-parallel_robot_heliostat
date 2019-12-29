%inputs needed are station co-ordinates longitude,latitude in degrees,min,sec ...
%altitude in kms and time of proposed sun visibility in year,month,day in Christian era,
%hr,min and sec in GMT (NOT BASED ON LOCAL TIME).
%Input the following
clear all
home
close all
clc
global OP t zc Rot rp

%% The following part of the code is used to mark point on the lab wall to show sun vector
% Bangalore Latitude
lad=12;
lam=58;
las=13;

% Bangalore Longitude
lod=77;
lom=33;
los=37;

latitude=lad+(lam+las/60)/60;  % in degrees
longitude=lod+(lom+los/60)/60; % in degrees

year=2016;
month=12;
day=20;
d=date(year,month,day);

time_incre=900; % in seconds or in increments of 15 minutes
LT=0:time_incre:24*3600;

[elevation, azimuth]=alpha_azimuthcalc(LT,d,longitude,latitude);
[xs ,ys, zs]=sunposition(elevation, azimuth);  % xs,ys,zs are coordinates of the sun
OS=[xs;ys;zs];

t_900 = 9/time_incre*3600+1;
t_915 = t_900 + 1;
t_930 = t_900 + 2;
t_945 = t_900 + 3;
t_1000 = 10/time_incre*3600+1;
t_1015 = t_1000 + 1;
t_1030 = t_1000 + 2;
t_1045 = t_1000 + 3;

t_1100 = 11/time_incre*3600+1;
t_1115 = t_1100 + 1;
t_1130 = t_1100 + 2;
t_1145 = t_1100 + 3;

t_1200 = 12/time_incre*3600+1;
t_1215 = t_1200 + 1;
t_1230 = t_1200 + 2;
t_1245 = t_1200 + 3;

t_100 = 13/time_incre*3600+1;
t_115 = t_100 + 1;
t_130 = t_100 + 2;
t_145 = t_100 + 3;

t_200 = 14/time_incre*3600+1;
t_215 = t_200 + 1;
t_230 = t_200 + 2;
t_245 = t_200 + 3;
t_300 = t_200 + 4;
t_315 = t_200 + 5;
t_330 = t_200 + 6;
t_345 = t_200 + 7;
t_400 = 16/time_incre*3600+1;
t_415 = t_400 + 1;
t_430 = t_400 + 2;
t_445 = t_400 + 3;
t_500 = 17/time_incre*3600+1;
t_515 = t_500 + 1;
t_530 = t_500 + 2;

%% Plot of azimuth and elevation angles
figure(1)
plot(LT(t_900:t_415)/3600,elevation((t_900:t_415)),'r','LineWidth',4)
hold on
plot(LT(t_900:t_415)/3600,azimuth((t_900:t_415)),'b','LineWidth',4)
xlabel('Time  \it{t / hr}')
ylabel('Sun vector angles  \it{\theta / (^0)}')
title('Azimuth and Elevation angles of the sun')
lgnd=legend('Elevation',' Azimuth');
grid on

%%
rb=0.400;
rp=0.250;
zc=1.64;
a = 0;  b = 0;  c = 6.72; % wrt TCS
OP = [a;b;c];  % wrt TCS; O is the origin of TCS which is R[10deg, Z] wrt gcs.
OO1 = [-14;3.45;0];
t = OO1;

% Rot = eye(3);
incl = 10;  % angle by which TCS and gcs are rotated
Rot_incl = [cosd(incl) -sind(incl) 0;
            sind(incl)  cosd(incl) 0;
                0           0      1];

O1R1=[rb;0;0];O1R2=[-.5*rb;sqrt(3)*rb/2;0];O1R3=[-.5*rb;-sqrt(3)*rb/2;0];
O1Rs=[O1R1 O1R2 O1R3];

GS1=[rp;0;0];   GS2=[-.5*rp;sqrt(3)*rp/2;0];   GS3=[-.5*rp;-sqrt(3)*rp/2;0];
GSs=[GS1 GS2 GS3];

gamma=0;   % the angle x_b(base coordinate system) makes with local-east
Rot=[ cosd(gamma) -sind(gamma) 0;
      sind(gamma)  cosd(gamma) 0;
          0            0       1]; 
% x0 = [0.2;0.61;0.5;-0.01;0.015];
x0 = [.7;-.6;0.1;-0.7;-0.01;-0.05];
% x0 = [-.7;-.6;-0.1;-0.7;0.01;-0.05];
% x0=[0.25; 0.6; 0.4; 0.35; 0.56];
% x0 = [1.7;.6;-1;.7;.1;.5;0.08;1]; 
options = optimoptions(@fsolve,'Algorithm','levenberg-marquardt');
% options = optimoptions('fsolve','Display','iter','TolX',1e-9);
% options = optimoptions('fsolve','Jacobian','off');

%%
% OS_d = [OS(:,t_700) OS(:,t_800) OS(:,t_900) OS(:,t_1000) OS(:,t_1100) OS(:,t_1200) OS(:,t_100) OS(:,t_200) OS(:,t_300) OS(:,t_400) OS(:,t_500)];% OS(:,t_600)];
OS_d = [OS(:,t_1000) OS(:,t_1015) OS(:,t_1030) OS(:,t_1045) ...
        OS(:,t_1100) OS(:,t_1115) OS(:,t_1130) OS(:,t_1145) ...
        OS(:,t_1200) OS(:,t_1215) OS(:,t_1230) OS(:,t_1245) ...
        OS(:,t_100)  OS(:,t_115) OS(:,t_130) OS(:,t_145) ...
        OS(:,t_200)  OS(:,t_215) OS(:,t_230) OS(:,t_245) ...
        OS(:,t_300)  OS(:,t_315) OS(:,t_330) OS(:,t_345) ...
        OS(:,t_400)  OS(:,t_415) OS(:,t_430) OS(:,t_445) ...
        OS(:,t_500)];

Act_req = [];
% tme = [7 8 9 10 11 12 1 2 3 4 5];
tme = [10 10.15 10.30 10.45 11.00 11.15 11.30 11.45 ...
       12.00 12.15 12.30 12.45 1.00 1.15 1.30 1.45 2.00 ...
       2.15 2.30 2.45 3.00 3.15 3.30 3.45 4.00 4.15 4.30 4.45 5.00];
       
   l10 = 1.645; l20 = 1.653;  l30 = 1.650;

for i=1:length(OS_d)
    fhandle=@(x)for_mainprog_3RPS(x,OS_d(:,i));
    [xval(:,i),fval(:,i)] = fsolve(fhandle,x0,options);
    x0=xval(:,i);
    n1=xval(1,i);n3=xval(2,i);o3=xval(3,i);xc=xval(4,i);yc=xval(5,i);
    B0t(:,i)=[xc;yc;zc];                        
    O1G=[xc;yc;zc];   % w.r.t base coordinate system (bcs)
    O1P_bcs=Rot'*(OP-t);      % O1P_bcs = O1P described in bcs
    g=O1P_bcs-O1G;
    GP=g/norm(g,2);
    GS=Rot'*OS_d(:,i);         % Sun vector in bcs
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
    
    Act_req = cat(1, Act_req, [tme(i), 1000*(l1(i)-l10) 1000*(l2(i)-l20) 1000*(l3(i)-l30)]);
end
% in 1 second the actuator moves 2.4 mm when the speed is 200(bits). 0-0 velocity and 255-max velocity
    Act_tme = [Act_req(:,1) Act_req(:,2)/2.4 Act_req(:,3)/2.4 Act_req(:,4)/2.4 ];
    L = [l1 l2 l3];
    
for i=1:length(Act_tme)-1
   act(i,:) = Act_tme(i+1,:) - Act_tme(i,:);
end
act(:,1) =  [];
act =  cat(1, Act_tme(1,2:end), act);
Act_time = [Act_tme(4,:); Act_tme(end,:); Act_tme(end,:)-Act_tme(4,:) ];


%% Arduino
global ard
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

[rw, cl] = size(act);
[sort_act, Ind_act] = sort(abs(act),2);
for i=1:rw
   tic
   while(toc <= sort_act(i,1))
      RunMotor1(act(i,1)) 
      RunMotor2(act(i,2)) 
      RunMotor3(act(i,3)) 
   end
   MOTORSTOP
%    fprintf(fileID, 'Elapsed Time is %8.3f \n',(toc))
%    pause(1)
%    fprintf(fileID, 'AfterPause,Elapsed Time is %4.2f \n',toc)
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
%      fprintf(fileID, 'Elapsed Time is %4.2f \n',toc)
   pause(910)
%      fprintf(fileID, 'AfterPause,Elapsed Time is %4.2f \n',toc)
end




