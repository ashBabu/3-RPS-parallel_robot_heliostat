clear all
close all
clc

% Here the base coordinate system and the global are same and located at
% the base of the Az-El heliostat. Hence O and O1 are the same points

%%
% ALL DIMENSIONS IN METERS
O1G = [0;0;1.58];  % height of the heliostat or vert dist btw bcs and mcs.
% Laser Coordintes.
xl = 6.38;  yl = -1.21;  zl = 2.56;
O1S = [xl;yl;zl];
% O1G + GS = O1S
GS = (O1S - O1G)/norm((O1S - O1G),2) ; % unit incident ray

% Coordintes of the Point to be tracked.
xp = -3.545; yp = 2.74;  zp = 3.10;
O1P = [xp; yp; zp];
GP = (O1P - O1G)/norm((O1P - O1G),2); % unit reflected ray
GN = (GP + GS)/norm((GP+GS),2);  % unit mirror normal
normal = GN;
%% To find out the azimuth and elevation angle of rotation

Motor_angle_Az_El = [];  Az_El_count = [];
N_azimuth= atan2d(normal(2),normal(1));
l=sqrt(normal(2)^2+normal(1)^2);
N_elevation = atand(normal(3)/l);  % ANGLE OF THE MIRRO NORMAL WRT HORIZONTAL.
                %IN THE AZ-EL HELIOSTAT MADE IN THE LAB, THE MOTOR ANGLE
                % SHOULD BE 90-N_elevation
Motor_angle_azimuth = N_azimuth  ;
Motor_angle_elevation = 90-N_elevation ;

counts_per_degree = 9091;
% 1 degree = 9091 counts of the motor
azi_count  = counts_per_degree*Motor_angle_azimuth;
ele_count  = counts_per_degree*Motor_angle_elevation;

Motor_angle_Az_El  = cat(1,Motor_angle_Az_El,[Motor_angle_azimuth Motor_angle_elevation])
Az_El_count = cat(1, Az_El_count, [azi_count -ele_count])  % elevation is multiplied by a -ve sign
                                        % bcoz the Az-El heliostat made
                                        % in the lab rotates in -ve
                                        % direction if a positive value
                                        % is given and vice versa. This
                                        % is bcoz of the direction of
                                        % the axis chosen.

%% GALIL
g = actxserver('galil');%set the variable g to the GalilTools COM wrapper
response = g.libraryVersion;%Retrieve the GalilTools library versions
disp(response);%display GalilTools library version
% g.address = '';%Open connections dialog box
g.address = '192.168.1.1';%Open connections dialog box
g.command('SH');
g.command('TPA=0');
g.command('TPB=0');

i=1;
str1 = ['''PAA=',num2str(Az_El_count(i,1)),''''];
str2 = ['''PAB=',num2str(Az_El_count(i,2)),''''];
g.command(eval(str1))
g.command(eval(str2))
g.command('BGA');
g.command('BGB');

% 
% g.command('PAA=0');
% g.command('PAB=0');
% g.command('BGA');
% g.command('BGB');
% 
% g.command('PAA=50000');
% g.command('PAB=500000');
% g.command('BGA');
% g.command('BGB');

g.command('PRA=-45455');
g.command('BGA');