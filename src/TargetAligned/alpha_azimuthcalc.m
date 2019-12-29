function [elevation, azimuth]=alpha_azimuthcalc(LT,d,longitude,latitude)
% the input to this function is
% LT=local time in seconds,
% d is the number of days since the start of the year
% longitude of the place in degrees, East is positive and west is negative
% latitude of the place in degrees, North is positive and south is negative
LSTM = 15*5.5 ; % (in degrees) local standard time meridian. 5.5 is the time differance between local time for india and GMT
B=360*(d-81)/365; % in degrees
%The equation of time (EoT) (in minutes) is an empirical equation that corrects for the eccentricity of the Earth's orbit and the Earth's axial tilt.
EOT=9.87*sind(2*B)-7.53*cosd(B)-1.5*sind(B); % in minutes
TC=4*(longitude-LSTM)+EOT; % Time Correction Factor (TC) in minutes
LST=LT./3600+TC/60; % in hours
HRA=15*(LST-12); % hour angle in degrees
delta=23.45*sind(B); % declination angle in degrees
elevation=asind(sind(delta)*sind(latitude)+cosd(delta)*cosd(latitude).*cosd(HRA));
a1=sind(delta)*cosd(latitude);
a2=cosd(delta)*sind(latitude);
a=acosd((a1-a2.*cosd(HRA))./cosd(elevation));
% for i=1:length(HRA)
%     if (HRA(i)<=0)
%     azimuth(i)=a(i);
% elseif(HRA(i)>0)
%     azimuth(i) =-a(i); %% VERY IMP: for some months like june azimuth = -a(1). otherwise time vs az-El graph will not be continuous
%     end
% end

for i=1:length(HRA)
        if (HRA(i)<=0 )
%         if (1<=d<=115 && 230<=d<=365)   HRA(i)>0 &&
               azimuth(i)=a(i);
        elseif( 116 <= d && d <= 229 && HRA(i)>0)
                azimuth(i) =-a(i); %% VERY IMP: for some months like june azimuth = -a(1). otherwise time vs az-El graph will not be continuous
        else%if(HRA(i)>0 && 1<=d<=115 && 230<=d<=365)
            azimuth(i) =360-a(i);
               
           end
    end
end
