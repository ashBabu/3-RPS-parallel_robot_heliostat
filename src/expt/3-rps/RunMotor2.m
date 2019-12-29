function RunMotor2(act_time)
global ard
% global fileID
%    fileID = fopen('EXP.txt','w');
   if act_time <0
%         disp('Motor 2 back')
%          fprintf(fileID,'Motor 2 back\n');
        % bwd
        ard.digitalWrite(32,0);
        ard.analogWrite(6,0);
        ard.digitalWrite(33,1);
        ard.analogWrite(7,200);
    else
%         disp('Motor 2 fwd')
%           fprintf(fileID,'Motor 2 fwd \n');
        % fwd
        ard.digitalWrite(33,0);
        ard.analogWrite(7,0);
        ard.digitalWrite(32,1);
        ard.analogWrite(6,200);
    end
