function RunMotor3(act_time)
global ard
% global fileID
%    fileID = fopen('EXP.txt','w');
   if act_time <0
%         disp('Motor 3 back')
%         fprintf(fileID,'Motor 3 back\n');
        % bwd
        ard.digitalWrite(34,0);
        ard.analogWrite(8,0);
        ard.digitalWrite(35,1);
        ard.analogWrite(9,200);
    else
%         disp('Motor 3 fwd')
%          fprintf(fileID,'Motor 3 fwd\n');
        % fwd
        ard.digitalWrite(35,0);
        ard.analogWrite(9,0);
        ard.digitalWrite(34,1);
        ard.analogWrite(8,200);
    end