function d=date(year,month,day)
if(rem(year,4)==0)
    if (month==1)
        d=day;
    elseif(month==2)
        d=day+31;
    elseif(month==3)
        d=day+60;
        elseif(month==4)
        d=day+91;
        elseif(month==5)
        d=day+121;
        elseif(month==6)
        d=day+152;
        elseif(month==7)
        d=day+182;
        elseif(month==8)
        d=day+213;
        elseif(month==9)
        d=day+244;
        elseif(month==10)
        d=day+274;
        elseif(month==11)
        d=day+305;
        elseif(month==12)
        d=day+335;
    end
end


if (rem(year,4)~=0)
    if (month==1)
        d=day;
    elseif(month==2)
        d=day+31;
    elseif(month==3)
        d=day+59;
        elseif(month==4)
        d=day+90;
        elseif(month==5)
        d=day+120;
        elseif(month==6)
        d=day+151;
        elseif(month==7)
        d=day+181;
        elseif(month==8)
        d=day+212;
        elseif(month==9)
        d=day+243;
        elseif(month==10)
        d=day+273;
        elseif(month==11)
        d=day+304;
        elseif(month==12)
        d=day+334;
    end
end


