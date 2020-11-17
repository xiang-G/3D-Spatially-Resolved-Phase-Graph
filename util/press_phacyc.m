function [phs1,phs2,phs3]=press_phacyc(nth,Nphaccyc)
switch Nphaccyc
    case 1
        phs1=0;phs2=90;phs3=90;
    case 2
        if nth==1; phs1=180;phs2=0;phs3=180;
        else;  phs1=0;phs2=0;phs3=0;end
    case 8
        % Hennig, JMR 96, 40-49, 1992
        if nth==1; phs1=90;phs2=0;phs3=90;end
        if nth==2; phs1=90;phs2=90;phs3=180;end
        if nth==3; phs1=90;phs2=180;phs3=270;end
        if nth==4; phs1=90;phs2=270;phs3=0;end
        if nth==5; phs1=270;phs2=0;phs3=0;end
        if nth==6; phs1=270;phs2=90;phs3=90;end
        if nth==7; phs1=270;phs2=180;phs3=180;end
        if nth==8; phs1=270;phs2=270;phs3=270;end
%     case 16
%         if nth<=4; phs1=wrapTo360(360/4*nth);phs2=0;phs3=0;
%         elseif nth<=8; phs1=wrapTo360(360/4*nth);phs2=0;phs3=180;
%         elseif nth<=12; phs1=wrapTo360(360/4*nth);phs2=180;phs3=0;
%         elseif nth<=16; phs1=wrapTo360(360/4*nth);phs2=180;phs3=180;end
    otherwise
        error('we do not support other Nphaccyc number for now');
end
end
