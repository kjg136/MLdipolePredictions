function [start,stop]=FindEvents(t,d,DT,RT)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Input
% t : Time 
% d : Dipole
% DT: Drop threshold (as % of mean dipole magnitude)
% RT: Recovery threshold (as % of mean dipole magnitude)
%
%      Output
% start: Time and dipole strength at start of each event
% stop : Time and dipole strength at end of each event
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start=[];
stop=[];

% flag = 1 during event, flag =0 otherwise
if abs(d(1))<DT
    flag=1;
    start=[start;t(1),d(1)];
else
    flag=0;
end

for ii=2:length(t)
    if flag==0 && abs(d(ii))<DT
        start=[start;t(ii),d(ii)];
        flag=1;
    elseif flag==0 && sign(d(ii-1))~=sign(d(ii))
        start=[start;t(ii),d(ii)];
        flag=1;
    elseif flag==1 && abs(d(ii))>=RT
        stop=[stop;t(ii),d(ii)];
        flag=0;
    end
end

        
        