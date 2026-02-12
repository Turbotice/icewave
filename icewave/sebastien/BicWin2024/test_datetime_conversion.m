drone_name = 'mesange';

Y = 2024;
M = 02;
D = 26;
H = 20; % local time of drone
MIN = 22;
S = 00;
MS = 509;
 
if strcmp(drone_name,'mesange') & Y == 2024
    TimeZone = 'Europe/Paris'; % mesange 
else 
    TimeZone = 'America/Montreal'; % bernache or Fulmar
end 

% initial time of recording
t0_UTC = datetime(Y,M,D,H,MIN,S,MS,'TimeZone',TimeZone); 
t0_UTC.TimeZone = 'UTC'; % converts time to UTC time 
t0_UTC.Format = 'yyyy-MM-dd HH:mm:ss.SSS';

t0_txt = char(t0_UTC);

t = milliseconds((1:10));
t0 = t0_UTC;
t0.Format = 'HH:mm:ss.SSS';
UTC_t = t0 + t;