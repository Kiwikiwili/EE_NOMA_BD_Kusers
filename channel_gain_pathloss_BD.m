%------ generate BD-users channel gains -------% 

%% input: distance (d), pathloss exponent(alpha)
%% output: channel gain (G)

function G = channel_gain_pathloss_BD(d,alpha)

%compute the channel gain according to the pathloss model
G=1./(d.^alpha);

end