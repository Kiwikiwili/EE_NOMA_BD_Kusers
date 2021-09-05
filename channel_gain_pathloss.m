%------ generate BS-BD and BS-users channel gains -------% 

%% input: distance (d), pathloss exponent(alpha), noise power(sigma)
%% output: channel gain (G), the arrangement of the channel gains  

function [G,I] = channel_gain_pathloss(d,alpha,sigma)

%compute the channel gain according to the pathloss model
H=1./(d.^alpha*sqrt(sigma));

%sort the channel gains (SIC ordering)
[G,I]=sort(H,'descend');

end