%time domain pulse for raised cosine, together with time vector to plot it against
%oversampling factor= how much faster than the symbol rate we sample at
%l=where to truncate response (multiple of symbol time) on each side of peak
%a = excess bandwidth
function [rc,time_axis] = raised_cosine(a,m,l)
l_os = floor(l*m); %number of samples on each side of peak
%time vector (in units of symbol interval) on one side of the peak
z = cumsum(ones(l_os,1))/m;
A= sin(pi*z)./(pi*z); %term 1
B= cos(pi*a*z); %term 2
C= 1 - (2*a*z).^2; %term 
zerotest = m/(2*a); %location of zero in denominator
%check whether any sample coincides with zero location
if (zerotest == floor(zerotest))
B(zerotest) = pi*a;
C(zerotest) = 4*a;
end
D = (A.*B)./C; %response to one side of peak
rc = [flipud(D);1;D]; %add in peak and other side
time_axis = [flipud(-z);0;z];
end