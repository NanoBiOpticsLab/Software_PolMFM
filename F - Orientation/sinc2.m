function y=sinc2(x)
%SINC Sin(x)/(x) function.

i=find(x==0);                                                              
x(i)= 1;      % From LS: don't need this is /0 warning is off                           
y = sin(x)./(x);                                                     
y(i) = 1;   