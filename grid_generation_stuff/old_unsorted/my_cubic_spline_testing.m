%% my_cubic_spline testing
clc; clear; close all;

t = linspace(0,2*pi,25);

f = @(x) sin(x);
df = @(x) cos(x);

s = my_cubic_spline(t,f(t),1.0e99,1.0e99);
s1 = my_linear_spline(t,f(t));


hold on;
fplot(f,[0,2*pi],'k')
plot(t,s.eval(t),'b')
plot(t,s.eval(t),'b--')

fplot(df,[0,2*pi],'r')
plot(t,s.deval(t),'g')
plot(t,s1.deval(t),'g--')
