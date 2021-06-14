clear;
clc;
close all

load data/eg3_continuous_0.mat;
t11=out.u_y.signal2.time;
y11=out.u_y.signal2.data;

load data/eg3_0001_0.mat;
t12=out.u_y.signal2.time;
y12=out.u_y.signal2.data;

load data/eg3_001_0.mat;
t13=out.u_y.signal2.time;
y13=out.u_y.signal2.data;


figure(1);
f1=figure(1);
f1.Position=[100 100 400 900];

plot(t12,y12,'r')
hold on;
plot(t11,y11,'b') % plot columns
% hold on;
% plot(t13,y13,'g')
xlabel('Time/s');
ylabel('Output');
grid on;
legend('Euler approximation with sampling period 0.001s','Continuous system')