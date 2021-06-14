clear;
clc;
clear all;
load data/eg3_0001_0.mat;
u=out.u_y.signal1.data;
y=out.u_y.signal2.data;


for Ny=1:8
    for Nu=1:8
        starting=9;
        % use NARMAX method to solve the linear regression problem
        [yp,pp]=pmatrix(y,u,Ny,Nu);
        [theta,ERR]=NARMAX(yp,pp);
        ERR_ALL{Ny,Nu}=ERR;
        SERR(Ny,Nu)=sum(ERR);
        THETA{Ny,Nu}=theta;
        % estimate y using computed theta y_pre=pp*theta
        y_pre=pp*theta';
        % mean square error of y_pre
        y_diff{Ny,Nu}=yp-y_pre;
        VPE(Ny,Nu)=var(y_diff{Ny,Nu});
        FPE(Ny,Nu)=(length(y)+(Ny+Nu+1))/(length(y)-(Ny+Nu+1))*VPE(Ny,Nu);
    end
end
[X, Y] = meshgrid(1:Ny, 1:Nu);

%PLOT SEER
f1=figure(1);
f1.Position=[100 100 400 900];
subplot(3,1,1);
surf(X,Y,SERR');
grid on
xlabel('X:y-lags');
ylabel('Y:u-lags');
zlabel('SEER');
subplot(3,1,2);
plot(SERR) % plot columns
xlabel('y-lags'); 
grid on;
subplot(3,1,3);
plot(SERR') % plot columns
xlabel('u-lags'); 
grid on;

%PLOT VPE
f2=figure(2);
f2.Position=[600 100 540 900];
subplot(3,1,1);
surf(X,Y,VPE');
grid on;
xlabel('X:y-lags');
ylabel('Y:u-lags');
zlabel('VPE');
subplot(3,1,2);
plot(VPE) % plot columns
xlabel('y-lags'); 
grid on;
subplot(3,1,3);
plot(VPE') % plot columns
xlabel('u-lags'); 
grid on;

%PLOT FPE
f3=figure(3);
f3.Position=[480 100 400 900];
subplot(3,1,1);
surf(X,Y,FPE');
grid on
xlabel('X:y-lags');
ylabel('Y:u-lags');
zlabel('FPE');
subplot(3,1,2);
plot(FPE) % plot columns
xlabel('y-lags'); 
grid on;
subplot(3,1,3);
plot(FPE') % plot columns
xlabel('u-lags'); 
grid on;

% construct y vector and p matrix.
function [yy,pp]=pmatrix(y,u,Ny,Nu)
Nyu=max(Ny,Nu)+1;%start point of y vector
if Nu==0
    pp(:,1)=ones(length(y)-Nyu+1,1);%constant column
    for ny=1:Ny
        pp(:,ny+1)=y(Nyu-ny:end-ny);
    end
    yy=y(Nyu:end);
    p=[yy,pp(:,1:Ny+1)];
else
    pp(:,1)=ones(length(y)-Nyu+1,1);%constant column
    for ny=1:Ny %
        pp(:,ny+1)=y(Nyu-ny:end-ny);
    end
    for nu=1:Nu
        pp(:,nu+1+Ny)=u(Nyu-nu:end-nu) ;
    end
    yy=y(Nyu:end);
   
end
end