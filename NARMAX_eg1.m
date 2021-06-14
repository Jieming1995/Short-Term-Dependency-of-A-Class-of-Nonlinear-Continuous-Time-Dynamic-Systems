clear;
clc;
close all
load data/eg1_0.01_0.0001.mat;
u=out.u_y.signal1.data;
y=out.u_y.signal2.data;


for Ny=1:100
    Nu=1;
        % use NARMAX method to solve the linear regression problem
        [yp,pp]=pmatrix(y,u,Ny,0);
        [theta,ERR]=NARMAX(yp,pp);
        ERR_ALL{Ny,Nu}=ERR;
        SERR(Ny,Nu)=sum(ERR);
        THETA{Ny,Nu}=theta;
        % estimate y using computed theta y_pre=pp*theta
        y_pre=pp*theta';
        % mean square error of y_pre
        y_diff{Ny,Nu}=yp-y_pre;
        VPE(Ny,Nu)=var(y_diff{Ny,Nu});
        FPE(Ny,Nu)=(10000+(Ny+Nu+1))/(10000-(Ny+Nu+1))*VPE(Ny,Nu);
    
end
%PLOT SEER
figure(1);
f1=figure(1);
f1.Position=[100 100 400 900];
plot(SERR,'-o') % plot columns
xlabel('y-lags');
ylabel('SEER');
grid on;


%PLOT VPE
figure(2);
f2=figure(2);
f2.Position=[600 100 540 900];
plot(VPE,'-o') % plot columns
xlabel('y-lags');
ylabel('VPE');
grid on;

%PLOT FPE
figure(3);
f3=figure(3);
f3.Position=[480 100 400 900];
plot(FPE,'-o') % plot columns
xlabel('y-lags');
ylabel('FPE');
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