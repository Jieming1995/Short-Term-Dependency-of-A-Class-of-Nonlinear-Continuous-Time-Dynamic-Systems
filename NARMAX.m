function [theta,err]=NARMAX(y,p)
% We are looking for theta such that y=p*theta

% w_norm:normalized candidate basis, every column in this matrix is orthoganal to all selected bases. 
% w_index: original index of the columnns in w_norm
% w0: selected normalized bases, orthogonal to each other
% g0: parameters of w0
% ell: index of w0 in the original matrix
%% initialization
M=size(p,2);
sigma = vecnorm(y);
w_norm=p;
w_index=1:M;
j=0;
while ~isempty(w_index)
    j=j+1;
    % find the closest column in w_norm to y
    cos_angle=(y'*w_norm./(sigma*vecnorm(w_norm))).^2;
    [ERR(j),l_w_norm]=max(cos_angle);
    ell(j)=w_index(l_w_norm);
    w0(:,j)=w_norm(:,l_w_norm);
    g0(j)=y'*w0(:,j)/(w0(:,j)'*w0(:,j));
    % update w_norm and w_index
    % 1. eliminate the chosen column
    % 2. normaliz the rest columns
    % 3. eleminate ~0 columns
    w_norm(:,l_w_norm)=[];
    w_index(:,l_w_norm)=[];
    w_norm=w_norm - w0(:,j)'*w_norm/(w0(:,j)'*w0(:,j)).*w0(:,j);
    zero_column_index=vecnorm(w_norm)<1e-6;
    w_norm(:,zero_column_index)=[];
    w_index(:,zero_column_index)=[];
end
%%% Compute A matrix
M0=length(g0);
A=eye(M0);
for j=2:M0
    A(1:j-1,j)=p(:,ell(j))'*w0(:,1:j-1)./(vecnorm(w0(:,1:j-1)).^2);
end
Theta=A\g0'; % every element corresponds to g0;
theta=zeros(1,M);
err=zeros(1,M);
theta(ell)=Theta;
err(ell)=ERR;
%SERR=sum(ERR);