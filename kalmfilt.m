%//////////////////////////////////////////////////////////////////////////
% Kalman Filter
% by M. R. Avendi
% The goal is to estimate signal S based on observations X
% we do not have actual data so we generate X based on Kalman model
% and later will estimate S from on X
% X[n]=S[n]+w[n],                   w[n] is N(0,sigma2_w) 
% S[n]= a S[n-1] + u[n],            u[n] is N(0, sigma2_u)

clc
clear all
close all

%% genrate observations according to Kalman model

N=100;% number of data points
a=.9; % state parameter
sigma2_u=1; % variance of state noise
sigma2_w=.9; % variance of observation noise 

S(1)=1; % initial state
X(1)=S(1)+sqrt(sigma2_u)*randn;
for i=2:N+1
    S(i)=a*S(i-1)+sqrt(sigma2_u)*randn;
    X(i)=S(i)+sqrt(sigma2_w^(i-1))*randn;
end

%% estimation of S based on observations X

% initialize
State  = zeros(1,N+1);
StateP =zeros(1,N+1); % prediced state
M =     zeros(1,N+1); % 
Mp=     zeros(1,N+1); % predicted covriance
K=      zeros(1,N+1); % Kalman gain

M(1)=1;
for t=2:N+1
    StateP(t)=a*State(t-1); % predic
    Mp(t)=(a^2)*M(t-1) + sigma2_u; % covariance of prediction
    K(t)=Mp(t)/(sigma2_w^(t-2)+Mp(t));% Kalman gain
    M(t)=(1-K(t))*Mp(t); 
    State(t)=StateP(t)+K(t)*(X(t-1)-StateP(t)); % update 
end
%%
subplot(1,3,1)
plot(K(2:N+1))
legend('Kalman Gain');
xlabel('time, n')
ylabel('K')

subplot(1,3,2)
plot(M(2:N+1),'r');
legend('Mean Square Error');
xlabel('time, n')
ylabel('MSE')

subplot(1,3,3)
plot(S(1:N),'b')
hold on
plot(State(2:N+1),'r')
legend('Actual signal','Estimated signal');
xlabel('time, n')
ylabel('signal')