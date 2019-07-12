% hcst_enviSensors_retrievePastData.m
%
%
%
% Grady Morrissey - July 10, 2019

% Read data from newest csv

% plot

d2 = webread('http://192.168.1.3/19700106.csv')
%Takes 367 seconds per file

figure(501)
plot()
xlabel('Time')
ylabel('Humidity')
set(gca,'FontSize',15)