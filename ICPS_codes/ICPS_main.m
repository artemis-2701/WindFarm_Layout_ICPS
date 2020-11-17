close all;
clear all;
clc;
tic;
%[bestSolCost, bestSol,windspeedmatrix,totalpow] = SA(200, 200, 0.1, 0.85, 4, 10)
%[bestSolCost, bestSol,windspeedmatrix,totalpow] =PSO(200,4,10,50)
[bestSol, bestSolCost,windspeedmatrix,totalpow]=GA(4,20,200,0.95,0.05,0.4,10)
toc;

%bestSolCost
%windspeedmatrix
%%windspeed is 12 m/s
%R= 20
%blades are 3
%beta = 2*pi*R/T*windspeed
%T is 2-2.5 sec
%for 3 blades beta(i.e tip speed ratio) is 4-5

positions = [0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 1 0];
count = 1;
for i=1:4
    for j=1:4
        if(bestSol(i,j)==1)
           positions(1,count) = i;
           positions(2,count) = j;
           count = count+1;
        end
    end
end
figure(1);
%plot(positions(1,:),positions(2,:),'r*','Markersize',12.5);
scatter(positions(1,:),positions(2,:),'filled');
grid on;

bestSol = reshape(bestSol.',1,[]); 
%disp(bestSol);
 
windspeedmatrix = reshape(windspeedmatrix.',1,[]);

u = [0 0 0 0 0 0 0 0 0 0];
u = u';
u1 = u;
k = u;
windspeed1 = [12 12 12 12 14 14 14 14 16 16 16 16 18 18 18 18];
count = 1;
N = 4;
for i=1:N*N
    if bestSol(i) == 1
        u(count) = windspeedmatrix(i);
        u1(count) = windspeed1(i);
        k(count) = i;
        count = count+1;
    end
end
%disp(u);
k = k';
disp(k);

theta = [0 0 0 0 0 0 0 0 0 0];
theta = theta';
beta = [0 0 0 0 0 0 0 0 0 0];
beta = beta';
%a = 2*pi*20/2.5;
a = 120;
%disp(a);
for i=1:10
    beta(i)=a/u(i);
end

%beta = a/u;
%disp(beta);
theta = sqrt(((beta-5.6)-0.17*exp(0.07*beta))/0.022);
%disp(theta);
theta = theta';
%disp(theta);
t = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
for i=1:10
    t(k(i)) = theta(i);
end
%disp(t);
t = reshape(t,[4,4]);   
%disp(t);
figure(2)
stem(k,theta);
title('blade pitch angle of each turbine');
xlabel('position of turbine');
ylabel('angle');

A = pi*20^2;
    rho = 1.2;
    Cp = 0.35;
    %Ng = 0.7;
    %Nb = 0.95;
    
    %pwr = 0.5 * rho * A * Cp * vel^3 * Ng * Nb;
    u = u.^3;
    pwr = 0.5 * rho * A * Cp * u ;
    %* Ng * Nb;
    pwr = pwr';
    disp(pwr);
    figure(3);
    stem(k,pwr);
    title('power of each turbine');
    xlabel('position of turbine');
    ylabel('power');
    
%{
u = u.^3;
cp = [1 1 1 1 1 1 1 1 1 1];
cp = 0.35*cp;
pow = 1/2*1.225*pi*20*20*cp*u;
disp(pow);
 
u1 = u1.^3;
pow = 1/2*1.225*pi*20*20*cp*u1;
disp(pow);
%}
