%% Project 3
% Name: Yuxuan Zhang
% Data: 2023/06/01 (last modification)

clear all; close all;

%%
% Loading the data of black hole
load('cluster1.mat');

% Setting the initial value of accela
a_x = zeros(1,length(hX)); a_y = zeros(1,length(hX));
a_x_sum = 0; a_y_sum = 0;
distance_sum = 0; distance_short = inf; distance_long = 0;
t = 0;
derta_t = 0.01;

% 
for interation = 1:1:10000
    disp(interation);
    x = -5 + 10*rand(1,1); x0 = x;
    v = 2 + 3*rand(1,1); v0 = v;
    theta = normrnd(pi/2,pi/12);
    y = -10;
    v_x = v*cos(theta);
    v_y = v*sin(theta);
    distance_sum = 0;
    t = 0;
    while y<=10
        %
        for i =1:1:length(hX)
            [a_x(i),a_y(i)] = accellerate(hX(i),hY(i),hM(i),x,y,1);
        end
        a_x_sum = 0;a_y_sum = 0;
        for i = 1:1:length(hX)
            a_x_sum = a_x(i)+a_x_sum;
            a_y_sum = a_y(i)+a_y_sum;
        end
        %
        a = norm([a_x_sum,a_y_sum],2);
        if a>4
            distance_sum = -1;
            break;
        end
        %
        if abs(x) > 10
            distance_sum = -1;
            break;
        end
        %
        [x,y,t,v_x,v_y,distance] = Eular(x,y,a_x_sum,a_y_sum,v_x,v_y,t,derta_t);
        distance_sum = distance_sum+distance;
    end
        if distance_short >= distance_sum && distance_sum > 0
            distance_short = distance_sum;
            x_short = x0;
            theta_short = theta;
            v_short = v0;
        end
        if distance_long <= distance_sum && distance_sum > 0
            distance_long = distance_sum;
            x_long = x0;
            theta_long = theta;
            v_long = v0;
        end
end

%% Start draw the picture of the 
figure(1);
scatter(hX,hY,'*');axis([-10 10 -10 10]);hold on;
xlabel('x = -10 to 10');
ylabel('y = -10 to 10')   


%% Draw the path of shortest distance
x =  x_short;
y = -10;
scatter(x,y,40,'ob','filled');hold on; 
v = v_short;
theta = theta_short;
v_x = v*cos(theta);
v_y = v*sin(theta);
while y<=10
    %
    for i =1:1:length(hX)
        [a_x(i),a_y(i)] = accellerate(hX(i),hY(i),hM(i),x,y,1);
    end
    a_x_sum = 0;a_y_sum = 0;
    for i = 1:1:length(hX)
        a_x_sum = a_x(i)+a_x_sum;
        a_y_sum = a_y(i)+a_y_sum;
    end
    a = norm([a_x_sum,a_y_sum],2);
    if a>4
        distance_sum = -1;
        break;
    end
    %
    if abs(x) > 10
        distance_sum = -1;
        break;
    end
    %
    [x,y,t,v_x,v_y,distance] = Eular(x,y,a_x_sum,a_y_sum,v_x,v_y,t,derta_t);
    distance_sum = distance_sum+distance;
    scatter(x,y,'.g');hold on; 
end
scatter(x,y,40,'or','filled');

%% Draw the path of longest distance
x =  x_long;
y = -10;
scatter(x,y,40,'ob','filled');hold on; 
v = v_long;
theta = theta_long;
v_x = v*cos(theta);
v_y = v*sin(theta);
while y<=10
    %
    for i =1:1:length(hX)
        [a_x(i),a_y(i)] = accellerate(hX(i),hY(i),hM(i),x,y,1);
    end
    a_x_sum = 0;a_y_sum = 0;
    for i = 1:1:length(hX)
        a_x_sum = a_x(i)+a_x_sum;
        a_y_sum = a_y(i)+a_y_sum;
    end
    a = norm([a_x_sum,a_y_sum],2);
    if a>4
        distance_sum = -1;
        break;
    end
    %
    if abs(x) > 10
        distance_sum = -1;
        break;
    end
    %
    [x,y,t,v_x,v_y,distance] = Eular(x,y,a_x_sum,a_y_sum,v_x,v_y,t,derta_t);
    distance_sum = distance_sum+distance;
    scatter(x,y,'.y');hold on; 
end
scatter(x,y,40,'or','filled');
hold off;
title(['The path of shoretst =',num2str(distance_short),'The path of longest =',num2str(distance_long)]);
%% functions
%
function [a_x,a_y] = accellerate(x_star,y_star,m_star,x_ship,y_ship,gravity)
    r = norm([x_star-x_ship,y_star-y_ship],2);
    direction = [x_star-x_ship,y_star-y_ship];
    a = (gravity*m_star/(r^3)).*direction;
    a_x = a(1);
    a_y = a(2);
end
%
function [x,y,t,v_x,v_y,distance] = Eular(x0,y0,a_x,a_y,v0_x,v0_y,t,derta_t)
    t = t+derta_t;
    v_x = v0_x + derta_t*a_x;
    v_y = v0_y + derta_t*a_y;
    x = x0 + derta_t*v_x;
    y = y0 + derta_t*v_y;
    distance = norm([x-x0,y-y0],2);
end