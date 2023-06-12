%% Project 3
% Name: Yuxuan Zhang
% Data: 2023/06/01 (last modification)

clear all; close all;

%% Finding the shortest and longest path
% Loading the data of black hole
load('cluster1.mat');
% Setting the initial value of acceleration (Two dimensions)
a_x = zeros(1,length(hX)); a_y = zeros(1,length(hX));
a_x_sum = 0; a_y_sum = 0;
% Setting the initial value of distance, where the shortest distance is infinity and the longest distance is 0
distance_sum = 0; distance_short = inf; distance_long = 0;
% Setting the initial value of time and derta_t which is the time slots (the increase of time)
t = 0;
derta_t = 0.01;
% Start the random interation we set it about 100000
for interation = 1:1:100000
    % Random x value of the first position of the ship. Store the first position into x0
    x = -5 + 10*rand(1,1); x0 = x;
    % Random v (velocity) value of the initial velocity of the ship. Store it into v0
    v = 2 + 3*rand(1,1); v0 = v;
    % Random angle of velocity. Which is normal distribution form pi/4 to 3pi/4 which is response to x array
    theta = normrnd(pi/2,pi/12);
    % The first position of y is from -10
    y = -10;
    % Using orthogonal decomposition
    v_x = v*cos(theta); v_y = v*sin(theta);
    % Every interation the sum of distance and time need to be clear (set 0)
    distance_sum = 0; t = 0;
    % If the y value is smaller than 10, start find path of the longest or shortest
    while y<=10
        % Calculate the accelaerate of the ship by the force of balck holes
        for i =1:1:length(hX)
            [a_x(i),a_y(i)] = accelerate(hX(i),hY(i),hM(i),x,y,1);
        end
        % clear the value of acceleration in x and y domain and calculate new value
        a_x_sum = 0;a_y_sum = 0;
        for i = 1:1:length(hX)
            a_x_sum = a_x(i)+a_x_sum;
            a_y_sum = a_y(i)+a_y_sum;
        end
        % Calculate the magnitude of acceleration
        a = norm([a_x_sum,a_y_sum],2);
        % If acceleration larger than 4 means failure
        if a>4
            distance_sum = -1;
            break;
        end
        % x distance value larger than 10 means failure
        if abs(x) > 10
            distance_sum = -1;
            break;
        end
        % Calculate the next hop of the new (x,y)
        [x,y,t,v_x,v_y,distance] = Eular_Integration(x,y,a_x_sum,a_y_sum,v_x,v_y,t,derta_t);
        distance_sum = distance_sum+distance;
    end
        % Update the short distance
        if distance_short >= distance_sum && distance_sum > 0
            distance_short = distance_sum;
            % Store the parameter of short distance
            x_short = x0;
            theta_short = theta;
            v_short = v0;
        end
        % Update the long distance
        if distance_long <= distance_sum && distance_sum > 0
            distance_long = distance_sum;
            % Store the parameter of long distance
            x_long = x0;
            theta_long = theta;
            v_long = v0;
        end
end

%% Start draw the scatter of the picture
figure(1);
scatter(hX,hY,'*');axis([-10.1 10.1 -10.1 10.1]);hold on;
xlabel('x = -10 to 10');
ylabel('y = -10 to 10')   

%% Draw the path of shortest distance
% The first position of the ship.
x =  x_short; y = -10;
% Start draw the path
scatter(x,y,40,'ob','filled');hold on; 
% The start velocity of the ship.
v = v_short;
% The start angle of the ship.
theta = theta_short;
v_x = v*cos(theta);
v_y = v*sin(theta);
while y<=10
    % Calculate the accelaerate of the ship by the force of balck holes
    for i =1:1:length(hX)
        [a_x(i),a_y(i)] = accelerate(hX(i),hY(i),hM(i),x,y,1);
    end
    % clear the value of acceleration in x and y domain and calculate new value
    a_x_sum = 0;a_y_sum = 0;
    for i = 1:1:length(hX)
        a_x_sum = a_x(i)+a_x_sum;
        a_y_sum = a_y(i)+a_y_sum;
    end
    % Calculate the magnitude of acceleration
    a = norm([a_x_sum,a_y_sum],2);
    % If acceleration larger than 4 means failure
    if a>4
        distance_sum = -1;
        break;
    end
    % x distance value larger than 10 means failure
    if abs(x) > 10
        distance_sum = -1;
        break;
    end
    % Calculate the next hop of the new (x,y)
    [x,y,t,v_x,v_y,distance] = Eular_Integration(x,y,a_x_sum,a_y_sum,v_x,v_y,t,derta_t);
    distance_sum = distance_sum+distance;
    % Draw the position of ship on the figure
    scatter(x,y,'.g');hold on; 
end
scatter(x,y,40,'or','filled');

%% Draw the path of longest distance
% The first position of the ship.
x =  x_long; y = -10;
% Start draw the path
scatter(x,y,40,'ob','filled');hold on; 
% The start velocity of the ship.
v = v_long;
% The start angle of the ship.
theta = theta_long;
v_x = v*cos(theta);
v_y = v*sin(theta);
while y<=10
% Calculate the accelaerate of the ship by the force of balck holes
    for i =1:1:length(hX)
        [a_x(i),a_y(i)] = accelerate(hX(i),hY(i),hM(i),x,y,1);
    end
    % clear the value of acceleration in x and y domain and calculate new value
    a_x_sum = 0;a_y_sum = 0;
    for i = 1:1:length(hX)
        a_x_sum = a_x(i)+a_x_sum;
        a_y_sum = a_y(i)+a_y_sum;
    end
    % Calculate the magnitude of acceleration
    a = norm([a_x_sum,a_y_sum],2);
    % If acceleration larger than 4 means failure
    if a>4
        distance_sum = -1;
        break;
    end
    % x distance value larger than 10 means failure
    if abs(x) > 10
        distance_sum = -1;
        break;
    end
    % Calculate the next hop of the new (x,y)
    [x,y,t,v_x,v_y,distance] = Eular_Integration(x,y,a_x_sum,a_y_sum,v_x,v_y,t,derta_t);
    distance_sum = distance_sum+distance;
    % Draw the position of ship on the figure
    scatter(x,y,'.y');hold on; 
end
scatter(x,y,40,'or','filled');
hold off;
% Show the detail value of distance of path
title(['The path of shoretst =',num2str(distance_short),'The path of longest =',num2str(distance_long)]);

%% functions
% Thie function named accelerate
% input:    x_blackhole is the x value of black hole position
%           y_blackhole is the y value of black hole position
%           m_blackhole is the mass value of black hole position 
%           x_ship is the x value of ship position
%           y_ship is the y value of ship position
%           gravity 
% output:   a_x is the accelerate in x domain
%           a_y is the accelerate in y domain
%
function [a_x,a_y] = accelerate(x_blackhole,y_blackhole,m_blackhole,x_ship,y_ship,gravity)
    % Distance between black hole and ship
    r = norm([x_blackhole-x_ship,y_blackhole-y_ship],2);
    % The vector with direction
    direction = [x_blackhole-x_ship,y_blackhole-y_ship];
    % Acceleration of the ship by F = ma
    a = (gravity*m_blackhole/(r^3)).*direction;
    a_x = a(1); a_y = a(2);
end
% This function named Eular_Intergration is to calculate the new position of the ship, distance, velocity and time
% input:    x0 is the x value of last position
%           y0 is the y value of last position
%           a_x is the accelerate in x domain
%           a_y is the accelerate in y domain
%           v0_x is the last value of velocity in x array
%           v0_y is the last value of velocity in y array
%           t is time following
%           derta_t time difference
% output:   x is the x value of now position
%           y is the y value of now position
%           t is time following
%           v_x is the now value of velocity in x array
%           v_y is the now value of velocity in y array
%           distance calculated by two position
function [x,y,t,v_x,v_y,distance] = Eular_Integration(x0,y0,a_x,a_y,v0_x,v0_y,t,derta_t)
    % Time increase
    t = t+derta_t;
    % Velocity increase by derta time
    v_x = v0_x + derta_t*a_x;
    v_y = v0_y + derta_t*a_y;
    % Distance increase by derta time
    x = x0 + derta_t*v_x;
    y = y0 + derta_t*v_y;
    % Calculating the distance between them
    distance = norm([x-x0,y-y0],2);
end