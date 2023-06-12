%%

%%
clear all; close all;
% Constant
load('cluster1.mat')
Upper = 10;
Left_MAX = -10;
Right_MAX = 10;

% The stars' position


% The ship first position normal distribution
x = 0;
y = -10;
v = normrnd(3.5,0.5,1,1);
theta = pi/4:0.1:3*pi/4;
v_x = v*cos(theta);
v_y = v*sin(theta);
a_x = zeros(1,length(hX));
a_y = zeros(1,length(hX));
a_x_sum = 0;
a_y_sum = 0;
distance_sum = 0;
distance_short = 99999;
distance_long = 0;
t = 0;
derta_t = 0.01;

% Simulate the trace of ship
for x0 = -5:0.01:5
    x0
    for theta = pi/4:0.01:3*pi/4
        for v = 2:1:5
        x = x0;
        %v = 2;
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
            v_short = v;
        end
        if distance_long <= distance_sum && distance_sum > 0
            distance_long = distance_sum;
            x_long = x0;
            theta_long = theta;
            v_long = v;
        end
        end
    end
end
%
figure(1);
scatter(hX,hY,'*');hold on;
%%
x =  x_short;
y = -10;
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
    scatter(x,y,'.b');axis([-10 10 -10 10]);hold on; 
end
%%
x =  x_long;
y = -10;
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
    end
    %
    if abs(x) > 10
        distance_sum = -1;
        break;
    end
    %
    [x,y,t,v_x,v_y,distance] = Eular(x,y,a_x_sum,a_y_sum,v_x,v_y,t,derta_t);
    distance_sum = distance_sum+distance;
    scatter(x,y,'.r');axis([-10 10 -10 10]);hold on; 
end
hold off;

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