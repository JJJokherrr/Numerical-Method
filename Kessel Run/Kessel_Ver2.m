%%

%%
% Constant
Upper = 10;
Left_MAX = -10;
Right_MAX = 10;

% The stars' position
figure(1);
X = [2];
Y = [3];
m = [0.2];
scatter(X,Y,'*');hold on;
axis([-10 10 -10 10]);
% The ship first position normal distribution
%  x = -5 + 10*rand(1);
x = 0;
y = -10;
%  v = 2 + 3*rand(1);
v = 2;
%  theta = pi/4 + pi/2*rand(1);
theta = pi/2;
v_x = v*cos(theta);
v_y = v*sin(theta);
a_x = zeros(length(X));
a_y = zeros(length(X));
a_x_sum = 0;
a_y_sum = 0;
t = 0;
derta_t = 0.1;

% Simulate the trace of ship

while y<=10 && abs(x)<=10
    %
    for(i = length(X))
        [a_x(i),a_y(i)] = accellerate(X(i),Y(i),m(i),x,y,1);
    end
    for(i = length(X))
        a_x_sum = a_x(i)+a_x_sum;
        a_y_sum = a_y(i)+a_y_sum;
    end
    a = sqrt(a_x_sum^2+a_y_sum^2);
    if a>4
        break;
    end
    %
    [x,y,t,v_x,v_y] = Eular(x,y,a_x_sum,a_y_sum,v_x,v_y,t,derta_t);
    scatter(x,y);axis([-10 10 -10 10]);hold on; 
end
hold off;
%% functions
%
function [a_x,a_y] = accellerate(x_star,y_star,m_star,x_ship,y_ship,gravity)
    r = sqrt((x_star-x_ship)^2+(y_star-y_ship)^2);
    a_x = (gravity.*m_star/(r^3))*(x_star-x_ship);
    a_y = (gravity.*m_star/(r^3))*(y_star-y_ship);
end
%
function [x,y,t,v_x,v_y] = Eular(x0,y0,a_x,a_y,v0_x,v0_y,t,derta_t)
    t = t+derta_t;
    v_x = v0_x + derta_t*a_x;
    v_y = v0_y + derta_t*a_y;
    x = x0 + derta_t*v_x;
    y = y0 + derta_t*v_y;
end