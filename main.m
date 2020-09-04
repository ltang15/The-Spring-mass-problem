%{Midterm: The Spring mass problem}%

% Description: Plotting two kinematic vectors: position and velocity over time, 
% and the total ennergy over time using three methods: 
% Forward Euler (explicit method), Backward Euler(implicit method,
% semi-implicit method and 4th order Runge Kutta method. 


clc, clear all, close all;

string = ['Select the Approximation Method:\n\n', ...
    '1. Forward Euler\n',...
    '2. Backward Euler\n',... 
    '3. Semi-implicit method\n',...
    '4. 4th order Runge Kutta method\n',...
    'Enter your selection\n'];


sel = input(string); % display the menu

%array holding constant values: first element is the mass, second element is the spring constant
cons = [5, 180]; 

m = cons(1);% mass of the object
k = cons(2);% spring constant

dt = 0.0005;


ti = 0; %initial time
tf = 10; %final time
nt = (tf-ti)/dt; %number of time steps
t = linspace (ti, tf, nt);% create an array of time

x = zeros(1, nt);% initialize the array of position 
v = zeros(1, nt);% initialize the array of  velocity
e = zeros(1, nt);% initialize the array of total energy

x(1) = 0.1;% initial position
e(1) = 1/2*(cons(2)*x(1)^2 + cons(1)*v(1)^2);%initial energy, there's only potential energy at the beginning

switch (sel)
    
    case 1
        %% Forward Euler
     
        for i = 1:nt-1
            [v(i+1),x(i+1)] = Euler (1, v(i), x(i), dt,cons);
            e(i+1) = 1/2*(k*x(i+1)^2 + m*v(i+1)^2); %calculate total enery at the next time step
        end 
      
        
    case 2
         %% Backward Euler 
          
        for i = 1:nt-1
             [v(i+1),x(i+1)] = Euler (2, v(i), x(i), dt,cons);
              e(i+1) = 1/2*(k*x(i+1)^2 + m*v(i+1)^2);
        end
        
    
    case 3
        %% Semi-implicit method
        
        for i=1:nt-1
            [v(i+1),x(i+1)] = Euler(3, v(i), x(i), dt,cons);
             e(i+1) = 1/2*(k*x(i+1)^2 + m*v(i+1)^2);
        end 
        
    case 4  
        %% 4th order Runge Kutta method
        
        for i=1:nt-1
            
            [v(i+1),x(i+1)] = Euler(4, v(i), x(i), dt,cons);
            e(i+1) = 1/2*(k*x(i+1)^2 + m*v(i+1)^2);
        end    
        
end        


%% Position and velocity vs time plot

figure(1);
pos = plot(t,x,'b',t,v,'r');
xlabel ('time (s)');
ylabel ('x(m), v(m/s)');
title ('Position(x) and velocity(v) over time');
xlim ([0 10]);
ylim ([1.1*min(v) 1.1*max(v)]);
set(pos, 'LineWidth', 2);
set(gcf, 'Position', [100 40 1000 600]);
set(gca, 'LineWidth', 2, 'FontSize', 15);
legend('x vs. t','v vs. t');
legend('Location','southeast')
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

%% Total energy over time plot

figure(2);
phase = plot (t, e);
xlabel ('t (s)');
ylabel ('E (J)');
title ('The total energy vs. time(t)');
xlim ([0 10]);
ylim ([0 1.1*max(e)]);
set(phase, 'LineWidth', 2);
set(gcf, 'Position', [200 40 900 600]);
set(gca, 'LineWidth', 2, 'FontSize', 15);
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';
yticks(0:0.1:1.2)

function [vkp1, xkp1] = Euler(method, vk, xk, dt,constants)
%EULER Calculate the position and velocity at the next time step
%   Method is either 1,2,3 or 4 corresponding to the explicit, implicit, 
%   semi-implicit and 4th order Runge Kutta method respectively. 
%   vk and xk are the velocity and position at the previous time step
%   dt is the time step; constant is an array holding 2 constants: mass of the
%   object(m) and the spring constant (k).



    m = constants(1);
    k = constants(2);
    a = (-k/m)*0.1; %initial acceleration is calculated based on x(1)= 0.1
    
    switch(method)
        
        case 1
         %% Forward Euler (explicit)
         
            vkp1 = vk + dt*(-k/m)*xk;
            xkp1 = xk + dt*vk;
            
        case 2         
         %% Backward Euler (implicit)
            
            a = (-k/m)*(xk + vk*dt);
            vkp1 = vk + dt*a;
            xkp1 = xk + dt*vkp1;
           
        case 3  
         %% Implicit method
         
            vkp1 = vk + dt*(-k/m)*xk;
            xkp1 = xk + dt*vkp1;
        
        case 4
          %% 4th order Runge Kutta method
          
           dvdt= @(x) -k*x/m;
           dxdt= @(v) v;
           
            
            c1v = dt*dvdt(xk);
            c2v = dt*dvdt(xk + 0.5*c1v);
            c3v = dt*dvdt(xk + 0.5*c2v);
            c4v = dt*dvdt(xk + c3v);
    
            c1x = dt*dxdt(vk);
            c2x = dt*dxdt(vk + 0.5*c1x);
            c3x = dt*dxdt(vk + 0.5*c2x);
            c4x = dt*dxdt(vk + c3x);
            
           vkp1 = vk + c1v/6 + c2v/3 + c3v/3 + c4v/6;
           xkp1 = xk + c1x/6 + c2x/3 + c3x/3 + c4x/6;
          
    end
    
end
