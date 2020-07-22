% Adapted for ChE453 from:
% matlab.cheme.cmu.edu/2011/08/09/phase-portraits-of-a-system-of-odes/
% 
% Inputs:
% s0: Initial susceptible population
% i0: Initial infected population
% tn: Time period (30-60 days for this disease)
%
% Enter inputs separated by commas. Example for CaseA:
% Epidemic_phase_portrait(100, 1, 0, 45)

function Epidemic_phase_portrait(s0, i0, r0, tn)

clc; close all;

totalPeople = s0 + i0 + r0;

iRate = 0.01;
rRate = 0.2;

range = totalPeople; % Confine our vector plot to a square region of + and - range.

sMin = 0; smax = range;
iMin = 0; imax = range; 
rMin = 0; rmax = range; 
npoints = 15;          % number of grid points to plot the portrait

s = linspace(sMin,smax,npoints);
i = linspace(iMin,imax,npoints);
r = linspace(rMin,rmax,npoints);

[x,y,z] = meshgrid(s,i,r);

%--------------------------------
% Step 2: Plot the vector field
%--------------------------------

% u and v are matrices of zeros with (npoints x npoints)elements
% We will store the RHS of our ODEs in these variables
u = zeros(size(x));
v = zeros(size(x));

% we can use a single loop over each element to compute the derivatives at
% each point (y1, y2)

t=0; % we want the derivatives at each point at the designated time
for i = 1:numel(x)
    Yprime = illRelations(t,[totalPeople - y(i) - z(i); y(i); z(i)]); 
    % x(i) is determined by totalPeople - y(i) - z(i) either than it self
    
    u(i) = Yprime(2);
    v(i) = Yprime(3);
end

    % Notice the use of the subplot function to plot 2 graphs (the other
    % comes shortly below) on the same window.  We use 'hold on' for the
    % first one since we will plot the solution of the ODE on the same
    % plot. The 'quiver' function plots a vector field.
    
    figure(1)
    subplot(2,1,2)
    
    hold on
    quiver(y,z,u,v,1,'r');
    xlabel('Infected')
    ylabel('Recovered')
    axis tight;
    
    figure(2)
    hold on
    quiver(y,z,u,v,1,'r');
    xlabel('Infected')
    ylabel('Recovered')
    axis tight;

    to =  0; % Initial time
     
     
    % Invoke the ODE solver here
    [ts,ys] = ode45(@illRelations,[to tn],[s0 i0 r0]);
    
    % Plot the solution of the ODE on the vector plot
    figure(1)
    subplot(2,1,2)
    plot(ys(:,2),ys(:,3),'g-')
    plot(ys(1,2),ys(1,3),'bo') % starting point
    plot(ys(end,2),ys(end,3),'ks') % ending point
    legend('derivative vector field','trajectory','start','end');
    hold off

     subplot(2,1,1)
     plot(ts,ys(:,2), ts,ys(:,3));
     xlabel('time')
     ylabel('Infected, Recovered')
     legend('Infected','Recovered');

    figure(2)
    plot(ys(:,2),ys(:,3),'g-')
    plot(ys(1,2),ys(1,3),'bo') % starting point
    plot(ys(end,2),ys(end,3),'ks') % ending point
    legend('derivative vector field','trajectory','start','end');
    hold off
  

    function dnumberOfPeople = illRelations(t,numberOfPeople)
        
        dnumberOfPeople = zeros(3,1);
        
        %Unpacking
        suspectible = numberOfPeople(1);
        infected = numberOfPeople(2);
        recovered = numberOfPeople(3);
        
        %Relations
        dSdt = - iRate*suspectible*infected; 
        dIdt = iRate*suspectible*infected - rRate*infected; 
        dRdt = rRate*infected;
        
        %Packing
        dnumberOfPeople(1) = dSdt;
        dnumberOfPeople(2) = dIdt;
        dnumberOfPeople(3) = dRdt;
        
    end
end