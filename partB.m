% Inputs: vacination percentage, death rate
% Input for Case A: (0,0)
% Input for Case B1: (0,0.2)
% Input for Case B2: (20,0.2), (50,0.2), (80,0.2)

function partB(vac,k3)
clc;

%Timeline limits
to = 0; 
tn = 45;

% Initial values for number of people
susceptible = 100;
infected = 1;
recovered = 0;
dead = 0;
totalPopulation = susceptible + infected + recovered;




% Infection and recovery rates
k1 = 0.01; 
k2 = 0.2;

% Vacination percentage
vacPer = vac/100;

% Re calculation after vaccination
vaccinated = susceptible * vacPer; % This is important to be before updating susceptible
susceptible = susceptible * ( 1 - vacPer);



% Actual solution implemented here
[tsoln,numberOfPeople] = ode45(@illRelations,[to tn],[susceptible infected recovered dead]);



%plot the soln

if k3 > 0   % with death
    plot(tsoln,numberOfPeople(:,1:2), tsoln,numberOfPeople(:,3), 'k', tsoln, numberOfPeople(:,4), 'dr', 'MarkerSize', 5)
    legend('Susceptible Individuals','Infected Individuals','Recovered Individuals','Dead Individuals');
    xlabel('Days')
    ylabel('Number of People')
    ylim([0 totalPopulation])
    
else % without death
    plot(tsoln,numberOfPeople(:,1:3))
    legend('Susceptible Individuals','Infected Individuals','Recovered Individuals');
    xlabel('Days')
    ylabel('Number of People')
    
end

% Change in sub-populations with time
%disp('Time(days)  Susceptible Infected Recovered  Dead')
%disp([tsoln numberOfPeople])

    function dnumberOfPeople = illRelations(t,numberOfPeople)
        dnumberOfPeople = zeros(4,1);
        
        % Unpacking
        susceptible = numberOfPeople(1);
        infected = numberOfPeople(2);
        recovered = numberOfPeople(3);
        dead = numberOfPeople(4);

        %Relations
        dSdt = - k1*infected*susceptible; 
        dIdt = k1*infected*susceptible - k2*infected - k3*infected; 
        dRdt = k2*infected;
        dDdt = k3*infected;
        
        %Packing
        dnumberOfPeople(1) = dSdt;
        dnumberOfPeople(2) = dIdt;
        dnumberOfPeople(3) = dRdt;
        dnumberOfPeople(4) = dDdt;
        
    end

Rmax = max(numberOfPeople(:,3));
[Rmax time] = max(numberOfPeople(:,3));

susceptibleList = numberOfPeople(:,1);
infectedList = numberOfPeople(:,2);
recoveredList = numberOfPeople(:,3);
deadList = numberOfPeople(:,4);

% Table data for applied case 

initialSusceptible = susceptibleList(1)   
finalSusceptible = susceptible                      
initialImmune = vaccinated 
finalImmune = recovered + vaccinated      
initialInfected = infectedList(1)
finalInfected = infected

finalDead = dead
finalSurvivor = susceptible+infected+recovered+vaccinated
finalRecovered = recovered

reproductionNumber = k1*totalPopulation/(k2+k3)  % Reproduction number of the disease
herdImmunityThresholdFraction=1-(1/reproductionNumber);         % Herd immunity treshold fraction

% Prints Herd immunity treshold for Case B2, vaccination with death.
if vac > 0 & k3 > 0

    treshH = herdImmunityThresholdFraction*totalPopulation
end

end