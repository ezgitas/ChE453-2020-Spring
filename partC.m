
function partC
    clc
    Days = 30;                     
    populationTotal = 101;
    
    % Infection and recovery rates
    k1 = 0.01; 
    k2 = 0.2;

    
    %%%%%%%%%%%%%%% Stochastic Model
    
    % creates a population matrix where everyone is susceptible (t<0)
    population = zeros(Days+1,populationTotal); 
    
    % randomly determines who is infected at t=0
    k = randi([1 populationTotal]);                  
  
    
    % In the population matrix, subpopulations defined with numbers:
    % 0: susceptible
    % 1: infected
    % 2: recovered
    
    % Determines the first infected person randomly
    population(1,k) = 1;                            

    for i=2:Days+1
        for j = 1:populationTotal
            for ip = 1:length(find(population(i-1,:)==1)) 
            % counts how many people were infected the day before and 
            % checks contact possibility of susceptible population
                c = randi([0 1]);                              
                 % contact probability: 50%
                 % c = 1:    
                 % population(i,j),the susceptible individual contacted 
                 % with an infected person
                    if c == 1
                        if population(i-1,j) == 0
                            s = randi(100,1);         % infection rate: 1%
                            if s == 1
                                population(i,j) = 1;
                            end
                        end
                    end
                end
                
                if population(i-1,j) == 1
                    % recovery rate: 20%
                    r = randi(5,1);
                    if r==5
                        population(i,j) = 2;
                    else
                        population(i,j) = 1;
                    end
                
                end
                
                if population(i-1,j) == 2
                    population(i,j) = 2;
                    % recovered people gain immunity and remain recovered
                end
        end
    end
    
    % Empty matrices created for daily subpopulation count  
    S = zeros(1,Days+1);    % Susceptible
    I = zeros(1,Days+1);    % Infected
    R = zeros(1,Days+1);    % Recovered
    
    % Daily subpopulation count
    for i = 1:Days+1
        
        S(1,i) = length(find(population(i,:)==0));  % Susceptible count
        I(1,i) = length(find(population(i,:)==1));  % Infected count
        R(1,i) = length(find(population(i,:)==2));  % Recovered count
        
    end
    
%%%%%%%%%%%%%%% Deterministic Model

% Deterministic Model Solution
[tsoln,pop] = ode45(@fun,[0 Days],[populationTotal-1 1 0]); 
% pop = [S I R]

time = [0:Days];
% Deterministic Model Plot
plot(tsoln,pop(:,1),'b:',tsoln,pop(:,2),'g:',tsoln,pop(:,3),'r:','LineWidth',2)          
hold on
% Stochastic Model Plot
plot(time,S,'b',time,I,'g',time,R,'r','LineWidth',1.5)                                      
legend({'Susceptible Individuals','Infected Individuals','Recovered Individuals','Susceptible Individuals','Infected Individuals','Recovered Individuals'},'FontSize', 10,'Orientation','vertical','NumColumns',2)
xlabel('Time (days)')
ylabel('Number of People')

    function dpop = fun(t,pop)   % Change in population with time function
        
        dpop = zeros(3,1);

        dpop(1) = - k1*pop(2)*pop(1);               % dSdt
        dpop(2) = k1*pop(2)*pop(1) - k2*pop(2);     % dIdt
        dpop(3) = k2*pop(2);                        % dRdt

    end
    
end