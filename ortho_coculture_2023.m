clear
clc
% Defining t range and intervals 
t0 = 0;
tf = 24; % hr
trange =linspace(t0, tf, 1000);

% Define model parameters as described in Table S6
umu = 0.88;% hr-1 
ump = 0.95;% hr-1 
Kso = 19.5; %mg/L 
Ksf = 0.5; %mg/L 
kdu = 0.2; %hr
kdp = 0.15; %hr
Yxsf = 0.4; % mg/mg 
Ysox = 0.125; %mg/mg
Ratio = 1;


%Set the initial conditions for the model as defined in Table S7
So0 = 0;%mg/L
Sf0 = 1000; %mg/L
Xp0 =2.5; %mg/L
Xu0 = Xp0*Ratio; %mg/L

%Store the initial conditions in a matrix for the solver
variable0 = [Xp0; Xu0; Sf0; So0]; 

%Solve ODE
[t,y]=ode45(@(t,variable)solveODEs(t,variable, Ksf, ump, umu, Yxsf, kdp, kdu, Kso, Ysox),trange,variable0); 
Xp = y(:,1);
Xu = y(:,2);
Sf = y(:,3);
So = y(:,4);      

%Export matrix to file 
M = [t,y];
writematrix(M, 'M_ortho_coculture.xls')

%plot of model output for model verification
figure(1)
clf
hold on 
plot(trange,Xp,'r-')
plot(trange,Xu,'g--')
plot(trange,Sf,'b:')
plot(trange,So,'b-')
hold off
xlim([0,tf])
legend('Xp','Xu','Sf','So')
xlabel('t (hr)')
ylabel('mg/L')
hold off

%Save plot graphic
saveas(gcf, 'ortho_coculture_plot', 'pdf')

function dsolveODEs=solveODEs(t,variable, Ksf, ump, umu, Yxsf, kdp, kdu, Kso, Ysox) 
     %Assign ys to the correct species
            Xp = variable(1); % Xp Cell Population
            Xu = variable(2); % Xu Cell Population
            Sf = variable(3); % Limiting nutrient in LB concentration
            So = variable(4); % OMeTyr concentration
     
     %The ODEs
    
            DSp= So/Kso;
            DSu= Sf/Ksf;
            if Sf < 0.0000001 
                Sf = 0;
            end

            %ODE for growth of the producer strain - Eqs 1-2 in Supplmentary Methods: 
            if real(Sf) > 0
                 ugp = ump*Sf / (Ksf + Sf);
                 unp= ugp - kdp;
                 dXpdt = Xp*unp;
            else
                 dXpdt = Xp*((ump*Sf / (Ksf))-kdp);
            end

            %ODE for production of OMeTyr - Eq 3 in Supplementary Methods
             dSodt = Xp*Ysox*(ump*Sf /(Ksf+Sf));

            %ODE for growth of the utilizer strain - Eq 4-5 in
            %Supplementary Methods
             if real(DSp) < real(DSu)
                 ugu = umu*So / (Kso + So);
                 unu= ugu - kdu;
                 dXudt = Xu*unu;
             else
                 ugu = umu*Sf / (Ksf + Sf);
                 unu= ugu - kdp;
                 dXudt = Xu*unu;
             end

             %ODE for consumption of the limiting nutrient for growth in LB
             %- Eq. 6 in the Supplementary Methods
             if real(Sf) > 0                 
                 dSfdt = -(dXpdt + dXudt)/Yxsf;
             else
                 dSfdt = 0;
             end
             
     %Set the solution matrix
             dsolveODEs = [dXpdt, dXudt, dSfdt, dSodt]';

 end