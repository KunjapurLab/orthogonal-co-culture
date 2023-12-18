clear
clc
% Defining t range and intervals 
t0 = 0;
tf = 48; % hr
trange =linspace(t0, tf, 1000);

% Defining model parameters as described in Table S6
umu = 0.88;% hr-1
umDEP = 0.8;% hr-1
Kso = 19.5; %mg/L 
Ksf = 0.5; %mg/L 
KsBipA = 0.3757*241.3/1000; %mg/L
kdp = 0.03; %hr
kdu = 0.05; %hr
YsoDEP = 0.0725; %mg/mg-h 
Yxsf = 0.4; % mg/mg

%Define initial conditions as described in Table S7
So0 = 0;%mg/L
Sf0 = 1000; %mg/L
Sb0 = 5*241.2/1000; %mg/L
Xp0 =2.5; %mg/L
Ratio = 1;
Xu0 = Xp0*Ratio; %mg/L


%Store the initial conditions in a matrix for the solver
variable0 = [Xp0; Xu0; Sf0; So0];  

%Solve ODE
[t,y]=ode45(@(t,variable)solveODEs(t,variable, Ksf, umDEP, umu, Yxsf, kdp, kdu, Kso, YsoDEP, KsBipA,Sb0),trange,variable0); 
Xp = y(:,1);
Xu = y(:,2);
Sf = y(:,3);
So = y(:,4);

%Export matrix to file 
M = [t,y];
writematrix(M, 'M_dualsynaux_coculture.xls')

%Figure output
figure(1)
clf
hold on
plot(trange,Xp,'k-')
plot(trange,Xu,'r--')
plot(trange,Sf,'b:')
plot(trange,So,'b-')
hold off
xlim([0,tf])
legend('DEPe5-p','N-adk.d6-u','Sf','So')
xlabel('t (hr)')
ylabel('mg/L')
    
%Save plot graphic
saveas(gcf, 'dualsynaux_coculture_plot', 'pdf')
    
%end

function dsolveODEs=solveODEs(t,variable, Ksf, ump, umu, Yxsf, kdp, kdu, Kso, Ysox, KsBipA, Sb) 
     %Assign ys to the correct species
            Xp = variable(1); % Xp Cell Population
            Xu = variable(2); % Xu Cell Population
            Sf = variable(3); % feed concentration
            So = variable(4); % OMeTyr concentration
     
     %The ODEs
    
            DSp= So/Kso;
            DSu= Sf/Ksf;
            DSb = Sb/KsBipA;

            if Sf < 0.0000001 
                Sf = 0;
            end

            
            %ODE for growth of the DEP.e5 producer strain - Eqs 7-8 in Supplmentary Methods:
            if real(DSb) < real(DSu)
                 ugp = ump*Sb / (KsBipA + Sb);
                 unp= ugp - kdp;
                 dXpdt = Xp*unp;
            else
                 ugp = ump*Sf / (Ksf + Sf);
                 unp= ugp - kdp;
                 dXpdt = Xp*unp;
            end

            %ODE for production of OMeTyr - Eq 9 in Supplementary Methods
            if real(DSb) < real(DSu)
                dSodt = Xp*Ysox*(ump*Sb / (KsBipA + Sb));
            else
                dSodt = Xp*Ysox*(ump*Sf / (Ksf + Sf));
            end

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