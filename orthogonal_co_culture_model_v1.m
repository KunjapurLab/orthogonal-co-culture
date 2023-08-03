clear
clc
% Defining t range and intervals 
t0 = 0;
tf = 24; % hr
trange =linspace(t0, tf, 1000);

% We also need to set the initial conditions
ump = 1.0;% hr-1 - best guess can keep fixed
umu = 0.9;% hr-1 - need to fix 
Ksf = 0.5; %mg/L
Yxs = 0.4; % mg/mg
kdp = 0.035; %hr
kdu = 0.05; %hr
Kso = 0.0998*195.22; %mg/L
qp = 0.05; %mg/mgh


So0 = 0;
Sf0 = 1000; 
Xp0 = 2.5; %mg/L
Ratio = 20;
Xu0 = Xp0*Ratio; %mg/L

%Store the initial conditions in a matrix for the solver
%for j = 1:3
    variable0 = [Xp0; Xu0; Sf0; So0];  
        [t,y]=ode45(@(t,variable)solveODEs(t,variable, Ksf, ump, umu, Yxs, kdp, kdu, Kso, qp),trange,variable0); 
        Xp = y(:,1);
        Xu = y(:,2);
        Sf = y(:,3);
        So = y(:,4);
    
        ratio = zeros(1000,1);
        for i = 1:length(ratio)
            ratio(i)=Xu(i)/Xp(i);
        end
    
    %output file information for plotting outside of MATLAB
    output_filename ='1-1-ratio-test1.csv';
    fid = fopen(output_filename, 'w');
    fprintf(fid, 'time,Mp,Mu,Sf,So,ratio\n');
    A = [trange' Xp Xu Sf So ratio];
    fprintf(fid, '%f\n',A);
    fclose(fid);
    
    %Figure output
    figure(1)
    clf
    hold on
    plot(trange,Xp,'r-')
    plot(trange,Xu,'g--')
    plot(trange,Sf,'b:')
    plot(trange,So,'b-')
    hold off
    xlim([0,tf])
    legend('Mp','Mu','Sf','So')
    xlabel('t (hr)')
    ylabel('mg/L')
    
    figure(2)
    clf
    hold on
    plot(trange, ratio, 'k-')
    hold off
    xlim([0,tf])
    legend('ratio Mu/Mp')
    xlabel('t(hr)')
    ylabel('Mu/Mp')
    
    disp(Xu(1000))
    
%end
    
    function dsolveODEs=solveODEs(t,variable, Ksf, ump, umu, Yxs, kdp, kdu, Kso, qp) 
     %Assign ys to the correct species
            Xp = variable(1); % Mp Cell Population
            Xu = variable(2); % Mu Cell Population
            Sf = variable(3); % feed concentration
            So = variable(4); % OMeTyr concentration
     
     %ODES  
            DSp= So/Kso;
            DSu= Sf/Ksf;
    %For producer biomass concentration:
            if Sf < .1 
                Sf = 0;
            end

            if real(Sf) > 0
                 ugp = ump*Sf / (Ksf + Sf);
                 unp= ugp - kdp;
                 dXpdt = Xp*unp;
            else
                 dXpdt = Xp*(-kdp);
            end

     %For utilizer biomass concentration:
             if real(DSp) < real(DSu)
                 ugu = umu*So / (Kso + So);
                 unu= ugu - kdu;
                 dXudt = Xu*unu;
             else
                 ugu = umu*Sf / (Ksf + Sf);
                 unu= ugu - kdp;
                 dXudt = Xu*unu;
             end

     %For feed substrate concentration:       
             if real(Sf) > 0                 
                 dSfdt = -(dXpdt + dXudt)/Yxs;
             else
                 dSfdt = 0;
             end

     %For OMeTyr concentration: 
             if real(Sf)> 0
                dSodt = Xp*qp;
             else
                 dSodt = qp;
             end

     %Set the solution matrix
             dsolveODEs = [dXpdt, dXudt, dSfdt, dSodt]';

 end