close all; 
clear all;

timestart= 0;
timeend=10;

a_1 = 0.35; 
a_3 = 0.25; 
a_2 = 1-a_1-a_3; 


k_1 = 6*10^4; 
k_2 = 10^9; 
k_3 = 3*10^8;

delta_1 = 10^-3;
delta_2 = 0.05; 
delta_3 = 0.462; 

lambda = 0.07; 
eta = 0.3954; 

a_ea = 2.4; 
gamma_sc = 10^(-4); 
gamma_ea = .1; 

xi_T = 0.35; 
xi_ea = 9.42*10^(-12); 

s = 1.4; 
d = 0.35; 
ell = 2/3; 
V= 1; 
k=10^9;

tau = 1; %Non-dimensionalization
Delta = 1;  %Non-dimensionalization
epsilon = 1; %Non-dimensionalization
sigma = 1;

par=[
    a_1; %1
a_3; %2
a_2; %3
k_1; %4
k_2; %5
k_3; %6
delta_1; %7
delta_2; %8
delta_3; %9
lambda; %10
eta; %11
a_ea; %12
gamma_sc; %13
gamma_ea; %14
xi_T; %15
xi_ea; %16
s; %17
d; %18
ell; %19
tau; %20
Delta; %21
epsilon;%22
sigma;%23
V;%24
k;%25
];
Data(:,1)=[3; 6; 9; 12; 15; 18; 21; 27; 30]; %time in days - from Wilson/Levy paper
Data(:,2)=[2.82; 9.83; 21.7; 25.08; 39.29; 60.12; 99.53; 169.38; 249.01]; %size of tumor (mm^2)



%non-dimensionalize data
Cellsize=0.01^2;%using average size of cell is 0.01 mm - this needs to be confirmed/investigated
ndData(:,1)=Data(:,1)/tau;
ndData(:,2)=(Data(:,2)/Cellsize)/sigma; 





%Initial conditions

InitCond(1)=0;%D_Sc=Y(1);
InitCond(2)=0;%D_Tumor=Y(2)
InitCond(3)=ndData(1,2);%Ea_Tumor=Y(3)



%the way MATLAB ode solvers work we can't really specify the times the
%solution is found. We will solve the ODEs for each segment of treatment
%time (using the previous segment of time as the initial condition)

TreatTimes=[0,1/48,14,14+1/48,21,21+1/48,47]/tau; %nondimensional (due to the division by tau)
%       treatment times, these are the times the treatment either turns on or off 
TreatValues=[2*10^5,0,2*10^5,0,2*10^5,0,0]/V; %amount of treatment given at specified times
% the first value is the amount given in between the first two treatment
% times, the second valued is the amount given between the second and third
% treatment times, the last value will not be used. these are
% nondimensional values because we are dividing by V

LL=1; %this is a "dummy" variable to allow the following to be concatinated

vb=0;

%parameter fitting
tstart=ndData(1,1);
tend=ndData(length(ndData),1);
%choosing which parameters to fit
kindex=[12]; %range: 1-25
kguessrange=[par(kindex)-2, par(kindex)+2];

%chooses m guesses b in range stated above for the parameter
m=10; %number of intial guesses want to try in fminsearch
 randkguess = kguessrange(1) + (kguessrange(2)-kguessrange(1)).*rand(m,1);

for i=1:length(randkguess)
    while randkguess(i)<0
        randkguess(i)=kguessrange(1) + (kguessrange(2)-kguessrange(1)).*rand;
    end
kguess=randkguess(i);


param = fminsearch(@(kguess) BC_fitdata(ndData, tstart, tend,InitCond, par,vb,kguess,kindex),kguess);



Error(i)=BC_fitdata(ndData, tstart, tend,InitCond, par,vb,param,kindex);

paramfit(i)=param;
end

index=find(Error==min(Error));
useparam=paramfit(index);

%we now have fit our parameter so we want the actual solution using that
%parameter
%[T,Xsoln]=ode45(@(t,x) SysDE(t,x,param),[tstart,tend], X0);
%%%%[Time,Y]=ode15s(@(y) fullDE(y,p,vb,param,kindex), [0,47/tau], InitCond);

%when incorporating treatment 
% for i=1:length(TreatTimes)-1
%    vb=TreatValues(i); %the treatment during this segment
% [Time,Y]=ode15s(@(t,y) fullDE(y,p,vb), [TreatTimes(i),TreatTimes(i+1)], InitCond); %solving the ODE for this segment
% InitCond=Y(length(Y),:); %setting the initial condition (the end result from this time segment) for the next segment
% FullY(LL:LL+length(Y)-1,:)=Y; %concatinating the solutions
% FullT(LL:LL+length(Time)-1)=Time; %concatinating the time
% LL=LL+length(Y); %updating where the next solution will go in the Full matrices
% end

    
 [Time,Y]=ode15s(@(t,y) BC_Model_fullDE(t,y,par,vb,useparam,kindex), [tstart,tend], InitCond);
 
 
 FullY=Y;
 FullT=Time;


figure
plot(FullT*tau,Delta*FullY(:,1))
xlabel('Time')
ylabel('D_Sc')
figure
plot(FullT*tau,sigma*FullY(:,2))
hold on
plot(Data(:,1),Data(:,2)/Cellsize,'*')
xlabel('Time')
ylabel('D_Tumor')
figure
plot(FullT*tau,epsilon*FullY(:,3))
xlabel('Time')
ylabel('E_a Tumor')