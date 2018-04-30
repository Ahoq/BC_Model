close all; 
clear all;

timestart= 0;
timeend=10;

a_1 = 0.35; 
a_2 = 0.65; 
a_3 = 1-a_1-a_2; 


k_1 = 6*10^4; 
k_2 = 10^12; 
k_3 = 3*10^8;

delta_1 = 0.0005;
delta_2 = 0.87222; 
delta_3 = 0.462; 

lambda = 1; 
eta = 1; 

a_ea = 2.4; 
gamma_sc = 10^(-9); 
gamma_ea = 10^(-9); 

xi_T = 0.35; 
xi_ea = 9.42*10^(-12); 

s = 1.4; 
d = 0.35; 
ell = 2/3; 
%V= 1; 
%k=10^9;
mu = 6.3;
alpha =1;
ab=0.3;
c1=100;
c2=300;
c3=300;
d=7*10^(-4);
delta0=10^(-5);
delta4=10^(-5);
f=0.62;
r=0.01;

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
mu; %20
alpha;%21
ab;%22
c1;%23
c2;%24
c3;%25
d;%26
delta0;%27
delta4;%28
f;%29
r;%30
tau; %31
Delta; %32
epsilon;%33
sigma;%34
];
Data(:,1)=[3; 6; 9; 12; 15; 18; 21; 27; 30]; %time in days - from Wilson/Levy paper
Data(:,2)=[2.82; 9.83; 21.7; 25.08; 39.29; 60.12; 99.53; 169.38; 249.01]; %size of tumor (mm^2)



%non-dimensionalize data
Cellsize=0.01^2;%using average size of cell is 0.01 mm - this needs to be confirmed/investigated
ndData(:,1)=Data(:,1)/tau;
ndData(:,2)=(Data(:,2)/Cellsize)/sigma; 





%Initial conditions

InitCond(1)=10^6*(1/79);%D_Sc=Y(1);
InitCond(2)=10^6*(1/79);%Ea_Tumor=Y(2)
InitCond(3)=10^6;%D_Tumor=Y(3)
InitCond(4)=10^4;%B
InitCond(5)=10^6*(20/79);%R



%the way MATLAB ode solvers work we can't really specify the times the
%solution is found. We will solve the ODEs for each segment of treatment
%time (using the previous segment of time as the initial condition)

%TreatTimes=[0,1/48,14,14+1/48,21,21+1/48,47]/tau; %nondimensional (due to the division by tau)
%       treatment times, these are the times the treatment either turns on or off 
%TreatValues=[2*10^5,0,2*10^5,0,2*10^5,0,0]/V; %amount of treatment given at specified times
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
kindex=17; %range: 1-34
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
plot(FullT,FullY(:,1))
xlabel('Time')
ylabel('D_Sc')

figure
plot(FullT,FullY(:,2))
xlabel('Time')
ylabel('Ea_Tumor')

figure
plot(FullT,FullY(:,3))
%hold on
%plot(Data(:,1),Data(:,2)/Cellsize,'*')
xlabel('Time')
ylabel('Tumor')

figure
plot(FullT,FullY(:,4))
xlabel('Time')
ylabel('B')

figure
plot(FullT,FullY(:,5))
xlabel('Time')
ylabel('R')