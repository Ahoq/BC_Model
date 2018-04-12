function E= BC_fitdata(ndData, tstart, tend,Y0,par,vb,kguess,kindex)

%par(kindex)=k;

[T,Ysoln]=ode15s(@(t,Y) BC_Model_fullDE(t,Y,par,vb,kguess,kindex),[tstart,tend], Y0);
%[T,Xsoln]=ode45(@(t,C) SysDE(t,C,k),[tstart,tend], X0);


%we can't specify the values for t in the ode solver but we only have
%specific values for t in the actual data. We want to be able to compare
%data to a value in the simulation so we find the time point in the
%simulation CLOSEST to that in the data
for j = 1:length(ndData(:,1))
    timepoint = ndData(j,1);
    I=find(T>= timepoint);
    SolnVal(j,1) = Ysoln(I(1),3);
end;
E = norm((SolnVal - ndData(:,2)), 2); %finding the total error





