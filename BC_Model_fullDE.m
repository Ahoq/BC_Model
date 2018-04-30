function dY= BC_Model_fullDE(t,Y,par,vb,kguess,kindex)

par(kindex)=kguess;

a_1=par(1); %1
a_3=par(2); %2
a_2=par(3); %3
k_1=par(4); %4
k_2=par(5); %5
k_3=par(6); %6
delta_1=par(7); %7
delta_2=par(8); %8
delta_3=par(9); %9
lambda=par(10); %10
eta=par(11); %11
a_ea=par(12); %12
gamma_sc=par(13); %13
gamma_ea=par(14); %14
xi_T=par(15); %15
xi_ea=par(16); %16
s=par(17); %17
d=par(18); %18
ell=par(19); %19
mu = par(20);
alpha=par(21);
ab=par(22);
c1=par(23);
c2=par(24);
c3=par(25);
d=par(26);
delta0=par(27);
delta4=par(28);
f=par(29);
r=par(30);

tau=par(31); %20
Delta=par(32); %21
epsilon=par(33);%22
sigma=par(34);%23
%V=par(24);%24
%k=par(25);%25






%%%%%%%%%%nondimensional parameters

p(1)=lambda*(a_1-a_3);
p(2)=delta_1;
p(3)=lambda;
p(4)=lambda*(a_2+2*(a_3));
p(5)=delta_2;
p(6)=gamma_sc;
p(7)=(2*lambda)/k_1;
p(8)=a_ea;
p(9)=1/k_3;
p(10)=delta_3;
p(11)=xi_ea;
p(12)=gamma_ea;
p(13)=eta;
p(14)=ell;
p(15)=s;
p(16)=d;
p(17)=1/k_2;
p(18)=k_1;
p(19)=f;
p(20)=c3;
p(21)=r;
p(22)=delta0;
p(23)=c1;
p(24)=ab;
p(25)=c2;
p(26)=d;
p(27)=delta4;

dY=zeros(5,1); 

D_Sc=Y(1);
Ea_Tumor=Y(2);
D_Tumor=Y(3);
B = Y(4);
R = Y(5);



%vaccine treatment: vb is set in the BC_Model_Main_Tumor 
vt=0;%we haven't looked at putting vaccine straight into the tumor yet


if ((D_Tumor<=10^(-9))||(Ea_Tumor<=10^(-9)))
    ScriptD=0;
else
    ScriptD=p(16)*((Ea_Tumor/D_Tumor)^(p(14))/(p(15)+(Ea_Tumor/D_Tumor)^(p(14))));
end



dY(1)= p(1)*D_Sc-(p(6)*Ea_Tumor*D_Sc)-(p(2)*D_Sc);%-(p(3)*(D_Sc)^2)/p(18);
%dY(2) = p(13)*D_Tumor*(1-(D_Tumor*p(17)))+(p(4)*D_Sc)-(p(5)*D_Tumor)-(ScriptD*D_Tumor)+p(7)*((D_Sc)^2);
dY(2) = 5 * Ea_Tumor + p(19)* Ea_Tumor/(1+p(20)*D_Tumor*B)-p(21)*Ea_Tumor-p(10)*Ea_Tumor-p(22)*R*Ea_Tumor;
%dY(3)= p(8)*Ea_Tumor*(1-(Ea_Tumor*p(9)))-p(10)*Ea_Tumor-(p(11)*D_Tumor+p(12)*D_Sc)*Ea_Tumor;
    
dY(3)= p(13)*D_Tumor*(1-D_Tumor*p(17))+ p(4)*D_Sc-p(5)*D_Tumor-ScriptD*D_Tumor-p(22)*Ea_Tumor*D_Tumor/(1+p(23)*B)+(2*p(3)*((D_Sc)^2))/p(18);
dY(4)= (p(24)*((D_Tumor)^2))/(p(25)+(D_Tumor)^2)-p(26)*B;
dY(5)= p(21)*Ea_Tumor-p(27)*R;
