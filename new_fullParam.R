#Parameters for the full system
#From Parameter Results in Dropbox->ReufCancerVaccines->MATLABCode->EmekParameters_LeeChoLee
alpha = 1 #cell
a_EaT=0.46 #1/day
a_D=0.23  # 1/day
a_Em=0.01  #1/day
a_EaS=0.12  #1/day
c=9.4*10^(-12)  #1/(cell.day)
d=0.35 #1/day
k=1.0*10^8 #cells
ell=0.6  #unitless
mu_B=27 #1/day
mu_BB=6.3 #1/day
mu_BS=2.8 #1/day
mu_BSE=0.0014 #1/day
#mu_BSE = 7.33*10^(-4)
#mu_SB= 0.012 #1/day
mu_TB=0.13 #1/day
mustar_SB=11 #1/day reduced spleen to blood
mu_normal=0.040 #1/day normal spleen to blood
m=1000 #cells/day
#m=2.4*10^4
q=200 #cells 
r=0.32 #1/day  
#r=0.3954
r_am=0.01  #1/day
s=2.3  #unitless  
theta_shut=0.12 # cells

#values from the non-dimensionalizations
Delta = 10^4 #cells
epsilon = 10^6 #cells
sigma = k #cells
V = (1/48) #cells/day - eqn 4 in DePillis
tau = 1/mu_BSE #in Matlab Code

#these are new parameters
gamma=14  # 1/day 
lambda=4/3 #1/day 
lambda_m=0.15 #1/day,
beta_d=(0.17+0.4)/2 #1/day (death of proliferating cells - from Rad 2012 (chapter), average of delta_A)
beta_E=13.37 #1/day #beta=beta_E+beta_d
beta=beta_E+beta_d #1/day 
omega=13 #unitless: cell/cell 


#these are variables in entire system, values taken from equilibrium values for REU writeup 
# D_blood=4
# Ea_blood=4*10^(5)
# Em_blood=14000

#vaccine schedule
dur1 = V/tau #duration of first dose in days
dur2 = V/tau #duration of second dose in days
dur3 = V/tau # duration of third dose in days


timestart = 0/tau #time of start
timedose1 = 0/tau #time from start of first dose
timedose2 = 2/tau #time from start of second dose
timedose3 = 7/tau #time from start of third dose
timedose4 = 11/tau #time from start of third dose
timeend = 12/tau #Last day

DD = 5*10^5 # dose
dose1 = DD/(V*Delta)*tau # amount of dose 1 
dose2 = DD/(V*Delta)*tau # amount of dose 2 
dose3 = DD/(V*Delta)*tau # amount of dose 3 
dose4 = DD/(V*Delta)*tau # amount of dose 4 

#Create a (2n+1) by 2 matrix where n = number of doses.  Currently n=3
vaccine = t(matrix(c(timestart, timestart+dur1, timedose2, timedose2+dur2, timedose3,
                     timedose3+dur3,timedose4, timedose4+dur3, timeend,
                     dose1, 0, dose2, 0, dose3, 0, dose4,0,0), 
                   nrow=9, ncol=2))
