

timestart = 0
timeend = 50

a_1 = 0.35 #Prob Symmetric Division of CSC 
a_2 = 0.65 # Prob of Asymmetric Division of a CSC
a_3 = 1-a_1-a_2 #Prob of Symmetric Differentiation (or commitment) of a CSC into 2 NSCC (non-stem cancerous cell)
#a_2 = 1-a_1-a_3 #Prob of Asymmetric Division of CSC


k_1 = 6*10^4 #Carrying Capacity of CSC  (cells)
k_2 = 10^(12) #Carrying Capacity of Tumor Cells (cells)- changed it! 
k_3 = 3*10^8 # Carrying Capacity of CD8a+ (cells)

delta_1 = 0#0.0005 #Natural Death Rate of CSC (1/time)
delta_2 = 0.87222#0.1 #Natural Death Rate of Tumor Cells (1/time)
delta_3 = 0.462 # Natural Death Rate of CD8a+ (1/time)

lambda = 1 #0.07 #Proliferation Rate of a CSC (1/days)
eta = 1 #0.3954 # Proliferation Rate ofor Tumor Cells (1/time)

a_ea = 2.4 #Recruitment Rate for CD8a+ (1/time)
gamma_sc = 10^(-9) # Interactions of CSCs and E_a leading to CSC death (1/(cell*time))
gamma_ea = 10^(-9)#3*10^(-3)#.1 #Interactions of CSCs and E_a leading to E_a death (1/(cell*time))

xi_T = 0.35 #Interaction of Tumor and E_a leading to Tumor cell death (1/(cell*time))
xi_ea = 9.42*10^(-12) # Interaction of Tumor and E_a leading to E_a death (1/(cell*time))

s = 1.4 #Value of (E_a/T)^l necessary for half-maximal activated CTL toxicity
d = 0.35 # Maximal fraction tumor kills by CTLs (1/day)
ell = 2/3 #Immune strength scaleing exponent


mu = 6.3 #from Immunokinetics paper
alpha = 1

# All of the paremeters below come from the Wilson-Levy paper
ab = 0.3
c1 = 100
c2 = 300
c3 = 300
d =7*10^(-4)
delta0 = 10^(-5) 
delta4 = 10^(-5) 
f =0.62
r = 0.01

  


tau = 1 #Non-dimensionalization
Delta = 1  #Non-dimensionalization
epsilon = 1 #Non-dimensionalization
