
  ODE<- function(t,y,parms){
    with(as.list(y),{
      if((Tumor<=10^(-9))|(E_a <=10^(-9)) ){
        ScriptD=0}
      else{
        ScriptD = d*(E_a/Tumor)^(ell)/(s+(E_a/Tumor)^(ell))
      }
      
      
      dS_cdt = lambda*(a_1-a_3)*S_c-gamma_sc*E_a*S_c-delta_1*S_c#-(lambda*(S_c)^2)/k_1
      
      # I wrote the following in my notes: (Wilson-Levy they also include the last mixed term??)
       dE_adt = 5*E_a+f*E_a*Tumor/(1+c3*Tumor*B)-r*E_a-delta_3*E_a-delta0*R*E_a
      
      #dE_adt = mu*Tumor*E_a/(alpha+Tumor)-delta_3*E_a-(xi_ea*Tumor+gamma_ea*S_c)*E_a-delta0*R*E_a
      
      #The Ea equation below has the recruitment term in it that involves T and B. 
      #dE_adt = f*E_a*Tumor/(1+c3*Tumor*B)-(xi_ea*Tumor+gamma_ea*S_c)*E_a#-delta4*E_a #a_ea*E_a*(1-E_a/k_3)-delta_3*E_a-(xi_ea*Tumor+gamma_ea*S_c)*E_a-r*E_a-delta0*R*E_a
      #dE_adt = a_ea*E_a*(1-E_a/k_3)-r*E_a-delta0*R*E_a-delta_3*E_a-(xi_ea*Tumor+gamma_ea*S_c)*E_a
      
      #The Ea equation below has the logistic growth
      #dE_adt = a_ea*E_a*(1-E_a/k_3)-delta_3*E_a-(xi_ea*Tumor+gamma_ea*S_c)*E_a-r*E_a-delta0*R*E_a
      dTumordt = eta*Tumor*(1-Tumor/k_2)+lambda*(a_2+2*a_3)*S_c-delta_2*Tumor-ScriptD*Tumor-delta0*E_a*Tumor/(1+c1*B)+2*lambda*(S_c)^2/k_1
      dBdt = ab*Tumor^2/(c2+Tumor^2)-d*B
      dRdt = r*E_a-delta4*R
      
      
   
      list(c(dS_cdt, dE_adt, dTumordt, dBdt, dRdt))
    })
  }
  

  