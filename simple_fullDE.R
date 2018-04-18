#Non-dimensionalized Full System with our new spleen equations
simple_fullDE<- function(t,y,vb){
  with(as.list(y),{
    if((Tumor<=10^(-9))|(Ea_spleen <=10^(-9)) ){
        ScriptD=0}
    else{
        ScriptD = (Ea_tumor/Tumor)^(p31)/(p28+(Ea_tumor/Tumor)^(p31))
    }
    
    vt=0 #currently we have v_blood = 0 if we want to change that, we should do the same as vt
   
    
    dD_blooddt = -p1*D_blood+p2*D_tumor+p3*vb
    dEa_blooddt = (p4+p24*(1/(1+D_spleen/p27)))*Ea_spleen-p5*Ea_blood
    dEm_blooddt = (p4+p24*1/(1+D_spleen/p27))*Em_spleen-p5*Em_blood
    
    dD_spleendt = p6*D_blood-p7*D_spleen+(p10)*W
    dEa_spleendt = p8*Ea_blood-(p24*1/(1+D_spleen/p27))-p14*Ea_spleen+p15*P
    dEm_spleendt = p16*Ea_spleen+p8*Em_blood-(p17+p4+p24*1/(1+D_spleen/p27))*Em_spleen-p26*Em_spleen*D_spleen
    dPdt = p11*W+p12*Em_spleen*D_spleen-p13*P
    dWdt = p9*D_spleen-p10*W
    
    dEa_tumordt = p5*1/(p29+Tumor)*Tumor*Ea_blood-p18*Ea_tumor-p19*Ea_tumor*Tumor
    dD_tumordt = p25*Tumor/(p30+Tumor)-p23*D_tumor+p3*vt
    dTumordt = p20*Tumor-p21*(Tumor)^2-p22*ScriptD*Tumor
    
    list(c(dD_blooddt, dEa_blooddt, dEm_blooddt, dD_spleendt, dWdt, dPdt, dEa_spleendt, dEm_spleendt, dEa_tumordt, dTumordt, dD_tumordt))
  })
}