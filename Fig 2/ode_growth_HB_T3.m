function dy =  ode_growth_HB_T3(t,y,par,flist)
    
       
    kmet = par(2);
    uR_aa =  par(3); 
    NR = par(4);
    NP = par(5);
    Mf = flist{1};
    Vf = flist{2};
    ktrans = flist{3};
    T3convert = flist{4};
    kmeteff = flist{5};
      
    
    dy = NaN*ones(3,1);
    %1 - 8: R - P - aa 
    if T3convert(y) > 0
    dy(1) = 1/NR*uR_aa*y(1)*ktrans(y);
    dy(2) = 1/NP*(1-uR_aa)*y(1)*ktrans(y);
    dy(3) = y(2)*kmeteff(y)-y(1)*ktrans(y); %-0.03*Mf(y);    
    end
    
end