
####### Thermalverhalten ##########


function MetalSubstrat(E,Ge,deltaT)
    ###### Cylindrical via ############
 

    ## Materials

    #Substrate
    E_s, v_s, alpha_s =E[1,:]
    
    #Metall
    E_m, v_m, alpha_m =E[2,:]
   
    
    ## Equations

    ## Verformung 
    u_m,u_s, u_s2 =symbols("u_m u_s u_s2")

    ## Spannung
    p=symbols("p")


    ## Geometrie

     r_s, r_m, l  = Ge 

    ## Load


        ## Kompatibilität
        E_m=Rini.E_Aluminium(deltaT)
        eq1 = u_s - ((1-v_s)/E_s)*((-p*r_m^2)/(r_s^2-r_m^2))*r_m+(((1+v_s)/E_s)*(r_m^2*r_s^2)/r_m)*(p/(r_s^2-r_m^2)) - r_m*alpha_s*deltaT
        eq2 = u_m - ((1-v_m)/E_m)*((+p*r_m^2)/(r_m^2))*r_m - r_m*alpha_m*deltaT
        eq3 = u_s  - u_m
        eq4 = u_s2 - ((1-v_s)/E_s)*(-p*r_m^2)/(r_s^2-r_m^2)*r_s+(((1+v_s)/E_s)*(r_m^2*r_s^2)/r_s)*(p/(r_s^2-r_m^2)) - r_s*alpha_s*deltaT

        ## solve
        sol=solve([eq1,eq2,eq3,eq4],
                [u_s,u_s2,u_m,p])
        


        return [sol[p], sol[u_m], sol[u_s2]]
    
end


####### Fatigue ##########

function StaticComposite(E,Ge,deltav)

    E3=E[1,:] #Silicon
    E2=E[2,:]
    E1=E[3,:]
    Es=E[4,:]

    t3,l3,w3 = Ge[1,:]
    t2,l2,w2 = Ge[2,:]
    t1,l1,w1 = Ge[3,:]
    ts,ls,ws = Ge[4,:]

    # Fatigue-Test
    ####### Parylene ###########

    ## Geometrie 
    d1=t1 #mm
    d2=t2 #mm
    d3=t3 #mm
    wp=w1 #mm

    ## Messaufbau Konfiguration

    ds=ts #mm
    w=ws #mm
    I=w*ds^3/12
    l=70

    ## Krafte und Torque
    F=deltav*48*Es*I/l^3
    M_max=(F/2)*l/2
    x_p=l/2
    M=(x_p/l)*M_max #N.mm

    ## Neutral axis berechnung

    c = symbols("c")
    eq=(1/2)*(w*Es*((ds-c)^2-(c)^2)+wp*E1*((ds-c+d1)^2-(ds-c)^2)+wp*E2*((ds-c+d1+d2)^2-(ds-c+d1)^2)+wp*E3*((ds-c+d1+d2+d3)^2-(ds-c+d1+d2)^2))
    sol=solve(eq,c)
    c_val=sol[c]
    ## Belastung

    m=symbols("m")
    eqb=M-m*(1/3)*(w*Es*((ds-c_val)^3-(c_val)^3)+wp*E1*((ds-c_val+d1)^3-(ds-c_val)^3)+wp*E2*((ds-c_val+d1+d2)^3-(ds-c_val+d1)^3)+wp*E3*((ds-c_val+d1+d2+d3)^3-(ds-c_val+d1+d2)^3)+w*Es*((c_val)^3))
    solb=solve(eqb,m)
    m_val=solb[m]

    sigma3=m_val*(ds-c_val+d1+d2+d3)*E3[]
    sigma2=m_val*(ds-c_val+d1+d2)*E2[]
    sigma1=m_val*(ds-c_val+d1)*E1[]
    sigmas=m_val*(c_val)*Es[]

    sigma=[sigma3, sigma2, sigma1, sigmas]
   

    return sigma

end