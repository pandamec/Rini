
####### Thermalverhalten ##########


function StaticBeam(E,Ge,deltav)

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