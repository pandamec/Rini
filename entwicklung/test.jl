using SymPy

x=symbols("x")  # Define symbolic variable

# Define the function to integrate
f = x^2 + 3x + 2

# Perform symbolic integration
integral_f = integrate(f, x)

println("Integral of f: ", integral_f)

using SymPy
function StaticBeam2(E,Ge,deltav)


    Ek=E[4,:]
    Es=E[5,:]

    tk,lk,wk = Ge[4,:]
    ts,ls,ws = Ge[5,:]

    # Fatigue-Test
    ####### Parylene ###########

    ## Geometrie 

    dk=tk
    wp=wk #mm

    ## Messaufbau Konfiguration

    ds=ts #mm
    w=ws #mm
    I=w*ds^3/12
    l=65

    ## Krafte und Torque
    F=(deltav*48*Es*I)/l^3
    M_max=(F/2)*l/2
    x_p=19
    M=(2*x_p/l)*M_max #N.mm
     print(F)
    ## Neutral axis berechnung

    c = symbols("c")
    eq=(w*Es*((ds-c)^2-(c)^2)+wp*Ek*((ds-c+dk)^2-(ds-c)^2))
    sol=solve(eq,c)
    c_val=sol[c]
    print(c_val)
    ## Belastung

    

    return c_val

end

E= [150000; 2800; 2800; 200000; 200000] #First layer Silicon
Ge=[0.32 20 10 ; 0.01 20 10; 0.01 20 10; 0.5 20 15; 0.5 65 15 ] 

v=0.5

c = StaticBeam2(E,Ge,v)