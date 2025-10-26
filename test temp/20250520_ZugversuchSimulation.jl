
using Pkg
Pkg.add("Plots")

using Plots


function TensileStress(E,delta)
    Sigma=E*delta
    return Sigma    

end


E_min=2.2*1e9
E_max=3.2*1e9
w=15*1e-3
t=15*1e-6
l=30*1e-3

As=w*t
delta=range(0,0.025,20)

Force_min=[]
Force_max=[]
for i in delta

    sigma=TensileStress(E_min,i)
    f=sigma*As
    push!(Force_min,f)
    sigma=TensileStress(E_max,i)
    f=sigma*As
    push!(Force_max,f)
end

scatter((delta.+1)*l*1000, Force_min,yticks=0.0:1:20, label="E_min 2.2GPa", xlabel="Zwischen Klemmung Lange (mm)", ylabel="Kraft(N)", title="Zugversuch Simulation")
scatter!((delta.+1)*l*1000, Force_max, label="E_max 3.2 GPa")
hline!([40*1e6*As], label="SigmaMax = 40MPa", color=:red, linestyle=:dash)