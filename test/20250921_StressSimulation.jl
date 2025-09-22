
using CairoMakie

w=2
l=0.05
s=w*l
Fy=[2.9566, 2.9581, 2.9752, 2.9881]./s
Fx=[1.5134, 1.516,1.5444, 1.5642]./s
N=[0 , 10000, 20000, 30000]
a=[130.6, 133.4, 167.82, 196.35]
a=[
    50,
    100,
    150,
    200,
    250,
    300,
    400,
    500,
    600,
    700,
    800,
    900,
    1000
]


Fx=[1.3476220409334019,
1.480923753028037,
1.530591614486184,
1.5665710215107538,
1.5966542534879409,
1.6236491426243447,
1.6763581680716015,
1.7279966218629852,
1.7794670916046016,
1.8320273816934787,
1.8820294358883984,
1.9291122432914563,
1.9774649338796735]./s

Fy=[2.8981567431474105,
2.9405367891304195,
2.9666839047567919,
2.989724988816306,
3.0115103101124987,
3.0327638013986871,
3.0795919359661639,
3.131694701500237,
3.1908060158602893,
3.2629877664148808,
3.34939130442217,
3.4574740896932781,
3.591245740884915]./s

font = 32
fig = Figure(resolution = (1000, 600))

ax = Axis(fig[1, 1],
    xlabel = L"a (μm)",
    ylabel = L"σ (MPa)",
    xlabelsize = font,
    ylabelsize = font,
    xticklabelsize = font - 5,
    yticklabelsize = font - 5,
    xgridvisible = false,   # kein vertikales Grid
    ygridvisible = false,
    yticks=5:5:40,
    xticks=50:150:1000
)

#ax.xticks = 0:5000:20000

CairoMakie.scatter!(ax, a, Fx,markersize=15,marker=:diamond)
CairoMakie.lines!(ax, a, Fx, linewidth=3,label="Stress in X-axis")
CairoMakie.scatter!(ax, a, Fy,markersize=15, marker=:rect)
CairoMakie.lines!(ax, a, Fy, linewidth=3,  label="Stress in Y-axis")
axislegend(ax; position=:lt,  labelsize=20)

save("stress.svg",fig)

fig