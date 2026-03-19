using JSON
using PyPlot, LaTeXStrings

rcParams = PyPlot.PyDict( PyPlot.matplotlib."rcParams" )
rcParams["font.size"] = 13.5
rcParams["text.usetex"] = true
rcParams["axes.linewidth"] = 1

d = JSON.parse(read( "data-consensus-times.json", String))

fig, ax = subplots()
fig.set_figheight(4)
fig.set_figwidth(6)

ax.tick_params( axis= "both",direction = "in", top = true, 
               right = true, length=7, width=1)
ax.tick_params( axis= "both", which="minor", direction = "in", 
               top = true, right = true, length= 3.5, width= 1)

vlines( 0.805, 8, 100, "gray", "dashed")
vlines( 0.954, 8, 100, "gray", "dashed")
vlines( 1.01, 8, 100, "gray", "dashed")
vlines( 1.04, 8, 100, "gray", "dashed")

ax.semilogy( d["Dx"], d["11"], lw=2, label=L"$N=11$" )
ax.semilogy( d["Dx"], d["15"], lw=2, label=L"$N=15$" )
ax.semilogy( d["Dx"], d["21"], lw=2, label=L"$N=21$" )
ax.semilogy( d["Dx"], d["31"], lw=2, label=L"$N=31$" )
ax.legend()

ax.set_ylim(8,100)
ax.set_yticks([10,100], labels=[L"$10^1$",L"$10^2$"])
    
ax.set_xlabel(L"$D_x$",fontsize="17")
ax.set_ylabel(L"$t_c$",fontsize="17")
tight_layout()     
