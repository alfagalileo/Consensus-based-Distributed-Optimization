using OrdinaryDiffEq, DiffEqCallbacks
using Graphs
using IterTools
using SparseArrays, BlockArrays, LinearAlgebra
using PyPlot
using JSON
using DelimitedFiles
using PyCall

matplotlib = pyimport("matplotlib")


#=
Dynamical solution to the problem:

\begin{align}
\vec{v} &= \dot\vec{y}
\dot{\vec{v}} &= -(\mathcal{L} - \mathbb{I})\vec{v} - \mathcal{L}^2 \vec{y} - \alpha(\mathbf{q}) \tilde{\mathbb{1}}
\end{align}
=#


function create_graph_n_nodes( n )
    # all connected graph
    g1 = Graphs.path_graph( n )
    for (i,j) in IterTools.subsets( 1:n, Val{2}() )
        add_edge!( g1, i, j )
    end

    # medium layer
    g2 = Graphs.path_graph(n)
    for i in 1:(n-1)
        add_edge!( g2, i, i+1 )
    end
    add_edge!( g2, n, 1 )

    add_edge!(g2, 1, 3)
    add_edge!(g2, 1, 6)
    add_edge!(g2, 2, 5)
    add_edge!(g2, 4, 6)

    # circular graph
    g3 = Graphs.path_graph(n)
    for i in 1:(n-1)
        add_edge!( g3, i, i+1 )
    end
    add_edge!( g3, n, 1 )

    return [g1, g2, g3]
end


function create_supralaplacian_op( n, d, dx )
    g = create_graph_n_nodes( n )

    L = []
    for ( i, gi ) in enumerate(g)
        push!( L, Graphs.LinAlg.laplacian_matrix(gi) )
    end
    
    sizeL = size( L[1], 1 )
    Ones = Matrix( I, sizeL, sizeL )
    SL = BlockArray{Float64}( zeros(3*sizeL, 3*sizeL), fill(sizeL, 3), fill(sizeL, 3)  )

    SL[Block(1,1)] = d[1]*L[1] + dx*Ones
    SL[Block(2,2)] = d[2]*L[2] + dx*Ones
    SL[Block(3,3)] = d[3]*L[3] + dx*Ones

    SL[Block(1,2)] = -dx*Ones
    SL[Block(2,1)] = -dx*Ones
    SL[Block(3,2)] = -dx*Ones
    SL[Block(2,3)] = -dx*Ones

    return SL
end


function plot_solution( sol::ODESolution, m::Int64 )
    s_ = size( sol.u, 1 )
    n_ = size( sol.u[1], 1 )
    z = Array{Float64}( undef, s_ , n_ )

    for i in 1:s_
        z[i,:] = sol.u[i][:]
    end

    rcParams = PyPlot.PyDict( PyPlot.matplotlib."rcParams" )
    rcParams["font.size"] = 14
    rcParams["text.usetex"] = true
    rcParams["axes.linewidth"] = 1

    fig, ax = subplots()
    fig.set_figheight(4)
    fig.set_figwidth(5.5)

    ax.tick_params( which= "both",direction = "in", top = true,
                    right = true, length=6, width=1 )
    ax.tick_params( which="minor", direction = "in",
                    top = true, right = true, length= 3, width= 1 )
    
    cstyle = ["b", "y", "r"]
    mstyle = ["^", "o", "s"]
    labi = [L"$\mathcal{G}$", L"$\mathcal{K}$",L"$\mathcal{T}$"]

    for (i,parti) in enumerate( partition(1:m, 7) )
        k = 0
        for j in parti
            lab = (k==0) ? labi[i] : ""
                ax.plot( sol.t, z[:,j], color = cstyle[i] )
            k += 1
        end
    end

    for (i,parti) in enumerate( partition(1:m, 7) )
        k = 0
        for j in parti
            lab = (k==0) ? labi[i] : ""
                ax.plot( sol.t[1:6:end], z[:,j][1:6:end], linestyle="", marker= mstyle[i], markersize=5,
                        color = cstyle[i], label=lab )
            k += 1
        end
        ax.legend()
    end

    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(5))
    ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(2.5))

    ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(5))
    ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(2.5))

    ax.set_ylabel(L"$q$",fontsize="17")
    ax.set_xlabel(L"$t$",fontsize="17")
    tight_layout()
end


function rhs!( du, u, p, t )
    n, ϕ, Pd, Gd, L, Id, q0, idx = p
    q = u[1:n]
    v = u[n+1:end]
    du[1:n] = v
    du[n+1:end] = -( L +  Id )*v - L^2*q - ( sum(q[i] for i in idx[1] ) + sum((1-ϕ)*q[i] for i in idx[2] ) + sum(q[i] for i in idx[3]) - Pd - Gd )*q0
end


function condition( u, t, integrator )
    ϵ = 0.0
    n = length(u)
    for (i,j) in subsets(1:n,2)
        ϵ += abs(u[i] - u[j])
    end
    @info "t = $t \t ϵ = $ϵ"
    return ϵ <= 1e-3
end


function integrate_system( n, d, dx; to::Float64 = 5000.0 )
    ϕ = 0.7
    Pd = 100.0
    Gd = 100.0

    L = create_supralaplacian_op( n, d, dx )
    n = size(L,1)
    idx = collect( partition( 1:n, Int64(n/3) ) )
    Id = Matrix( I, n, n )

    q0 = ones(n)
    q0[ idx[2][1]:idx[2][end] ] .= 1 - ϕ

    p = [ n, ϕ, Pd, Gd, L, Id, q0, idx ]
    tspan = ( 0.0, to )
    u0 = vcat( 1:1:(n*1), fill(0.0,n) )
    
    prob = ODEProblem{true}( rhs!, u0, tspan, p )

    affect!( integrator ) = terminate!( integrator )
    cb = DiscreteCallback( condition, affect! )

    sol = solve( prob, AutoTsit5(Rosenbrock23()), callback=cb )
    plot_solution( sol, n )

    return sol.t[end]
end

function main()
    Dx = LinRange(0.1,1.0,100)
    Times = Dict()
    for nnodes in [11]
        Tc = Dict()
        Threads.@threads for dx in Dx
            Tc[dx] = integrate_system( nnodes, [1.0,1.0,1.0], dx; to = 1020.0 )
        end
        
        tc = sort( collect(Tc), by = x->x[1] )
        Times["$nnodes"] = map( x->x[2], tc )
    end
    Times["Dx"] = Dx

    open("data_threelayers.json","w") do f 
        write(f, JSON.json(Times) )
    end
end

integrate_system( 7, [0.2,0.8,0.2], 0.6; to = 30.0 )
# main()
