using OrdinaryDiffEq, DiffEqCallbacks
using Graphs
using IterTools
using SparseArrays, BlockArrays, LinearAlgebra
using PyPlot
using JSON
using DelimitedFiles

#=
Dynamical solution to the problem:

\begin{align}
\vec{v} &= \dot\vec{y}
\dot{\vec{v}} &= -(\mathcal{L}- \mathbb{I})\vec{v} - \mathcal{L}^2 \vec{y}
\end{align}
=#

ζ(t) = 1.0/(t + 10000 )

function create_graph_n_nodes( n )
    g1 = Graphs.path_graph(n)
    for i in 1:(n-1)
        add_edge!( g1, i, i+1 )
    end
    add_edge!( g1, 1, n )

    g2 = Graphs.path_graph( n )
    for (i,j) in IterTools.subsets( 1:n, Val{2}() )
        add_edge!( g2, i, j )
    end

    return g1, g2
end


function create_supralaplacian_op( n, d1, d2, dx )
    g1, g2  = create_graph_n_nodes( n )
    L1 = Graphs.LinAlg.laplacian_matrix(g1)
    L2 = Graphs.LinAlg.laplacian_matrix(g2)
    
    sizeL = size( L1, 1 )
    Ones = Matrix( I, sizeL, sizeL )
    SL = BlockArray{Float64}( undef_blocks, [sizeL,sizeL], [sizeL,sizeL]  )

    SL[Block(1,1)] = d1*L1 + dx*Ones
    SL[Block(1,2)] = -dx*Ones
    SL[Block(2,1)] = -dx*Ones
    SL[Block(2,2)] = d2*L2 + dx*Ones
    return Array(SL)
end


function plot_solution( sol::ODESolution, m::Int64 )
    s_ = size( sol.u, 1 )
    n_ = size( sol.u[1], 1 )
    z = Array{Float64}( undef, s_ , n_ )

    for i in 1:s_
        z[i,:] = sol.u[i][:]
    end

    #rcParams = PyPlot.PyDict( PyPlot.matplotlib."rcParams" )
    #rcParams["font.size"] = 13.5
    #rcParams["text.usetex"] = true
    #rcParams["axes.linewidth"] = 1

    fig, ax = subplots()
    fig.set_figheight(4)
    fig.set_figwidth(6)

    ax.tick_params( axis= "both",direction = "in", top = true, 
                    right = true, length=7, width=1 )
    ax.tick_params( axis= "both", which="minor", direction = "in", 
                    top = true, right = true, length= 3.5, width= 1 )
    
    for i in 1:m
        style = (i <= Int(m/2) ) ? "dashed" : "solid"
        ax.plot( sol.t, z[:,i], linestyle = style )
    end

    #ax.set_ylabel(L"$y$",fontsize="17")
    #ax.set_xlabel(L"$t$",fontsize="17")
    tight_layout()
end


function integrate_system( n, d1, d2, dx; to::Float64 = 5000.0 )
    L = create_supralaplacian_op( n, d1, d2, dx )
    n = size(L,1)
    Identity = Matrix(I,n,n)

    function rhs!( du, u, p, t)
        du[1:n] = u[n+1:end]
        du[n+1:end] = -( L +  Identity )*u[n+1:end] - L^2*u[1:n]
    end
    
    function condition( u, t, integrator ) 
        s = 0.0
        for (i,j) in subsets(1:n,2)
            s += abs(u[i] - u[j])
        end
        @show s
        return s <= 1e-2
    end

    affect!( integrator ) = terminate!( integrator )
    cb = DiscreteCallback( condition, affect! )

    p = []
    tspan = ( 0.0, to )
    u0 = vcat( 1:1:(n*1), fill(0.0,n) )
    
    prob = ODEProblem{true}( rhs!, u0, tspan, p )
    sol = solve( prob, RK4(), callback=cb )
    plot_solution( sol, n )

    return sol.t[end]
end


function integrate_system_gd( n, d1, d2, dx )
    L = create_supralaplacian_op( n, d1, d2, dx )
    n = size(L,1)
    Identity = Matrix(I,n,n)

    function rhs!( du, u, p, t)
        du[1:n] = -( L +  ζ(t)*Identity )*u[1:n] - ζ(t)*[1:n...]
    end
    
    function condition( u, t, integrator ) 
        s = 0.0
        for (i,j) in subsets( 1:n, 2 )
            s += abs(u[i] - u[j])
        end
        @show t,s
        return s <= 1e-3
    end

    affect!( integrator ) = terminate!(integrator)
    cb = DiscreteCallback( condition, affect! )

    p = ()
    tspan = ( 0.0, 500.0 )
    u0 = 2:2:(2*n)
    
    prob = ODEProblem{true}( rhs!, u0, tspan, p )
    sol = solve( prob, RK4(), callback=cb )
    plot_solution( sol, n )

    return sol.t[end]
end


function main()
    Dx = LinRange(0.5,1.5,1000)
    Times = Dict()

    for nnodes in [ 11, 15, 21, 31 ]
        Tc = Dict()
        Threads.@threads for dx in Dx
            Tc[dx] = integrate_system( nnodes, 1, 1, dx )
        end
        
        tc = sort( collect(Tc), by = x->x[1] )
        Times["$nnodes"] = map( x->x[2], tc )
    end
    Times["Dx"] = Dx

#     open("data.json","w") do f
#         write(f, JSON.json(Times) )
#     end
end

integrate_system( 5, 1.0, 1.0, 0.1; to = 130.0 )
# integrate_system( 11, 1.0, 1.0, 1.0; to = 6.0 )
