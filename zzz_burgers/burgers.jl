cd(@__DIR__) #using this directory as working directory
using Pkg; Pkg.activate("."); Pkg.instantiate() # activate virtual environment


"""
Numerical solution
    - Time integration using Runge-Kutta third order
    - 5th-order Compact WENO scheme for spatial terms
    - Dirichlet Boundary Conditions
    

    xlim :tuple: min value of x and max value of x
    dx :float: step value of x
    tlim :tuple: min value of t and max value of t
    dt :float: step value of t
    ic :function: initial condition function
    u :array: solution field, dims: (nx, ns)
    f :function: flux, function of 'u'
    ns :int, "all": number of solutions to output 

#Julia-Tip: array start at index 1, not 0
#Julia-Tip: functions can mutate the arguments
    use ! to mark them explicitly
#Julia-Tip: tuples are immutable
#Julia-Tip: map functions to arrays using '.'
    Ex: foo.(x)
"""
function numerical_weno_d!(xlim, dx, tlim, dt, ic, u, f, ns)
    # nx is the number of x points
    nx = Int64(abs((xlim[2]-xlim[1])/dx))
    nt = Int64(abs((tlim[2]-tlim[1])/dt))

    x = LinRange(xlim[1], xlim[2], nx) # grid x
    un = Array{Float64}(undef, nx) # solution u for each time
    u_aux = Array{Float64}(undef, nx) # auxiliary array for RK3 integration
    r = Array{Float64}(undef, nx)

    # we can output 'u' for each time ("all") or at a certain frequency
    if ns == "all" || typeof(ns) == String
        ns = nt
    end
    freq = Int64(nt/ns)

    # store the solution at t=0
    u[:, 1] = ic.(x)
    

end
