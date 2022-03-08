cd(@__DIR__) #using this directory as working directory
using Pkg; Pkg.activate("."); Pkg.instantiate() # activate virtual environment


"""
Numerical solution
    - Time integration using Runge-Kutta third order
    - 5th-order Compact WENO scheme for spatial terms
    - Periodic Boundary Conditions
    

    xlim :tuple: min value of x and max value of x
    dx :float: step value of x
    tlim :tuple: min value of t and max value of t
    dt :float: step value of t
    icf :function: initial condition function
    ff :function: flux, function of 'u'
    ns :int, "all": number of solutions to output
    
    Modifies:
        u :array: solution field, dims: (nx, ns)

#Julia-Tip: array start at index 1, not 0
#Julia-Tip: functions can mutate the arguments
    use ! to mark them explicitly
#Julia-Tip: tuples are immutable
#Julia-Tip: map functions to arrays using '.'
    Ex: foo.(x)
"""
function numerical_weno_p!(xlim, dx, tlim, dt, icf, ff, ns, u)
    # nx is the number of x points
    nx = Int64(abs((xlim[2]-xlim[1])/dx))
    # nt is the number of iterations
    nt = Int64(abs((tlim[2]-tlim[1])/dt))

    x = LinRange(xlim[1], xlim[2], nx) # grid x
    un = Array{Float64}(undef, nx) # solution u for each time
    u_aux = Array{Float64}(undef, nx) # auxiliary array for RK3 integration
    r = Array{Float64}(undef, nx)

    # we can output 'u' for each iteration ("all") or at a certain frequency
    if ns == "all" || typeof(ns) == String
        ns = nt
    end
    freq = Int64(nt/ns)

    # store the solution at t=0
    u[:, 1] = icf.(x)
    
    for n = 1:nt # iterations (time step)

        # rhs?

end


"""
Right-hand side terms

'The conservative form of the inviscid Burgers' eq allows us to use
the flux splitting method and WENO reconstruction to compute the flux
at the interface.

    #1 We first compute the flux at the nodal points.
    #2 Then split it into positive and negative components. 
        This process is called as Lax-Friedrichs flux splitting.
        f+- = 0.5*(f +- a*u)
        a := absolute value of the flux Jacobian (wave speed)
'

    nx :int: number of x points
    dx :float: step value of x
    u :array: solution field, dims: (nx, ns)
    ff :function: flux, function of 'u'

    Modifies:
        r :array: right-hand side term

"""
function rhs!(nx, dx, u, ff, r)
    #1 We first compute the flux at the nodal points.
    f = ff.(u)
    #2 Lax-Friedrichs flux splitting.
    fP = 0.5*()
end

"""
Computes the absolute value of the flux Jacobian

'The wavespeed is the absolute value of the flux Jacobian.
We chose the values of `a` as the maximum value of u_i over the
stencil used for WENO-5 reconstruction, i.e.,
    a = max(|u_i-2|,|u_i-1|,|u_i|,|u_i+1|,|u_i+2|)
'

Periodic Boundary Conditions!! (1,2,n-1,n) 3 possible ways:
    1. Hardcoded: as simple as that, poor generalization
    2. Circular Array: not currently implemented (negative index,zero-index)
    3. Modular indexing: zero-index doesn't exist

    nx :int: number of x points
    u :array: solution field, dims: (nx, ns)

    Modifies:
        a :array: alpha, wave speed

"""
function wavespeed(nx, u, a)
    u_aux = Array{Float64}(undef, nx+4)
    u_aux[3:nx+2] = u
    u_aux[1] = u_aux[nx+1]
    u_aux[2] = u_aux[nx+2]
    u_aux[nx+3] = u_aux[3]
    u_aux[nx+4] = u_aux[4]

    for i = 3:nx+2
        a[i] = max(abs(u_aux[i-2]), abs(u_aux[i-1]),abs(u_aux[i]),abs(u_aux[i+1]),abs(u_aux[i+2]))
    end
    # for i = 3:nx-2
    #     a[i] = max(abs(u[i-2]), abs(u[i-1]),abs(u[i]),abs(u[i+1]),abs(u[i+2]))
    # end
    # # periodicity
	# i = 1
	# ps[i] = max(abs(u[n-1]), abs(u[n]), abs(u[i]), abs(u[i+1]), abs(u[i+2]))
	# i = 2
	# ps[i] = max(abs(u[n]), abs(u[i-1]), abs(u[i]), abs(u[i+1]), abs(u[i+2]))
	# i = n-1
	# ps[i] = max(abs(u[i-2]), abs(u[i-1]), abs(u[i]), abs(u[i+1]), abs(u[1]))
	# i = n
	# ps[i] = max(abs(u[i-2]), abs(u[i-1]), abs(u[i]), abs(u[1]), abs(u[2]))
end