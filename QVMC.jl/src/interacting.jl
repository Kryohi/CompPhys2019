using LinearAlgebra, Statistics, DataFrames, CSV, ProgressMeter
using StaticArrays # statically sized arrays, with faster linear algebra methods
include("common_math.jl")
include("operators.jl")

workdir = ".."
any(x->x=="Data", readdir(workdir)) || mkdir(string(workdir,"/Data"))
any(x->x=="Plots", readdir(workdir)) || mkdir(string(workdir,"/Plots"))

const h2m = 0.5

## Interacting case, with Coulomb interaction modeled by the Jastrow function

macro wfname(N)
    return Symbol(:wf_,N)
end
macro enname(N)
    return Symbol(:en_,N)
end

# real basis set (0,0 solution plus real and imaginary parts of 0,1 solution)
φ1(x,y) = exp(-ω*(x^2+y^2)/2)
φ2(x,y) = x*exp(-ω*(x^2+y^2)/2)
φ3(x,y) = y*exp(-ω*(x^2+y^2)/2)
φ = @SVector [φ1, φ2, φ3]

# Jastrow function, where b_parallel and b_antiparallel are variational parameters
function jas(R, up2down, b_p, b_a)
    N = size(R,1)
    expon = 0.
    spin = [ones(up2down,1); -ones(N-up2down,1)]
    for j=2:N
        for i=1:j-1
            r = sqrt((R[i,1]-R[j,1])^2 + (R[i,2]-R[j,2])^2)
            if spin[i]*spin[j] > 0
                expon += r/(3*(1+b_p*r))
            else
                expon += r/(1+b_a*r)
            end
        end
    end
    return exp(sum(expon))
end

function V_C(R,N)
    VC = 0.
    for j=2:N
        for i=1:j-1
            r = sqrt((R[i,1]-R[j,1])^2 + (R[i,2]-R[j,2])^2)
            VC += 1/r
        end
    end
    return VC
end

# In the interacting case the laplacian of the total wavefunction also contains two
# additional terms coming from the Jastrow factor
function T_L(R, SDu, SDd)
    T_SD = -0.5*sum(lap_psi(R,SDu,SDd))
    T_J =
    T_mixed =
    return T_SD + T_J + T_mixed
end


# Performs a Variational Monte Carlo sampling
function VMC_i(N, ω, rmax, Δ, b_p, b_a; MAX_ITER=10^6, phi_stride=10, energy_stride=10)
    # Initial configuration
    R = rand(N,2)
    R = (R.-0.5).*(rmax/2)
    Ψ_old = 1e-60

    # Slater determinants initialization
    SDu = @MMatrix zeros(ceil(Int,N/2), ceil(Int,N/2))
    SDd = @MMatrix zeros(floor(Int,N/2), floor(Int,N/2))

    # saved data
    wf = DataFrame(iteration=Int64[], x=Float64[], y=Float64[], psi=Float64[], jas=Float64[])
    en = DataFrame(iteration=Int64[], V=Float64[], T=Float64[], TJF=Float64[])

    up2down = ceil(Int,N/2) # where to switch to the spin down electrons
    println("The system has $up2down electrons with spin up and $(N-up2down) with spin down")
    mean_acc = 0

    println("Starting VMC sampling with N=$N")

    for n = 1:MAX_ITER
        # new trial positions
        R_ = R .+ Δ.*(rand(N,2).-0.5)
        R_ .= R_ .- rmax.*round.(R_./rmax) #PBC

        # construction of the Slater matrices
        SDu = [φ[j](R_[i,1],R_[i,2]) for j=1:up2down, i=1:up2down]
        SDd = [φ[j](R_[i,1],R_[i,2]) for j=1:up2down-N%2, i=up2down+1:N]

        # wavefunction in R_
        jastrow = jas(R_,up2down,b_p,b_a)
        Ψ = det(SDu)*det(SDd)*jastrow

        # save data on the wavefunction sampling
        if n%phi_stride == 0
            for i=1:N
                push!(wf, [n; R_[i,1]; R_[i,2]; Ψ; jastrow])
            end
        end

        # Energy
        if n%energy_stride == 0
            VL = 0.5*ω^2*sum(R_[:,1].^2 .+ R_[:,2].^2) + V_C(R_,N)
            TL = T_L(R_,SDu,SDd)
            TJF = 0.#T_JF(R_,SDu,SDd)
            push!(en, [n; VL; TL; TJF])
        end

        # acceptance ratio
        acc = Ψ*Ψ / (Ψ_old*Ψ_old)

        if acc > rand()
            R = R_
            Ψ_old = Ψ
            mean_acc += 1
        end
    end

    mean_acc /= MAX_ITER
    println("Finished VMC sampling for N=$N, with average acceptance ratio $mean_acc")
    return wf, en
end

# Gradient descent for the minimization of energy with
# respect to the b↑↑ and b↑↓ parameters of the Jastrow function
function energy_min(N, ω, rmax, Δ, γ, Δb; K_grid=10, tol=1e-5)

    b_grid = LinRange(0.1, 0.6, K_grid)
    E_grid = zeros(K_grid, K_grid)
    ΔE_grid = zeros(K_grid, K_grid)

    @showprogress 1 "exploring grid..." for i=1:K_grid, j=1:K_grid
        @info "Sampling at b_p=$(b_grid[i]), b_a=$(b_grid[j])"
        wf, en = VMC_i(N, ω, rmax, Δ, b_grid[i], b_grid[j], MAX_ITER=1*10^6)
        E_grid[i,j] = mean(en.T + en.V)
        ΔE_grid[i,j] = std(en.T + en.V)
        @info "Finished with E = $(round(E_grid[i,j],digits=4)) ± $(round(ΔE_grid[i,j],digits=4))"
        println(" ")

        #CSV.write("$workdir/Data/psi_grid_$(i)_$(j)_$N.csv", wf)
        CSV.write("$workdir/Data/energy_grid_$(i)_$(j)_$N.csv", en)
    end

    return b_grid, E_grid, ΔE_grid
end



rmax = 6.
Δ = 0.5
ω = 1
γ = 0.1
Δb = 0.01
N = 3

b_grid, E_grid, ΔE_grid = energy_min(N, ω, rmax, Δ, γ, Δb; tol=1e-5)
