using LinearAlgebra, Printf, DataFrames, CSV
using StaticArrays # statically sized arrays, with faster linear algebra methods
include("common_math.jl")
include("operators.jl")

workdir = ".."
any(x->x=="Data", readdir(workdir)) || mkdir(string(workdir,"/Data"))
any(x->x=="Plots", readdir(workdir)) || mkdir(string(workdir,"/Plots"))

const h2m = 0.5

function T_L(R, SDu, SDd)
    return -0.5*sum(lap_psi(R,SDu,SDd))
end

function T_JF(R, SDu, SDd)
    gg = grad_psi(R,SDu,SDd)
    return 0.25*sum(gg.*gg) - 0.25*sum(lap_psi(R,SDu,SDd))
end

# real basis set (0,0 solution plus real and imaginary parts of 0,1 solution)
φ1(x::Float64,y::Float64) = exp(-ω*(x^2+y^2)/2)
φ2(x::Float64,y::Float64) = x*exp(-ω*(x^2+y^2)/2)
φ3(x::Float64,y::Float64) = y*exp(-ω*(x^2+y^2)/2)
φ = @SVector [φ1, φ2, φ3]


function VMC_ni(N::Int64, ω, rmax::Float64, Δ::Float64; MAX_ITER=10^6, phi_stride=10, energy_stride=10)
    # Initial configuration
    R = rand(N,2)
    R = (R.-0.5).*(rmax/2)
    R_ = zeros(size(R))
    Ψ, Ψ_old = 1e-60, 1e-60

    # Slater matrices initialization
    SDu = zeros(ceil(Int,N/2), ceil(Int,N/2))
    SDd = zeros(floor(Int,N/2), floor(Int,N/2))

    # saved data
    wf = DataFrame(iteration=Int64[], x=Float32[], y=Float32[], psi=Float32[])
    en = DataFrame(iteration=Int64[], V=Float64[], T=Float64[], TJF=Float64[])

    up2down = ceil(Int,N/2) # where to switch to the spin down electrons
    mean_acc = 0

    @info "Starting VMC sampling with N=$N"

    for n = 1:MAX_ITER
        # new trial positions
        R_ = R .+ Δ.*(rand(N,2).-0.5)
        R_ = @. R_ - rmax*round(R_/rmax) #PBC

        # construction of the Slater matrices
        SDu = [φ[j](R_[i,1],R_[i,2]) for j=1:up2down, i=1:up2down]
        SDd = [φ[j](R_[i,1],R_[i,2]) for j=1:up2down-N%2, i=up2down+1:N]

        # wavefunction in R_
        Ψ = det(SDu)*det(SDd)

        # save data on the wavefunction sampling
        if n%phi_stride == 0
            for i=1:N
                push!(wf, [n; Float32(R_[i,1]); Float32(R_[i,2]); Float32(Ψ)])
            end
        end

        # Energy
        if n%energy_stride == 0
            VL = 0.5*ω^2*sum(R_[:,1].^2 .+ R_[:,2].^2)
            TL = T_L(R_,SDu,SDd)
            TJF = T_JF(R_,SDu,SDd)
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
    @info "Finished VMC sampling for N=$N, with average acceptance ratio $mean_acc\n"
    return wf, en
end

macro wfvar(N)
    return Symbol(:wf_,N)
end
macro envar(N)
    return Symbol(:en_,N)
end

rmax = 6.
Δ = 0.5 # should be probably optimized
ω = 1

Threads.@threads for N = [2,3,4,5,6]
    @time @wfvar(N), @envar(N) = VMC_ni(N, ω, rmax, Δ, MAX_ITER=15*10^6)
    CSV.write("$workdir/Data/psi_ni_$N.csv", @wfvar(N))
    CSV.write("$workdir/Data/energy_ni_$N.csv", @envar(N))
end


## Testing stuff
N = 5
wf, en = VMC_ni(N, ω, rmax, Δ, MAX_ITER=1*10^6)
Juno.@profiler VMC_ni(N, ω, rmax, Δ, MAX_ITER=2*10^6)
CSV.write("$workdir/Data/psi_ni_$N.csv", wf)
CSV.write("$workdir/Data/energy_ni_$N.csv", en)
