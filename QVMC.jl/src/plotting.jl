using DataFrames, CSV, Plots, LaTeXStrings, ColorSchemes, Statistics
using ThreadTools
using Plots.PlotMeasures
workdir = ".."

N = 6
wf = DataFrame(CSV.File("$workdir/Data/psi_ni_$N.csv"))
en = DataFrame(CSV.File("$workdir/Data/energy_ni_$N.csv"))

mean(wf.y), mean(wf.x)

function binn(wf, nb)
    pdf = wf.psi .^ 2
    binned_psi = zeros(Float64, nb, nb)
    binned_psi2 = zeros(Float64, nb, nb)
    max_r = ceil(maximum(abs.([wf.x; wf.y])))/2
    X = LinRange(-max_r, max_r, nb)
    @inbounds for i = 1:nb-1, j = 1:nb-1
        idxs = (wf.x .> X[i]) .* (wf.x .< X[i+1]) .* (wf.y .> X[j]) .* (wf.y .< X[j+1])
        binned_psi[i,j] = sum(wf.psi[idxs])
        binned_psi2[i,j] = sum(pdf[idxs])
    end
    return X, binned_psi, binned_psi2
end

X, binned_psi, binned_psi2 = binn(wf, 50)
p1 = Plots.contourf(X, X, binned_psi2, aspect_ratio=1)


## Comparison plots

wf = [DataFrame(CSV.File("$workdir/Data/psi_ni_$n.csv")) for n = [2,3,4,5,6]]
nb = 50 # number of bins for the contour plot

data = map(n->binn(wf[n],nb), [1,2,3,4,5])

pyplot(size=(1200, 830))
fnt = "DejaVu Sans"
default(titlefont=Plots.font(fnt,24), guidefont=Plots.font(fnt,18), tickfont=Plots.font(fnt,12), legendfont=Plots.font(fnt,14))

# Grafico con prob densities
p1 = Plots.contourf(data[1][1], data[1][1], data[1][3], aspect_ratio=1, legend=:none);
p2 = Plots.contourf(data[2][1], data[2][1], data[2][3], aspect_ratio=1, legend=:none);
p3 = Plots.contourf(data[3][1], data[3][1], data[3][3], aspect_ratio=1, legend=:none);
p4 = Plots.contourf(data[4][1], data[4][1], data[4][3], aspect_ratio=1, legend=:none);
p5 = Plots.contourf(data[5][1], data[5][1], data[5][3], aspect_ratio=1, legend=:none);
Plots.plot(p1,p2,p3,p4,p5, layout=5, title=["N=2" "N=3" "N=4" "N=5" "N=6"],
title_location=:center, left_margin=[0mm 0mm], bottom_margin=16px, xrotation=60, reuse=false, legend=:none)

savefig("$workdir/Plots/comparison_psi2.pdf")

# Grafico con wavefunctions NOTE sembrerebbe essere inutile/dannoso
p1 = Plots.contourf(data[1][1], data[1][1], data[1][2], aspect_ratio=1, legend=:none);
p2 = Plots.contourf(data[2][1], data[2][1], data[2][2], aspect_ratio=1, legend=:none);
p3 = Plots.contourf(data[3][1], data[3][1], data[3][2], aspect_ratio=1, legend=:none);
p4 = Plots.contourf(data[4][1], data[4][1], data[4][2], aspect_ratio=1, legend=:none);
p5 = Plots.contourf(data[5][1], data[5][1], data[5][2], aspect_ratio=1, legend=:none);
Plots.plot(p1,p2,p3,p4,p5, layout=5, title=["N=2" "N=3" "N=4" "N=5" "N=6"],
title_location=:center, left_margin=[0mm 0mm], bottom_margin=16px, xrotation=60, reuse=false, legend=:none)

savefig("$workdir/Plots/comparison_psi.pdf")


## Energy

mean(en.T[2:end] .+ en.V[2:end])
std(en.T[2:end] .+ en.V[2:end])
std(en.V[2:end])
std(en.T[2:end])
mean(en.TJF[2:end] .+ en.V[2:end])
std(en.TJF[2:end])

plot(en.V[2:end] .+ en.T[2:end])


## Energy comparison
ens = [DataFrame(CSV.File("$workdir/Data/energy_ni_$n.csv")) for n = [2,3,4,5,6]]

avg_tot = [mean(ens[n].T .+ ens[n].V) for n in [1,2,3,4,5]]
std_tot = [std(ens[n].T .+ ens[n].V) for n in [1,2,3,4,5]]




## Interacting case
p1 = Plots.contourf(b_grid, b_grid, E_grid, aspect_ratio=1, legend=:none)





## Autocorrelation function, to choose the stride parameters
using FFTW
function fft_acf(H::Array{Float64,1}, k_max::Int)

    Z = H .- mean(H)
    fvi = rfft(Z)
    acf = fvi .* conj.(fvi)
    acf = ifft(acf)
    acf = real.(acf)
    C_H = acf[1:k_max]

    return C_H./C_H[1]
end

acf_E = fft_acf(en.T .+ en.V, length(en.T)÷4)
plot(acf_E[1:1000])
