using Statistics, Plots, DataFrames, CSVFiles, LaTeXStrings
gr()

## CSV import
cd("./Data_scattering")
df_E = DataFrame(load("eigenvalues.csv"))
df_u = DataFrame(load("eigenvectors.csv"))

## Harmonic oscillator
l = 0
X = LinRange(0, 12, length(df_u.y00)-2)
u = convert(Matrix, df_u[3:end, (df_u[1,:] .== l) .& (df_u[2,:] .== 0)])
P0 = Plots.plot(X, u, label=string("u",l,0), show=true)

for n ∈ 0:0
    u = convert(Matrix, df_u[3:end, (df_u[1,:] .== l) .& (df_u[2,:] .== n)])
    Plots.plot!(X, u, label=string("u",l,n))
    gui()
end
u = convert(Matrix, df_u[3:end, (df_u[1,:] .== 1) .& (df_u[2,:] .== 0)])
Plots.plot!(X, u, label=string("u",1,0))
u = convert(Matrix, df_u[3:end, (df_u[1,:] .== 2) .& (df_u[2,:] .== 0)])
Plots.plot!(X, u, label=string("u",2,0))
u = convert(Matrix, df_u[3:end, (df_u[1,:] .== 3) .& (df_u[2,:] .== 0)])
Plots.plot!(X, u, label=string("u",3,0))


file = string("./eigenstates3D_l",l,".pdf")
savefig(P0,file)


## SCATTERING
cd("../Data_scattering")
df_u = DataFrame(load("wavefunction_l6.csv"))
N = size(df_u,1)
X = LinRange(0,10,N)
S = plot(X, [zeros(900,1);df_u.Veff[901:end]], label="Veff") # the munber is xmin
plot(X, df_u.y, label="u(x)")

file = string("./scatterwave_l",6,".pdf")
savefig(S,file)


# total cross-section
df_tcs = DataFrame(load("total_csection.csv"))
tcs = plot(df_tcs.E, df_tcs.sigma_tot, label=L"\sigma_{tot}")
savefig(tcs,"./crosssectionfine.pdf")

# phase-shifts
df_pss = DataFrame(load("phase_shifts.csv"))
pss = plot(df_pss.E, df_pss.l0, label=L"\delta_{0}")
plot!(df_pss.E, df_pss.l1, label=L"\delta_{1}")
plot!(df_pss.E, df_pss.l2, label=L"\delta_{2}")
plot!(df_pss.E, df_pss.l3, label=L"\delta_{3}")
plot!(df_pss.E, df_pss.l4, label=L"\delta_{4}")
plot!(df_pss.E, df_pss.l5, label=L"\delta_{5}")
plot!(df_pss.E, df_pss.l6, label=L"\delta_{6}")
savefig(pss,"./phase_shifts.pdf")

pssc = plot(df_pss.E, sin.(df_pss.l0).^2, label=L"\delta_{0}")
plot!(df_pss.E, sin.(df_pss.l1).^2 .*3, label=L"\delta_{1}")
plot!(df_pss.E, sin.(df_pss.l2).^2 .*5, label=L"\delta_{2}")
plot!(df_pss.E, sin.(df_pss.l3).^2 .*7, label=L"\delta_{3}")
plot!(df_pss.E, sin.(df_pss.l4).^2 .*9, label=L"\delta_{4}")
plot!(df_pss.E, sin.(df_pss.l5).^2 .*11, label=L"\delta_{5}")
plot!(df_pss.E, sin.(df_pss.l6).^2 .*13, label=L"\delta_{6}")
savefig(pssc,"./phase_shifts_contributions.pdf")
