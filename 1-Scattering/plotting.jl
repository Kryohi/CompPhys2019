using Statistics, PyCall, Plots, DataFrames, CSVFiles, QuadGK
gr()

## CSV import
cd("./Data_scattering")
df_E = DataFrame(load("eigenvalues.csv"))
df_u = DataFrame(load("eigenvectors.csv"))

## Plots
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
df_u = DataFrame(load("wavefunction_l3.csv"))
N = size(df_u,1)
X = LinRange(0,7,N)
S = plot(X, [zeros(900,1);df_u.Veff[901:end]], label="Veff") # the munber is xmin
plot!(X, df_u.y, label="u(x)")

file = string("./scatterwave_l",3,".pdf")
savefig(S,file)
