using Statistics, PyCall, Plots, DataFrames, CSVFiles, QuadGK
pyplot()

## CSV import
#cd("./Data")
df_E = DataFrame(load("eigenvalues.csv"))
df_u = DataFrame(load("eigenvectors.csv"))

## Plots
l = 0
X = LinRange(0, 12, length(df_u.y00)-2)
u = convert(Matrix, df_u[3:end, (df_u[1,:] .== l) .& (df_u[2,:] .== 0)])
P0 = Plots.plot(X, u, label=string("u",l,0), show=true)

for n ∈ 1:5
    u = convert(Matrix, df_u[3:end, (df_u[1,:] .== l) .& (df_u[2,:] .== n)])
    Plots.plot!(X, u, label=string("u",l,n))
end


file = string("./eigenstates3Dho_l",l,".pdf")
savefig(P0,file)
