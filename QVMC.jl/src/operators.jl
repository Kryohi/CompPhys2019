## Operators

function grad_orb(x::Float64,y::Float64)
    grad = @SMatrix [ φ1(x,y)*(-ω*x) φ1(x,y)*(-ω*y) ;
                      φ1(x,y)*(1-ω*x^2) φ1(x,y)*(-ω*x*y) ;
                      φ1(x,y)*(-ω*x*y) φ1(x,y)*(1-ω*y^2) ]
    return grad # returns J*2 matrix
end

function lap_orb(x::Float64,y::Float64)
    r2 = x^2 + y^2
    lap = [ φ1(x,y) * (ω^2*r2 -2ω),
            φ2(x,y) * (ω^2*r2 -4ω),
            φ3(x,y) * (ω^2*r2 -4ω)]
    return lap # returns array of J elements
end

# returns the gradient of all the particle wavefunctions (N*2 Matrix)
function grad_psi(R, SDu, SDd)
    g = zeros(size(R))
    Nu, Nd = size(SDu,1), size(SDd,1)
    SDu_1, SDd_1  = inv(SDu), inv(SDd)
    for i=1:Nu
        g_orb = grad_orb(R[i,1], R[i,2])
        g_ = zeros(2,1)
        for j=1:Nu
            g_ .+= g_orb[j,:] .* SDu_1[i,j]
        end
        g[i,:] = g_
    end
    for i=Nu+1:size(R,1)
        g_orb = grad_orb(R[i,1], R[i,2])
        g_ = zeros(2,1)
        for j=1:Nd
            g_ .+= g_orb[j,:] .* SDd_1[i-Nu,j]
        end
        g[i,:] = g_
    end
    return g
end

# returns the laplacian of all the particle wavefunctions
function lap_psi(R, SDu, SDd)
    l = zeros(size(R,1),1)
    Nu, Nd = size(SDu,1), size(SDd,1) # number of electrons with spin
    SDu_1, SDd_1  = inv(SDu), inv(SDd)
    for i=1:Nu
        lo = lap_orb(R[i,1], R[i,2])
        for j=1:Nu
            l[i] += lo[j] * SDu_1[i,j]
        end
    end
    for i=Nu+1:size(R,1)
        lo = lap_orb(R[i,1], R[i,2])
        for j=1:Nd
            l[i] += lo[j] * SDd_1[i-Nu,j] # check indices
        end
    end
    return l
end


## Operators for the interacting case (concerning the Jastrow function)

function grad_J(R, up2down, b_p, b_a)
    g = zeros(size(R))
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

end





#rad_solution(r::Vector,ω,m,k)  = ω^(abs(m)/2) .* r.^(abs(m)) .* exp.(-ω.*r.^2 ./2) .* L(ω.*r.^2,abs(m),k)
#phi_nonint(r,k,m,ω) = rad_solution(r,ω,m,k) .* ang_solution(m)
