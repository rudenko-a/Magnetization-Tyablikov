#!/vol/tcm30/rudenko/julia-1.7.3/bin/julia
using Printf
using Plots
using LinearAlgebra
#using LoopVectorization
const μB=5.7883818060e-5 # eV/T
const nkx=50; const nky=nkx
const J = [10*0.001,0*0.001]; const D=J[1]/10
const S=3/2
const kB=8.617333262e-5
const Niter=200
#const z=[4,4]
const δ=1e-4
const latt="triangular"
const SIA="SSWT" #HP0, HP1, 

#TO-DO:

#println(Threads.nthreads())


function sumbz(a,b)::Float64
    if a+b >= 1.0
        return modf(a+b)[1] # frac part of a number
    elseif a+b < 0.0
        return modf(a+b+1)[1]
    else
        return a+b
    end
end
@printf("%10.12f\n",sumbz(0,0.999))


function loop()
N=0.0

for kx=1:nkx , ky=1:nky
for qx=1:nkx , qy=1:nky
#
    px = sumbz(kx/nkx, qx/nkx)
    py = sumbz(ky/nky, qy/nky)
    N += exp(px)/exp(-py)
#
end
end
N
@printf("%10.12f\n",N)
end

function n0(β,ε)
    n0 = (exp(β*ε) - 1)^(-1)
end

function Sz(S,ϕ)
    if !isinteger(S*2)
        println("WARNING! S is neither integer nor half-integer")
    end
    if S==1/2
        num = S
        den = 1 + 2ϕ
    elseif S==1
        num = S*(1 + 2ϕ)
        den = 1 + 3ϕ + 3ϕ^2
     #   num = S
     #   den = 1 + 2ϕ
    else
    num = (S - ϕ)*(1 + ϕ)^(2S+1) + (1 + S + ϕ)*ϕ^(2S+1)
    den = (1 + ϕ)^(2S+1) - ϕ^(2S+1)
    end
    return num/den
end

function γk(x,y)
    gamma = Float64[]
    if latt == "square"
        γ1 = 2*( cos(2π*x) + cos(2π*y) )  
        γ2 = 4*( cos(2π*x)*cos(2π*y) )
    end
    if latt == "triangular"
        γ1 = 2*(cos(2π*x) + 2*cos(2π*x/2)*cos(2π*√3y/2) )
        γ2 = 2*(cos(2π*√3y) + 2*cos(2π*3x/2)*cos(2π*√3y/2) )
    end
 #   if latt == "honeycomb"
 #       gamma = Complex{Float64}[]
 #       γ1 = 1 + 2*cos(2π*x/2)*exp(im*2π√3ky/2)
 #       γ2 = 2*(cos(2π*x) + 2*cos(2π*x/2)*cos(2π*√3y/2) )
 #   end
    push!(gamma,γ1)
    push!(gamma,γ2)
    return gamma
end

function magnetization(T,method,Bz)
    if (method != "Tyablikov") && 
       (method != "Callen") &&
       (method != "SSWT")
        error("Unknown method")
    end
#    @printf("%s \n",method)
    β=1/kB/T
    σ = S
    σ_prev = 0.0
    a = 1
    q = zeros(size(J,1))
    final = false
global dispersion = zeros(nkx,nkx)
    P0 = 0.0
    for iter = 1:Niter
        P0i = 0.0
        Q0 = 0.0
        qi = zeros(size(J,1))
 #       (1 - 2*P0/S) < 0 && (Δ=0) #avoids negative Ωk. single-line 'if' ???
  #     (1 - 2*P0/S) < 0 && (P0 = S/2) #avoids negative Ωk. single-line 'if'
        #
        if     method=="Tyablikov"
            ν = zeros(size(J,1))
        elseif method=="SSWT"
            ν = q
        elseif method=="Callen"
            a = (σ/S)^1 # or (σ/S)^3
            ν = q.*(σ/S)*a
        else
            error("Unknown method for magnon dispersion")
        end
        #
        if     SIA=="AC"
             SzSz = S*(S+1) - σ*(1 + 2*P0)
             Φ = 2*σ*(1 - (S*(S+1) - SzSz) / (2*S^2))
             ΩSIA = D*Φ
        elseif SIA=="HP0"
             ΩSIA = (2S-1)*D
        elseif SIA=="HP1"
             ΩSIA = (2S-1)*D*(1-2*P0/S)
        elseif SIA=="SSWT"
             ΩSIA = (2S-1)*D*(σ/S)^2
        else
             error("Unknown method for SIA decoupling")
        end
        J0 = γk(0,0).*J
        for ikx = 1 : nkx, iky = 1 : nky
            kx = ikx - nkx/2 ; ky = iky - nky/2
            if kx==0&&ky==00 
                continue
            end
            γ = γk(kx/nkx,ky/nky)
            Jk = γ.*J
            Ωk = (σ .+ ν)⋅(J0 - Jk) + ΩSIA     + 2*μB*Bz  #
            P0i += n0(β,Ωk)  / (nkx*nky) 
            qi += n0(β,Ωk).*γ./γk(0,0)  / (nkx*nky) 
    #        if (final) 
    #            dispersion[ikx,iky] = Ωk
    #        end
            if abs(Ωk) < 1/β/1e3
   #             println(Ωk," ",1/β/1e3," ", n0(β,Ωk)," ",P0i)
            end
        end
        #
    P0 = P0i
    Q0 = 1.0 + 2*P0
    q = qi
    method == "SSWT" ? σ = S - P0i : σ = Sz(S,P0)  #σ = S / Q0
    #
        if (abs(σ-σ_prev)<δ)&&(iter>3)
            (!final) ? final=true : break
        elseif ((σ<=0)||(σ>S))#||(P0>S*0.9)) 
            σ=0
            break 
        end 
#     @printf("%5i %12.10f %8.6f %6i \n", iter, σ,  ΩSIA, T)
    σ_prev = σ
    end
    return σ
end

@time begin

    println(γk(0,0))
    println("M=",magnetization(10,"Tyablikov",0.0))
#    println("dispersion(1,1)=",dispersion[10,10])
    
    figure = plot(xlabel="Temperature, T (K)", ylabel = "Magnetization, σ", title="M(T), square lattice")
    for method ∈ ["Tyablikov", "SSWT", "Callen"]
        println(method)
        plot!(figure,1:10:500, T -> magnetization(T,method,0.0), label=method, m=:circle,ls=:solid)
    end
    hline!(figure,[0.0,S], color=:black, ls=:dash, label=false)
    figure
    
        
 #   figure2 = plot(xlabel="Temperature, T (K)", ylabel = "Magnetization, σ", title="M(T), square lattice")
 #   for method ∈ ["Tyablikov", "SSWT", "Callen"]
 #       println(method)
 #       plot!(figure2,1:10:300, T -> magnetization(T,method,1.0), label=method, m=:circle,ls=:solid)
 #   end
 #   figure2
    

    #      loop()

end











