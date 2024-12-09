module Born_Calderon_3d_Radial
export EigDtN, eigDtN, FqB, Fourier, Integral
using ArbNumerics, SpecialFunctions, Jacobi, FFTW


######### Auxiliary functions used in EigDtN


function f1(r, a, l)
    if a == 0
        return r^(l + 1)
    elseif a > 0
        return sqrt(r) * besseli(l + 1 // 2, sqrt(a) * r)
    else
        return sqrt(r) * besselj(l + 1 // 2, sqrt(-a) * r)
    end
end


function df1(r, a, l)
    if a == 0
        return (l + 1) * (r^l)
    elseif a > 0
        return ((l + 1) / sqrt(r)) * besseli(l + 1 // 2, sqrt(a) * r) + sqrt(a * r) * besseli(l + 3 // 2, sqrt(a) * r)
    else
        return ((l + 1) / sqrt(r)) * besselj(l + 1 // 2, sqrt(-a) * r) - sqrt(-a * r) * besselj(l + 3 // 2, sqrt(-a) * r)
    end
end


function f2(r, a, l)
    if a == 0
        return r^(-l)
    elseif a > 0
        return sqrt(r) * besselk(l + 1 // 2, sqrt(a) * r)
    else
        return sqrt(r) * bessely(l + 1 // 2, sqrt(-a) * r)
    end
end


function df2(r, a, l)
    if a == 0
        return (-l) * (r^(-l - 1))
    elseif a > 0
        return ((l + 1) / sqrt(r)) * besselk(l + 1 // 2, sqrt(a) * r) - sqrt(a * r) * besselk(l + 3 // 2, sqrt(a) * r)
    else
        return ((l + 1) / sqrt(r)) * bessely(l + 1 // 2, sqrt(-a) * r) - sqrt(-a * r) * bessely(l + 3 // 2, sqrt(-a) * r)
    end
end


#########


"""
    EigDtN(f, l, N)

Compute ``\\lambda_l[f,0]`` by splitting the interval [0,1] into N equally sized subintervals. 
"""
function EigDtN(f, l, N)

    l = ArbReal(l)

    R = ArbReal.([i // N for i in 0:N])

    F = f.(R)

    λ = df1(R[2], F[1], l) / f1(R[2], F[1], l)

    for i in 2:N
        x = df2(R[i], F[i], l) - (f2(R[i], F[i], l) * λ)
        y = -df1(R[i], F[i], l) + (f1(R[i], F[i], l) * λ)
        λ = ((x * df1(R[i+1], F[i], l)) + (y * df2(R[i+1], F[i], l))) / ((x * f1(R[i+1], F[i], l)) + (y * f2(R[i+1], F[i], l)))
    end

    return λ - 1
end


"""
    eigDtN(κ, l)

Compute ``\\lambda_l[0,\\kappa]``. 
"""
function eigDtN(κ, l)

    l = ArbReal(l)

    if κ == 0
        return l
    elseif κ > 0
        return l - (sqrt(κ) * besselj(l + 3 // 2, sqrt(κ)) / besselj(l + 1 // 2, sqrt(κ)))
    else
        return l + (sqrt(-κ) * besseli(l + 3 // 2, sqrt(-κ)) / besseli(l + 1 // 2, sqrt(-κ)))
    end
end


"""
    FqB(κ, Eig, ξ)

Compute ``\\mathcal{F}q_{\\kappa}^\\mathrm{B}(\\xi)`` with as many terms as eigenvalues in Eig.
"""
function FqB(κ, Eig, ξ)

    L = ArbReal.(0:(length(Eig)-1))

    if κ == 0
        dependence_l = ((Eig - eigDtN.(κ, L)) .* ((-1) .^ L)) ./ (gamma.(L .+ 3 // 2) .* factorial.(L))
    elseif κ > 0
        dependence_l = (2 * gamma(ArbReal(3 // 2)) * (Eig - eigDtN.(κ, L)) .* (2 * L .+ 1) .* (besselj.(L .+ 1 // 2, sqrt(κ)) .^ 2) / sqrt(κ))
    else
        dependence_l = (2 * gamma(ArbReal(3 // 2)) * ((-1) .^ L) .* (Eig - eigDtN.(κ, L)) .* (2 * L .+ 1) .* (besseli.(L .+ 1 // 2, sqrt(-κ)) .^ 2) / sqrt(-κ))
    end

    if κ == 0
        dependence_ξ_l = (ξ / 2) .^ (2 * L)
    else
        dependence_ξ_l = legendre.(1 - (ξ^2) / (2 * κ), L)
    end

    return (2 * ArbReal(pi)^(3 // 2)) * sum(dependence_l .* dependence_ξ_l)
end


"""
    Fourier(g, N, dr)

Compute the Fourier transform of a radial function g by sampling N points spaced by dr.
Return a vector of points in Fourier space and the Fourier transform at such points. 
"""
function Fourier(g, N, dr)

    R = (1:N) * dr

    aux = Float64.(g.(R) .* R)
    GR = [-reverse(aux); 0; aux]

    FGΞ = -(2 * pi * dr) * imag(fftshift(fft(ifftshift(GR))))[N+2:end]
    Ξ = (1:N) * ((2 * pi) / ((2 * N + 1) * dr))
    FG = FGΞ ./ Ξ

    return Ξ, FG
end

"""
    Integral(R, F)

Compute the integral of a radial function sampled at R with values F. 
"""
function Integral(R, F)

    return (4 * pi / length(R)) * sum((R .^ 2) .* F)
end


end # module