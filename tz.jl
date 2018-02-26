using Interpolations
using DifferentialEquations
#import GR
#GR.opengks()
#GR.openws(10, "", 210)
#GR.activatews(10)
import Plots

# Needed Packages
if false
    Pkg.add("Interpolations")
    Pkg.add("DifferentialEquations")
    Pkg.add("Plots")
    Pkg.add("GR")
end


# Versions check
if false
    println("computing versions...")
    versions = Pkg.installed()
    
    assert(VERSION == v"0.6.2")
    assert(versions["Interpolations"] == v"0.7.3")
    assert(versions["DifferentialEquations"] == v"4.1.0")
    assert(versions["Plots"] == v"0.15.0")
    assert(versions["GR"] == v"0.26.0")
end

#=

Cells and boundaries indexing for NZ=2, NT=3

  Z  
  |
  +--07--+--08--+--09--+
  |      |      |      |
 14 (4) 15 (5) 16 (6) 17
  |      |      |      |
  +--04--+--05--+--06--+
  |      |      |      |
 10 (1) 11 (2) 12 (3) 13
  |      |      |      |
  +--01--+--02--+--03--+-- T

iexm[6] == 6
iexp[6] == 9
ihym[6] == 16
ihyp[6] == 17

NN == 17 # nombre de variables
NE == 8  # nombre d'équations (1) (2) (3) (4) (5) (6) et (1-2-5-4) (2-3-6-5)
Nfree == 9 # conditions aux limites 09 08 07 14 10 01 02 03 + 04

=#

"""
ZT problem structure

nz : number of cells in z direction
nt : number of cells in t direction
nn : number of cell-borders
ne : number of equations
nf : number of free (ie boundary conditions) variables (nn == ne + nf)
dz : cell height (z direction)
dt : cell width (t direction)
# ibloc : nz×nt cells indices
iexm : nz×nt south-cell-borders indexing iz,it -> icb
iexp : nz×nt north-cell-borders indexing iz,it -> icb
ihym : nz×nt west-cell-borders indexing  iz,it -> icb
ihyp : nz×nt east-cell-borders indexing  iz,it -> icb
ivv : vector of variable-cell-borders indices 1:ne -> icb
ivc : vector of constant-cell-borders indices 1:nf -> icb
ivvi : variable-cell-borders indexing iz,it -> icb
ivci : constant-cell-borders indexing iz,it -> icb
eqmatrix (internal) : ne×nn matrix eqmatrix * cell_borders_vector = 0
eqmatrix is splitted in A and B matrices :
    eqmatrix * cell_borders_vector =
        A * variable_cell_borders_vector - B * constant_cell_border_vector,
        giving A * vv = B * vc
A : ne×ne matrix
B : ne×nf matrix
"""
struct ZT
    nz::Int64
    nt::Int64
    nn::Int64
    ne::Int64
    nf::Int64
    dz::Float64
    dt::Float64
    # main rectangular cell (ibloc is unused in fact)
    # ibloc::Array{Int64,2}
    # indices of Ex and Hy nodes, minus ou plus relative to main rectangular cell
    iexm::Array{Int64,2}
    iexp::Array{Int64,2}
    ihym::Array{Int64,2}
    ihyp::Array{Int64,2}
    ivv::Array{Int64,1}
    ivc::Array{Int64,1}
    ivvi::SparseVector{Int64,Int64} # not so sparse!
    ivci::SparseVector{Int64,Int64} # really sparse
    A::SparseMatrixCSC{Float64,Int64}
    B::SparseMatrixCSC{Float64,Int64}
    # eqmatrix::SparseMatrixCSC{Float64,Int64} # just for check

    function ZT(ZMIN::Float64, ZMAX::Float64, TMIN::Float64, TMAX::Float64, NZ::Int64, NT::Int64, boundaries::Symbol)
        DZ = (ZMAX-ZMIN)/NZ
        DT = (TMAX-TMIN)/NT
        ibloc = Matrix{Int64}(NZ,NT)
        iexm = Matrix{Int64}(NZ,NT)
        iexp = Matrix{Int64}(NZ,NT)
        ihym = Matrix{Int64}(NZ,NT)
        ihyp = Matrix{Int64}(NZ,NT)
        
        for it in 1:NT
            for iz in 1:NZ
                ibloc[iz, it] = NT*(iz-1) + it
                iexm[iz, it] = NT*(iz-1) + it
                iexp[iz, it] = NT*(iz) + it
            end
        end
        assert(maximum(iexm) == NZ*NT)
        assert(maximum(iexp) == (NZ+1)*NT)
        for it in 1:NT
            for iz in 1:NZ
                ihym[iz, it] = (NZ+1)*NT + (NT+1)*(iz-1) + it
                ihyp[iz, it] = (NZ+1)*NT + (NT+1)*(iz-1) + it+1
            end
        end
        NN = NT*(NZ+1) + (NT+1)*NZ
        assert(maximum(ihyp) == NN)
        
        NE = NZ*NT + (NZ-1)*(NT-1)
        NFree = NN - NE
        assert(NFree == 2*NT + 2*NZ - 1)
        
        vI = Vector{Int64}(4*NE)
        vJ = Vector{Int64}(4*NE)
        vV = Vector{Float64}(4*NE)
        
        ieq = 0
        ii = 0
        for it in 1:NT
            for iz in 1:NZ
                ieq += 1
                ii += 1
                vI[ii] = ieq
                vJ[ii] = iexp[iz, it]
                vV[ii] = DT
                ii+=1
                vI[ii] = ieq
                vJ[ii] = iexm[iz, it]
                vV[ii] = -DT
                ii+=1
                vI[ii] = ieq
                vJ[ii] = ihyp[iz, it]
                vV[ii] = DZ
                ii+=1
                vI[ii] = ieq
                vJ[ii] = ihym[iz, it]
                vV[ii] = -DZ
            end
        end
        
        for it in 1:NT-1
            for iz in 1:NZ-1
                ieq += 1
                ii += 1
                vI[ii] = ieq
                vJ[ii] = ihyp[iz+1, it]
                vV[ii] = DT
                ii += 1
                vI[ii] = ieq
                vJ[ii] = ihyp[iz, it]
                vV[ii] = -DT
                ii += 1
                vI[ii] = ieq
                vJ[ii] = iexp[iz, it+1]
                vV[ii] = DZ
                ii += 1
                vI[ii] = ieq
                vJ[ii] = iexp[iz, it]
                vV[ii] = -DZ
            end
        end
        
        eqmatrix = sparse(vI, vJ, vV)
        
        ivv = Vector{Int64}(NE)
        ivc = Vector{Int64}(NFree)

        if boundaries == :FDTD
            ii = 0
            # bottom boundary
            for it in 1:NT
                ii += 1
                ivc[ii] = iexm[1, it]
            end
            # left 2nd layer
            for iz in 1:NZ-1
                ii += 1
                ivc[ii] = iexp[iz, 1]
            end
            # top boundary
            for it in 1:NT
                ii += 1
                ivc[ii] = iexp[NZ, it]
            end
            # left 1st layer
            for iz in 1:NZ
                ii += 1
                ivc[ii] = ihym[iz, 1]
            end
            ii = 0
            for it in 2:NT
                for iz in 1:NZ-1
                    ii += 1
                    ivv[ii] = iexp[iz, it]
                end
            end
            for it in 1:NT
                for iz in 1:NZ
                    ii += 1
                    ivv[ii] = ihyp[iz, it]
                end
            end
        elseif boundaries == :QUASICLOSED
            ii = 0
            # bottom boundary
            for it in 1:NT
                ii += 1
                ivc[ii] = iexm[1, it]
            end
            # top boundary
            for it in 1:NT
                ii += 1
                ivc[ii] = iexp[NZ, it]
            end
            # left boundary
            for iz in 1:NZ
                ii += 1
                ivc[ii] = ihym[iz, 1]
            end
            # right boundary
            for iz in 1:NZ-1
                ii += 1
                ivc[ii] = ihyp[iz, NT]
            end
            
            ii = 0
            for it in 1:NT
                for iz in 1:NZ-1
                    ii += 1
                    ivv[ii] = iexp[iz, it]
                end
            end
            for it in 1:NT-1
                for iz in 1:NZ
                    ii += 1
                    ivv[ii] = ihyp[iz, it]
                end
            end
            ii += 1
            ivv[ii] = ihyp[NZ, NT]   
        end
        
        ivci = sparsevec(ivc, 1:length(ivc), NN)
        ivvi = sparsevec(ivv, 1:length(ivv), NN)
        
        A = eqmatrix[:,ivv]
        B = -eqmatrix[:,ivc]

        new(NZ, NT, NN, NE, NFree, DZ, DT, iexm, iexp, ihym, ihyp, ivv, ivc, ivvi, ivci, A, B)
    end
end


"""
Two contrapropagative waves vacuum spacetime filling, with wavelength = 1
    zmin: minimal z coordinate
    zmax: maximal z coordinate
    tmin: minimal t coordinate
    tmax: maximal t coordinate
    nz: number of cells in z direction
    nt: number of cells in t direction

Returns a TZ structure, and electric and magnetic fields on center of cell boundaries
"""
function make_contra(zmin, zmax, tmin, tmax, nz, nt, boundaries)
    
    zt = ZT(zmin, zmax, tmin, tmax, nz, nt, boundaries)

    kzZ = 2*pi*zt.dz # wavelength = 1 -> kz = 2 pi
    # dispersion given by (cos(kz*Z) - 1)*DT^2 + (cos(kt*T) -1)*DZ^2 = 0
    ktT = acos(1 - (1 - cos(kzZ))*(zt.dt/zt.dz)^2)

    exa = Array{Float64,2}(zt.nz+1, zt.nt)
    hya = Array{Float64,2}(zt.nz, zt.nt+1)

    # contrapropagative waves spacetime filling, Ex boundaries
    for it in 1:zt.nt
        for iz in 1:1:zt.nz+1
            kzz = kzZ*(iz-1-zt.nz÷2)
            ktt = ktT*(it-0.5)
            # maximum field is 2*pi, maximum flux is 1
            exa[iz, it] = (cos(kzz - ktt) + cos(-kzz - ktt))*pi
        end
    end
    # contrapropagative waves spacetime filling, Hy boundaries
    for it in 1:zt.nt+1
        for iz in 1:1:zt.nz
            kzz = kzZ*(iz-0.5-zt.nz÷2)
            ktt = ktT*(it-1)
            hya[iz, it] = (cos(kzz - ktt) - cos(-kzz - ktt))*pi
        end
    end

    return (zt, exa, hya)
end
    
"""
Initialize a two contrapropagative waves structure and fields in vacuum,
solve the TZ problem with defined field on one layer bottom boundary,
one layer top boundary and two-layer left boundary,
and compare computed and solved fields.
This is an FDTD situation, where DT<DZ condition is required
"""
function check_fdtd_contra(zmin, zmax, tmin, tmax, nz, nt)
    (zt, exa, hya) = make_contra(zmin, zmax, tmin, tmax, nz, nt, :FDTD)

    # boundary conditions vector
    vc = zeros(Float64, zt.nf)
    # field matrices; the boundary components are pushed yet for further comparison
    exac = Array{Float64,2}(zt.nz+1, zt.nt)
    hyac = Array{Float64,2}(zt.nz, zt.nt+1)
    for it in 1:zt.nt
        vc[zt.ivci[zt.iexm[1, it]]] = exa[1, it]           # bottom (zmin) boundary
        exac[1, it]                 = exa[1, it]
        vc[zt.ivci[zt.iexp[zt.nz, it]]] = exa[zt.nz+1, it] # top (zmax) boundary
        exac[zt.nz+1, it]               = exa[zt.nz+1, it]
    end
    for iz in 1:zt.nz
        vc[zt.ivci[zt.ihym[iz, 1]]] = hya[iz, 1]           # left (tmin) boundary, 1st layer
        hyac[iz, 1]                 = hya[iz, 1]
    end
    for iz in 1:zt.nz-1
        vc[zt.ivci[zt.iexp[iz, 1]]] = exa[iz+1, 1]         # left (tmin) boundary, 2nd layer
        exac[iz+1, 1]               = exa[iz+1, 1]
    end
    
    b = zt.B*vc
    vv = zt.A\b

    # Ex matrix filling with solved terms 
    for it in 2:zt.nt
        for iz in 1:zt.nz-1
            exac[iz+1, it] = vv[zt.ivvi[zt.iexp[iz, it]]]
        end
    end
    
    # Hy matrix filling with solved terms 
    for it in 1:zt.nt
        for iz in 1:zt.nz
            hyac[iz, it+1] = vv[zt.ivvi[zt.ihyp[iz, it]]]
        end
    end

    # comparison between computed and solved fields
    return(max(maximum(abs.(hya - hyac)),maximum(abs.(exa - exac))))
end


"""
Initialize a two contrapropagative waves structure and fields in vacuum,
solve the TZ problem with defined field on one layer bottom boundary,
one layer top boundary
one layer left boundary
one leyer right bondary (with one open segment)
and compare computed and solved fields.
"""
function check_quasiclosed_contra(zmin, zmax, tmin, tmax, nz, nt)
    (zt, exa, hya) = make_contra(zmin, zmax, tmin, tmax, nz, nt, :QUASICLOSED)

    # boundary conditions vector
    vc = zeros(Float64, zt.nf)
    # field matrices; the boundary components are pushed yet for further comparison
    exac = Array{Float64,2}(zt.nz+1, zt.nt)
    hyac = Array{Float64,2}(zt.nz, zt.nt+1)
    for it in 1:zt.nt
        vc[zt.ivci[zt.iexm[1, it]]] = exa[1, it]           # bottom (zmin) boundary
        exac[1, it]                 = exa[1, it]
        vc[zt.ivci[zt.iexp[zt.nz, it]]] = exa[zt.nz+1, it] # top (zmax) boundary
        exac[zt.nz+1, it]               = exa[zt.nz+1, it]
    end
    for iz in 1:zt.nz
        vc[zt.ivci[zt.ihym[iz, 1]]] = hya[iz, 1]           # left (tmin) boundary
        hyac[iz, 1]                 = hya[iz, 1]
    end
    for iz in 1:zt.nz-1
        vc[zt.ivci[zt.ihyp[iz, zt.nt]]] = hya[iz, zt.nt+1] # right (tmax) boundary
        hyac[iz, zt.nt+1]               = hya[iz, zt.nt+1]
    end
    
    b = zt.B*vc
    vv = zt.A\b

    # Ex matrix filling with solved terms 
    for it in 1:zt.nt
        for iz in 1:zt.nz-1
            exac[iz+1, it] = vv[zt.ivvi[zt.iexp[iz, it]]]
        end
    end
    
    # Hy matrix filling with solved terms 
    for it in 1:zt.nt-1
        for iz in 1:zt.nz
            hyac[iz, it+1] = vv[zt.ivvi[zt.ihyp[iz, it]]]
        end
    end

    hyac[zt.nz, zt.nt+1] = vv[zt.ivvi[zt.ihyp[zt.nz, zt.nt]]]

    # comparison between computed and solved fields
    return(max(maximum(abs.(hya - hyac)),maximum(abs.(exa - exac))))
end



# checks FDTD condition
println("DT=DZ       ", check_fdtd_contra(-1.0, 1.0, 0.0, 1.0, 200, 100))
println("DT<DZ       ", check_fdtd_contra(-1.0, 1.0, 0.0, 1.0, 200, 101))
println("DT>DZ(bad)  ", check_fdtd_contra(-1.0, 1.0, 0.0, 1.0, 200, 99))

# pas de condition pour le contour clos, mais il peut y avoir des modes (champ à conditions aux limites nulles)
println("DT=DZ(mode) ", check_quasiclosed_contra(-1.0, 1.0, 0.0, 1.0, 200, 100))
println("DT<DZ       ", check_quasiclosed_contra(-1.0, 1.0, 0.0, 1.0, 200, 101))
println("DT>DZ       ", check_quasiclosed_contra(-1.0, 1.0, 0.0, 1.0, 200, 99))

(zt, exa, hya) = make_contra(-1.0, 1.0, 0.0, 1.0, 200, 100, :FDTD)

# pair of ex and hy knots pseudo-vectors z and t coordinates
exknots = (zt.dt*(0.5+(0:zt.nt-1)), -1.0 + zt.dz*(0:zt.nz))
hyknots = (zt.dt*(0:zt.nt), -1.0 + zt.dz*(0.5+(0:zt.nz-1)))
# Ex and Yy interpolators
ex = interpolate(exknots, exa.', Gridded(Linear()))
hy = interpolate(hyknots, hya.', Gridded(Linear()))
#=

GKCL opengks GKOP closegks GKCL
GKOP openws WSOP closews GKOP
WSOP activatews WSAC deactivatews WSOP
WSAC createseg SGOP closeseg WSAC

=#

tt = linspace(0,1.0,150)
Plots.plot(tt, ex[tt, 0], label="z=0", title="\$E_x(t)\$", xlabel="\$t\$", ylabel="\$E_x\$")
Plots.plot!(tt, ex[tt, 0.24], label="z=0.24")
Plots.plot!(tt, ex[tt, 0.5], label="z=0.5")
Plots.savefig("figures/plot1.png")

Plots.plot(tt, ex[tt, 0.05], label="\$E_x\$", title="\$E_x(z = 0.05, t), H_y(z = 0.05, t)\$", xlabel="\$t\$", ylabel="Field")
Plots.plot!(tt, hy[tt, 0.05], label="\$H_y\$")
Plots.savefig("figures/plot2.png")


#=
d Phie = Hy dt - Dx dz 
d Phim = -Ex dt + By dz

d Phie = 0 for (dt, dx) = (Dx, Hy)
d Phim = 0 for (dt, dx) = (By, Ex)

Phie line (with flux = 0 on Minkowski perpendicular direction), with dl Eudlidean length = dl*(Hy, Dx)/sqrt(Hy^2 + Dx^2)
Phim line (with flux = 0 on Minkowski perpendicular direction), with dl Eudlidean length = dl*(Ex, By)/sqrt(Ex^2 + By^2)

Flux on these lines

d Phie = dl*(Hy^2 - Dx^2)/sqrt(Hy^2 + Dx^2)
d Phim = dl*(-Ex^2 - By^2)/sqrt(Ex^2 + By^2)

=#

function phie_vector(t, z)
    vt = ex[t, z]
    vz = hy[t, z]
    return [vt, vz]
end

function phim_vector(t, z)
    vt = hy[t, z]
    vz = ex[t, z]
    return [vt, vz]
end

"""
Returns dt/dl, dz/dl, dphi/dl, and d(euclidean length)/dl
"""
function phie(u, sign, t)
    vt = sign*ex[u[1], u[2]]
    vz = sign*hy[u[1], u[2]]
    v = hypot(vt, vz)
    vtn = vt/v
    vzn = vz/v
    return [vtn, vzn, vzn^2 - vtn^2, v]
end

function phim(u, sign, t)
    vt = sign*hy[u[1], u[2]]
    vz = sign*ex[u[1], u[2]]
    v = hypot(vt, vz)
    vtn = vt/v
    vzn = vz/v
    return [vtn, vzn, vzn^2 - vtn^2, v]
end

function didi(u0, tspan, phif, sign)
    prob = ODEProblem(phif, u0, tspan, sign)
    sol = solve(prob, reltol=1e-8, abstol=1e-8)

    t   = Vector{Float64}(length(sol.u))
    z   = Vector{Float64}(length(sol.u))
    phi = Vector{Float64}(length(sol.u))
    s   = Vector{Float64}(length(sol.u))
    
    for (i,v) in enumerate(sol.u)
        t[i]   = v[1]
        z[i]   = v[2]
        phi[i] = v[3]
        s[i]   = v[4]
    end
    return (t, z, phi, s)
end


(t, z, phi, s) = didi([0.65, 0.0, 0, 0], (0.0,1.0), phim, 1)
Plots.scatter(t, z, aspect_ratio=:equal, title="\$\\Phi_m\$", xlabel="\$t\$", ylabel="\$z\$", label="")
Plots.plot!(t, z, label="")
(t, z, phi, s) = didi([0.6, 0.0, 0, 0], (0.0,1.0), phim, 1)
Plots.plot!(t, z, label="")
(t, z, phi, s) = didi([0.5, 0.0, 0, 0], (0.0,1.0), phim, 1)
Plots.plot!(t, z, label="")
(t, z, phi, s) = didi([0.55, 0.0, 0, 0], (0.0,1.0), phim, 1)
Plots.plot!(t, z, label="")
(t, z, phi, s) = didi([0.51, 0.0, 0, 0], (0.0,1.0), phim, 1)
Plots.plot!(t, z, label="")
(t, z, phi, s) = didi([0.501, 0.0, 0, 0], (0.0,1.0), phim, 1)
Plots.scatter!(t, z, label="")
Plots.plot!(t, z, label="")
(t, z, phi, s) = didi([0.501, 0.0, 0, 0], (0.0,1.0), phim, 1)
Plots.scatter!(t, z, label="")
Plots.plot!(t, z, label="")
Plots.savefig("figures/plot3.png")


(t, z, phi, s) = didi([0.001, 0.0, 0, 0], (0.0,2.0), phim, 1)
Plots.scatter(t, z, aspect_ratio=:equal, title="\$\\Phi_m\$", xlabel="\$t\$", ylabel="\$z\$", label="")
Plots.plot!(t, z, label="")

(t, z, phi, s) = didi([0.001, 0.5, 0, 0], (0.0,2.0), phim, 1)
Plots.scatter(t, z, aspect_ratio=:equal, title="\$\\Phi_m\$", xlabel="\$t\$", ylabel="\$z\$", label="")
Plots.plot!(t, z, label="")
(t, z, phi, s) = didi([0.001, 0.5, 0, 0], (0.0,2.0), phie, -1)
Plots.scatter!(t, z)
Plots.plot!(t, z, label="")

(t, z, phi, s) = didi([0.7, 0.2, 0, 0], (0.0,2.0), phim, 1)
Plots.plot(t, z, aspect_ratio=:equal, title="Flux lines", xlabel="\$t\$", ylabel="\$z\$", label="\$\\Phi_m\$")
(t, z, phi, s) = didi([0.7, 0.2, 0, 0], (0.0,2.0), phie, -1)
Plots.plot!(t, z, label="\$\\Phi_e\$")
Plots.savefig("figures/flux_lines.png")

Plots.plot(s, phi)
Plots.savefig("figures/plot4.png")

(t, z, phi, s) = didi([0.65, 0.0, 0, 0], (0.0,1.0), phim, 1)
Plots.plot(s, phi)
(t, z, phi, s) = didi([0.60, 0.0, 0, 0], (0.0,1.0), phim, 1)
Plots.plot!(s, phi)
(t, z, phi, s) = didi([0.55, 0.0, 0, 0], (0.0,1.0), phim, 1)
Plots.plot!(s, phi)
(t, z, phi, s) = didi([0.5001, 0.0, 0, 0], (0.0,1.0), phim, 1)
Plots.plot!(s, phi)
Plots.scatter!(s, phi)
(t, z, phi, s) = didi([0.500001, 0.0, 0, 0], (0.0,1.0), phim, 1)
Plots.plot!(s, phi)
Plots.savefig("figures/plot5.png")

