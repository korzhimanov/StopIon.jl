module StopIon
using Unitful
using SpecialFunctions: erf

export Particle, bethe_bloch, fixed_background

const mc² = Unitful.me*Unitful.c0*Unitful.c0
const e⁴_mc² = (Unitful.q^2/(4π*Unitful.ɛ0))^2/mc²
const mc³_e² = mc²*1.0u"c"*4π*Unitful.ɛ0/(Unitful.q^2)
const mc²_ħ = mc²/Unitful.ħ

@derived_dimension Concentration inv(Unitful.𝐋^3)
@derived_dimension ParticleDensityFlux inv(Unitful.𝐋^2)

const shell_correction_coefficient = 0.5

struct Particle
    z::Float64
    m::Unitful.Mass
end

function γ²(ɛ::Unitful.Energy, p::Particle)
    return (1 + ɛ/(p.m*1u"c^2"))^2
end

function β²(ɛ::Unitful.Energy, p::Particle)
    return (γ²(ɛ,p)-1)/γ²(ɛ,p)
end

function MIP(z) :: Unitful.Energy
    return 13.6u"eV"*z
end

function EffectiveMIP(Z, Z_eff) :: Unitful.Energy
    return 13.6u"eV"*Z^2/(Z-Z_eff)
end

function z_p(ɛ::Unitful.Energy, p::Particle)
    return p.z*(1 - 1.034*exp(-137.04*sqrt(β²(ɛ, p))*p.z^(-0.688)))
end

function G(ɛ::Unitful.Energy, p::Particle, Te::Unitful.Energy)
    x = sqrt(mc²*β²(ɛ, p)*0.5/Te)
    return erf(ustrip(uconvert(Unitful.NoUnits, x))) - 2*x/sqrt(π)*exp(-x*x)
end

function L(ɛ::Unitful.Energy, p::Particle, target::Particle)
    return log(2*mc²/MIP(target.z)*(γ²(ɛ,p)-1)) - β²(ɛ,p) - shell_correction_coefficient/target.z
end

function L_bound(ɛ::Unitful.Energy, p::Particle, target::Particle, Z_eff)
    return log(2*mc²/EffectiveMIP(target.z, Z_eff)*(γ²(ɛ,p)-1)) - β²(ɛ,p) - shell_correction_coefficient/(target.z - Z_eff)
end

function L_free(ɛ::Unitful.Energy, p::Particle, N_a::Concentration, Te::Unitful.Energy, Z_eff)
    ω_p = sqrt(Z_eff*N_a*Unitful.q^2/(Unitful.ɛ0*Unitful.me))
    u = sqrt(2*Te/mc² + β²(ɛ,p))
    b_min = max(
        z_p(ɛ,p)/(mc³_e²*u^2),
        0.5/(mc²_ħ*u)
    )
    #b_min = z_p(ɛ,p)/(mc³_e²*u^2)
    η = (b_min*ω_p)/sqrt(β²(ɛ,p))
    return 0.764*log((1.0+η)/η)
end

function bethe_bloch(
    projectile::Particle,
    target::Particle,
    density::Unitful.Density,
    energy_samplings::Vector{EnergyT},
    energy_distribution_function,
    x_samplings::Vector{LengthT}) where EnergyT<:Unitful.Energy where LengthT<:Unitful.Length
    
    @assert length(energy_samplings) == length(energy_distribution_function) "energy samplings and energy distribution function must be of the same length!"
    
    KNaZ = 4π*(projectile.z)^2*e⁴_mc²*density/target.m*target.z
    
    x_diff = diff(x_samplings)
    deposition = fill(0.0u"MeV", length(x_diff))
    
    for (ɛ0, f) in zip(energy_samplings, energy_distribution_function)
        ɛ = ɛ0
        for (i, (_, Δx)) in enumerate(zip(deposition, x_diff))
            if (ɛ < 0u"MeV")
                break
            end
            if (L(ɛ, projectile, target) < 0)
                break
            end
            Δɛ = KNaZ/β²(ɛ, projectile)*L(ɛ, projectile, target)*Δx
            deposition[i] += Δɛ*f
            ɛ -= Δɛ
        end
    end
    
    return deposition
end

function fixed_background(
    projectile::Particle,
    target::Particle,
    density::Unitful.Density,
    electron_temperature::Unitful.Energy,
    effective_charge,
    energy_samplings::Vector{EnergyT},
    energy_distribution_function,
    x_samplings::Vector{LengthT}) where EnergyT<:Unitful.Energy where LengthT<:Unitful.Length
    
    @assert length(energy_samplings) == length(energy_distribution_function) "energy samplings and energy distribution function must be of the same length!"
    @assert effective_charge ≤ target.z "effective charge must not be greater than target charge!"
    
    KNa = 4π*e⁴_mc²*density/target.m
    KNaZ = KNa*(target.z-effective_charge)
    KNaZ_eff = KNa*effective_charge
    
    x_diff = diff(x_samplings)
    deposition = fill(0.0u"MeV", length(x_diff))
    
    for (ɛ0, f) in zip(energy_samplings, energy_distribution_function)
        ɛ = ɛ0
        for (i, (_, Δx)) in enumerate(zip(deposition, x_diff))
            if (ɛ < 0u"MeV")
                break
            end
            Δɛ_free = KNaZ_eff*z_p(ɛ,projectile)^2/β²(ɛ, projectile)*G(ɛ, projectile, electron_temperature)*L_free(ɛ, projectile, density/target.m, electron_temperature, effective_charge)*Δx
            Δɛ_bound = KNaZ*z_p(ɛ,projectile)^2/β²(ɛ, projectile)*L_bound(ɛ, projectile, target, effective_charge)*Δx
            Δɛ = Δɛ_free+Δɛ_bound
            deposition[i] += Δɛ*f
            ɛ -= Δɛ
        end
    end
    
    return deposition
end

function self_consistent(
    projectile::Particle,
    particle_density_flux::ParticleDensityFlux,
    target::Particle,
    density::Unitful.Density,
    energy_samplings::Vector{EnergyT},
    energy_distribution_function,
    x_samplings::Vector{LengthT}) where EnergyT<:Unitful.Energy where LengthT<:Unitful.Length
    
    @assert length(energy_samplings) == length(energy_distribution_function) "energy samplings and energy distribution function must be of the same length!"
    
    Na = density/target.m
    KNa = 4π*e⁴_mc²*Na
    
    x_diff = diff(x_samplings)
    deposition = fill(0.0u"MJ/g", length(x_diff))
    T_e = fill(0.01u"eV", length(x_diff))
    Z_eff = fill(0.01, length(x_diff))
    
    for (ɛ0, f) in zip(reverse(energy_samplings), reverse(energy_distribution_function))
        ɛ = ɛ0
        for (i, Δx) in enumerate(x_diff)
            if (ɛ < 0u"MeV")
                break
            end
            
            KNaZ = KNa*(target.z-Z_eff[i])
            KNaZ_eff = KNa*Z_eff[i]
            
            Δɛ = 0u"eV"
            
            Δɛ_free = KNaZ_eff*z_p(ɛ,projectile)^2/β²(ɛ, projectile)*G(ɛ, projectile, T_e[i])*L_free(ɛ, projectile, density/target.m, T_e[i], Z_eff[i])*Δx
            if Δɛ_free > 0u"MeV"
                T_e[i] += f*particle_density_flux*Δɛ_free/(Na*Δx*Z_eff[i])/1.5
                Δɛ += Δɛ_free
            end
            
            Z_eff_old = Z_eff[i]
            Δɛ_bound = KNaZ*z_p(ɛ,projectile)^2/β²(ɛ, projectile)*L_bound(ɛ, projectile, target, Z_eff[i])*Δx
            if (Δɛ_bound > 0u"MeV") && (Z_eff[i] < target.z)
                ionization_energy = f*particle_density_flux*Δɛ_bound/(Na*Δx)
                while ionization_energy > 0u"eV"
                    ionization_energy_limit = (floor(Z_eff[i]+1)-Z_eff[i])*MIP(floor(target.z-Z_eff[i]))*(target.z/(target.z-Z_eff[i]))^2
                    if ionization_energy > ionization_energy_limit
                        Z_eff[i] = floor(Z_eff[i]+1)
                        ionization_energy -= ionization_energy_limit
                    else
                        Z_eff[i] += ionization_energy/ionization_energy_limit * (floor(Z_eff[i]+1)-Z_eff[i])
                        break
                    end
                end
                
                T_e[i] *= Z_eff_old/Z_eff[i]
                
                Δɛ += Δɛ_bound
            end
            
            deposition[i] += Δɛ*f*particle_density_flux/(density*Δx)
            ɛ -= Δɛ
        end
    end
    
    return deposition, T_e, Z_eff
end

end
