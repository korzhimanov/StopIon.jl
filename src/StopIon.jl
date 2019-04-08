module StopIon
using Unitful

export Particle, bethe_bloch

const mc² = Unitful.me*Unitful.c0*Unitful.c0

struct Particle
    z::Float64
    m::Unitful.Mass
end

function MIP(z) :: Unitful.Energy
    return 11.4u"eV"*z
end

function γ²(ɛ::Unitful.Energy, p::Particle)
    return (1 + ɛ/(p.m*1u"c^2"))^2
end

function β²(ɛ::Unitful.Energy, p::Particle)
    return (γ²(ɛ,p)-1)/γ²(ɛ,p)
end

const shell_correction_coefficient = 0.5

function L(ɛ::Unitful.Energy, p::Particle, target::Particle)
    return log(2*mc²/MIP(target.z)*(γ²(ɛ,p)-1)) - β²(ɛ,p) - shell_correction_coefficient/target.z
end

const e⁴_mc² = (Unitful.q^2/(4π*Unitful.ɛ0))^2/mc²

function bethe_bloch(
    projectile::Particle,
    target::Particle,
    density::Unitful.Density,
    energy_samplings::Vector{EnergyT},
    energy_distribution_function,
    x_samplings::Vector{LengthT}) where EnergyT<:Unitful.Energy where LengthT<:Unitful.Length
    
    @assert length(energy_samplings) == length(energy_distribution_function) "energy samplings and energy distribution function must be of the same length!"
    
    KNa = 4π*(projectile.z)^2*e⁴_mc²*density/target.m
    
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
            Δɛ = KNa*target.z/β²(ɛ, projectile)*L(ɛ, projectile, target)*Δx
            deposition[i] += Δɛ*f
            ɛ -= Δɛ
        end
    end
    
    return deposition
end

end
