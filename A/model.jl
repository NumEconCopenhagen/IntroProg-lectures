module Model

using Random

# model solution
c_star(w::Real, e::Real, θ::NamedTuple) = θ.γ*(1 - θ.τ)*w + θ.γ*e

l_star(w::Real, e::Real, θ::NamedTuple) = 1 - θ.γ + (1 - θ.γ)/(1 - θ.τ)*(e/w)

# objective function
function obj_fun(θ::NamedTuple,
                 w::AbstractVector,
                 mom_data::Tuple{Float64,Float64,Float64},
                 moments_fun::Function,
                 rng::AbstractRNG=MersenneTwister(893245) # Allow possing in an RNG to use when drawing random numbers
                 )

    S = 100 # number of simulation draws for each observed individual
    n = length(w)    # number of individuals
    w = repeat(w, S) # stack observed wages S times
    e = θ.σ*randn(rng, n*S)

    con = c_star.(w, e, Ref(θ))
    lab = l_star.(w, e, Ref(θ))

    # Calculate moments based on simulated data for the value of θ
    mom_sim = moments_fun(w, con, lab)

    # Calculate objective function as squared difference
    Q = sum(abs2(data - sim) for (data, sim) in zip(mom_data, mom_sim))

    return Q
end

end
