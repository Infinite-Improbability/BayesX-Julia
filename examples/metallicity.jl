using BayesJ

function convert_to_anders(elements::Vector{String}, mass_fractions::Vector{T})::Vector{T} where {T<:AbstractFloat}
    """
    Convert mass fractions abundances to values relative to Anders and Grevesse (1989) solar values
    """

    # Validate inputs
    if !("H" in elements)
        throw(DomainError("Hydrogen not found in elements"))
    end

    # load in the Anders values
    # Convert mass fractions abundances to values relative to Anders and Grevesse (1989) solar values
    element_names = [
        "H", "He", "C", "N", "O", "Ne", "Na", "Mg", "Al", "Si", "S", "Ar", "Ca", "Fe", "Ni"
    ]

    # Number fraction w.r.t hydrogen using A&G values
    anders_dict = Dict(i[1] => i[2] for i in zip(element_names, BayesJ.ander_Ni_per_NH))
    element_mass = Dict(i[1] => i[2] for i in zip(element_names, BayesJ.nucleon_total))

    mf_dict::Dict{String,T} = Dict(i[1] => i[2] for i in zip(elements, mass_fractions))

    if "Other" in elements
        delete!(mf_dict, "Other")
    end

    # Calculate mass fractions for A&G
    anders_fi_per_fH = Dict{String,T}()
    for el in anders_dict
        name::String, Ni_per_NH::T = el
        anders_fi_per_fH[name] = Ni_per_NH * element_mass[name]
    end
    anders_fi = Dict{String,T}()
    anders_sum = 1 / sum(i.second for i in anders_fi_per_fH)
    for el in anders_fi_per_fH
        name::String, fi_per_fH::T = el
        anders_fi[name] = fi_per_fH / anders_sum
    end

    # For elements missing from input allocate the leftover mass
    # Scaled by the Anders mass proportional to all missing elements
    # This fails to capture correlations between elements
    # but I don't have anything better right now
    leftover_mass_fraction = 1.0 - sum(i.second for i in mf_dict)
    missing_elements = Vector{String}()
    for el in anders_dict
        name::String, _ = el
        if !(name in elements)
            push!(missing_elements, name)
        end
    end
    # Fraction of anders mass composed by missing elements
    missing_fraction_anders = sum(anders_fi[i] for i in missing_elements)
    for i in missing_elements
        mf_dict[i] = leftover_mass_fraction * anders_fi[i] / missing_fraction_anders
    end

    @assert length(mf_dict) == length(anders_dict)
    @assert sum(i.second for i in mf_dict) ≈ 1.0

    # Calculate relative abundances
    relative_dict = Dict{String,T}()
    for el in mf_dict
        name::String, fi::T = el
        # Mi = fi * MT
        # Ni = Mi / Ai = fi * MT / Ai

        # Ni / NH = (fi * MT / Ai) / (fH * MT / AH)
        #         = (fi * MT * AH) / (fH * MT * Ai)
        #         = (fi * AH) / (fH * Ai)
        # Note that AH = 1
        # Ni / NH = fi / fH / Ai

        # Mass fraction / mass fraction of hydrogen
        fi_per_fH = fi / mf_dict["H"]

        # Number fraction w.r.t. hydrogen
        Ni_per_NH = fi_per_fH / element_mass[name]

        # Number fractions w.r.t hydrogen relative to Anders and Grevesse valyes
        Ni_per_NH_rel_anders = Ni_per_NH / anders_dict[name]
        relative_dict[name] = Ni_per_NH_rel_anders
    end

    return [relative_dict[i] for i in element_names]
end
