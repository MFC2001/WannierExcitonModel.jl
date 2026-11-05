function (bse::BSEspinblind_wiTDA)(::Val{:elegap}, q::AbstractVector{<:Real}, wave::AbstractVector{<:Number})
    kgrid = bse.kgrid
    bandk = bse.bandk
    q = ReducedCoordinates(q)
    bandkq = BAND(map(k -> k + q, kgrid), bse.eleH)

    Nvck = length(bse.vck)
    Eg = sum(Base.OneTo(Nvck)) do ivck
        (v, c, ik) = bse.vck[ivck]
        abs2(wave[i]) * (bandkq[c, ik] - bandk[ik].values[v] + bse.scissor)
    end
    return Eg
end
