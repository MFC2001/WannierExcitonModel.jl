export spin_ik2vck
function spin_ik2vck(bse::AbstractBSE, q::ReducedCoordinates{<:Real}, upidx::AbstractVector{<:Integer}; calc_eleband::Bool = true)
	# 1. 目前先假定，wannier基均具有确定的自旋，故自旋算符实际上与动量无关。
	# 2. 要求输入关于一种自旋的wannier基编号，下述可直接进行波函数分量求和。

	norb = numorb(bse.TB)
	upidx = Int.(collect(upidx))
	dnidx = setdiff(1:norb, upidx)

	bandk = bse.bandk
	if calc_eleband
		bandkq = BAND(map(k -> k + q, bse.kgrid), bse.TB; vector = true)
		phase_normalize!(bandkq, :real_sum)
	else
		bandkq = bse.bandkq
	end

	vckmap = bse.vckmap
	Nvck = length(bse.vckmap)
	operator_vck = Matrix{ComplexF64}(undef, Nvck, Nvck)
	for vck2 in Base.OneTo(Nvck)
		(v2, c2, k2) = vckmap[vck2]
		ψv2 = view(bandk[k2].vectors, :, v2)
		ψc2 = view(bandkq[k2].vectors, :, c2)
		for vck1 in Base.OneTo(Nvck)
			(v1, c1, k1) = vckmap[vck1]
			ψv1 = view(bandk[k1].vectors, :, v1)
			ψc1 = view(bandkq[k1].vectors, :, c1)
			c = conj(ψc1) .* ψc2
			v = conj(ψv1) .* ψv2
			operator_vck[vck1, vck2] = (sum(c[upidx]) - sum(c[dnidx]) + sum(v[dnidx]) - sum(v[upidx])) / 2
		end
	end
	return Hermitian(operator_vck)
end



function operator_ik2vck(bse::AbstractBSE, q::ReducedCoordinates{<:Real}, operator_ckq, operator_vk; calc_eleband::Bool = true)

	bandk = bse.bandk
	if calc_eleband
		bandkq = BAND(map(k -> k + q, bse.kgrid), bse.TB; vector = true)
		phase_normalize!(bandkq, :real_sum)
	else
		bandkq = bse.bandkq
	end

	vckmap = bse.vckmap
	Nvck = length(bse.vckmap)
	operator_vck = Matrix{ComplexF64}(undef, Nvck, Nvck)
	for vck2 in Base.OneTo(Nvck)
		(v2, c2, k2) = vckmap[vck2]
		for vck1 in Base.OneTo(Nvck)
			(v1, c1, k1) = vckmap[vck1]
			# operator_vck[vck1, vck2] =
			#     _operator_expectation(view(bandkq[k1].vectors, :, c1), operator_ckq, view(bandkq[k2].vectors, :, c2)) + 
			#     _operator_expectation(view(bandk[k1].vectors, :, v1), operator_vk, view(bandk[k2].vectors, :, v2))
			# 	bandkq[k1].vectors[:, c1] ⋅ (operator_ckq * bandkq[k2].vectors[:, c2]) -
			# 	bandk[k1].vectors[:, v1] ⋅ (operator_vk * bandk[k2].vectors[:, v2])
		end
	end
	return operator_vck
end
function _operator_expectation(ψ_left::AbstractVector{<:Number}, operator::AbstractMatrix{<:Number}, ψ_right::AbstractVector{<:Number})
	return ψ_left ⋅ (operator * ψ_right)
end
