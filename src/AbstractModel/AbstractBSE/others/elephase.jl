
#From Rohlfing and Louie, we need to control the arbitrary phase of electron wavefunction.

function _sum_wave_is_real!(band::Eigen)
	_sum_wave_is_real!(band.vectors)
end
function _sum_wave_is_real!(vectors::AbstractMatrix{<:Complex})
	for i in axes(vectors, 2)
		_sum_wave_is_real!(view(vectors, :, i))
	end
	return vectors
end
function _sum_wave_is_real!(vector::AbstractVector{<:Complex})
	t = sum(vector)
	if abs(imag(t)) > 1e-10 #Here this judge make sure abs(t) is't zero.
		e_neg_iθ = conj(t) / abs(t)
		vector .*= e_neg_iθ
	end
	return vector
end
function _sum_wave_is_real!(vectors::AbstractMatrix{<:Real})
	return vectors
end
function _sum_wave_is_real!(vector::AbstractVector{<:Real})
	return vector
end

function _sum_wave_is_real(vector::VT)::VT where {T <: Number, VT <: AbstractVector{T}}
	t = sum(vector)
	if abs(imag(t)) > 1e-10 #Here this judge make sure abs(t) is't zero.
		e_neg_iθ = conj(t) / abs(t)
		return vector .* e_neg_iθ
	else
		return vector
	end
end
