"""
	soft/ smooth minimum/ mellowmin
	created: 2023, March
	author©: Alois Pichler
"""

#	╭────────────────────────────────────────────────────────────────
#	│	LogSumExp minimum function with the mean
#	│	extend the min function from Base
function SoftMin(x; λ= 0.0)
	m= minimum(x); σ= Vector{Float64}(undef, size(x))
	if λ == 0
		for i in eachindex(x)
			if x[i] ≤ m; σ[i]= 1.0; else; σ[i]= 0.0; end
		end
		Eσ= sum(σ)/ length(σ)
	else	#	LogSumExp
		for (i,xi) in enumerate(x)
			σ[i]= exp((m.- xi)/ λ)
		end
		Eσ= sum(σ)/ length(σ); m-= λ* log(Eσ)
	end
	return (sMin= m, σ= σ/ Eσ, σ0= σ)
end

function SoftMin2(x, p; λ= 0.0)
	length(p)== length(x) || error("SoftMinimum.jl: x and p have different length.")
	@show m= minimum(x); σ= Vector{Float64}(undef, size(x))
	if λ == 0
		for i in eachindex(x)
			if x[i] ≤ m; σ[i]= 1.0 else σ[i]= 0.0 end
		end
		Eσ= p'* σ; Eσ > 0 || (@show x, p; error("SoftMinimum.jl: denominator ≤ 0."))
	else	#	LogSumExp
		for i in eachindex(x)
			σ[i]= exp((m.- x[i])/ λ)	# note, that 0 ≤ σ[i] ≤ 1
		end
		Eσ= p'* σ; Eσ > 0 || (@show x, p; error("SoftMinimum.jl: denominator ≤ 0."))
		m-= λ* log(Eσ)
	end
	return (sMin= m, σ= σ/ Eσ, σ0= σ)
end

function SoftMin(x, p; λ= 0.0)
	length(p)== length(x) || error("SoftMinimum.jl: x and p have different length.")
	m= Inf; σ= Vector{Float64}(undef, size(x))
	for (i,xi) in enumerate(x)	#	compute the minimum
		if xi ≤ m && p[i] ≠ 0.0; m= xi; end
	end
	if λ == 0
		for i in eachindex(x)
			if x[i] ≤ m; σ[i]= 1.0 else σ[i]= 0.0 end
		end
		Eσ= p'* σ; Eσ > 0 || error("SoftMinimum.jl: denominator ≤ 0.")
	else	#	LogSumExp
		for (i,xi) in enumerate(x)
			if p[i] > 0.0
				σ[i]= exp((m.- xi)/ λ)	# note, that 0 ≤ σ[i] ≤ 1
			else
				if xi ≤ m
					σ[i]= 1.0
				else
					σ[i]= 0.0
				end
			end
		end
		Eσ= p'* σ; Eσ > 0 || error("SoftMinimum.jl: denominator ≤ 0.")
		m-= λ* log(Eσ)
	end
	return (sMin= m, σ= σ/ Eσ, σ0= σ)
end
