"""
	Quantization with Entropic Regularization: with explicit, manual derivative
	created: 2023
	author©: Alois Pichler
"""

using LinearAlgebra, Distributions#, ForwardDiff #, Revise
using Gnuplot; Gnuplot.options.gpviewer= true	# use external viewer
include("ontoSimplex.jl")
include("SoftMinimum.jl")


function myObjective(y, ξ, pSupp; rWasserstein=2, λ= 0.0)
	pNorm= 2
	(d,m)= size(y)
	distVector= Vector{Float64}(undef, m)
	for j= 1:m
		distVector[j]= norm(y[:,j]- ξ, pNorm)
	end
	sMin= SoftMin(distVector.^rWasserstein, pSupp, λ= λ)	#	softmin function
	qGrad= Array{Float64}(undef, size(y))
	for j= 1:m
		for i= 1:d
			qGrad[i,j]= begin
						pSupp[j] * sMin.σ[j] *
								rWasserstein* distVector[j]^ (rWasserstein- pNorm) *
								abs(y[i,j]- ξ[i])^ (pNorm- 1) *
								sign(y[i,j]- ξ[i])
	end; end; end
	return (distance= sMin.sMin, σ= sMin.σ, qGrad= qGrad, σ0= sMin.σ0)
end

a(k)= 1/ (10+ k)^0.65    #   learning rate


#	╭────────────────────────────────────────────────────────────────
#	│	testscript
printstyled("\n\t ═════ ", basename(@__FILE__), " ═════\n"; bold= true, color= 3)

#𝒟 = MvNormal([0.0, 0.0], [1.0 -0.0;-0.0 1.0])
#𝒟 = MvNormal([0.0, 0.0], [3.0 1.0;1.0 3.0])
#𝒟 = Product(fill(Uniform(0, 1), 2))

m= 16  # Number of quantizers

d= 2   # Dimensions

qSupp= rand(𝒟, m)  # starting quantizers

pSupp= ones(m); pSupp/= sum(pSupp)  # starting probabilities 

λ= 0.1; dWasserstein= 0.0
sampleSize= 1000000
plotTime= -Inf
for k= 1:sampleSize
    ξNew= rand(𝒟 )
	sMin= myObjective_1d(qSupp, ξNew, pSupp/sum(pSupp); rWasserstein= 2, λ= λ)
	global dWasserstein+= sMin.distance # add to Wasserstein distance
	global qSupp-= a(k)* sMin.qGrad		# modify all locations
	global pSupp= pSupp*k/(k+1) + pSupp.* sMin.σ/ (k+1)	# modify all probabilities
	ontoSimplexOrthogonal!(pSupp) #ontoSimplex!(pSupp)
	if plotTime+ 1 ≤ time() || k == sampleSize     # update every  second
	global plotTime= time()
	@gp "reset"
	@gp "set title 'dWasserstein= $(round(dWasserstein/ k, digits=5)), $k samples'" "set key outside right top"
	@gp:-  qSupp[1,:] qSupp[2,:] max.(.5/m, 5*sqrt.(pSupp)/maximum(sqrt.(pSupp))) "using 1:2:3 with points pt 6 ps var lc 11 title 'quantizers'"
	end
end
@show pSupp # optimal probabilities
@show dWasserstein/ sampleSize
@show qSupp # optimal quantizers
