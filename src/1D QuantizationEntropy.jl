"""
	Quantization with Entropic Regularization: with explicit, manual derivative
	created: 2023
	author©: Rajmadan Lakshmanan and Alois Pichler
"""

using LinearAlgebra, Distributions#, ForwardDiff #, Revise
using Gnuplot; Gnuplot.options.gpviewer= true	# use external viewer
include("ontoSimplex.jl")
include("SoftMinimum.jl")


function cdfPoints(xSupp::Vector, p::Vector, Δ= nothing)
    if isnothing(Δ); Δ= max(0.1, 0.05*(maximum(xSupp) - minimum(xSupp))); end
    n= length(xSupp); xTmp= zeros(2+ 2n); pTmp= zeros(2+ 2n)
    perm= sortperm(xSupp)
    for i= 0:n
        if i==0
            xTmp[1]= minimum(xSupp)- Δ
            xTmp[2]= xSupp[perm[1]]
            pTmp[1]= pTmp[2]= 0
        elseif i< n
            xTmp[2i+1]= xSupp[perm[i]]
            xTmp[2i+2]= xSupp[perm[i+1]]
            pTmp[2i+1]= pTmp[2i+2]= pTmp[2i]+ p[perm[i]]
        else
            xTmp[2n+1]= xSupp[perm[n]]
            xTmp[2n+2]= maximum(xSupp)+ Δ
            pTmp[2n+1]= pTmp[2n+2]= 1   # pTmp[2n]+ p[n]
        end
    end
    return xTmp, pTmp
end

function myObjective(y, ξ, pSupp; rWasserstein=2, λ= 0.0)
	pNorm= 2
	m = length(y)
	distVector= Vector{Float64}(undef, m)
	for j= 1:m
		distVector[j]= norm(y[j]- ξ, pNorm)
	end
	sMin= SoftMin(distVector.^rWasserstein, pSupp, λ= λ)	#	softmin function
	qGrad= Vector{Float64}(undef, length(y))
	for j= 1:m
		qGrad[j]= begin
						pSupp[j] * sMin.σ[j] *
								rWasserstein* distVector[j]^ (rWasserstein- pNorm) *
								abs(y[j]- ξ)^ (pNorm- 1) *
								sign(y[j]- ξ)
	end; end
	return (distance= sMin.sMin, σ= sMin.σ, qGrad= qGrad, σ0= sMin.σ0)
end

a(k)= 1/ (10+ k)^0.65    #   learning rate


#	╭────────────────────────────────────────────────────────────────
#	│	testscript
printstyled("\n\t ═════ ", basename(@__FILE__), " ═════\n"; bold= true, color= 3)

###  Problem of Interest

#𝒟 = Exponential(1.0)
#𝒟 = Normal(0,1)
𝒟 = Uniform(0.0, 1.0)
		

m= 8    # No of quantizers
qSupp= rand(𝒟, 8)  # starting quantizers

pSupp= ones(m); pSupp/= sum(pSupp) # starting probabilities 

λ= 0.1; dWasserstein= 0.0
sampleSize= 1000000
plotTime= -Inf
plotSupp= LinRange(-0.1, 1.1, 100); plotP= zeros(length(plotSupp))
for k= 1:sampleSize
    ξNew= rand(𝒟 )
	sMin= myObjective(qSupp, ξNew, pSupp/sum(pSupp); rWasserstein= 2, λ= λ)
	global plotP+= ξNew .≤ plotSupp
	global dWasserstein+= sMin.distance # add to Wasserstein distance
	global qSupp-= a(k)* sMin.qGrad		# modify all locations
	global pSupp= pSupp*k/(k+1) + pSupp.* sMin.σ/ (k+1)	# modify all probabilities
	ontoSimplexOrthogonal!(pSupp) #ontoSimplex!(pSupp)
	if plotTime + 1 ≤ time() || k == sampleSize     # update every  second
        global plotTime= time()
        (plotRange, yCdf)= cdfPoints(qSupp, pSupp)
        @gp "reset"
        @gp "set title 'dWasserstein= $(round(dWasserstein/ k, digits=4)), $k samples'" "set key outside right top"
        @gp :- plotRange yCdf "with lines notitle lc '#0918e8'" 
        plotRange= collect(LinRange(minimum(plotRange), maximum(plotRange), 100))
        #@gp :- plotRange cdf(𝒟 , plotRange) "with lines title 'distribution'"
        @gp :- plotSupp plotP/ plotP[end] "with lines notitle lc '#000000'" 
        #@gp key = "outside center bottom"
	@gp:- "set grid front"
         end
end
@show pSupp # optimal probabilities
@show dWasserstein/ sampleSize
@show qSupp # optimal quantizers

