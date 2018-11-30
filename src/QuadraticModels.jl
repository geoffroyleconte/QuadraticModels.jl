module QuadraticModels

using FastClosures
using LinearOperators
using NLPModels
using Requires

using LinearAlgebra
using SparseArrays

import NLPModels:
    objgrad, objgrad!, obj,grad, grad!,
    hess_coord, hess, hess_op, hprod, hprod!,
    cons, cons!,
    jac_coord, jac, jac_op, jprod, jprod!, jtprod, jtprod!

export AbstractQuadraticModel, QuadraticModel
export objgrad, objgrad!, obj,grad, grad!,
       hess_coord, hess, hess_op, hprod, hprod!,
       cons, cons!,
       jac_coord, jac, jac_op, jprod, jprod!, jtprod, jtprod!

include("qpmodel.jl")

function __init__()
    @require QPSReader = "758ba83c-e923-11e8-036a-3f93d0cb3d0c" include("qps.jl")
end

end # module
