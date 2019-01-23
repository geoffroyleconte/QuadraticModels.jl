using QPSReader

function QuadraticModel(qps::QPSData; kwargs...)
    QuadraticModel(qps.c, qps.Q, opHermitian(qps.Q), qps.A,
                   qps.lcon, qps.ucon, qps.lvar, qps.uvar, c0=qps.c0; kwargs...)
end

QuadraticModel(filename::AbstractString) = QuadraticModel(readqps(filename))
