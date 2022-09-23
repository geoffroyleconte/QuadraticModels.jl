struct FreeLinearDoubletonColumn{T, S} <: PresolveOperation{T, S}
  j::Int
  i::Int
  k::Int
  aij::T
  akj::T
  yi::T
  conival::T
end

function free_linear_doubleton_columns!(
  operations::Vector{PresolveOperation{T, S}},
  hcols::Vector{Col{T}},
  arows::Vector{Row{T}},
  acols::Vector{Col{T}},
  c::AbstractVector{T},
  c0::T,
  lcon::AbstractVector{T},
  ucon::AbstractVector{T},
  lvar::AbstractVector{T},
  uvar::AbstractVector{T},
  ly::AbstractVector{T},
  uy::AbstractVector{T},
  nvar,
  row_cnt,
  col_cnt,
  kept_rows,
  kept_cols,
  s,
) where {T, S}
  free_dsc_pass = false
  c0_offset = zero(T)
  for j = 1:nvar
    (kept_cols[j] && (col_cnt[j] == 2) && isempty(hcols[j].nzind) &&
      (lvar[j] == -T(Inf)) && (uvar[j] == T(Inf))) || continue
    # check infinity bounds and no hessian contribution
    colj = acols[j]
    i = 0
    aij = zero(T)
    k = 0
    akj = zero(T)
    c = 0
    for (l, alj) in zip(colj.nzind, colj.nzval) # get aij and akj
      kept_rows[l] || continue
      c += 1
      if c == 1
        i = l
        aij = alj
      elseif c == 2
        k = l
        akj = alj
      end
    end
    if lcon[k] == ucon[k] && lcon[i] != ucon[i] # swap i and k. i must be an equality constraint
      l = i
      i = k
      k = l
      alj = aij
      aij = akj
      akj = alj
    end
    lcon[i] == ucon[i] || continue # should not exit?
    s[k] == 0 || continue
    s[k] = i # concatenate row i to row k if fills-in

    arowk = rows[k]
    arowi = rows[i]

    # suppose aij is now removed. column j is a singleton with only element akj.
    # remove variable j from constraint
    yi = c[j] / aij
    nzcj = c[j] != zero(T)
    nzcj && (c0_offset += yi * ucon[i]) # update c0
    # remove arowi
    for (l, ail) in zip(rowi.nzind, rowi.nzval)
      # add all col elements except the col j2 == j
      (kept_cols[l] && l != j) || continue
      nzcj && (c[l] -= yi * ail) # update c if c[j] != 0
      col_cnt[l] -= 1
    end
    conival = ucon[i] # constant for postsolve

    # check aij not too small?
    c_k = 1
    c_i = 1
    # update arowk
    row_cnt_k_offset = 0
    for l in 1:nvar
      kept_cols[l] || continue
      if l == arowk.nzind[c_k]
        akl = arowk.nzval[c_k]
        c_k += 1
      else
        akl = zero(T) 
      end
      if l == arowi.nzind[c_i]
        ail = arowi.nzval[c_i]
        c_i += 1
      else
        ail = zero(T)
      end

      # akl2 = akl - akj * ail / aij
      # if ail == 0, then akl is not changed
      # if ail != 0, we store - akj * ail / aij in ail
      # then, at each akl call, we have to add ail
      if ail != zero(T)
        arowi.nzval[c_i - 1] = - akj * ail / aij # ail2 = - akj * ail / aij (use this for postsolve)
        akl == zero(T) && (row_cnt_k_offset += 1)
      end
    end
    lvar[k] -= akj * lvar[i] / aij
    uvar[k] -= akj * uvar[i] / aij
    if sign(aij) == sign(akj)
      ly[i] = (c[j] - akj * uy[k]) / aij
      uy = (c[j] - akj * ly[k]) / aij
    else
      ly[i] = (c[j] - akj * ly[k]) / aij
      uy = (c[j] - akj * uy[k]) / aij
    end
    if abs(akj) > eps(T)
      if sign(aij) == sign(akj)
        ly[i] = (c[j] - aij * uy[k]) / akj
        uy[i] = (c[j] - aij * ly[k]) / akj
      else
        ly[k] = (c[j] - aij * ly[k]) / akj
        uy[k] = (c[j] - aij * uy[k]) / akj
      end
    end

    #####

    push!(operations, FreeLinearDoubletonColumn{T, S}(j, i, k, aij, akj, yi, conival))

    kept_cols[j] = false # col j removed
    col_cnt[j] = -1
    kept_rows[i] = false # row i removed
    row_cnt[i] = -1
    row_cnt[k] -= 1 # row k updated and lost akj
    qmp.free_dsc_pass = true
  end
  qmp.c0 += c0_offset
end

function postsolve!(
  sol::QMSolution,
  operation::FreeLinearDoubletonColumn{T, S},
  psd::PresolvedData{T, S},
) where {T, S}
  x = sol.x
  kept_rows, kept_cols = psd.kept_rows, psd.kept_cols
  i, j, k = operation.i, operation.j, operations.k
  aij, akj = operation.aij, operation.akj
  arowi = psd.arows[i]
  arowk = psd.arows[k]
  kept_rows[i] = true
  kept_cols[j] = true
  # x[j] = (coival - Σₖ Aik x[k]) / aij , where k ≂̸ j
  x[j] = operation.conival

  c_k = 1
  c_i = 1
  for l in 1:nvar
    (kept_cols[l] && l != j) || continue
    if l == arowk.nzind[c_k]
      akl2 = arowk.nzval[c_k]
      c_k += 1
    else
      akl2 = zero(T) 
    end
    if l == arowi.nzind[c_i]
      ail2 = arowi.nzval[c_i]
      c_i += 1
    else
      ail2 = zero(T)
    end 

    if ail != zero(T)
      ail = -ail2 * aij / akj # ail2 = - akj * ail / aij (use this for postsolve)
    end

    x[j] -= ail * x[l]
  end
  x[j] /= aij
  sol.s_l[j] = zero(T)
  sol.s_u[j] = zero(T)
  sol.y[i] = (c[j] - akj * y[k]) / aij
end

function pass_arowi_cat(arows::Vector{Row{T}}, s, i, kept_rows, kept_cols, nvar) where {T}
  arowi = arows[i]
  for j=1:nvar
    kept_col[j] || continue
    aij = zero(T)
    while s[i] != 0
      
    end

end