struct FreeLinearDoubletonColumn{T, S} <: PresolveOperation{T, S}
  i::Int
  j::Int
  aij::T
  arowi::Row{T}
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
    for (l, alj) in zip(colj.nzind, colj.nzval)
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
    if lcon[k] == ucon[k] && lcon[i] != ucon[i] # swap i and k
      l = i
      i = k
      k = l
      alj = aij
      aij = akj
      akj = alj
    end
    lcon[i] == ucon[i] || continue

    arowk = rows[k]
    arowi = rows[i]
    # remove variable j from constraint
    yi = c[j] / aij
    nzcj = c[j] != zero(T)
    nzcj && (c0_offset += yi * ucon[i]) # update c0
    conival = lcon[i]
    # new row to store for postsolve:
    nb_elem_i = row_cnt[i] - 1
    rowi = arows[i] # i-th row
    rowi2 = Row(zeros(Int, nb_elem_i), zeros(T, nb_elem_i))
    c_i = 1
    for k = 1:length(rowi.nzind)
      j2 = rowi.nzind[k]
      # add all col elements to rowi2 except the col j2 == j
      if kept_cols[j2] && j2 != j
        rowi2.nzind[c_i] = j2
        aij2 = rowi.nzval[k]
        rowi2.nzval[c_i] = aij2
        nzcj && (c[j2] -= yi * aij2) # update c if c[j] != 0
        col_cnt[j2] -= 1
        c_i += 1
      end
    end
    push!(operations, FreeLinearDoubletonColumn{T, S}(i, j, aij, rowi2, yi, conival))


    # check aij not too small?
    c_k = 1
    c_i = 1
    for l in 1:nvar
      kept_cols[l] || continue
      if l == arowk.nzind[c_k]
        akl = arowk.nzval[c_k]
        c_k += 1
      else
        akl = zero(T) 
      end
      if l == arowi.nzind[c_i]
        ail = arowi.nzvil[c_i]
        c_i += 1
      else
        ail = zero(T)
      end

      akl2 = akl - akj * ail / aij
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
  end
end