function doubleton_rows!(
  arows::Vector{Row{T}},
  lcon::S,
  ucon::S,
  lvar::S,
  uvar::S,
  ncon,
  kept_rows,
  kept_cols,
  row_cnt,
) where {T, S}
  doubl_row_pass = false
  for i=1:ncon
    (kept_rows[i] && row_cnt[i] == 2 && lcon[i] == ucon[i]) || continue
    doubl_row_pass = true
    lconi = lcon[i]
    rowi = arows[i]
    k = 0
    aik = zero(T)
    j = 0
    aij = zero(T)
    c = 0
    for (l, ail) in zip(rowi.nzind, rowi.nzval)
      kept_cols[l] || continue
      c += 1
      if c == 1
        k = l
        aik = ail
      elseif c == 2
        j = l
        aij = ail
      end
    end
    if sign(aij) == sign(aik)
      lvar[k] = max(lvar[k], (lconi - aij * uvar[j]) / aik)
      uvar[k] = min(uvar[k], (lconi - aij * lvar[j]) / aik)
    else
      lvar[k] = max(lvar[k], (lconi - aij * lvar[j]) / aik)
      uvar[k] = min(uvar[k], (lconi - aij * uvar[j]) / aik)
    end
    lvar[j] = -T(Inf)
    uvar[j] = T(Inf)
  end
  return doubl_row_pass
end