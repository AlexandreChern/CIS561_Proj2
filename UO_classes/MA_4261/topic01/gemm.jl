function mygemm_pji!(C, A, B)
  ## Get the matrix sizes
  m, n = size(C)
  _, k = size(A)

  ## Check the matrix sizes
  @assert size(C) == (m, n)
  @assert size(A) == (m, k)
  @assert size(B) == (k, n)

  @inbounds begin
    for p = 1:k
      for j = 1:n
        for i = 1:m
          C[i, j] += A[i, p] * B[p, j]
        end
      end
    end
  end

  return C
end

function mygemm_jpi!(C, A, B)
  ## Get the matrix sizes
  m, n = size(C)
  _, k = size(A)

  ## Check the matrix sizes
  @assert size(C) == (m, n)
  @assert size(A) == (m, k)
  @assert size(B) == (k, n)

  @inbounds begin
    for j = 1:n
      for p = 1:k
        for i = 1:m
          C[i, j] += A[i, p] * B[p, j]
        end
      end
    end
  end

  return C
end

function mygemm_jip_pji!(C, A, B, mb, nb, kb)
  ## Get the matrix sizes
  m, n = size(C)
  _, k = size(A)

  ## Check the matrix sizes
  @assert size(C) == (m, n)
  @assert size(A) == (m, k)
  @assert size(B) == (k, n)

  @inbounds for J = 1:mb:m
    jrng = J:min(n, J+nb-1)
    for I = 1:nb:n
      irng = I:min(m, I+mb-1)
      for P = 1:kb:k
        prng = P:min(k, P+kb-1)
        @views mygemm_pji!(C[irng, jrng], A[irng, prng], B[prng, jrng])
      end
    end
  end
end

function mygemm_ji_pji!(C, A, B, mb, nb)
  ## Get the matrix sizes
  m, n = size(C)
  _, k = size(A)

  ## Check the matrix sizes
  @assert size(C) == (m, n)
  @assert size(A) == (m, k)
  @assert size(B) == (k, n)

  @inbounds for J = 1:mb:m
    J2 = min(n, J+nb-1)
    for I = 1:mb:m
      I2 = min(m, I+mb-1)
      @views mygemm_pji!(C[I:I2, J:J2], A[I:I2, :], B[:, J:J2])
    end
  end
end

# mr and nr need to be known at compile time!
# h/t to https://github.com/Sacha0/GemmDemo.jl for this realization
function mygemm_jpi_ji_pji!(C, A, B, mb, nb, kb, ::Val{mr}, ::Val{nr}) where {mr, nr}
  ## Get the matrix sizes
  m, n = size(C)
  _, k = size(A)

  ## Check the matrix sizes
  @assert size(C) == (m, n)
  @assert size(A) == (m, k)
  @assert size(B) == (k, n)

  @inbounds begin
    # Cache blocking
    for jb1 = 1:nb:n
      jb2 = min(n, jb1+nb-1)
      for pb1 = 1:kb:k
        pb2 = min(k, pb1+kb-1)
        for ib1 = 1:mb:m
          ib2 = min(n, ib1+mb-1)
          # Register Blocking
          for jr = jb1:nr:jb2
            for ir = ib1:mr:ib2
              #=
              @views mygemm_pji!(C[ir:ir+mr-1, jr:jr+nr-1],
                                 A[ir:ir+mr-1, pb1:pb2],
                                 B[pb1:pb2, jr:jr+nr-1])
              =#
              mygemm_pji!(C, A, B, ir, jr, pb1, pb2, Val(mr), Val(nr))
            end
          end
        end
      end
    end
  end
  return C
end

using SIMD
@inline function mygemm_pji!(C, A, B, ir, jr, p1, p2, ::Val{4}, ::Val{4})
  T = eltype(C)
  m = size(C, 1)
  ptrC = pointer(C, ir + (jr-1) * m)
  c_1 = vload(Vec{4, T}, ptrC + 0m*sizeof(T))
  c_2 = vload(Vec{4, T}, ptrC + 1m*sizeof(T))
  c_3 = vload(Vec{4, T}, ptrC + 2m*sizeof(T))
  c_4 = vload(Vec{4, T}, ptrC + 3m*sizeof(T))

  ptrA = pointer(A, ir)
  @inbounds for p = p1:p2
    a_p = vload(Vec{4, T}, ptrA + (p-1)*m*sizeof(T))

    b_pj = Vec{4,T}(B[p, jr+0])
    c_1 = fma(b_pj, a_p, c_1)

    b_pj = Vec{4,T}(B[p, jr+1])
    c_2 = fma(b_pj, a_p, c_2)

    b_pj = Vec{4,T}(B[p, jr+2])
    c_3 = fma(b_pj, a_p, c_3)

    b_pj = Vec{4,T}(B[p, jr+3])
    c_4 = fma(b_pj, a_p, c_4)
  end
  vstore(c_1, ptrC + 0m*sizeof(T))
  vstore(c_2, ptrC + 1m*sizeof(T))
  vstore(c_3, ptrC + 2m*sizeof(T))
  vstore(c_4, ptrC + 3m*sizeof(T))

  return C
end

@inline function mygemm_pji!(C, A, B, ir, jr, p1, p2, ::Val{mr}, ::Val{4}) where mr
  @assert mod(mr, 4) == 0 && mr <= 12

  T = eltype(C)
  m = size(C, 1)

  ptrC = pointer(C, ir + (jr-1) * m)
  mr > 0 && (c_1a = vload(Vec{4, T}, ptrC + (0m + 0)*sizeof(T)))
  mr > 4 && (c_1b = vload(Vec{4, T}, ptrC + (0m + 4)*sizeof(T)))
  mr > 8 && (c_1c = vload(Vec{4, T}, ptrC + (0m + 8)*sizeof(T)))
  mr > 0 && (c_2a = vload(Vec{4, T}, ptrC + (1m + 0)*sizeof(T)))
  mr > 4 && (c_2b = vload(Vec{4, T}, ptrC + (1m + 4)*sizeof(T)))
  mr > 8 && (c_2c = vload(Vec{4, T}, ptrC + (1m + 8)*sizeof(T)))
  mr > 0 && (c_3a = vload(Vec{4, T}, ptrC + (2m + 0)*sizeof(T)))
  mr > 4 && (c_3b = vload(Vec{4, T}, ptrC + (2m + 4)*sizeof(T)))
  mr > 8 && (c_3c = vload(Vec{4, T}, ptrC + (2m + 8)*sizeof(T)))
  mr > 0 && (c_4a = vload(Vec{4, T}, ptrC + (3m + 0)*sizeof(T)))
  mr > 4 && (c_4b = vload(Vec{4, T}, ptrC + (3m + 4)*sizeof(T)))
  mr > 8 && (c_4c = vload(Vec{4, T}, ptrC + (3m + 8)*sizeof(T)))

  ptrA = pointer(A, ir)
  @inbounds for p = p1:p2
    mr > 0 && (a_pa = vload(Vec{4, T}, ptrA + ((p-1)*m + 0)*sizeof(T)))
    mr > 4 && (a_pb = vload(Vec{4, T}, ptrA + ((p-1)*m + 4)*sizeof(T)))
    mr > 8 && (a_pc = vload(Vec{4, T}, ptrA + ((p-1)*m + 8)*sizeof(T)))

    b_pj = Vec{4,T}(B[p, jr+0])
    mr > 0 && (c_1a = fma(b_pj, a_pa, c_1a))
    mr > 4 && (c_1b = fma(b_pj, a_pb, c_1b))
    mr > 8 && (c_1c = fma(b_pj, a_pc, c_1c))

    b_pj = Vec{4,T}(B[p, jr+1])
    mr > 0 && (c_2a = fma(b_pj, a_pa, c_2a))
    mr > 4 && (c_2b = fma(b_pj, a_pb, c_2b))
    mr > 8 && (c_2c = fma(b_pj, a_pc, c_2c))

    b_pj = Vec{4,T}(B[p, jr+2])
    mr > 0 && (c_3a = fma(b_pj, a_pa, c_3a))
    mr > 4 && (c_3b = fma(b_pj, a_pb, c_3b))
    mr > 8 && (c_3c = fma(b_pj, a_pc, c_3c))

    b_pj = Vec{4,T}(B[p, jr+3])
    mr > 0 && (c_4a = fma(b_pj, a_pa, c_4a))
    mr > 4 && (c_4b = fma(b_pj, a_pb, c_4b))
    mr > 8 && (c_4c = fma(b_pj, a_pc, c_4c))
  end

  mr > 0 && (vstore(c_1a, ptrC + (0m + 0)*sizeof(T)))
  mr > 4 && (vstore(c_1b, ptrC + (0m + 4)*sizeof(T)))
  mr > 8 && (vstore(c_1c, ptrC + (0m + 8)*sizeof(T)))
  mr > 0 && (vstore(c_2a, ptrC + (1m + 0)*sizeof(T)))
  mr > 4 && (vstore(c_2b, ptrC + (1m + 4)*sizeof(T)))
  mr > 8 && (vstore(c_2c, ptrC + (1m + 8)*sizeof(T)))
  mr > 0 && (vstore(c_3a, ptrC + (2m + 0)*sizeof(T)))
  mr > 4 && (vstore(c_3b, ptrC + (2m + 4)*sizeof(T)))
  mr > 8 && (vstore(c_3c, ptrC + (2m + 8)*sizeof(T)))
  mr > 0 && (vstore(c_4a, ptrC + (3m + 0)*sizeof(T)))
  mr > 4 && (vstore(c_4b, ptrC + (3m + 4)*sizeof(T)))
  mr > 8 && (vstore(c_4c, ptrC + (3m + 8)*sizeof(T)))

  return C
end

function mygemm_jpi_ji_pji_packed!(C, A, B, mb, nb, kb, ::Val{mr}, ::Val{nr}) where {mr, nr}
  ## Get the matrix sizes
  m, n = size(C)
  _, k = size(A)

  ## Check the matrix sizes
  @assert size(C) == (m, n)
  @assert size(A) == (m, k)
  @assert size(B) == (k, n)

  Bt = similar(B, nr, kb, div(nb, nr))
  At = similar(A, mr, kb, div(mb, mr))

  @inbounds begin
    # Cache blocking
    for jb1 = 1:nb:n
      jb2 = min(n, jb1+nb-1)
      for pb1 = 1:kb:k
        pb2 = min(k, pb1+kb-1)
        for (jcol, j) in enumerate(jb1:nr:jb2)
          Bt[:, 1:pb2-pb1+1, jcol] = B[pb1:pb2, j:j+nr-1]'
        end
        for ib1 = 1:mb:m
          ib2 = min(n, ib1+mb-1)
          for (icol, i) in enumerate(ib1:mr:ib2)
            At[:, 1:pb2-pb1+1, icol] = A[i:i+mr-1, pb1:pb2]
          end
          # Register Blocking
          for (jcol, jr) in enumerate(jb1:nr:jb2)
            for (icol, ir) in enumerate(ib1:mr:ib2)
              mygemm_pji_packed!(C, At, Bt, ir, jr, icol, jcol, pb2-pb1+1, kb, Val(mr), Val(nr))
            end
          end
        end
      end
    end
  end
  return C
end

@inline function mygemm_pji_packed!(C, At, Bt, ir, jr, icol, jcol, pend, kb, ::Val{4}, ::Val{4})
  T = eltype(C)
  m = size(C, 1)
  ptrC = pointer(C, ir + (jr-1) * m)
  c_1 = vload(Vec{4, T}, ptrC + 0m*sizeof(T))
  c_2 = vload(Vec{4, T}, ptrC + 1m*sizeof(T))
  c_3 = vload(Vec{4, T}, ptrC + 2m*sizeof(T))
  c_4 = vload(Vec{4, T}, ptrC + 3m*sizeof(T))

  mr = 4
  ptrAt = pointer(At, 1 + (icol-1) * mr * kb)
  @inbounds for p = 1:pend
    a_p = vload(Vec{4, T}, ptrAt + (p-1)*mr*sizeof(T))

    b_pj = Vec{4,T}(Bt[1, p, jcol])
    c_1 = fma(b_pj, a_p, c_1)

    b_pj = Vec{4,T}(Bt[2, p, jcol])
    c_2 = fma(b_pj, a_p, c_2)

    b_pj = Vec{4,T}(Bt[3, p, jcol])
    c_3 = fma(b_pj, a_p, c_3)

    b_pj = Vec{4,T}(Bt[4, p, jcol])
    c_4 = fma(b_pj, a_p, c_4)
  end
  vstore(c_1, ptrC + 0m*sizeof(T))
  vstore(c_2, ptrC + 1m*sizeof(T))
  vstore(c_3, ptrC + 2m*sizeof(T))
  vstore(c_4, ptrC + 3m*sizeof(T))

  return C
end

@inline function mygemm_pji_packed!(C, At, Bt, ir, jr, icol, jcol, pend, kb, ::Val{mr}, ::Val{4}) where mr
  @assert mod(mr, 4) == 0 && mr <= 12

  T = eltype(C)
  m = size(C, 1)

  ptrC = pointer(C, ir + (jr-1) * m)
  mr > 0 && (c_1a = vload(Vec{4, T}, ptrC + (0m + 0)*sizeof(T)))
  mr > 4 && (c_1b = vload(Vec{4, T}, ptrC + (0m + 4)*sizeof(T)))
  mr > 8 && (c_1c = vload(Vec{4, T}, ptrC + (0m + 8)*sizeof(T)))
  mr > 0 && (c_2a = vload(Vec{4, T}, ptrC + (1m + 0)*sizeof(T)))
  mr > 4 && (c_2b = vload(Vec{4, T}, ptrC + (1m + 4)*sizeof(T)))
  mr > 8 && (c_2c = vload(Vec{4, T}, ptrC + (1m + 8)*sizeof(T)))
  mr > 0 && (c_3a = vload(Vec{4, T}, ptrC + (2m + 0)*sizeof(T)))
  mr > 4 && (c_3b = vload(Vec{4, T}, ptrC + (2m + 4)*sizeof(T)))
  mr > 8 && (c_3c = vload(Vec{4, T}, ptrC + (2m + 8)*sizeof(T)))
  mr > 0 && (c_4a = vload(Vec{4, T}, ptrC + (3m + 0)*sizeof(T)))
  mr > 4 && (c_4b = vload(Vec{4, T}, ptrC + (3m + 4)*sizeof(T)))
  mr > 8 && (c_4c = vload(Vec{4, T}, ptrC + (3m + 8)*sizeof(T)))

  ptrAt = pointer(At, 1 + (icol-1) * mr * kb)
  @inbounds for p = 1:pend
    mr > 0 && (a_pa = vload(Vec{4, T}, ptrAt + ((p-1)*mr + 0)*sizeof(T)))
    mr > 4 && (a_pb = vload(Vec{4, T}, ptrAt + ((p-1)*mr + 4)*sizeof(T)))
    mr > 8 && (a_pc = vload(Vec{4, T}, ptrAt + ((p-1)*mr + 8)*sizeof(T)))

    b_pj = Vec{4,T}(Bt[1, p, jcol])
    mr > 0 && (c_1a = fma(b_pj, a_pa, c_1a))
    mr > 4 && (c_1b = fma(b_pj, a_pb, c_1b))
    mr > 8 && (c_1c = fma(b_pj, a_pc, c_1c))

    b_pj = Vec{4,T}(Bt[2, p, jcol])
    mr > 0 && (c_2a = fma(b_pj, a_pa, c_2a))
    mr > 4 && (c_2b = fma(b_pj, a_pb, c_2b))
    mr > 8 && (c_2c = fma(b_pj, a_pc, c_2c))

    b_pj = Vec{4,T}(Bt[3, p, jcol])
    mr > 0 && (c_3a = fma(b_pj, a_pa, c_3a))
    mr > 4 && (c_3b = fma(b_pj, a_pb, c_3b))
    mr > 8 && (c_3c = fma(b_pj, a_pc, c_3c))

    b_pj = Vec{4,T}(Bt[4, p, jcol])
    mr > 0 && (c_4a = fma(b_pj, a_pa, c_4a))
    mr > 4 && (c_4b = fma(b_pj, a_pb, c_4b))
    mr > 8 && (c_4c = fma(b_pj, a_pc, c_4c))
  end

  mr > 0 && (vstore(c_1a, ptrC + (0m + 0)*sizeof(T)))
  mr > 4 && (vstore(c_1b, ptrC + (0m + 4)*sizeof(T)))
  mr > 8 && (vstore(c_1c, ptrC + (0m + 8)*sizeof(T)))
  mr > 0 && (vstore(c_2a, ptrC + (1m + 0)*sizeof(T)))
  mr > 4 && (vstore(c_2b, ptrC + (1m + 4)*sizeof(T)))
  mr > 8 && (vstore(c_2c, ptrC + (1m + 8)*sizeof(T)))
  mr > 0 && (vstore(c_3a, ptrC + (2m + 0)*sizeof(T)))
  mr > 4 && (vstore(c_3b, ptrC + (2m + 4)*sizeof(T)))
  mr > 8 && (vstore(c_3c, ptrC + (2m + 8)*sizeof(T)))
  mr > 0 && (vstore(c_4a, ptrC + (3m + 0)*sizeof(T)))
  mr > 4 && (vstore(c_4b, ptrC + (3m + 4)*sizeof(T)))
  mr > 8 && (vstore(c_4c, ptrC + (3m + 8)*sizeof(T)))

  return C
end
