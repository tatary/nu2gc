SUBROUTINE star_formation
  use global_var
  implicit none
  DOUBLE PRECISION :: ab, taueff, fbh
  DOUBLE PRECISION :: abinv, tmp1, tmp2
  INTEGER, PARAMETER :: b_or_q = 2
  DOUBLE PRECISION, PARAMETER :: EPS = 1.d-20 

  fbh = param%fbh(b_or_q)
  IF(Mc0 > EPS) THEN
     Zc0 = MZc0 / Mc0
     IF(Zc0 > 1.d0) THEN
       print *, 'quiescent'
       print '(3G11.4)', MZc0, Mc0, Zc0
     ENDIF

     ab     = ssp%alp(b_or_q) + beta + fbh

     abinv  = 1.d0 / ab
     taueff = taust * abinv
  
     tmp1   = tlife / taueff
     tmp1   = exp(-tmp1)
     dMcold = Mc0 * (1.d0 - tmp1)
  
     tmp2   = dMcold * abinv
     dMstar = ssp%alp(b_or_q) * tmp2; dMhot = beta * tmp2
     dMbh_norm = tmp2 * fbh
     tmp2   = ssp%p(b_or_q) * tlife / taust
     Zc     = min(1.d0, Zc0 + tmp2)
     dMZh   = beta * abinv &
                   * ((ssp%p(b_or_q) * abinv + Zc) * dMcold &
                   - (Zc - Zc0) * Mc0) ! eq.(20) of Nagashima+05
  ELSE
     Zc0 = 0.d0; dMstar = 0.d0; dMcold = 0.d0; dMhot = 0.d0; dMbh_norm = 0.d0
     Zc  = 0.d0; dMZh   = 0.d0
  ENDIF
END SUBROUTINE star_formation
!!$==========================================================================
SUBROUTINE star_formation_burst(zp1)
  use global_var
  implicit none
  DOUBLE PRECISION, INTENT(IN) :: zp1
  DOUBLE PRECISION :: ab
  DOUBLE PRECISION :: abinv, tmp, fbh
  INTEGER, PARAMETER :: b_or_q = 1
  DOUBLE PRECISION, PARAMETER :: EPS = 1.d-20 
!!$  DOUBLE PRECISION :: diff, diff_max = 0.d0, diff_min = 100.d0

  fbh = param%fbh(b_or_q) !* zp1
  IF(Mc0 > EPS) THEN
     Zc0 = MZc0 / Mc0
     IF(Zc0 > 1.d0) THEN
       print *, 'burst'
       print '(3G11.4)', MZc0, Mc0, Zc0
     ENDIF

     ab     = ssp%alp(b_or_q) + beta + fbh
     abinv  = 1.d0 / ab
     dMcold = Mc0
     tmp    = dMcold * abinv
     dMstar = tmp * ssp%alp(b_or_q); dMhot = tmp * beta
   
     dMbh = tmp * fbh
  
     Zc   = 0.d0
     dMZh = beta * abinv * dMcold * (ssp%p(b_or_q) * abinv + Zc0)
            ! eq.(20) of Nagashima+05 w/ dMcold = Mc0
!!$        diff = (ssp%p(b_or_q) + Zc0) / (ssp%p(b_or_q) * abinv + Zc0) - 1.d0
!!$        IF(diff > diff_max) diff_max = diff
!!$        IF(diff < diff_min) diff_min = diff
!!$        print '(4G10.3)', abinv, diff, diff_min, diff_max  
  ELSE
     Zc0 = 0.d0; dMstar = 0.d0; dMcold = 0.d0; dMhot = 0.d0
     Zc  = 0.d0; dMZh = 0.d0
     dMbh = 0.d0
  ENDIF
END SUBROUTINE star_formation_burst
