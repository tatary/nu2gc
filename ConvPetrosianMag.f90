!!$ ifort -o conv2pet.out ConvPetrosianMag.f90
  PROGRAM ConvPetrosianMag
    implicit none

    INTEGER :: i_in, i_out, i_tab, ier, i_nbody = 1
    CHARACTER(LEN=50) :: fname_tab = '/home/nugc/Pet2TotRatio.dat'
    CHARACTER(LEN=50) :: fbase(2) = &
         (/'/home/kobayasi/nugc/data/LAE/2013Sep/z0p0LAE100Mpc', &
           '/home/kobayasi/nugc/data/LAE/2013Sep/z0p0LAE400Mpc'/)
    DOUBLE PRECISION :: weight(2) = (/70.d0**-3.d0, 280.d0**-3.d0/) 
                        ! 70^-3 [h^3/Mpc^3] or 280^-3 [h^3/Mpc^3]
    CHARACTER(LEN=500) :: fname_in, fname_out
    CHARACTER(LEN=10000) :: buf

    TYPE Pet2TotRatio
       INTEGER :: nbin_r, nbin_bt
       DOUBLE PRECISION, ALLOCATABLE :: rb2rd(:) ! Re(bulge)/Re(disk) ratio
       DOUBLE PRECISION, ALLOCATABLE :: BT(:)    ! L(bulge)/L(total) ratio
       DOUBLE PRECISION, ALLOCATABLE :: LP2Ltot(:,:) ! L(Petrosian)/L(total) ratio
       DOUBLE PRECISION, ALLOCATABLE :: magP2tot(:,:) ! mag(total)+magP2tot --> mag(Petrosian)
       DOUBLE PRECISION, ALLOCATABLE :: rp2rb(:,:) ! R(Petrosian)/Re(bulge) ratio
       DOUBLE PRECISION, ALLOCATABLE :: rp2rd(:,:) ! R(Petrosian)/Re(disk) ratio
    END TYPE Pet2TotRatio
    TYPE(Pet2TotRatio) :: p
    INTEGER :: ncolumn, i_r, i_bt, i_rm1, i_btm1
    DOUBLE PRECISION :: x

    INTEGER, PARAMETER :: N_COLUMN   = 33 ! # of double precision column
    INTEGER, PARAMETER :: N_COLUMN_I = 10 ! # of integer column
    INTEGER, PARAMETER :: NWAVE = 6
    TYPE LuminosityFunction
       INTEGER :: Nbin, N1, N2
       DOUBLE PRECISION :: step, base, invstep
       DOUBLE PRECISION, ALLOCATABLE :: bin(:), n(:,:,:)
    END TYPE LuminosityFunction
    TYPE(LuminosityFunction) :: lf(2), lf_d(2) ! 1:original, 2:Petrosian
    INTEGER :: i_lf, ic, iloop, mor, mor_d, iband, ibin, i, i1, i2
    INTEGER(KIND=8), DIMENSION(:) :: ixin(N_COLUMN_I)
    DOUBLE PRECISION, DIMENSION(:) :: xin(N_COLUMN)
    DOUBLE PRECISION :: LB_b,  LB_d,  Ltot,   BT,   BTm1,   rb, rd, rb2rd
    DOUBLE PRECISION :: LBd_b, LBd_d, Ltot_d, BT_d, BTm1_d
    DOUBLE PRECISION, PARAMETER :: EPS = 1.d-10
    DOUBLE PRECISION :: LP2Ltot,   magP2tot,   rp2rb,   rp2rd,   rp,   mag(NWAVE)
    DOUBLE PRECISION :: LP2Ltot_d, magP2tot_d, rp2rb_d, rp2rd_d, rp_d, mag_d(NWAVE)
    DOUBLE PRECISION :: LinearInterp ! function


    ! --------------------------------------------------------------------------
    !  read the data table in which conversion from total to Petrosian
    !      quantities (i.e., luminosity, radius, and concentration) is written
    ! --------------------------------------------------------------------------
    i_tab = 20
    ! --- count the number of column in the data table
    !      and allocate the structure for the table
    open(i_tab, file=trim(fname_tab), status='old', iostat=ier); call &
            CheckIerror(ier, '# fail to open file='//trim(fname_tab))
    read(i_tab, *); read(i_tab, *) ! header
    ier = 0; ncolumn = 0
    DO WHILE(ier == 0)
       read(i_tab, *, iostat=ier)
       ncolumn = ncolumn + 1
    ENDDO
    close(i_tab)
!!$    print '(A,I10, A,F10.5)', '# ncolumn = ', ncolumn, ': ncolumn^0.5 = ', sqrt(dble(ncolumn))
    p%nbin_r = sqrt(dble(ncolumn)) + 1; p%nbin_bt = sqrt(dble(ncolumn))
!!$    print '(2(A,I10))', '# nbin_r = ', p%nbin_r, ', nbin_bt = ', p%nbin_bt
    allocate(p%rb2rd(p%nbin_r), stat=ier); call CheckIerror(ier, &
         '# allocation fault: p%rb2rd')
    allocate(p%BT(p%nbin_bt), stat=ier); call CheckIerror(ier, &
         '# allocation fault: p%BT')
    allocate(p%LP2Ltot(p%nbin_r, p%nbin_bt), stat=ier); call CheckIerror(ier, &
         '# allocation fault: p%LP2Ltot')
    allocate(p%magP2tot(p%nbin_r, p%nbin_bt), stat=ier); call CheckIerror(ier, &
         '# allocation fault: p%magP2tot')
    allocate(p%rp2rb(p%nbin_r, p%nbin_bt), stat=ier); call CheckIerror(ier, &
         '# allocation fault: p%rp2rb')
    allocate(p%rp2rd(p%nbin_r, p%nbin_bt), stat=ier); call CheckIerror(ier, &
         '# allocation fault: p%rp2rd')
    ! --- substitute the quantities in the table into the structure
    open(i_tab, file=trim(fname_tab), status='old', iostat=ier); call &
            CheckIerror(ier, '# fail to open file='//trim(fname_tab))
    read(i_tab, *); read(i_tab, *) ! header
    ier = 0; i_r = 1; i_bt = 1
    DO WHILE(ier == 0)
!!$       print '(2I5)', i_r, i_bt
       read(i_tab, '(A)', iostat=ier) buf
       IF(ier == 0) &
          read(buf, *) p%rb2rd(i_r), p%BT(i_bt), p%rp2rb(i_r, i_bt), p%LP2Ltot(i_r, i_bt), &
               p%rp2rd(i_r, i_bt), x, p%magP2tot(i_r, i_bt), x, x
       i_bt = i_bt + 1
       IF(i_bt > p%nbin_bt) THEN
          i_bt = 1
          i_r = i_r + 1
       ENDIF
    ENDDO
    close(i_tab)


    ! --------------------------------------------------------------------------
    !   read the data in which nuGC outputs are written (default = *1.dat),
    !     convert total mag. into Petrosian mag. according to the B/T ratio
    !     and rb/rd ratio of each galaxy, and calculate luminosity function
    ! --------------------------------------------------------------------------
    ! --- allocate the structure for luminosity functions
    DO i_lf = 1, 2 ! i_lf = 1: original (total mag.), 2: Petrosian mag.
       lf(i_lf)%Nbin = 200; lf(i_lf)%N1 = NWAVE; lf(i_lf)%N2 = 4; lf(i_lf)%base = -40.25d0
       lf(i_lf)%step = 0.5d0; lf(i_lf)%invstep = 1.d0 / lf(i_lf)%step
       allocate(lf(i_lf)%bin(lf(i_lf)%Nbin), stat=ier); call CheckIerror(ier, &
            '# allocation fault: lf(i_lf)%bin')
       allocate(lf(i_lf)%n(lf(i_lf)%N1, lf(i_lf)%N2, lf(i_lf)%Nbin), stat=ier); call &
            CheckIerror(ier, '# allocation fault: lf(i_lf)%n')
       DO ibin = 1, lf(i_lf)%Nbin
          lf(i_lf)%bin(ibin) = lf(i_lf)%step * dble(ibin) + lf(i_lf)%base
          DO i1 = 1, lf(i_lf)%N1
             DO i2 = 1, lf(i_lf)%N2
                lf(i_lf)%n(i1, i2, ibin) = 0.d0
             ENDDO
          ENDDO
       ENDDO
       lf_d(i_lf)%Nbin = 200; lf_d(i_lf)%N1 = NWAVE; lf_d(i_lf)%N2 = 4; lf_d(i_lf)%base = -40.25d0
       lf_d(i_lf)%step = 0.5d0; lf_d(i_lf)%invstep = 1.d0 / lf_d(i_lf)%step
       allocate(lf_d(i_lf)%bin(lf_d(i_lf)%Nbin), stat=ier); call CheckIerror(ier, &
            '# allocation fault: lf_d(i_lf)%bin')
       allocate(lf_d(i_lf)%n(lf_d(i_lf)%N1, lf_d(i_lf)%N2, lf_d(i_lf)%Nbin), stat=ier); call &
            CheckIerror(ier, '# allocation fault: lf_d(i_lf)%n')
       DO ibin = 1, lf_d(i_lf)%Nbin
          lf_d(i_lf)%bin(ibin) = lf_d(i_lf)%step * dble(ibin) + lf_d(i_lf)%base
          DO i1 = 1, lf_d(i_lf)%N1
             DO i2 = 1, lf_d(i_lf)%N2
                lf_d(i_lf)%n(i1, i2, ibin) = 0.d0
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    ! --- set pointers for the input and output files
    i_in = 20; i_out = 21
    fname_in  = trim(fbase(i_nbody))//'1.dat'
    fname_out = trim(fbase(i_nbody))//'1_Pet.dat'
    open(i_in, file=trim(fname_in), status='old', iostat=ier); call &
            CheckIerror(ier, '# fail to open file='//trim(fname_in))
    open(i_out, file=trim(fname_out), status='unknown', iostat=ier); call &
            CheckIerror(ier, '# fail to open file='//trim(fname_out))
    ! --- main loop
    ier = 0; iloop = 1
    DO WHILE(ier == 0)
       read(i_in, '(A)', iostat=ier) buf
       IF(iloop <= 5) THEN ! for headers
          write(i_out, '(A, \)') trim(buf)
          IF(iloop == 2) THEN
             write(i_out, '(A)') ' (44,45)R_Pet(w/o,w/ dust)[kpc/h] (46,47)L_B,Petrosian'//&
                  '(w/o,w/dust) (48,49)Delta(mag[Petrosian])(w/o,w/ dust)'
          ELSE
             write(i_out, *)
          ENDIF
          iloop = iloop + 1
       ELSE ! for main data
          write(i_out, '(A, \)') trim(buf)
          IF(ier == 0) THEN ! in the case that input file does not reach its end
             read(buf, *) (ixin(ic), ic=1, N_COLUMN_I), (xin(ic), ic=1, N_COLUMN)
             ! --- substitute the read data into the adequate quantities
             mor = ixin(3); mor_d = ixin(4)
             LB_b = xin(15); LBd_b = xin(16); LB_d = xin(17); LBd_d = xin(18)
             Ltot = LB_b + LB_d; Ltot_d = LBd_b + LBd_d
             BT = LB_b / (LB_b + LB_d); BT_d = LBd_b / (LBd_b + LBd_d)
             BTm1 = BT - 1.d0; BTm1_d = BT_d - 1.d0
             IF(BTm1   < 0.d0) BTm1   = -BTm1
             IF(BTm1_d < 0.d0) BTm1_d = -BTm1_d
             rb = xin(11) ! [kpc/h]
             rd = xin(13) ! [kpc/h]
             DO iband = 1, NWAVE
                i = 2*(iband-1) + 19
                mag(iband) = xin(i); mag_d(iband) = xin(i+1)
             ENDDO
             ! --- conversion into Petrosian quantities for w/o dust
             IF(BT < EPS) THEN ! pure disk
                i_r = 1; i_bt = 1
                LP2Ltot = p%LP2Ltot(i_r, i_bt); magP2tot = p%magP2tot(i_r, i_bt)
                rp2rb = 0.d0; rp2rd = p%rp2rd(i_r, i_bt); rp = rp2rd * rd
             ELSEIF(BTm1 < EPS) THEN ! pure bulge
                rd = 0.d0
                LP2Ltot = 0.8165d0; magP2tot = -2.5d0 * log10(LP2Ltot)
                rp2rb = 1.714d0; rp2rd = 0.d0; rp = rp2rb * rb
             ELSE
                rb2rd = rb / rd; i_r = 1; i_bt = 1
                IF(rb2rd > p%rb2rd(p%nbin_r)) THEN
                   i_r = p%nbin_r; print '(A, G10.3)', '[rb2rd(1)] ', rb2rd
                ELSE
                   DO WHILE(rb2rd > p%rb2rd(i_r))
                      i_r = i_r + 1
                   ENDDO ! p%rb2rd(i_r-1) < rb2rd <= p%rb2rd(i_r)
                ENDIF
                IF(BT > p%BT(p%nbin_bt)) THEN
                   i_bt = p%nbin_bt!; print '(A, G10.3)', '[BT] ', BTm1
                ELSE
                   DO WHILE(BT > p%BT(i_bt))
                      i_bt = i_bt + 1
                   ENDDO ! p%BT(i_bt-1) < BT <= p%BT(i_bt)
                ENDIF
                IF(i_r == 1) THEN
                   IF(i_bt == 1) THEN
                      LP2Ltot = p%LP2Ltot(i_bt, i_r); magP2tot = p%magP2tot(i_r, i_bt)
                      rp2rb = 0.d0; rp2rd = p%rp2rd(i_r, i_bt); rp = rp2rd * rd
                   ELSE
                      i_btm1 = i_bt - 1
                      LP2Ltot = LinearInterp(BT, p%BT(i_btm1), p%BT(i_bt), &
                           p%LP2Ltot(i_r, i_btm1), p%LP2Ltot(i_r, i_bt))
                      magP2tot = LinearInterp(BT, p%BT(i_btm1), p%BT(i_bt), &
                           p%magP2tot(i_r, i_btm1), p%magP2tot(i_r, i_bt))
                      rp2rb = LinearInterp(BT, p%BT(i_btm1), p%BT(i_bt), &
                           p%rp2rb(i_r, i_btm1), p%rp2rb(i_r, i_bt))
                      rp2rd = LinearInterp(BT, p%BT(i_btm1), p%BT(i_bt), &
                           p%rp2rd(i_r, i_btm1), p%rp2rd(i_r, i_bt))
                      rp = rp2rb * rb
                   ENDIF
                ELSE
                   IF(i_bt == 1) THEN
                      i_rm1 = i_r - 1
                      LP2Ltot = LinearInterp(rb2rd, p%rb2rd(i_rm1), p%rb2rd(i_r), &
                           p%LP2Ltot(i_rm1, i_bt), p%LP2Ltot(i_r, i_bt))
                      magP2tot = LinearInterp(rb2rd, p%rb2rd(i_rm1), p%rb2rd(i_r), &
                           p%magP2tot(i_rm1, i_bt), p%magP2tot(i_r, i_bt))
                      rp2rb = LinearInterp(rb2rd, p%rb2rd(i_rm1), p%rb2rd(i_r), &
                           p%rp2rb(i_rm1, i_bt), p%rp2rb(i_r, i_bt))
                      rp2rd = LinearInterp(rb2rd, p%rb2rd(i_rm1), p%rb2rd(i_r), &
                           p%rp2rd(i_rm1, i_bt), p%rp2rd(i_r, i_bt))
                   ELSE
                      call TwoDimensLinearInterp(0, BT, rb2rd, i_bt, i_r, p, LP2Ltot)
                      call TwoDimensLinearInterp(1, BT, rb2rd, i_bt, i_r, p, magP2tot)
                      call TwoDimensLinearInterp(2, BT, rb2rd, i_bt, i_r, p, rp2rb)
                      call TwoDimensLinearInterp(3, BT, rb2rd, i_bt, i_r, p, rp2rd)
                   ENDIF
                   rp = rp2rb * rb
                ENDIF
             ENDIF
             ! --- conversion into Petrosian quantities for w/ dust
             IF(BT_d < EPS) THEN ! pure disk
                i_r = 1; i_bt = 1
                LP2Ltot_d = p%LP2Ltot(i_r, i_bt); magP2tot_d = p%magP2tot(i_r, i_bt)
                rp2rb_d = 0.d0; rp2rd_d = p%rp2rd(i_r, i_bt); rp_d = rp2rd_d * rd
             ELSEIF(BTm1_d < EPS) THEN ! pure bulge
                rd = 0.d0
                LP2Ltot_d = 0.8165d0; magP2tot_d = -2.5d0 * log10(LP2Ltot_d)
                rp2rb_d = 1.714d0; rp2rd_d = 0.d0; rp_d = rp2rb_d * rb
             ELSE
                rb2rd = rb / rd; i_r = 1; i_bt = 1
                IF(rb2rd > p%rb2rd(p%nbin_r)) THEN
                   i_r = p%nbin_r; print '(A, G10.3)', '[rb2rd(2)] ', rb2rd
                ELSE
                   DO WHILE(rb2rd > p%rb2rd(i_r))
                      i_r = i_r + 1
                   ENDDO ! p%rb2rd(i_r-1) < rb2rd <= p%rb2rd(i_r)
                ENDIF
                IF(BT_d > p%BT(p%nbin_bt)) THEN
                   i_bt = p%nbin_bt!; print '(A, G10.3)', '[BT_d] ', BTm1_d
                ELSE
                   DO WHILE(BT_d > p%BT(i_bt))
                      i_bt = i_bt + 1
                   ENDDO ! p%BT(i_bt-1) < BT <= p%BT(i_bt)
                ENDIF
                IF(i_r == 1) THEN
                   IF(i_bt == 1) THEN
                      LP2Ltot_d = p%LP2Ltot(i_bt, i_r); magP2tot_d = p%magP2tot(i_r, i_bt)
                      rp2rb_d = 0.d0; rp2rd_d = p%rp2rd(i_r, i_bt); rp_d = rp2rd_d * rd
                   ELSE
                      i_btm1 = i_bt - 1
                      LP2Ltot_d = LinearInterp(BT, p%BT(i_btm1), p%BT(i_bt), &
                           p%LP2Ltot(i_r, i_btm1), p%LP2Ltot(i_r, i_bt))
                      magP2tot_d = LinearInterp(BT, p%BT(i_btm1), p%BT(i_bt), &
                           p%magP2tot(i_r, i_btm1), p%magP2tot(i_r, i_bt))
                      rp2rb_d = LinearInterp(BT, p%BT(i_btm1), p%BT(i_bt), &
                           p%rp2rb(i_r, i_btm1), p%rp2rb(i_r, i_bt))
                      rp2rd_d = LinearInterp(BT, p%BT(i_btm1), p%BT(i_bt), &
                           p%rp2rd(i_r, i_btm1), p%rp2rd(i_r, i_bt))
                      rp_d = rp2rb_d * rb
                   ENDIF
                ELSE
                   IF(i_bt == 1) THEN
                      i_rm1 = i_r - 1
                      LP2Ltot_d = LinearInterp(rb2rd, p%rb2rd(i_rm1), p%rb2rd(i_r), &
                           p%LP2Ltot(i_rm1, i_bt), p%LP2Ltot(i_r, i_bt))
                      magP2tot_d = LinearInterp(rb2rd, p%rb2rd(i_rm1), p%rb2rd(i_r), &
                           p%magP2tot(i_rm1, i_bt), p%magP2tot(i_r, i_bt))
                      rp2rb_d = LinearInterp(rb2rd, p%rb2rd(i_rm1), p%rb2rd(i_r), &
                           p%rp2rb(i_rm1, i_bt), p%rp2rb(i_r, i_bt))
                      rp2rd_d = LinearInterp(rb2rd, p%rb2rd(i_rm1), p%rb2rd(i_r), &
                           p%rp2rd(i_rm1, i_bt), p%rp2rd(i_r, i_bt))
                   ELSE
                      call TwoDimensLinearInterp(0, BT_d, rb2rd, i_bt, i_r, p, LP2Ltot_d)
                      call TwoDimensLinearInterp(1, BT_d, rb2rd, i_bt, i_r, p, magP2tot_d)
                      call TwoDimensLinearInterp(2, BT_d, rb2rd, i_bt, i_r, p, rp2rb_d)
                      call TwoDimensLinearInterp(3, BT_d, rb2rd, i_bt, i_r, p, rp2rd_d)
                   ENDIF
                   rp_d = rp2rb_d * rb
                ENDIF
             ENDIF
             ! --- write the Petrosian quantities into the output file
             write(i_out, '(4G13.5, 2F11.4)') &
                  rp, rp_d, LP2Ltot * Ltot, LP2Ltot_d * Ltot_d, magP2tot, magP2tot_d
             ! --- calculate the contributon to luminosity function
             DO iband = 1, NWAVE
                ! --- intrinsic LF w/ total mag.
                ibin = int(anint((mag(iband) - lf(1)%base) * lf(1)%invstep))
                IF(ibin <= lf(1)%Nbin .and. ibin >= 1) &
                     lf(1)%n(iband, mor_d, ibin) = lf(1)%n(iband, mor_d, ibin) + weight(i_nbody)
                ! --- intrinsic LF w/ Petrosian mag.
                ibin = int(anint((mag(iband)+magP2tot - lf(2)%base) * lf(2)%invstep))
                IF(ibin <= lf(2)%Nbin .and. ibin >= 1) &
                     lf(2)%n(iband, mor_d, ibin) = lf(2)%n(iband, mor_d, ibin) + weight(i_nbody)
                ! --- observable LF w/ total mag.
                ibin = int(anint((mag_d(iband) - lf_d(1)%base) * lf_d(1)%invstep))
                IF(ibin <= lf_d(1)%Nbin .and. ibin >= 1) &
                     lf_d(1)%n(iband, mor_d, ibin) = lf_d(1)%n(iband, mor_d, ibin) + weight(i_nbody)
                ! --- observable LF w/ Petrosian mag.
                ibin = int(anint((mag_d(iband)+magP2tot_d - lf_d(2)%base) * lf_d(2)%invstep))
                IF(ibin <= lf_d(2)%Nbin .and. ibin >= 1) &
                     lf_d(2)%n(iband, mor_d, ibin) = lf_d(2)%n(iband, mor_d, ibin) + weight(i_nbody)
             ENDDO
          ENDIF ! IF(ier == 0) ELSE
       ENDIF ! IF(iloop <= 5) ELSE
    ENDDO ! DO WHILE(ier == 0)
    close(i_in); close(i_out)
    deallocate(p%rb2rd,    stat=ier); call CheckIerror(ier, '# deallocation fault: p%rb2rd')
    deallocate(p%BT,       stat=ier); call CheckIerror(ier, '# deallocation fault: p%BT')
    deallocate(p%LP2Ltot,  stat=ier); call CheckIerror(ier, '# deallocation fault: p%LP2Ltot')
    deallocate(p%magP2tot, stat=ier); call CheckIerror(ier, '# deallocation fault: p%magP2tot')
    deallocate(p%rp2rb,    stat=ier); call CheckIerror(ier, '# deallocation fault: p%rp2rb')
    deallocate(p%rp2rd,    stat=ier); call CheckIerror(ier, '# deallocation fault: p%rp2rd')

    ! --------------------------------------------------------------------------
    !   write the calculated luminosity functions
    ! --------------------------------------------------------------------------
    fname_out = trim(fbase(i_nbody))//'3_Pet.dat'
    open(i_out, file=trim(fname_out), status='unknown', iostat=ier); call &
            CheckIerror(ier, '# fail to open file='//trim(fname_out))
    write(i_out, '(A)') '# morphology-dependent LFs in each band'
    write(i_out, '(A)') '# (1)M-5log(h)[mag] (2-7,20-25:38-43,56-61)E, '//&
         '(8-13,26-31:44-49,62-67)S0, (14-19,32-37:50-55,68-73)S (w/o,w/ dust)(total:Petrosian)'
    write(i_out, '(A)') '# --- total magnitude'
    write(i_out, '(A)') '# (2,  8, 14[w/o dust]: 20, 26, 32[w/ dust])B'
    write(i_out, '(A)') '# (3,  9, 15[w/o dust]: 21, 27, 33[w/ dust])V'
    write(i_out, '(A)') '# (4, 10, 16[w/o dust]: 22, 28, 34[w/ dust])Kp'
    write(i_out, '(A)') '# (5, 11, 17[w/o dust]: 23, 29, 35[w/ dust])B_rest'
    write(i_out, '(A)') '# (6, 12, 18[w/o dust]: 24, 30, 36[w/ dust])V_rest'
    write(i_out, '(A)') '# (7, 13, 19[w/o dust]: 25, 31, 37[w/ dust])Kp_rest'
    write(i_out, '(A)') '# --- Petrosian magnitude'
    write(i_out, '(A)') '# (38, 44, 50[w/o dust]: 56, 62, 68[w/ dust])B'
    write(i_out, '(A)') '# (39, 45, 51[w/o dust]: 57, 63, 69[w/ dust])V'
    write(i_out, '(A)') '# (40, 46, 52[w/o dust]: 58, 64, 70[w/ dust])Kp'
    write(i_out, '(A)') '# (41, 47, 53[w/o dust]: 59, 65, 71[w/ dust])B_rest'
    write(i_out, '(A)') '# (42, 48, 54[w/o dust]: 60, 66, 72[w/ dust])V_rest'
    write(i_out, '(A)') '# (43, 49, 55[w/o dust]: 61, 67, 73[w/ dust])Kp_rest'
    DO ibin = 1, lf(1)%Nbin
       write(i_out, '(F13.6, \)') lf(1)%bin(ibin)
       DO mor = 1, 3 ! intrinsic LF for E,S0,S in total mag.
          DO iband = 1, NWAVE
             IF(lf(1)%n(iband, mor, ibin) < EPS) THEN
                write(i_out, '(2X, E13.6, \)') 1.d-10
             ELSE
                write(i_out, '(2X, E13.6, \)') lf(1)%n(iband, mor, ibin) * lf(1)%invstep
             ENDIF
          ENDDO
       ENDDO
       DO mor = 1, 3 ! observable LF for E,S0,S in total mag.
          DO iband = 1, NWAVE
             IF(lf_d(1)%n(iband, mor, ibin) < EPS) THEN
                write(i_out, '(2X, E13.6, \)') 1.d-10
             ELSE
                write(i_out, '(2X, E13.6, \)') lf_d(1)%n(iband, mor, ibin) * lf_d(1)%invstep
             ENDIF
          ENDDO
       ENDDO
       DO mor = 1, 3 ! intrinsic LF for E,S0,S in Petrosian mag.
          DO iband = 1, NWAVE
             IF(lf(2)%n(iband, mor, ibin) < EPS) THEN
                write(i_out, '(2X, E13.6, \)') 1.d-10
             ELSE
                write(i_out, '(2X, E13.6, \)') lf(2)%n(iband, mor, ibin) * lf(2)%invstep
             ENDIF
          ENDDO
       ENDDO
       DO mor = 1, 3 ! observable LF for E,S0,S in Petrosian mag.
          DO iband = 1, NWAVE
             IF(lf_d(2)%n(iband, mor, ibin) < EPS) THEN
                write(i_out, '(2X, E13.6, \)') 1.d-10
             ELSE
                write(i_out, '(2X, E13.6, \)') lf_d(2)%n(iband, mor, ibin) * lf_d(2)%invstep
             ENDIF
          ENDDO
       ENDDO
       write(i_out, *)
    ENDDO
    close(i_out)
    DO i_lf = 1, 2
       deallocate(lf(i_lf)%bin,   stat=ier); call CheckIerror(ier, &
            '# deallocation fault: lf(i_lf)%bin')
       deallocate(lf(i_lf)%n,     stat=ier); call CheckIerror(ier, &
            '# deallocation fault: lf(i_lf)%n')
       deallocate(lf_d(i_lf)%bin, stat=ier); call CheckIerror(ier, &
            '# deallocation fault: lf_d(i_lf)%bin')
       deallocate(lf_d(i_lf)%n,   stat=ier); call CheckIerror(ier, &
            '# deallocation fault: lf_d(i_lf)%n')
    ENDDO
!!$============================================================================
  CONTAINS
!!$============================================================================
    SUBROUTINE CheckIerror(ier, cerr)
      implicit none
      INTEGER, INTENT(IN) :: ier
      CHARACTER(LEN=*), INTENT(IN) :: cerr

      IF(ier /= 0) THEN
         print '(A, I4)', trim(cerr)//', stat = ', ier; stop
      ENDIF
    END SUBROUTINE CheckIerror
!!$============================================================================
    SUBROUTINE TwoDimensLinearInterp(i_type, x, y, i_x, i_y, p, z)
      INTEGER, INTENT(IN) :: i_type, i_x, i_y
      DOUBLE PRECISION, INTENT(IN)   :: x, y
      TYPE(Pet2TotRatio), INTENT(IN) :: p
      DOUBLE PRECISION, INTENT(INOUT) :: z
      INTEGER :: i_xm1, i_ym1
      DOUBLE PRECISION :: x1, x2, y1, y2, z11, z21, z12, z22, z1x, z2x
      DOUBLE PRECISION :: LinearInterp

      i_xm1 = i_x - 1; i_ym1 = i_y - 1
      x1 = p%BT(i_xm1); x2 = p%BT(i_x); y1 = p%rb2rd(i_ym1); y2 = p%rb2rd(i_y)
      IF(i_type == 0) THEN ! LP2Ltot
         z11 = p%LP2Ltot(i_ym1, i_xm1); z21 = p%LP2Ltot(i_y, i_xm1)
         z12 = p%LP2Ltot(i_ym1, i_x);   z22 = p%LP2Ltot(i_y, i_x)
      ELSEIF(i_type == 1) THEN ! magP2tot
         z11 = p%magP2tot(i_ym1, i_xm1); z21 = p%magP2tot(i_y, i_xm1)
         z12 = p%magP2tot(i_ym1, i_x);   z22 = p%magP2tot(i_y, i_x)
      ELSEIF(i_type == 2) THEN ! rp2rb
         z11 = p%rp2rb(i_ym1, i_xm1); z21 = p%rp2rb(i_y, i_xm1)
         z12 = p%rp2rb(i_ym1, i_x);   z22 = p%rp2rb(i_y, i_x)
      ELSEIF(i_type == 3) THEN ! rp2rd
         z11 = p%rp2rd(i_ym1, i_xm1); z21 = p%rp2rd(i_y, i_xm1)
         z12 = p%rp2rd(i_ym1, i_x);   z22 = p%rp2rd(i_y, i_x)
      ENDIF
      z1x = LinearInterp(x, x1, x2, z11, z21); z2x = LinearInterp(x, x1, x2, z12, z22)
      z = LinearInterp(y, y1, y2, z1x, z2x)
    END SUBROUTINE TwoDimensLinearInterp
!!$============================================================================
  END PROGRAM ConvPetrosianMag
!!$============================================================================
    DOUBLE PRECISION FUNCTION LinearInterp(x, x1, x2, y1, y2) RESULT(y)
      DOUBLE PRECISION, INTENT(IN) :: x, x1, x2, y1, y2

      y = ((y1 - y2) / (x1 - x2)) * (x - x1) + y1
    END FUNCTION LinearInterp
!!$============================================================================
