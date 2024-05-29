module MCMCrelated
   use global_var
   implicit none

   include 'mpif.h'

   INTEGER :: iMCMC, itnum, j, time ! tmp variables
   DOUBLE PRECISION :: mag_obs_K(4,50), phi_obs_K(4,50, 3), phi_obs_err_K(4,50, 3)
   DOUBLE PRECISION :: mag_obs_r(50), phi_obs_r(50, 3), phi_obs_err_r(50, 3)
   DOUBLE PRECISION :: MHI_obs(25), phi_obs_MHI(25), phi_obs_err_MHI(25)
   DOUBLE PRECISION :: mag_obs_q(3,50), phi_obs_q(3,50), phi_obs_err_q(3,50)
   DOUBLE PRECISION :: MB_obs(242), ML_obs_med(50), ML_obs_75(50), ML_obs_25(50)
   DOUBLE PRECISION :: ML_obs, ML_model
   DOUBLE PRECISION :: obs_tmp
   INTEGER :: num_in_bin
   DOUBLE PRECISION :: p_tmp, chisqtot
   DOUBLE PRECISION :: chisq(NObsProp) 
   DOUBLE PRECISION :: chisq_pre, chisq_tmp, phi_tmp, phi_err_tmp, &
        lf_j, lf_j1, mf_j, mf_j1, tmp_rand, d, dd, d_tmp
   DOUBLE PRECISION :: alpst_tmp, Vst_tmp, tau0st_tmp, eps_SF_tmp 
   DOUBLE PRECISION :: emax_tmp, emin_tmp, Zch_tmp, taud_th_tmp
   DOUBLE PRECISION :: Vhot1_tmp, Vhot2_tmp, alphot1_tmp, alphot2_tmp, &
                       alp_ret_tmp
   DOUBLE PRECISION :: fmerge_tmp, fmajor_tmp, Krem_tmp, fdiss_tmp, &
                       Mh0_tmp, alpha_halo_tmp
   DOUBLE PRECISION :: Vcut_tmp, kappa_tmp, eta_tmp, alp_bow_tmp, eps_bow_tmp
   DOUBLE PRECISION :: em_tmp, fdi_tmp
   DOUBLE PRECISION :: fbh1_tmp, fbh2_tmp, Mbhseed_tmp, tacc0_tmp, &
                       tacc1_tmp, alpha_ad_tmp, gamma_ad_tmp, Edd_max_tmp
   REAL :: ran1
   INTEGER :: idiv, idiv_dex, key, param_bounds_key
   DOUBLE PRECISION :: idiv_inv, param_ID
   DOUBLE PRECISION :: CSFR_model, CSFR_obs
!!$==================================================
CONTAINS
!!$==================================================
  SUBROUTINE ReadObsData
    INTEGER :: ierr, ibin, ibinm
    DOUBLE PRECISION :: errm, errp 
    DOUBLE PRECISION :: dum1, dum2
    DOUBLE PRECISION, PARAMETER :: invh3 = 2.92d0 ! [1/(Mpc^3 dex)] >> [(h/Mpc)^3/dex] 
    DOUBLE PRECISION, PARAMETER :: logh5 = 0.77d0 ! x = mag-5logh = mag+logh5
    DOUBLE PRECISION, PARAMETER:: hconv = 1.1d0 ! for QSOLF conversion of the
                                                ! hubble parameter 
                                                ! 0.7021**3 / param%h**3
    CHARACTER(LEN=50) :: obsdir = '/home/shirakata/workSA/ObsData/MCMC/'

    IF(inode == 0) print *, "# Reading observational data..."

    ! === at z ~ 0 === !
    open(1, file = trim(obsdir)//'GAMA_KbandLF_Driver2012.dat', &
         status = 'old', iostat = ierr)
    call CheckIerr(ierr, 'K-band LF: fail to open file')
    DO ibin = 1, 27
       IF(ibin < 10) THEN
          read(1,*)
       ELSE
          ibinm = ibin - 9
          read(1, *) mag_obs_K(1,ibinm), &
               phi_obs_K(1,ibinm,1), phi_obs_err_K(1,ibinm,1), & ! all gal.
               phi_obs_K(1,ibinm,2), phi_obs_err_K(1,ibinm,2), & ! blue gal.
               phi_obs_K(1,ibinm,3), phi_obs_err_K(1,ibinm,3)    ! red gal.
       ENDIF
    ENDDO
    phi_obs_K(1,:,:)     = phi_obs_K(1,:,:)     * 2.d0
    phi_obs_err_K(1,:,:) = phi_obs_err_K(1,:,:) * 2.d0
    close(1)

    open(1, file = trim(obsdir)//'GAMA_rbandLF_Driver2012.dat', &
         status = 'old', iostat = ierr)
    call CheckIerr(ierr, 'r-band LF: fail to open file')
    DO ibin = 1, 29
       IF(ibin < 10) THEN
          read(1,*)
       ELSE
          ibinm = ibin - 9
          read(1, *) mag_obs_r(ibinm), &
               phi_obs_r(ibinm,1), phi_obs_err_r(ibinm,1), & ! all gal.
               phi_obs_r(ibinm,2), phi_obs_err_r(ibinm,2), & ! blue gal.
               phi_obs_r(ibinm,3), phi_obs_err_r(ibinm,3)    ! red gal.
       ENDIF
    ENDDO
    phi_obs_r(:,:)     = phi_obs_r(:,:)     * 2.d0
    phi_obs_err_r(:,:) = phi_obs_err_r(:,:) * 2.d0
    close(1)

    !!--- M_B vs M/L
    !IF(param%ML)THEN
    !   open(1, file = 'ML_obs_median.dat', status = 'old', iostat = ierr)
    !   call CheckIerr(ierr, 'ML_obs_median.dat: fail to open file')
    !   DO ibin = 1, 31
    !      read(1, *) MB_obs(ibin), ML_obs_med(ibin), ML_obs_75(ibin), ML_obs_25(ibin)
    !   ENDDO
    !   close(1)
    !ENDIF

    !--- HI mass function
    open(1, file = trim(obsdir)//'HIMF_Martin2010.dat', status = 'old', iostat = ierr)
    call CheckIerr(ierr, 'HIMF_Martin2010.dat: fail to open file')
    DO ibin = 1, 28
       IF(ibin < 5) THEN
          read(1,*)
       ELSE
          ibinm = ibin - 4
          read(1, *) MHI_obs(ibinm), phi_obs_MHI(ibinm), phi_obs_err_MHI(ibinm)
       ENDIF
    ENDDO
    close(1)

    ! === at z ~ 0.4 === !
    ! --- QSO luminosity function (added by Shirakata)
    open(1, file = trim(obsdir)//'Ueda14_z04.dat', &
         status = 'old', iostat = ierr)
    call CheckIerr(ierr, 'AGN LF: fail to open file')
    DO ibin = 1, 11
       IF(ibin < 5) THEN
          read(1,*)
       ELSE
          ibinm = ibin - 4
          read(1, *) mag_obs_q(1,ibinm), &
               phi_obs_q(1,ibinm), errp, errm
          phi_obs_err_q(1,ibinm) = errp - errm
       ENDIF
    ENDDO
    phi_obs_q(1,:)     = hconv * phi_obs_q(1,:)
    phi_obs_err_q(1,:) = hconv * phi_obs_err_q(1,:)
    close(1)

    ! === at z ~ 1 === !
    open(1, file = trim(obsdir)//'Cirasuolo10_z10.txt', &
         status = 'old', iostat = ierr)
    call CheckIerr(ierr, 'K-band LF (z~1): fail to open file')
    DO ibin = 1, 13
       IF(ibin < 4) THEN
          read(1,*)
       ELSE
          ibinm = ibin - 3
          read(1, *) mag_obs_K(2,ibinm), phi_obs_K(2,ibinm,1), &
                     dum1, errm, dum2, errp  
          mag_obs_K(2,ibinm) = mag_obs_K(2,ibinm) + logh5 
          phi_obs_K(2,ibinm,1)     = invh3 * 10.d0 ** phi_obs_K(2,ibinm,1)
          phi_obs_err_K(2,ibinm,1) = invh3 * (10.d0 ** errp - 10.d0 ** errm)
       ENDIF
    ENDDO
    close(1)

    ! --- QSO luminosity function (added by Shirakata)
    open(1, file = trim(obsdir)//'Ueda14_z10.dat', &
         status = 'old', iostat = ierr)
    call CheckIerr(ierr, 'AGN LF (z~1): fail to open file')
    DO ibin = 1, 11
       IF(ibin < 5) THEN
          read(1,*)
       ELSE
          ibinm = ibin - 4
          read(1, *) mag_obs_q(2,ibinm), &
               phi_obs_q(2,ibinm), errp, errm
          phi_obs_err_q(2,ibinm) = errp - errm
       ENDIF
    ENDDO
    phi_obs_q(2,:)     = hconv * phi_obs_q(2,:) 
    phi_obs_err_q(2,:) = hconv * phi_obs_err_q(2,:)
    close(1)

    ! === at z ~ 2 === !
    open(1, file = trim(obsdir)//'Cirasuolo10_z20.txt', &
         status = 'old', iostat = ierr)
    call CheckIerr(ierr, 'K-band LF (z~1): fail to open file')
    DO ibin = 1, 11
       IF(ibin < 4) THEN
          read(1,*)
       ELSE
          ibinm = ibin - 3
          read(1, *) mag_obs_K(3,ibinm), phi_obs_K(3,ibinm,1), &
                     dum1, errm, dum2, errp  
          mag_obs_K(3,ibinm)       = mag_obs_K(3,ibinm) + logh5
          phi_obs_K(3,ibinm,1)     = invh3 * 10.d0 ** phi_obs_K(3,ibinm,1)
          phi_obs_err_K(3,ibinm,1) = invh3 * (10.d0 ** errp - 10.d0 ** errm)
       ENDIF
    ENDDO
    close(1)

    ! --- QSO luminosity function (added by Shirakata)
    open(1, file = trim(obsdir)//'Ueda14_z20.dat', &
         status = 'old', iostat = ierr)
    call CheckIerr(ierr, 'AGN LF: fail to open file')
    DO ibin = 1, 10
       IF(ibin < 5) THEN
          read(1,*)
       ELSE
          ibinm = ibin - 4
          read(1, *) mag_obs_q(3,ibinm), &
               phi_obs_q(3,ibinm), errp, errm
          phi_obs_err_q(3,ibinm) = errp - errm
       ENDIF
    ENDDO
    phi_obs_q(3,:)     = hconv * phi_obs_q(3,:)
    phi_obs_err_q(3,:) = hconv * phi_obs_err_q(3,:)
    close(1)

    ! === at z ~ 3 === !
    open(1, file = trim(obsdir)//'Cirasuolo10_z30.txt', &
         status = 'old', iostat = ierr)
    call CheckIerr(ierr, 'K-band LF (z~1): fail to open file')
    DO ibin = 1, 10
       IF(ibin < 4) THEN
          read(1,*)
       ELSE
          ibinm = ibin - 3
          read(1, *) mag_obs_K(4,ibinm), phi_obs_K(4,ibinm,1), &
                     dum1, errm, dum2, errp  
          mag_obs_K(4,ibinm)       = mag_obs_K(4,ibinm) + logh5
          phi_obs_K(4,ibinm,1)     = invh3 * 10.d0 ** phi_obs_K(4,ibinm,1)
          phi_obs_err_K(4,ibinm,1) = invh3 *(10.d0 ** errp - 10.d0 ** errm)
       ENDIF
    ENDDO
    close(1)
  END SUBROUTINE ReadObsData
!!$=============================================================
  SUBROUTINE InitParams
    ! << STAR FORMATION>>
    IF(param%SFmodel == 1 .or. param%SFmodel == 2) THEN
       IF(mcmc%alpst) THEN
          call random_number(d)
          param%alpst = mcmc%alpst_min &
                        + (mcmc%alpst_max - mcmc%alpst_min) * d
       ENDIF
       IF(mcmc%Vst) THEN
          call random_number(d)
          param%Vst = mcmc%Vst_min + (mcmc%Vst_max - mcmc%Vst_min) * d
       ENDIF
       call random_number(d)
       IF(param%SFmodel == 1 .and. mcmc%tau0st) THEN ! CSF
          param%tau0st = (mcmc%tau0st_min + (mcmc%tau0st_max &
                            - mcmc%tau0st_min) * d) / param%th
       ELSEIF(param%SFmodel == 2 .and. mcmc%eps_SF) THEN ! DSF
          param%eps_SF = (mcmc%eps_SF_min &
                         + (mcmc%eps_SF_max - mcmc%eps_SF_min) * d)
       ENDIF
    ELSE IF(param%SFmodel == 3) THEN
       IF(mcmc%emax) THEN
          call random_number(d)
          param%emax = mcmc%emax_min + (mcmc%emax_max - mcmc%emax_min) * d
       ENDIF
       IF(mcmc%emin) THEN
          call random_number(d)
          param%emin = mcmc%emin_min + (mcmc%emin_max - mcmc%emin_min) * d
          param%emin = 10.d0 ** param%emin
       ENDIF
       IF(mcmc%Zch) THEN
          call random_number(d)
          param%Zch = mcmc%Zch_min + (mcmc%Zch_max - mcmc%Zch_min) * d
       ENDIF
       IF(mcmc%taud_th) THEN
          call random_number(d)
          param%taud_th = mcmc%taud_th_min &
                          + (mcmc%taud_th_max - mcmc%taud_th_min) * d
       ENDIF
    ENDIF

    ! << SN FEEDBACK AND RETURNED GAS >>
    IF(mcmc%Vhot(2)) THEN
       call random_number(d)
       param%Vhot(2) = mcmc%Vhot_min(2) &
                       + (mcmc%Vhot_max(2) - mcmc%Vhot_min(2)) * d
    ENDIF
    IF(mcmc%Vhot(1)) THEN
       IF(mcmc%Vhot_min(1) == mcmc%Vhot_min(2) &
          .and. mcmc%Vhot_max(1) == mcmc%Vhot_max(2)) THEN
          param%Vhot(1) = param%Vhot(2)
       ELSE
          call random_number(d)
          param%Vhot(1) = mcmc%Vhot_min(1) &
                          + (mcmc%Vhot_max(1) - mcmc%Vhot_min(1)) * d
       ENDIF
    ENDIF
    IF(mcmc%alphot(2)) THEN
       call random_number(d)
       param%alphot(2) = mcmc%alphot_min(2) &
                         + (mcmc%alphot_max(2) - mcmc%alphot_min(2)) * d
    ENDIF
    IF(mcmc%alphot(1)) THEN
       IF(mcmc%alphot_min(1) == mcmc%alphot_min(2) &
          .and. mcmc%alphot_max(1) == mcmc%alphot_max(2)) THEN
          param%alphot(1) = param%alphot(2)
       ELSE
          call random_number(d)
          param%alphot(1) = mcmc%alphot_min(1) &
                            + (mcmc%alphot_max(1) - mcmc%alphot_min(1)) * d
       ENDIF
    ENDIF
    IF(mcmc%alp_ret) THEN
       call random_number(d)
       param%alp_ret = mcmc%alp_ret_min + (mcmc%alp_ret_max - mcmc%alp_ret_min) * d
    ENDIF

    ! << MERGERS OF GALAXIES >> 
    IF(mcmc%fmerge) THEN
       call random_number(d)
       param%fmerge = mcmc%fmerge_min + (mcmc%fmerge_max - mcmc%fmerge_min) * d
    ENDIF
    IF(mcmc%fmajor) THEN
       call random_number(d)
       param%fmajor = mcmc%fmajor_min + (mcmc%fmajor_max - mcmc%fmajor_min) * d
    ENDIF
    IF(mcmc%Krem) THEN
       call random_number(d)
       param%Krem = mcmc%Krem_min + (mcmc%Krem_max - mcmc%Krem_min) * d
    ENDIF
    IF(mcmc%fdiss) THEN
       call random_number(d)
       param%fdiss = mcmc%fdiss_min + (mcmc%fdiss_max - mcmc%fdiss_min) * d
    ENDIF
    IF(mcmc%Mh0) THEN
       call random_number(d)
       param%Mh0 = mcmc%Mh0_min + (mcmc%Mh0_max - mcmc%Mh0_min) * d
    ENDIF
    IF(mcmc%alpha_halo) THEN
       call random_number(d)
       param%alpha_halo = mcmc%alpha_halo_min &
                          + (mcmc%alpha_halo_max - mcmc%alpha_halo_min) * d
    ENDIF

    ! << AGN FEEDBACKS >>
    IF(param%AGNFB_key == 1 .and. mcmc%Vcut) THEN
       call random_number(d)
       param%Vcut = mcmc%Vcut_min + (mcmc%Vcut_max - mcmc%Vcut_min) * d
       param%Vcut = 10.d0 ** param%Vcut
    ELSEIF(param%AGNFB_key == 2) THEN
       IF(mcmc%kappa_croton) THEN 
          call random_number(d)
          param%kappa_croton = mcmc%kappa_croton_min &
                               + (mcmc%kappa_croton_max - mcmc%kappa_croton_min) * d
       ENDIF
       IF(mcmc%eta_croton) THEN
          call random_number(d)
          param%eta_croton = mcmc%eta_croton_min &
                             + (mcmc%eta_croton_max - mcmc%eta_croton_min) * d
       ENDIF
    ELSEIF(param%AGNFB_key == 3) THEN
       IF(mcmc%alp_bower) THEN
          call random_number(d)
          param%alp_bower = mcmc%alp_bower_min &
                            + (mcmc%alp_bower_max - mcmc%alp_bower_min) * d
       ENDIF
       IF(mcmc%eps_bower) THEN
          call random_number(d)
          param%eps_bower = mcmc%eps_bower_min &
                            + (mcmc%eps_bower_max - mcmc%eps_bower_min) * d
          param%eps_bower = 10.d0**param%eps_bower
       ENDIF
    ENDIF

    ! << DISC INSTABILITIES >>
    IF(mcmc%em) THEN
       call random_number(d)
       param%em = mcmc%em_min + (mcmc%em_max - mcmc%em_min) * d
    ENDIF
    IF(mcmc%fdi) THEN
       call random_number(d)
       param%fdi = mcmc%fdi_min + (mcmc%fdi_max - mcmc%fdi_min) * d
    ENDIF

    ! << SMBHS AND AGNS >>
    IF(mcmc%fbh(1)) THEN
       call random_number(d)
       param%fbh(1) = mcmc%fbh_min(1) + (mcmc%fbh_max(1) - mcmc%fbh_min(1)) * d 
    ENDIF
    IF(mcmc%fbh(2)) THEN
       call random_number(d)
       param%fbh(2) = mcmc%fbh_min(2) + (mcmc%fbh_max(2) - mcmc%fbh_min(2)) * d 
    ENDIF
    IF(mcmc%Mbhseed) THEN
       call random_number(d)
       param%Mbhseed = mcmc%Mbhseed_min &
                       + (mcmc%Mbhseed_max - mcmc%Mbhseed_min) * d
       param%Mbhseed = 10.d0**param%Mbhseed
    ENDIF
    IF(mcmc%tacc_0) THEN
       call random_number(d)
       param%tacc_0 = mcmc%tacc_0_min &
                      + (mcmc%tacc_0_max - mcmc%tacc_0_min) * d 
    ENDIF
    IF(mcmc%tacc_1) THEN
       call random_number(d)
       param%tacc_1 = mcmc%tacc_1_min &
                      + (mcmc%tacc_1_max - mcmc%tacc_1_min) * d 
    ENDIF
    IF(mcmc%alpha_ad) THEN
       call random_number(d)
       param%alpha_ad = mcmc%alpha_ad_min &
                      + (mcmc%alpha_ad_max - mcmc%alpha_ad_min) * d 
    ENDIF
    IF(mcmc%gamma_ad) THEN
       call random_number(d)
       param%gamma_ad = mcmc%gamma_ad_min &
                      + (mcmc%gamma_ad_max - mcmc%gamma_ad_min) * d 
    ENDIF
    IF(mcmc%Edd_max) THEN
       call random_number(d)
       param%Edd_max = mcmc%Edd_max_min &
                      + (mcmc%Edd_max_max - mcmc%Edd_max_min) * d 
    ENDIF
  END SUBROUTINE InitParams
!!$===========================================================================
  SUBROUTINE SynthesizeParams

    call MPI_Bcast(param_bounds_key, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)

    !--- Synthesize Params. ---
    IF(param%SFmodel == 1 .or. param%SFmodel == 2) THEN
       call MPI_Bcast(param%alpst,     1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)
       call MPI_Bcast(param%Vst,   1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)
       IF(param%SFmodel == 1) THEN !CSF
          call MPI_Bcast(param%tau0st, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)
       ELSEIF(param%SFmodel == 2) THEN !DSF
          call MPI_Bcast(param%eps_SF, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)
       ENDIF 
    ELSE IF(param%SFmodel == 3) THEN
       call MPI_Bcast(param%emin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)
       call MPI_Bcast(param%emax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)
       call MPI_Bcast(param%Zch,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)
       call MPI_Bcast(param%taud_th, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)
    ENDIF

    call MPI_Bcast(param%Vhot(1),   1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)
    call MPI_Bcast(param%Vhot(2),   1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)
    call MPI_Bcast(param%alphot(1), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)
    call MPI_Bcast(param%alphot(2), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)
    call MPI_Bcast(param%alp_ret,   1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)

    call MPI_Bcast(param%fmerge,    1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)
    call MPI_Bcast(param%fmajor,    1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)
    call MPI_Bcast(param%Krem,    1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)
    call MPI_Bcast(param%fdiss,    1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)
    call MPI_Bcast(param%Mh0,    1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)
    call MPI_Bcast(param%alpha_halo,    1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)

    IF(param%AGNFB_key == 1) THEN
       call MPI_Bcast(param%Vcut, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)
    ELSE IF(param%AGNFB_key == 2) THEN
       call MPI_Bcast(param%kappa_croton, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)
       call MPI_Bcast(param%eta_croton,   1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)
    ELSE IF(param%AGNFB_key == 3) THEN
       call MPI_Bcast(param%alp_bower, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)
       call MPI_Bcast(param%eps_bower, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)
    ENDIF

    call MPI_Bcast(param%em,        1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)
    call MPI_Bcast(param%fdi,       1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)

    call MPI_Bcast(param%fbh(1),        1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)
    call MPI_Bcast(param%fbh(2),        1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)
    call MPI_Bcast(param%Mbhseed,   1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)
    call MPI_Bcast(param%tacc_0,        1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)
    call MPI_Bcast(param%tacc_1,        1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)
    call MPI_Bcast(param%alpha_ad,        1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)
    call MPI_Bcast(param%gamma_ad,        1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)
    call MPI_Bcast(param%Edd_max,        1, MPI_DOUBLE, 0, MPI_COMM_WORLD, istat)
  END SUBROUTINE SynthesizeParams
!!!$===========================================================================
  SUBROUTINE RenewalParams
    IF(inode == 0) print '(A)', ' Renewal parameters...'
    itnum = 0

    ! << STAR FORMATION>>
    IF(param%SFmodel == 1 .or. param%SFmodel == 2) THEN
       IF(mcmc%alpst) THEN
          param%alpst = alpst_tmp+gauss_range(mcmc%alpst_min, mcmc%alpst_max)
          DO WHILE (param%alpst < mcmc%alpst_min .or. mcmc%alpst_max < param%alpst)
             param%alpst = alpst_tmp+gauss_range(mcmc%alpst_min, mcmc%alpst_max)
             itnum = itnum + 1
             IF(itnum > 1000) THEN
                print *, 'RenewalParams: itnum exceeds 1000 at alpst ', &
                         param%alpst 
                stop
             ENDIF
          ENDDO
       ENDIF
       itnum = 0
       IF(mcmc%Vst) THEN
          param%Vst = Vst_tmp + gauss_range(mcmc%Vst_min, mcmc%Vst_max)
          DO WHILE (param%Vst < mcmc%Vst_min .or. mcmc%Vst_max < param%Vst)
             param%Vst = Vst_tmp + gauss_range(mcmc%Vst_min, mcmc%Vst_max)
             itnum = itnum + 1
             IF(itnum > 1000) THEN
                print *, 'RenewalParams: itnum exceeds 1000 at Vst ', &
                         param%Vst 
                stop
             ENDIF
          ENDDO
       ENDIF
       IF(param%SFmodel == 1) THEN ! CSF
          itnum = 0
          IF(mcmc%tau0st) THEN
             param%tau0st = tau0st_tmp * param%th &
                            + gauss_range(mcmc%tau0st_min, mcmc%tau0st_max)
             DO WHILE (param%tau0st < mcmc%tau0st_min .or. mcmc%tau0st_max < param%tau0st)
                param%tau0st = tau0st_tmp * param%th &
                               + gauss_range(mcmc%tau0st_min, mcmc%tau0st_max)
                itnum = itnum + 1
                IF(itnum > 1000) THEN
                   print *, 'RenewalParams: itnum exceeds 1000 at tau0st ', &
                            param%tau0st 
                   stop
                ENDIF
             ENDDO
             param%tau0st = param%tau0st / param%th
          ENDIF
       ELSEIF(param%SFmodel == 2) THEN ! DSF
          itnum = 0
          IF(mcmc%eps_SF) THEN
             param%eps_SF = eps_SF_tmp & 
                            + gauss_range(mcmc%eps_SF_min, mcmc%eps_SF_max)
             DO WHILE (param%eps_SF < mcmc%eps_SF_min .or. mcmc%eps_SF_max < param%eps_SF)
                param%eps_SF = eps_SF_tmp & 
                               + gauss_range(mcmc%eps_SF_min, mcmc%eps_SF_max)
                itnum = itnum + 1
                IF(itnum > 1000) THEN
                   print *, 'RenewalParams: itnum exceeds 1000 at eps_SF ', &
                            param%eps_SF 
                   stop
                ENDIF
             ENDDO
          ENDIF
       ENDIF
    ELSE IF(param%SFmodel == 3) THEN
       itnum = 0
       IF(mcmc%emin) THEN
          param%emin = log10(emin_tmp) &
                       + gauss_range(mcmc%emin_min, mcmc%emin_max)
          DO WHILE(param%emin < mcmc%emin_min .or. mcmc%emin_max < param%emin)
             param%emin = log10(emin_tmp) &
                          + gauss_range(mcmc%emin_min, mcmc%emin_max)
             itnum = itnum + 1
             IF(itnum > 1000) THEN
                print *, 'RenewalParams: itnum exceeds 1000 at emin ', &
                         param%emin
                stop
             ENDIF
          ENDDO
          param%emin = 10.0**param%emin
       ENDIF
       itnum = 0
       IF(mcmc%emax) THEN
          param%emax = emax_tmp + gauss_range(mcmc%emax_min, mcmc%emax_max)
          DO WHILE(param%emax < mcmc%emax_min .or. mcmc%emax_max < param%emax)
             param%emax = emax_tmp + gauss_range(mcmc%emax_min, mcmc%emax_max)
             itnum = itnum + 1
             IF(itnum > 1000) THEN
                print *, 'RenewalParams: itnum exceeds 1000 at emax ', &
                         param%emax
                stop
             ENDIF
          ENDDO
       ENDIF
       itnum = 0
       IF(mcmc%Zch) THEN
          param%Zch = Zch_tmp + gauss_range(mcmc%Zch_min, mcmc%Zch_max)
          DO WHILE(param%Zch < mcmc%Zch_min .or. mcmc%Zch_max < param%Zch)
             param%Zch = Zch_tmp + gauss_range(mcmc%Zch_min, mcmc%Zch_max)
             itnum = itnum + 1
             IF(itnum > 1000) THEN
                print *, 'RenewalParams: itnum exceeds 1000 at Zch ', &
                         param%Zch
                stop
             ENDIF
          ENDDO
       ENDIF
       itnum = 0
       IF(mcmc%taud_th) THEN
          param%taud_th = taud_th_tmp &
             + gauss_range(mcmc%taud_th_min, mcmc%taud_th_max)
          DO WHILE(param%taud_th < mcmc%taud_th_min .or. mcmc%taud_th_max < param%taud_th)
             param%taud_th = taud_th_tmp &
                + gauss_range(mcmc%taud_th_min, mcmc%taud_th_max)
             itnum = itnum + 1
             IF(itnum > 1000) THEN
                print *, 'RenewalParams: itnum exceeds 1000 at taud_th ', &
                         param%taud_th
                stop
             ENDIF
          ENDDO
       ENDIF
    ENDIF
    ! << SN FEEDBACK AND RETURNED GAS >> 
    itnum = 0
    IF(mcmc%Vhot(2)) THEN
       param%Vhot(2) = Vhot2_tmp &
          + gauss_range(mcmc%Vhot_min(2), mcmc%Vhot_max(2))
       DO WHILE (param%Vhot(2) < mcmc%Vhot_min(2) .or. mcmc%Vhot_max(2) < param%Vhot(2))
          param%Vhot(2) = Vhot2_tmp &
             + gauss_range(mcmc%Vhot_min(2), mcmc%Vhot_max(2))
          itnum = itnum + 1
          IF(itnum > 1000) THEN
             print *, 'RenewalParams: itnum exceeds 1000 at Vhot(2) ', &
                      param%Vhot(2)
             stop
          ENDIF
       ENDDO
    ENDIF
    itnum = 0
    IF(mcmc%Vhot(1)) THEN
       IF(mcmc%Vhot_min(1) == mcmc%Vhot_min(2) &
          .and. mcmc%Vhot_max(1) == mcmc%Vhot_max(2)) THEN
          param%Vhot(1) = param%Vhot(2) 
          ! use the same Vhot for normal SF gal. and SB gal.
       ELSE
          param%Vhot(1) = Vhot1_tmp &
             + gauss_range(mcmc%Vhot_min(1), mcmc%Vhot_max(1))
          DO WHILE (param%Vhot(1) < mcmc%Vhot_min(1) .or. mcmc%Vhot_max(1) < param%Vhot(1))
             param%Vhot(1) = Vhot1_tmp &
                + gauss_range(mcmc%Vhot_min(1), mcmc%Vhot_max(1))
             itnum = itnum + 1
             IF(itnum > 1000) THEN
                print *, 'RenewalParams: itnum exceeds 1000 at Vhot(1) ', &
                         param%Vhot(1)
                stop
             ENDIF
          ENDDO
       ENDIF
    ENDIF
    itnum = 0
    IF(mcmc%alphot(2)) THEN
       param%alphot(2) = alphot2_tmp &
          + gauss_range(mcmc%alphot_min(2), mcmc%alphot_max(2))
       DO WHILE (param%alphot(2) < mcmc%alphot_min(2) .or. mcmc%alphot_max(2) < param%alphot(2))
          param%alphot(2) = alphot2_tmp &
             + gauss_range(mcmc%alphot_min(2), mcmc%alphot_max(2))
          itnum = itnum + 1
          IF(itnum > 1000) THEN
             print *, 'RenewalParams: itnum exceeds 1000 at alphot(2) ', &
                      param%alphot(2)
             stop
          ENDIF
       ENDDO
    ENDIF
    itnum = 0
    IF(mcmc%alphot(1)) THEN
       IF(mcmc%alphot_min(1) == mcmc%alphot_min(2) &
          .and. mcmc%alphot_max(1) == mcmc%alphot_max(2)) THEN
          param%alphot(1) = param%alphot(2) 
          ! use the same alphot for normal SF gal. and SB gal.
       ELSE
          param%alphot(1) = alphot1_tmp &
             + gauss_range(mcmc%alphot_min(1), mcmc%alphot_max(1))
          DO WHILE (param%alphot(1) < mcmc%alphot_min(1) .or. mcmc%alphot_max(1) < param%alphot(1))
             param%alphot(1) = alphot1_tmp &
                + gauss_range(mcmc%alphot_min(1), mcmc%alphot_max(1))
             itnum = itnum + 1
             IF(itnum > 1000) THEN
                print *, 'RenewalParams: itnum exceeds 1000 at alphot(1) ', &
                         param%alphot(1)
                stop
             ENDIF
          ENDDO
       ENDIF
    ENDIF
    itnum = 0
    IF(mcmc%alp_ret) THEN
       param%alp_ret = alp_ret_tmp + gauss_range(mcmc%alp_ret_min, mcmc%alp_ret_max)
       DO WHILE(param%alp_ret < mcmc%alp_ret_min .or. mcmc%alp_ret_max < param%alp_ret)
          param%alp_ret = alp_ret_tmp + gauss_range(mcmc%alp_ret_min, mcmc%alp_ret_max)
          itnum = itnum + 1
          IF(itnum > 1000) THEN
             print *, 'RenewalParams: itnum exceeds 1000 at alp_ret ', &
                      param%alp_ret
             stop
          ENDIF
       ENDDO
    ENDIF
    ! << MERGERS OF GALAXIES >>
    itnum = 0
    IF(mcmc%fmerge) THEN
       param%fmerge = fmerge_tmp + gauss_range(mcmc%fmerge_min,mcmc%fmerge_max)
       DO WHILE (param%fmerge < mcmc%fmerge_min .or. mcmc%fmerge_max < param%fmerge)
          param%fmerge = fmerge_tmp + gauss_range(mcmc%fmerge_min,mcmc%fmerge_max)
          itnum = itnum + 1
          IF(itnum > 1000) THEN
             print *, 'RenewalParams: itnum exceeds 1000 at fmerge ', &
                      param%fmerge
             stop
          ENDIF
       ENDDO
    ENDIF
    itnum = 0
    IF(mcmc%fmajor) THEN
       param%fmajor = fmajor_tmp + gauss_range(mcmc%fmajor_min,mcmc%fmajor_max)
       DO WHILE (param%fmajor < mcmc%fmajor_min .or. mcmc%fmajor_max < param%fmajor)
          param%fmajor = fmajor_tmp + gauss_range(mcmc%fmajor_min,mcmc%fmajor_max)
          itnum = itnum + 1
          IF(itnum > 1000) THEN
             print *, 'RenewalParams: itnum exceeds 1000 at fmajor ', &
                      param%fmajor
             stop
          ENDIF
       ENDDO
    ENDIF
    itnum = 0
    IF(mcmc%Krem) THEN
       param%Krem = Krem_tmp + gauss_range(mcmc%Krem_min,mcmc%Krem_max)
       DO WHILE (param%Krem < mcmc%Krem_min .or. mcmc%Krem_max < param%Krem)
          param%Krem = Krem_tmp + gauss_range(mcmc%Krem_min,mcmc%Krem_max)
          itnum = itnum + 1
          IF(itnum > 1000) THEN
             print *, 'RenewalParams: itnum exceeds 1000 at Krem ', &
                      param%Krem
             stop
          ENDIF
       ENDDO
    ENDIF
    itnum = 0
    IF(mcmc%fdiss) THEN
       param%fdiss = fdiss_tmp + gauss_range(mcmc%fdiss_min,mcmc%fdiss_max)
       DO WHILE (param%fdiss < mcmc%fdiss_min .or. mcmc%fdiss_max < param%fdiss)
          param%fdiss = fdiss_tmp + gauss_range(mcmc%fdiss_min,mcmc%fdiss_max)
          itnum = itnum + 1
          IF(itnum > 1000) THEN
             print *, 'RenewalParams: itnum exceeds 1000 at fdiss ', &
                      param%fdiss
             stop
          ENDIF
       ENDDO
    ENDIF
    itnum = 0
    IF(mcmc%Mh0) THEN
       param%Mh0 = Mh0_tmp + gauss_range(mcmc%Mh0_min,mcmc%Mh0_max)
       DO WHILE (param%Mh0 < mcmc%Mh0_min .or. mcmc%Mh0_max < param%Mh0)
          param%Mh0 = Mh0_tmp + gauss_range(mcmc%Mh0_min,mcmc%Mh0_max)
          itnum = itnum + 1
          IF(itnum > 1000) THEN
             print *, 'RenewalParams: itnum exceeds 1000 at Mh0 ', &
                      param%Mh0
             stop
          ENDIF
       ENDDO
    ENDIF
    itnum = 0
    IF(mcmc%alpha_halo) THEN
       param%alpha_halo = alpha_halo_tmp + gauss_range(mcmc%alpha_halo_min,mcmc%alpha_halo_max)
       DO WHILE (param%alpha_halo < mcmc%alpha_halo_min .or. mcmc%alpha_halo_max < param%alpha_halo)
          param%alpha_halo = alpha_halo_tmp + gauss_range(mcmc%alpha_halo_min,mcmc%alpha_halo_max)
          itnum = itnum + 1
          IF(itnum > 1000) THEN
             print *, 'RenewalParams: itnum exceeds 1000 at alpha_halo ', &
                      param%alpha_halo
             stop
          ENDIF
       ENDDO
    ENDIF
    ! << AGN FEEDBACKS >>
    IF(param%AGNFB_key == 1 .and. mcmc%Vcut) THEN
       param%Vcut = log10(Vcut_tmp) + gauss_range(mcmc%Vcut_min, mcmc%Vcut_max)
       DO WHILE (param%Vcut < mcmc%Vcut_min .or. mcmc%Vcut_max < param%Vcut)
          param%Vcut = log10(Vcut_tmp) + gauss_range(mcmc%Vcut_min, mcmc%Vcut_max)
          itnum = itnum + 1
          IF(itnum > 1000) THEN
             print *, 'RenewalParams: itnum exceeds 1000 at Vcut ', &
                      param%Vcut
             stop
          ENDIF
       ENDDO
       param%Vcut = 10.d0 ** param%Vcut
    ELSEIF(param%AGNFB_key == 2) THEN
       itnum = 0
       IF(mcmc%kappa_croton) THEN
          param%kappa_croton = kappa_tmp &
                               + gauss_range(mcmc%kappa_croton_min, &
                                 mcmc%kappa_croton_max)
          DO WHILE (param%kappa_croton < mcmc%kappa_croton_min &
                    .or. mcmc%kappa_croton_max < param%kappa_croton)
             param%kappa_croton = kappa_tmp &
                                  + gauss_range(mcmc%kappa_croton_min, &
                                    mcmc%kappa_croton_max)
             itnum = itnum + 1
             IF(itnum > 1000) THEN
                print *, 'RenewalParams: itnum exceeds 1000 at kappa_croton ', &
                         param%kappa_croton
                stop
             ENDIF
          ENDDO
       ENDIF
       itnum = 0
       IF(mcmc%eta_croton) THEN
          param%eta_croton = eta_tmp &
                               + gauss_range(mcmc%eta_croton_min, &
                                 mcmc%eta_croton_max)
          DO WHILE (param%eta_croton < mcmc%eta_croton_min &
                    .or. mcmc%eta_croton_max < param%eta_croton)
             param%eta_croton = eta_tmp &
                                  + gauss_range(mcmc%eta_croton_min, &
                                    mcmc%eta_croton_max)
             itnum = itnum + 1
             IF(itnum > 1000) THEN
                print *, 'RenewalParams: itnum exceeds 1000 at eta_croton ', &
                         param%eta_croton
                stop
             ENDIF
          ENDDO
       ENDIF
    ELSE IF(param%AGNFB_key == 3) THEN
       itnum = 0
       IF(mcmc%alp_bower) THEN
          param%alp_bower = alp_bow_tmp &
                            + gauss_range(mcmc%alp_bower_min, mcmc%alp_bower_max)
          DO WHILE (param%alp_bower < mcmc%alp_bower_min &
                    .or. mcmc%alp_bower_max < param%alp_bower)
             param%alp_bower = alp_bow_tmp &
                               + gauss_range(mcmc%alp_bower_min, mcmc%alp_bower_max)
             itnum = itnum + 1
             IF(itnum > 1000) THEN
                print *, 'RenewalParams: itnum exceeds 1000 at alp_bower ', &
                         param%alp_bower
                stop
             ENDIF
          ENDDO
       ENDIF
       itnum = 0
       IF(mcmc%eps_bower) THEN
          param%eps_bower = log10(eps_bow_tmp) &
                            + gauss_range(mcmc%eps_bower_min, mcmc%eps_bower_max)
          DO WHILE (param%eps_bower < mcmc%eps_bower_min &
                    .or. mcmc%eps_bower_max < param%eps_bower)
             param%eps_bower = log10(eps_bow_tmp) &
                               + gauss_range(mcmc%eps_bower_min, mcmc%eps_bower_max)
             itnum = itnum + 1
             IF(itnum > 1000) THEN
                print *, 'RenewalParams: itnum exceeds 1000 at eps_bower ', &
                         param%eps_bower
                stop
             ENDIF
          ENDDO
          param%eps_bower = 10.d0**param%eps_bower
       ENDIF
    ENDIF
    ! << DISC INSTABILITIES >>
    itnum = 0
    IF(mcmc%em) THEN
       param%em = em_tmp + gauss_range(mcmc%em_min,mcmc%em_max)
       DO WHILE (param%em < mcmc%em_min .or. mcmc%em_max < param%em)
          param%em = em_tmp + gauss_range(mcmc%em_min,mcmc%em_max)
          itnum = itnum + 1
          IF(itnum > 1000) THEN
             print *, 'RenewalParams: itnum exceeds 1000 at em ', &
                      param%em
             stop
          ENDIF
       ENDDO
    ENDIF
    itnum = 0
    IF(mcmc%fdi) THEN
       param%fdi = fdi_tmp + gauss_range(mcmc%fdi_min,mcmc%fdi_max)
       DO WHILE (param%fdi < mcmc%fdi_min .or. mcmc%fdi_max < param%fdi)
          param%fdi = fdi_tmp + gauss_range(mcmc%fdi_min,mcmc%fdi_max)
          itnum = itnum + 1
          IF(itnum > 1000) THEN
             print *, 'RenewalParams: itnum exceeds 1000 at fdi ', &
                      param%fdi
             stop
          ENDIF
       ENDDO
    ENDIF
    ! << SMBHS AND AGNS >>
    itnum = 0
    IF(mcmc%fbh(1)) THEN
       param%fbh(1) = fbh1_tmp + gauss_range(mcmc%fbh_min(1),mcmc%fbh_max(1))
       DO WHILE (param%fbh(1) < mcmc%fbh_min(1) .or. mcmc%fbh_max(1) < param%fbh(1))
          param%fbh(1) = fbh1_tmp + gauss_range(mcmc%fbh_min(1),mcmc%fbh_max(1))
          itnum = itnum + 1
          IF(itnum > 1000) THEN
             print *, 'RenewalParams: itnum exceeds 1000 at fbh(1) ', &
                      param%fbh(1)
             stop
          ENDIF
       ENDDO
    ENDIF
    itnum = 0
    IF(mcmc%fbh(2)) THEN
       param%fbh(2) = fbh2_tmp + gauss_range(mcmc%fbh_min(2),mcmc%fbh_max(2))
       DO WHILE (param%fbh(2) < mcmc%fbh_min(2) .or. mcmc%fbh_max(2) < param%fbh(2))
          param%fbh(2) = fbh2_tmp + gauss_range(mcmc%fbh_min(2),mcmc%fbh_max(2))
          itnum = itnum + 1
          IF(itnum > 1000) THEN
             print *, 'RenewalParams: itnum exceeds 1000 at fbh(2) ', &
                      param%fbh(2)
             stop
          ENDIF
       ENDDO
    ENDIF
    itnum = 0
    IF(mcmc%Mbhseed) THEN
       param%Mbhseed = log10(Mbhseed_tmp) + gauss_range(mcmc%Mbhseed_min,mcmc%Mbhseed_max)
       DO WHILE (param%Mbhseed < mcmc%Mbhseed_min .or. mcmc%Mbhseed_max < param%Mbhseed)
          param%Mbhseed = log10(Mbhseed_tmp) + gauss_range(mcmc%Mbhseed_min,mcmc%Mbhseed_max)
          itnum = itnum + 1
          IF(itnum > 1000) THEN
             print *, 'RenewalParams: itnum exceeds 1000 at Mbhseed ', &
                      param%Mbhseed
             stop
          ENDIF
       ENDDO
       param%Mbhseed = 10.d0 ** param%Mbhseed
    ENDIF
    itnum = 0
    IF(mcmc%tacc_0) THEN
       param%tacc_0 = tacc0_tmp + gauss_range(mcmc%tacc_0_min,mcmc%tacc_0_max)
       DO WHILE (param%tacc_0 < mcmc%tacc_0_min .or. mcmc%tacc_0_max < param%tacc_0)
          param%tacc_0 = tacc0_tmp + gauss_range(mcmc%tacc_0_min,mcmc%tacc_0_max)
          itnum = itnum + 1
          IF(itnum > 1000) THEN
             print *, 'RenewalParams: itnum exceeds 1000 at tacc_0 ', &
                      param%tacc_0
             stop
          ENDIF
       ENDDO
    ENDIF
    itnum = 0
    IF(mcmc%tacc_1) THEN
       param%tacc_1 = tacc1_tmp + gauss_range(mcmc%tacc_1_min,mcmc%tacc_1_max)
       DO WHILE (param%tacc_1 < mcmc%tacc_1_min .or. mcmc%tacc_1_max < param%tacc_1)
          param%tacc_1 = tacc1_tmp + gauss_range(mcmc%tacc_1_min,mcmc%tacc_1_max)
          itnum = itnum + 1
          IF(itnum > 1000) THEN
             print *, 'RenewalParams: itnum exceeds 1000 at tacc_1 ', &
                      param%tacc_1
             stop
          ENDIF
       ENDDO
    ENDIF
    itnum = 0
    IF(mcmc%alpha_ad) THEN
       param%alpha_ad = alpha_ad_tmp &
                        + gauss_range(mcmc%alpha_ad_min,mcmc%alpha_ad_max)
       DO WHILE (param%alpha_ad < mcmc%alpha_ad_min &
                  .or. mcmc%alpha_ad_max < param%alpha_ad)
          param%alpha_ad = alpha_ad_tmp &
                           + gauss_range(mcmc%alpha_ad_min,mcmc%alpha_ad_max)
          itnum = itnum + 1
          IF(itnum > 1000) THEN
             print *, 'RenewalParams: itnum exceeds 1000 at alpha_ad ', &
                      param%alpha_ad
             stop
          ENDIF
       ENDDO
    ENDIF
    itnum = 0
    IF(mcmc%gamma_ad) THEN
       param%gamma_ad = gamma_ad_tmp &
                        + gauss_range(mcmc%gamma_ad_min,mcmc%gamma_ad_max)
       DO WHILE (param%gamma_ad < mcmc%gamma_ad_min &
                  .or. mcmc%gamma_ad_max < param%gamma_ad)
          param%gamma_ad = gamma_ad_tmp &
                           + gauss_range(mcmc%gamma_ad_min,mcmc%gamma_ad_max)
          itnum = itnum + 1
          IF(itnum > 1000) THEN
             print *, 'RenewalParams: itnum exceeds 1000 at alpha_ad ', &
                      param%alpha_ad
             stop
          ENDIF
       ENDDO
    ENDIF
    itnum = 0
    IF(mcmc%Edd_max) THEN
       param%Edd_max = Edd_max_tmp &
                        + gauss_range(mcmc%Edd_max_min,mcmc%Edd_max_max)
       DO WHILE (param%Edd_max < mcmc%Edd_max_min &
                  .or. mcmc%Edd_max_max < param%Edd_max)
          param%Edd_max = Edd_max_tmp &
                           + gauss_range(mcmc%Edd_max_min,mcmc%Edd_max_max)
          itnum = itnum + 1
          IF(itnum > 1000) THEN
             print *, 'RenewalParams: itnum exceeds 1000 at Edd_max ', &
                      param%Edd_max
             stop
          ENDIF
       ENDDO
    ENDIF
    IF(inode == 0) print '(A)', '--done.'
  END SUBROUTINE RenewalParams
!!$=========================================================================
!  SUBROUTINE RenewalParams_mean
!    DOUBLE PRECISION :: means(9), sigmas(9)
!    !! 1: alpst 2: eps_SF 3: alphot 4: Vhot
!    !! 5: fbh 6: fburst 7: alp_gas 8,9: AGNFB related
!    IF(inode == 0) print '(A)', 'Renewal parameters...'
!
!    IF(param%AGNFB_key == 1) THEN
!       means  = [-1.141, 0.216, 3.075, 131.755, 0.005, 150.0, 0.5, 0.00, 0.00]
!       sigmas = [0.184,  0.018, 0.139,   4.135, 0.001,    5.000, 0.01, 0.00, 0.00]
!    ELSE IF(param%AGNFB_key == 2) THEN
!      ! not fixed now (Shirakata 2016/07/11)
!       means  = [-1.141, 0.216, 3.075, 131.755, 9.183, -0.552, 5.263, 0.00, 0.00]
!       sigmas = [0.184,  0.018, 0.139,   4.135, 0.454,  0.363, 0.405, 0.00, 0.00]
!    ELSE IF(param%AGNFB_key == 3) THEN
!       means  = [-1.141, 0.216, 3.075, 131.755, 9.183, -0.552, 5.263, 0.00, 0.00]
!       sigmas = [0.184,  0.018, 0.139,   4.135, 0.454,  0.363, 0.405, 0.00, 0.00]
!    ENDIF    
!
!    param%alpst = means(1)+gauss(sigmas(1))
!
!    param%eps_SF = means(2) + gauss(sigmas(2))
! 
!    param%alphot(2) = means(3) + gauss(sigmas(3))
!    param%alphot(1) = param%alphot(2) ! use the same alphot for normal SF gal.
!                                      !   and SB gal.
!    param%Vhot(2) = means(4) + gauss(sigmas(4))
!    param%Vhot(1) = param%Vhot(2) ! use the same Vhot for normal SF gal. and SB gal.
!    !param%fbh = means(5) + gauss(sigmas(5))
!    !param%fburst = means(7) + gauss(sigmas(7))
!
!    IF(param%AGNFB_key == 1) THEN
!      param%Vcut = means(8) + gauss(sigmas(8))
!!      param%alp_cut = means(9) + gauss(sigmas(9))
!    ELSE IF(param%AGNFB_key == 2) THEN
!      param%kappa_croton = means(8) + gauss(sigmas(8))
!      param%eta_croton = means(9) + gauss(sigmas(9))
!    ELSE IF(param%AGNFB_key == 3) THEN 
!      !param%alp_bower = means(8) + gauss(sigmas(8))
!
!      !param%eps_bower = means(9) + gauss(sigmas(9))
!      !param%eps_bower = 10.d0**param%eps_bower
!    ENDIF
!  END SUBROUTINE RenewalParams_mean
!!$=========================================================================
  SUBROUTINE CalcChisq(inv_V, dof)
    DOUBLE PRECISION, INTENT(IN) :: inv_V
    INTEGER, INTENT(INOUT) :: dof
    INTEGER :: ii, k, dof_tmp
    DOUBLE PRECISION :: mag_data, mass_data, invstep, chisq0, tmp
    INTEGER :: nz, itnum
    DOUBLE PRECISION :: Square ! function

    !!! calc chi-sq !!!
    chisqtot = 0.d0; chisq(:) = 0.d0
    chisq0 = Square(0.08d0)
    phi_tmp = 0.d0;  phi_err_tmp = 0.d0
    dof_tmp = 0
    itnum = 0

    invstep = 1.d0 / stepLF

    nz = 1 ! z ~ 0
    !--- K-band LF
    DO ii = 2, 16
       IF(phi_obs_K(1,ii,1) <= 0.d0) phi_obs_K(1,ii,1) = 1.d-20
       j = (mag_obs_K(1,ii) - baseLF) / stepLF
       mag_data = baseLF + stepLF * dble(j)
       lf_j1 = (lf(nz)%n(4,1,j+1) + lf(nz)%n(4,2,j+1) + lf(nz)%n(4,3,j+1)) * invstep
       lf_j  = (lf(nz)%n(4,1,j)   + lf(nz)%n(4,2,j)   + lf(nz)%n(4,3,j))   * invstep
       phi_tmp = (lf_j1 - lf_j) * invstep * (mag_obs_K(1,ii) - mag_data) + lf_j
       IF(phi_tmp <= inv_V) phi_tmp = inv_V
       IF(phi_obs_K(1,ii,1) > inv_V * invstep) THEN 
          dof_tmp = dof_tmp + 1
          chisq(1) = chisq(1) &
               + Square(phi_obs_K(1,ii,1) - phi_tmp) &
                 / (Square(phi_obs_err_K(1,ii,1)) &
                    + Square(sqrt(phi_tmp*inv_V*invstep)) &
                    + Square(0.2*phi_obs_K(1,ii,1)))
       ENDIF
    ENDDO

    !--- r-band LF
    DO ii = 2, 18
       IF(phi_obs_r(ii,1) <= 0.d0) phi_obs_r(ii,1) = 1.d-20
       j = (mag_obs_r(ii) - baseLF) / stepLF
       mag_data = baseLF + stepLF * dble(j)
       lf_j1 = (lf(nz)%n(3,1,j+1) + lf(nz)%n(3,2,j+1) + lf(nz)%n(3,3,j+1)) * invstep
       lf_j  = (lf(nz)%n(3,1,j)   + lf(nz)%n(3,2,j)   + lf(nz)%n(3,3,j))   * invstep
       phi_tmp = (lf_j1 - lf_j) * invstep * (mag_obs_r(ii) - mag_data) + lf_j
       IF(phi_tmp <= inv_V) phi_tmp = inv_V
       IF(phi_obs_r(ii,1) > inv_V * invstep) THEN 
          dof_tmp = dof_tmp + 1
          chisq(2) = chisq(2) &
               + Square(phi_obs_r(ii,1) - phi_tmp) &
                 / (Square(phi_obs_err_r(ii,1)) &
                    + Square(sqrt(phi_tmp*inv_V*invstep)) &
                    + Square(0.2*phi_obs_r(ii,1)))
       ENDIF
    ENDDO

    invstep = 1.d0 / stepMF
    !--- HI mass function
     DO ii = 8, 23
        IF(MHI_obs(ii) <= 0.d0) MHI_obs(ii) = 1.d-20
        j = (MHI_obs(ii) - 0.d0) * invstep
        mass_data = stepMF * dble(j)
        mf_j1 = (mf(nz)%n(10,1,j+1) + mf(nz)%n(10,2,j+1) + mf(nz)%n(10,3,j+1)) * invstep
        mf_j  = (mf(nz)%n(10,1,j)   + mf(nz)%n(10,2,j)   + mf(nz)%n(10,3,j))   * invstep
        phi_tmp = (mf_j1 - mf_j) / stepMF * (MHI_obs(ii) - mass_data) + mf_j
        IF(phi_tmp <= inv_V) phi_tmp = inv_V
        IF(phi_obs_MHI(ii) > inv_V * invstep) THEN 
           dof_tmp = dof_tmp + 1
           chisq(3) = chisq(3) &
                + Square(phi_obs_MHI(ii) - phi_tmp) &
                  / (Square(phi_obs_err_MHI(ii)) &
                     + Square(sqrt(phi_tmp*inv_V*invstep)) &
                     + Square(0.2*phi_obs_MHI(ii)))
        ENDIF
     ENDDO

     DO j = 30, 50 ! log(Mbulge[Msun]) = [9,13]
     ! --- Fitting function: log(Mbh/Msun) = 1.16 log(Mbulge/Msun) - 4.12
     ! --- Based on Eq. 11 in Kormendy & Ho (2013).
     ! --- Intrinsic scatter is 0.29 dex. (Shirakata 2018/06/19)
       IF(MbhMbulge_xn_all(1,4,j) > 0.d0) THEN
          dof_tmp = dof_tmp + 1
          obs_tmp = 1.16d0 * MbhMbulge_bin(j) - 4.12d0
          chisq(4) = chisq(4) &
                   + Square(obs_tmp - MbhMbulge_xn_all(1,4,j)) &
                      / (Square(0.29d0) &
                         + Square(MbhMbulge_xxn_all(1,4,j)) &
                         + Square(0.2*obs_tmp))
       ENDIF
     ENDDO

    !!--- M_B vs M/L
    ! IF(param%ML == .true.) THEN
    !   DO ii = 9, 15
    !      ML_model = ML_med(1,ii)
    !      IF(ML_model < 1.d-20) ML_model = 1.d-20
    !      ML_model = log10(ML_model)
    !      ML_obs = log10(ML_obs_med(ii))
    !      chisq(5) = chisq(5) + Square(ML_model - ML_obs) / Square(0.3d0)
    !   ENDDO
    ! ENDIF

     DO j = 4, 8 ! log(Vd) = [1.8,2.6]
     ! --- Fitting function: log(Rd [kpc]) = 1.10 log(Vd [km/s]) - 1.93
     ! --- Based on Courteau et al. (2007).
     ! --- Intrinsic scatter is 0.18 dex. (Shirakata 2018/06/19)
       IF(DiskScale_xn_all(1,3,j) > 0.d0) THEN
          dof_tmp = dof_tmp + 1
          obs_tmp = 1.1d0 * DiskScale_bin(j) - 1.93d0
          chisq(6) = chisq(6) &
                   + Square(obs_tmp - DiskScale_xn_all(1,3,j)) &
                      / (Square(0.18d0) &
                         + Square(DiskScale_xxn_all(1,3,j)) &
                         + Square(0.2*obs_tmp))
       ENDIF
     ENDDO

     DO j = 12, 24 ! Mk-5logh = [-24,-18]
     ! --- Fitting function: 2log(Vd [km/s]) + log(Rd [kpc/h]) = 
     ! ---- -0.48 * (MK - 5logh) - 5.23 
     ! --- Based on Forbes et al. (2008) and The fitting is done by Shirakata
     ! --- Intrinsic scatter is 0.250 dex. (Shirakata 2018/06/19)
       IF(SphScale_xn_all(1,3,j) > 0.d0) THEN ! E
          dof_tmp = dof_tmp + 1
          obs_tmp = -0.48d0 * SphScale_bin(j) - 5.23d0
          chisq(7) = chisq(7) &
                   + Square(obs_tmp - SphScale_xn_all(1,3,j)) &
                      / (Square(0.25d0) &
                         + Square(SphScale_xxn_all(1,3,j)) &
                         + Square(0.2*obs_tmp))
       ENDIF
     ENDDO

    !!--- cosmic SFR density
     DO ii = 9, 49, 5
       CSFR_model = log10(CSFR_all(2,ii) * 2.92d0 / 1.d+6 * param%munit / param%th_yr)
       tmp = CSFR_all(1,ii) - 1.d0
!       CSFR_obs = log10((0.01d0 + 0.1d0 * tmp) &
!                        / (1.d0 + (tmp / 3.3d0)**5.3d0) * param%h)
       CSFR_obs = log10((0.01d0 + 0.18d0 * tmp) &
                        / (1.d0 + (tmp / 3.3d0)**4.3d0) * param%h)
       IF (tmp < 4.d0) THEN
          dof_tmp = dof_tmp + 1
          chisq(8) = chisq(8) + Square(CSFR_model - CSFR_obs) / Square(0.1d0)
       ENDIF
     ENDDO

    ! --- z ~ 0.4
    nz = 8 
    !--- hard X-ray QLF
    invstep = 1.d0 / stepLF
    DO ii = 2, 7
       IF(phi_obs_q(1,ii) <= 0.d0) phi_obs_q(1,ii) = 1.d-20
       j = (mag_obs_q(1,ii) - baseLF) / stepLF
       mag_data = baseLF + stepLF * dble(j)
       lf_j1 = lf_q(nz)%n(4,1,j+1) * invstep
       lf_j  = lf_q(nz)%n(4,1,j)   * invstep
       phi_tmp = (lf_j1 - lf_j) / stepLF * (mag_obs_q(1,ii) - mag_data) + lf_j
       IF(phi_tmp <= inv_V) phi_tmp = inv_V
       IF(phi_obs_q(1,ii) > inv_V * invstep) THEN 
          dof_tmp = dof_tmp + 1
          chisq(9) = chisq(9) &
               + Square(phi_obs_q(1,ii) - phi_tmp) &
                 / (Square(phi_obs_err_q(1,ii)) &
                    + Square(sqrt(phi_tmp*inv_V*invstep)) &
                    + Square(0.2*phi_obs_q(1,ii)))
       ENDIF
    ENDDO
    ! --- z ~ 1.0
    nz = 17
    !--- K-band LF
    DO ii = 2, 11
       IF(phi_obs_K(2,ii,1) <= 0.d0) phi_obs_K(2,ii,1) = 1.d-20
       j = (mag_obs_K(2,ii) - baseLF) / stepLF
       mag_data = baseLF + stepLF * dble(j)
       lf_j1 = (lf(nz)%n(4,1,j+1) + lf(nz)%n(4,2,j+1) + lf(nz)%n(4,3,j+1)) * invstep
       lf_j  = (lf(nz)%n(4,1,j)   + lf(nz)%n(4,2,j)   + lf(nz)%n(4,3,j))   * invstep
       phi_tmp = (lf_j1 - lf_j) * invstep * (mag_obs_K(2,ii) - mag_data) + lf_j
       IF(phi_tmp <= inv_V) phi_tmp = inv_V
       IF(phi_obs_K(2,ii,1) > inv_V * invstep) THEN 
          dof_tmp = dof_tmp + 1
          chisq(10) = chisq(10) &
               + Square(phi_obs_K(2,ii,1) - phi_tmp) &
                 / (Square(phi_obs_err_K(2,ii,1)) &
                    + Square(sqrt(phi_tmp*inv_V*invstep)) &
                    + Square(0.2*phi_obs_K(2,ii,1)))
       ENDIF
    ENDDO

    !--- hard X-ray QLF
    invstep = 1.d0 / stepLF
    DO ii = 2, 8
       IF(phi_obs_q(2,ii) <= 0.d0) phi_obs_q(2,ii) = 1.d-20
       j = (mag_obs_q(2,ii) - baseLF) / stepLF
       mag_data = baseLF + stepLF * dble(j)
       lf_j1 = lf_q(nz)%n(4,1,j+1) * invstep
       lf_j  = lf_q(nz)%n(4,1,j)   * invstep
       phi_tmp = (lf_j1 - lf_j) / stepLF * (mag_obs_q(1,ii) - mag_data) + lf_j
       IF(phi_tmp <= inv_V) phi_tmp = inv_V
       IF(phi_obs_q(2,ii) > inv_V * invstep) THEN 
          dof_tmp = dof_tmp + 1
          chisq(11) = chisq(11) &
               + Square(phi_obs_q(2,ii) - phi_tmp) &
                 / (Square(phi_obs_err_q(2,ii)) &
                    + Square(sqrt(phi_tmp*inv_V*invstep)) &
                    + Square(0.2*phi_obs_q(2,ii)))
       ENDIF
    ENDDO

    ! --- z ~ 2.0
    nz = 26
    !--- K-band LF
    DO ii = 1, 9
       IF(phi_obs_K(3,ii,1) <= 0.d0) phi_obs_K(3,ii,1) = 1.d-20
       j = (mag_obs_K(3,ii) - baseLF) / stepLF
       mag_data = baseLF + stepLF * dble(j)
       lf_j1 = (lf(nz)%n(4,1,j+1) + lf(nz)%n(4,2,j+1) + lf(nz)%n(4,3,j+1)) * invstep
       lf_j  = (lf(nz)%n(4,1,j)   + lf(nz)%n(4,2,j)   + lf(nz)%n(4,3,j))   * invstep
       phi_tmp = (lf_j1 - lf_j) * invstep * (mag_obs_K(3,ii) - mag_data) + lf_j
       IF(phi_tmp <= inv_V) phi_tmp = inv_V
       IF(phi_obs_K(3,ii,1) > inv_V * invstep) THEN 
          dof_tmp = dof_tmp + 1
          chisq(12) = chisq(12) &
               + Square(phi_obs_K(3,ii,1) - phi_tmp) &
                 / (Square(phi_obs_err_K(3,ii,1)) &
                    + Square(sqrt(phi_tmp*inv_V*invstep)) &
                    + Square(0.2*phi_obs_K(3,ii,1)))
       ENDIF
    ENDDO

    !--- hard X-ray QLF
    invstep = 1.d0 / stepLF
    DO ii = 1, 7
       IF(phi_obs_q(3,ii) <= 0.d0) phi_obs_q(3,ii) = 1.d-20
       j = (mag_obs_q(3,ii) - baseLF) / stepLF
       mag_data = baseLF + stepLF * dble(j)
       lf_j1 = lf_q(nz)%n(4,1,j+1) * invstep
       lf_j  = lf_q(nz)%n(4,1,j)   * invstep
       phi_tmp = (lf_j1 - lf_j) / stepLF * (mag_obs_q(1,ii) - mag_data) + lf_j
       IF(phi_tmp <= inv_V) phi_tmp = inv_V
       IF(phi_obs_q(3,ii) > inv_V * invstep) THEN 
          dof_tmp = dof_tmp + 1
          chisq(13) = chisq(13) &
               + Square(phi_obs_q(3,ii) - phi_tmp) &
                 / (Square(phi_obs_err_q(3,ii)) &
                    + Square(sqrt(phi_tmp*inv_V*invstep)) &
                    + Square(0.2*phi_obs_q(3,ii)))
       ENDIF
    ENDDO

    ! --- z ~ 3.0
    nz = 32
    !--- K-band LF
    DO ii = 1, 8
       IF(phi_obs_K(4,ii,1) <= 0.d0) phi_obs_K(4,ii,1) = 1.d-20
       j = (mag_obs_K(4,ii) - baseLF) / stepLF
       mag_data = baseLF + stepLF * dble(j)
       lf_j1 = (lf(nz)%n(4,1,j+1) + lf(nz)%n(4,2,j+1) + lf(nz)%n(4,3,j+1)) * invstep
       lf_j  = (lf(nz)%n(4,1,j)   + lf(nz)%n(4,2,j)   + lf(nz)%n(4,3,j))   * invstep
       phi_tmp = (lf_j1 - lf_j) * invstep * (mag_obs_K(4,ii) - mag_data) + lf_j
       IF(phi_tmp <= inv_V) phi_tmp = inv_V
       IF(phi_obs_K(4,ii,1) > inv_V * invstep) THEN 
          dof_tmp = dof_tmp + 1
          chisq(14) = chisq(14) &
               + Square(phi_obs_K(4,ii,1) - phi_tmp) &
                 / (Square(phi_obs_err_K(4,ii,1)) &
                    + Square(sqrt(phi_tmp*inv_V*invstep)) &
                    + Square(0.2*phi_obs_K(4,ii,1)))
       ENDIF
    ENDDO

    DO k=1, NObsProp
       chisqtot = chisqtot + chisq(k)
    ENDDO
    IF(inode == 0) THEN 
       print *, "chisq at z=0 (1)KLF (2)rLF, (3)HIMF, (4) CSFR (5)Mbh-Mbulge, (6)Vd-Rd, (7)Bulge FP)"
       print *, "chisq at z>0 (8-10)KLF z = 1, 2, 3 (11-13) AGNLF z = 0.4, 1, 2"
       print *, chisq(1),chisq(2),chisq(3),chisq(8)
       print *, chisq(4),chisq(6),chisq(7), chisq(10)
       print *, chisq(12),chisq(14),chisq(9),chisq(11),chisq(13)
       print *, "chisq : ", chisqtot
       IF(dof == 0) dof = dof_tmp
    ENDIF
  END SUBROUTINE CalcChisq
!!$============================================================================
  REAL FUNCTION gauss(dev)
    DOUBLE PRECISION :: dev, d1, d2

    call random_number(d1)
    call random_number(d2)
    IF(d1 == 0.d0) d1 = 1.d-10
    gauss = dev * sqrt(-2.d0 * log(d1)) * sin(2.d0 * 3.141592d0 * d2)
  END FUNCTION gauss
!!$============================================================================
  REAL FUNCTION gauss_range(pmin, pmax)
    DOUBLE PRECISION :: pmin, pmax, dev, d1, d2

    call random_number(d1)
    call random_number(d2)
    dev = (pmax - pmin) * 1.d-2 ! 1% of the parameter range
                                ! Shirakata (2018/Jun/26)
    IF(d1 == 0.d0) d1 = 1.d-10
    gauss_range = dev * sqrt(-2.d0 * log(d1)) * sin(2.d0 * 3.141592d0 * d2)
  END FUNCTION gauss_range
!$$=============================================================================
END module MCMCrelated
