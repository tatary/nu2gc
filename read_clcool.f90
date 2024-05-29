!!$ ========================================================================
! --- << NOTE >> --- !
! When you change the input cooling file, please check 
! read_clcool.f90, moduleCLCOOL.f90, and global_var.f90
! whether the array structure should change.
! (Added by Shirakata 18/Aug/20)
!!$ ========================================================================
SUBROUTINE ReadCLCoolingData 
  use global_var; use CLCOOLrelated
  implicit none

  INTEGER :: iz, imet, inh, item, ivar
  INTEGER :: index_met, index_nh
  DOUBLE PRECISION :: ztmp
  DOUBLE PRECISION :: read_data(5)

  !! Definition of reading data 
  CHARACTER(LEN=60)  :: data_dir = 'inpdata/coolfn_TO/'
  CHARACTER(LEN=500) :: openfile
  CHARACTER(LEN=100) :: filename
  INTEGER :: uni = 20
!!$ ================================================
  DO iz = 1, N_CLCOOL_RED_BIN
    ztmp = cool_redshift(iz)
    DO imet = 1, N_CLCOOL_MET_BIN
      index_met = imet-1
      DO inh = 1, N_CLCOOL_DEN_BIN
        index_nh = inh-1
        write(filename, '("rahmati_z", i2.2, ".", i2.2, "_Z", i2.2, "_nH", i2.2, ".dat")') &
             int(ztmp), int(100*(ztmp-int(ztmp))), index_met, index_nh
        openfile = trim(data_dir)//trim(filename)
        open(unit=uni, file=openfile, action='read',form='unformatted', access='stream', status='old')

        item = 1
        DO
          read(unit=uni, end=100) read_data
          DO ivar = 1, N_CLCOOL_VAR_BIN
             clcool_table(iz,imet,inh,item,ivar) = read_data(ivar)
          ENDDO
          item = item + 1
          IF(item == paramCLCOOL%N_VALUE) goto 100
        ENDDO
           

100     continue
        close(uni)
      ENDDO ! NH
    ENDDO ! Z
  ENDDO ! redshift

  ! Collisional ionizing equilibrium before reionization
  DO imet = 1, N_CLCOOL_MET_BIN
     index_met = imet - 1
     DO inh = 1, N_CLCOOL_DEN_BIN
        index_nh = inh - 1
        write(filename, '("rahmati_CIE_Z", i2.2, "_nH", i2.2, ".dat")') &
             index_met, index_nh
        openfile = trim(data_dir)//trim(filename)
        open(unit=uni, file=openfile, action='read',form='unformatted', access='stream', status='old')

        item = 1
        DO
          read(unit=uni, end=110) read_data
          DO ivar = 1, N_CLCOOL_VAR_BIN
             clcool_table(N_CLCOOL_RED_BIN + 1,imet,inh,item,ivar) = read_data(ivar)
          ENDDO
          item = item + 1
          IF(item == paramCLCOOL%N_VALUE) goto 110
        ENDDO
110     continue
        close(uni)
     ENDDO
  ENDDO

  call CalcDensBin(cool_log_NH); call CalcTemBin(cool_log_T)
END SUBROUTINE ReadCLCoolingData 
!!$ ========================================================================
SUBROUTINE LoadCLCoolingRedshift(zp1,tmpclcoolz)
   use global_var; use CLCOOLrelated
   implicit none
   DOUBLE PRECISION, INTENT(in) :: zp1
   DOUBLE PRECISION, INTENT(out) :: tmpclcoolz(N_CLCOOL_MET_BIN,N_CLCOOL_DEN_BIN,N_CLCOOL_TEM_BIN,N_CLCOOL_VAR_BIN) 
   ! Z, NH, T, variables 
   INTEGER :: imet, inh, item, ivar
   INTEGER :: izred
   DOUBLE PRECISION :: zred
   DOUBLE PRECISION :: var1, var2
   DOUBLE PRECISION :: z1, z2, dz
   izred = 0
   zred = zp1 - 1.d0

   IF(zred > param%zreion) THEN
      DO imet = 1, N_CLCOOL_MET_BIN   
         DO inh = 1, N_CLCOOL_DEN_BIN
            DO item = 1, N_CLCOOL_TEM_BIN
               DO ivar = 1, N_CLCOOL_VAR_BIN
                  tmpclcoolz(imet, inh, item, ivar) &
                     = clcool_table(N_CLCOOL_RED_BIN+1, imet, inh, item, ivar) 
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ELSE
      IF(izred < 1 .or. izred == N_CLCOOL_RED_BIN + 1) THEN
         izred = CLCOOL_BinarySearch(cool_redshift, N_CLCOOL_RED_BIN, zred)
      ELSE IF(izred < N_CLCOOL_RED_BIN)  THEN
         IF(cool_redshift(izred) > zred .or. zred >= cool_redshift(izred + 1)) THEN 
            izred = CLCOOL_BinarySearch(cool_redshift, N_CLCOOL_RED_BIN, zred)
         ENDIF
      ENDIF

      DO imet = 1, N_CLCOOL_MET_BIN   
         DO inh = 1, N_CLCOOL_DEN_BIN
            DO item = 1, N_CLCOOL_TEM_BIN
               DO ivar = 1, N_CLCOOL_VAR_BIN
                  var1 = clcool_table(izred,imet,inh,item,ivar)
                  var2 = clcool_table(izred+1,imet,inh,item,ivar)
                  z1 = cool_redshift(izred); z2 = cool_redshift(izred+1)
                  dz = zred - z1
                  tmpclcoolz(imet, inh, item, ivar) &
                     = var1 + ((var2 - var1)/(z2 - z1)) * dz
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDIF
END SUBROUTINE LoadCLCoolingRedshift
!!$ ========================================================================
SUBROUTINE CalcHeatCoolRate(Temp, Metal, NH, zp1, tmpclcoolz, heatcoolrate)
   use global_var; use CLCOOLrelated
   implicit none

   DOUBLE PRECISION, INTENT(in) :: Temp, Metal, NH, zp1
   DOUBLE PRECISION, INTENT(IN) :: tmpclcoolz(N_CLCOOL_MET_BIN, N_CLCOOL_DEN_BIN,N_CLCOOL_TEM_BIN, N_CLCOOL_VAR_BIN)
   DOUBLE PRECISION, INTENT(OUT) :: heatcoolrate

   DOUBLE PRECISION :: tmpclcoolzmetnh(N_CLCOOL_TEM_BIN, N_CLCOOL_VAR_BIN)
   INTEGER :: itemp
   DOUBLE PRECISION :: zred
   DOUBLE PRECISION :: log_T
   DOUBLE PRECISION :: var1, var2
   DOUBLE PRECISION :: x1, x2, dx
   DOUBLE PRECISION :: log_lambda, log_heat
   DOUBLE PRECISION :: lambda_cool, ne

   DOUBLE PRECISION :: FindNe ! function 

   call InterpTables(Metal, NH, tmpclcoolz, tmpclcoolzmetnh)  

   log_T = log10(Temp); zred = zp1 - 1.d0

   IF(log_T < paramCLCOOL%Log_T_min) THEN
      log_lambda = tmpclcoolzmetnh(1,1)
      log_heat   = tmpclcoolzmetnh(1,2)
   ELSE IF(log_T > paramCLCOOL%Log_T_max) THEN
      log_lambda = tmpclcoolzmetnh(N_CLCOOL_TEM_BIN,1)
      log_heat   = tmpclcoolzmetnh(N_CLCOOL_TEM_BIN,2)
   ELSE
      itemp = int((log_T - paramCLCOOL%Log_T_min)/paramCLCOOL%D_Log_T)
      x1    = paramCLCOOL%Log_T_min + itemp * paramCLCOOL%D_Log_T
      x2    = paramCLCOOL%Log_T_min + (itemp+1) * paramCLCOOL%D_Log_T
      dx    = Log_T - x1
      var1  = tmpclcoolzmetnh(itemp,1)
      var2  = tmpclcoolzmetnh(itemp+1,1)
      log_lambda = var1 + ((var2 - var1) / (x2 - x1)) * dx
      var1  = tmpclcoolzmetnh(itemp,2)
      var2  = tmpclcoolzmetnh(itemp+1,2)
      log_heat = var1 + ((var2 - var1) / (x2 - x1)) * dx
   ENDIF
   lambda_cool = 10.d0 ** log_lambda

   IF(zred > param%zreion) THEN
      ne = FindNe(log_T, tmpclcoolzmetnh) 
      lambda_cool = lambda_cool + 5.65d-36 * (10.d0**zp1) * (Temp-2.728d0 * zp1) * ne / NH
   ENDIF

   heatcoolrate = log10(abs(10.d0 ** log_heat - lambda_cool))
END SUBROUTINE CalcHeatCoolRate 
!!$ ========================================================================
SUBROUTINE CalcDensBin(nh_table)
   use CLCOOLrelated

   DOUBLE PRECISION, INTENT(out) :: nh_table(N_CLCOOL_DEN_BIN)
   INTEGER :: i
   DO i = 1, N_CLCOOL_DEN_BIN
      nh_table(i) = paramCLCOOL%Log_NH_min + paramCLCOOL%D_Log_NH * i
   ENDDO
END SUBROUTINE CalcDensBin
!!$ ========================================================================
SUBROUTINE CalcTemBin(tem_table)
   use CLCOOLrelated

   DOUBLE PRECISION, INTENT(out) :: tem_table(N_CLCOOL_TEM_BIN)
   INTEGER :: i
   DO i = 1, N_CLCOOL_TEM_BIN
      tem_table(i) = paramCLCOOL%Log_T_min + paramCLCOOL%D_Log_T * i
   ENDDO
END SUBROUTINE CalcTemBin
!!$ ========================================================================
SUBROUTINE InterpTables(met_red, nh_red, tmpclcoolz, tmptable)
   use global_var; use CLCOOLrelated
   implicit none
   DOUBLE PRECISION, INTENT(in) :: met_red, nh_red
   DOUBLE PRECISION, INTENT(in) :: tmpclcoolz(N_CLCOOL_MET_BIN,N_CLCOOL_DEN_BIN,N_CLCOOL_TEM_BIN,N_CLCOOL_VAR_BIN)
   DOUBLE PRECISION, INTENT(out) :: tmptable(N_CLCOOL_TEM_BIN,N_CLCOOL_VAR_BIN)

   INTEGER :: imet, inh, item, ivar
   DOUBLE PRECISION :: Z, log_nh
   DOUBLE PRECISION :: tmpclcoolzmet(N_CLCOOL_DEN_BIN,N_CLCOOL_TEM_BIN,N_CLCOOL_VAR_BIN)
   DOUBLE PRECISION :: var1, var2
   DOUBLE PRECISION :: x1, x2, dx

   ! --- metallicity
   IF(met_red / const%Zsun < 10.d0 ** cool_metal(1)) THEN
      Z = 10.d0 ** cool_metal(1)
   ELSE
      Z = met_red / const%Zsun
   ENDIF

   IF(Z >= 10.d0**cool_metal(N_CLCOOL_MET_BIN)) THEN
      DO inh = 1, N_CLCOOL_DEN_BIN
         DO item = 1, N_CLCOOL_TEM_BIN
            DO ivar = 1, N_CLCOOL_VAR_BIN
               tmpclcoolzmet(inh,item,ivar) &
                  = tmpclcoolz(N_CLCOOL_MET_BIN,inh,item,ivar)
            ENDDO
         ENDDO
      ENDDO
   ELSE IF(Z >= 10.d0**cool_metal(1)) THEN
      imet = CLCOOL_BinarySearch(cool_metal, N_CLCOOL_MET_BIN, log10(Z))
      x1 = 10.d0 ** cool_metal(imet)
      x2 = 10.d0 ** cool_metal(imet+1)
      dx = Z - x1

      DO inh = 1, N_CLCOOL_DEN_BIN
         DO item = 1, N_CLCOOL_TEM_BIN
            DO ivar = 1, N_CLCOOL_VAR_BIN
               var1 = tmpclcoolz(imet,inh,item,ivar)
               var2 = tmpclcoolz(imet+1,inh,item,ivar)
               tmpclcoolzmet(inh,item,ivar) &
                  = var1 + ((var2 - var1) / (x2 - x1)) * dx
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   ! --- density
   log_nh = log10(nh_red)
   IF(log_nh < cool_log_NH(1)) THEN 
      DO item = 1, N_CLCOOL_TEM_BIN
         DO ivar = 1, N_CLCOOL_VAR_BIN
            tmptable(item,ivar) = tmpclcoolzmet(1,item,ivar)
         ENDDO
      ENDDO
   ELSE IF(log_nh >= cool_log_NH(N_CLCOOL_DEN_BIN)) THEN
      DO item = 1, N_CLCOOL_TEM_BIN
         DO ivar = 1, N_CLCOOL_VAR_BIN
            tmptable(item,ivar) = tmpclcoolzmet(N_CLCOOL_DEN_BIN,item,ivar)
         ENDDO
      ENDDO
   ELSE
      inh = CLCOOL_BinarySearch(cool_log_NH, N_CLCOOL_DEN_BIN, log_nh) 
      x1  = cool_log_NH(inh)
      x2  = cool_log_NH(inh+1)
      dx  = log_nh - x1
      DO item = 1, N_CLCOOL_TEM_BIN
         DO ivar = 1, N_CLCOOL_VAR_BIN
            var1 = tmpclcoolzmet(inh,item,ivar)
            var2 = tmpclcoolzmet(inh+1,item,ivar)
            tmptable(item,ivar) = var1 + ((var2 - var1) / (x2 - x1)) * dx
         ENDDO
      ENDDO
   ENDIF
END SUBROUTINE InterpTables
!!$ ========================================================================
DOUBLE PRECISION FUNCTION FindNe(log_T, tmpclcoolzmetnh) RESULT(Ne) 
   use global_var; use CLCOOLrelated
   implicit none

   DOUBLE PRECISION, INTENT(IN) :: log_T
   DOUBLE PRECISION, INTENT(IN) :: tmpclcoolzmetnh(N_CLCOOL_TEM_BIN,N_CLCOOL_VAR_BIN)

   DOUBLE PRECISION :: log_ne
   INTEGER itemp
   DOUBLE PRECISION :: var1, var2, x1, x2, dx

   IF(log_T < paramCLCOOL%Log_T_min) THEN
      log_ne = tmpclcoolzmetnh(1,5)
   ELSE IF(log_T > paramCLCOOL%Log_T_max) THEN
      log_ne = tmpclcoolzmetnh(N_CLCOOL_TEM_BIN,5)
   ELSE
      itemp = int((log_T - paramCLCOOL%Log_T_min)/paramCLCOOL%D_Log_T)
      x1    = paramCLCOOL%Log_T_min + itemp * paramCLCOOL%D_Log_T
      x2    = paramCLCOOL%Log_T_min + (itemp+1) * paramCLCOOL%D_Log_T
      dx    = Log_T - x1
      var1  = tmpclcoolzmetnh(itemp,5)
      var2  = tmpclcoolzmetnh(itemp+1,5)
      log_ne = var1 + ((var2 - var1) / (x2 - x1)) * dx
   ENDIF
   Ne = 10.d0 ** log_ne
END FUNCTION FindNe
!!$ ========================================================================
