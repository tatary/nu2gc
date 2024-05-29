!!$ ========================================================================
! --- << NOTE >> --- !
! When you change the input cooling file, please check 
! read_clcool.f90, moduleCLCOOL.f90, and global_var.f90
! whether the array structure should change.
! (Added by Shirakata 18/Aug/20)
!!$ ========================================================================
module CLCOOLrelated
  TYPE ParameterForCoolFNTO
     DOUBLE PRECISION :: Log_T_min = 1.d0
     DOUBLE PRECISION :: Log_T_max = 9.d0
     DOUBLE PRECISION :: D_Log_T = 0.25d0

     DOUBLE PRECISION :: Log_NH_min = -6.d0
     DOUBLE PRECISION :: Log_NH_max = 6.d0
     DOUBLE PRECISION :: D_Log_NH = 1.d0

     INTEGER :: N_VALUE = 34
  END TYPE ParameterForCoolFNTO
  TYPE(ParameterForCoolFNTO) :: paramCLCOOL

  CHARACTER(LEN=100) :: emsg
  INTEGER :: N_CLCOOL_RED_BIN = 30 
  INTEGER :: N_CLCOOL_MET_BIN = 8 
  INTEGER :: N_CLCOOL_DEN_BIN = 13 
  INTEGER :: N_CLCOOL_TEM_BIN = 33 
  INTEGER :: N_CLCOOL_VAR_BIN = 5 

  DOUBLE PRECISION :: cool_redshift(30) = (/0.d0, 0.25d0, 0.5d0, 0.75d0, &
    1.d0, 1.25d0, 1.5d0, 1.75d0, 2.d0, 2.25d0, 2.5d0, 2.75d0, 3.d0, 3.25d0, &
    3.5d0, 3.75d0, 4.d0, 4.5d0, 5.d0, 5.5d0, 6.d0, 7.d0, &
    8.d0, 9.d0, 10.d0, 11.d0, 12.d0, 13.d0, 14.d0, 15.d0/)
  DOUBLE PRECISION :: cool_metal(8) = (/-10.d0, -3.d0, -2.d0, -1.d0, &
    -0.5d0, 0.d0, 0.5d0, 1.d0/)
 DOUBLE PRECISION :: cool_log_NH(13), cool_log_T(33)

  !! array for the cooling function
  DOUBLE PRECISION, DIMENSION(31,8,13,33,5) :: clcool_table 
                    ! z (N_CLCOOL_RED_BIN+1), Z, NH, T, variables
                    ! variables: log10(cooling rate), log10(heating rate), log10(u erg/g), 
                    ! mean  molecular weitht, log10(ne)
  DOUBLE PRECISION, DIMENSION(8,13,33,5) :: clcool_table_z
!!$ ========================================================================
   CONTAINS
!!$ ========================================================================
   INTEGER FUNCTION CLCOOL_BinarySearch(table, table_size, read_value) RESULT(middle)
      implicit none
      DOUBLE PRECISION, INTENT(in) :: table(:)
      INTEGER, INTENT(in) :: table_size
      DOUBLE PRECISION, INTENT(in) :: read_value
      INTEGER :: first, last

      first = 0; last = table_size - 1; middle = (first + last) / 2

      DO WHILE (first < middle)
         IF(table(middle) < read_value) THEN
            first = middle
         ELSE IF (table(middle) > read_value) THEN
            last = middle
         ELSE
           exit 
         ENDIF
         middle = (first + last) / 2
      ENDDO
   END FUNCTION CLCOOL_BinarySearch 
!!$ ========================================================================
END module CLCOOLrelated 
