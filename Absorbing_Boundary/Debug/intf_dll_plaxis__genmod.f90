        !COMPILER-GENERATED INTERFACE MODULE: Thu Mar 16 09:18:31 2017
        MODULE INTF_DLL_PLAXIS__genmod
          INTERFACE 
            SUBROUTINE INTF_DLL_PLAXIS(DT,MATPROP,DEPSINC,SIGINOUT,STATV&
     &)
              REAL(KIND=8), INTENT(IN) :: DT
              REAL(KIND=8), INTENT(IN) :: MATPROP(26)
              REAL(KIND=8), INTENT(IN) :: DEPSINC(3)
              REAL(KIND=8), INTENT(INOUT) :: SIGINOUT(4)
              REAL(KIND=8), INTENT(INOUT) :: STATV(*)
            END SUBROUTINE INTF_DLL_PLAXIS
          END INTERFACE 
        END MODULE INTF_DLL_PLAXIS__genmod
