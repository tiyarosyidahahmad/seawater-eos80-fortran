program Projek_PDOK_MANTAPJIWA_BANYAKBANGET
    implicit none
    
integer :: pilihan
  
do    
    print *, "Algorithms for computation of fundamental properties of seawater"
    print *, "======================== MAIN MENU ============================="
    print *, "1. Conductivity to Salinity Conversion"
    print *, "2. Salinity to Conductivity Conversion"
    print *, "3. Specific Volume Anomaly and Density Anomaly of Seawater"
    print *, "4. Pressure to Depth Conversion"
    print *, "5. Freezing Point Temperature of Seawater"
    print *, "6. Specific Heat of Seawater"
    print *, "7. Adiabatic Lapse Rate"
    print *, "8. Potential Temperature"
    print *, "9. Sound Speed in Seawater"
    print *, "10. EXIT"
    print *, "================================================================"
    
    print *, "What are you trying to count?"
    read *, pilihan
    
    select case (pilihan)
        case (1)
            call ALGO1_Conductivity_to_Salinity
        case (2)
            call ALGO2_Salinity_to_Conductivity
        case (3)
            call ALGO3_EOS80
        case (4)
            call ALGO4_Pressure_to_Depth
        case (5)
            call ALGO5_Freezing_Point_Seawater
        case (6)
            call ALGO6_SPECIFIC_HEAT
        case (7)
            call ALGO7_ADIABATIC_LAPSE
        case (8)
            call ALGO8_POTENTIAL_TEMPERATURE
        case (9)
            call ALGO9_SOUND_SPEED
        case (10)
            print *, "==========================THANK YOU====================K-11====="
            print *, ' _   _   __  ___   ____  '
            print *, '| | | | |  \ /  | / __ \ '
            print *, '| |_| | | | \/| | | | | |'
            print *, '|  _  | | |   | | | |_| |'
            print *, '|_| |_| |_|   |_|  \___/ '

            exit
        case default
            print *, "invalid"
    end select
print *, ''
print *, 'Click ENTER to go back'
read (*,*)
end do

contains 
!/////////////////////////////PROGRAM 1/////////////////////////////////////!
SUBROUTINE ALGO1_Conductivity_to_Salinity
IMPLICIT NONE

! =====================================================
! Deklarasi variabel utama
! =====================================================
REAL :: R, t, p
REAL :: rt, rp, rtn
REAL :: S, S0, dS
REAL :: sqrtRtn

! =====================================================
! Konstanta rt(t)
! =====================================================
REAL, PARAMETER :: c0 = 0.6766097
REAL, PARAMETER :: c1 = 2.00564E-2
REAL, PARAMETER :: c2 = 1.104259E-4
REAL, PARAMETER :: c3 = -6.9698E-7
REAL, PARAMETER :: c4 = 1.0031E-9

! =====================================================
! Konstanta Rp
! =====================================================
REAL, PARAMETER :: e1 = 2.070E-5
REAL, PARAMETER :: e2 = -6.370E-10
REAL, PARAMETER :: e3 = 3.989E-15

REAL, PARAMETER :: d1 = 3.426E-2
REAL, PARAMETER :: d2 = 4.464E-4
REAL, PARAMETER :: d3 = 4.215E-1
REAL, PARAMETER :: d4 = -3.107E-3

! =====================================================
! Konstanta salinitas
! =====================================================
REAL, PARAMETER :: a0 = 0.0080
REAL, PARAMETER :: a1 = -0.1692
REAL, PARAMETER :: a2 = 25.3851
REAL, PARAMETER :: a3 = 14.0941
REAL, PARAMETER :: a4 = -7.0261
REAL, PARAMETER :: a5 = 2.7081

REAL, PARAMETER :: b0 = 0.0005
REAL, PARAMETER :: b1 = -0.0056
REAL, PARAMETER :: b2 = -0.0066
REAL, PARAMETER :: b3 = -0.0375
REAL, PARAMETER :: b4 = 0.0636
REAL, PARAMETER :: b5 = -0.0144

REAL, PARAMETER :: k  = 0.0162

! =====================================================
! Input dari user
! =====================================================
PRINT *, 'Masukkan conductivity ratio R:'
READ *, R

PRINT *, 'Masukkan suhu t (C):'
READ *, t

PRINT *, 'Masukkan tekanan p (dbar):'
READ *, p

! =====================================================
! Peringatan batas valid
! =====================================================
IF (t .LT. -2.0 .OR. t .GT. 35.0) THEN
    PRINT *, 'PERINGATAN: Suhu di luar range valid (-2 sampai 35 C)'
END IF

! =====================================================
! Hitung rt(t)
! =====================================================
rt = c0 + c1*t + c2*t*t + c3*t*t*t + c4*t*t*t*t

! =====================================================
! Hitung rp(R,t,p)
! =====================================================
rp = 1.0 + ( p*(e1 + e2*p + e3*p*p) ) &
     / ( 1.0 + d1*t + d2*t*t + (d3 + d4*t)*R )

! =====================================================
! Hitung conductivity ratio ternormalisasi
! =====================================================
rtn = R / (rp * rt)

! =====================================================
! Hitung salinitas
! =====================================================
sqrtRtn = SQRT(rtn)

S0 = a0 + a1*sqrtRtn + a2*rtn + a3*rtn*sqrtRtn &
     + a4*rtn*rtn + a5*rtn*rtn*sqrtRtn

dS = ((t - 15.0)/(1.0 + k*(t - 15.0))) * &
     ( b0 + b1*sqrtRtn + b2*rtn + b3*rtn*sqrtRtn &
     + b4*rtn*rtn + b5*rtn*rtn*sqrtRtn )

S = S0 + dS

! =====================================================
! Output
! =====================================================
PRINT *
PRINT *, '---------------------------------------------'
PRINT *, '   Conductivity to Salinity RESULT TABLE'
PRINT *, '---------------------------------------------'
PRINT '(A,F8.2)', ' Conductivity Ratio (R)        : ', R
PRINT '(A,F8.2)', ' Temperature (C)     : ', t
PRINT '(A,F8.2)', ' Pressure (dbar)      : ', p
PRINT *, '---------------------------------------------'
PRINT '(A,F12.3)', ' rp        : ', rp
PRINT '(A,F12.3)', ' rt        : ', rt
PRINT '(A,F12.3)', ' rtn       : ', rtn
PRINT *, '---------------------------------------------'
PRINT '(A,F12.3)', ' Salinity   : ', S
PRINT *, '---------------------------------------------'

END SUBROUTINE ALGO1_Conductivity_to_Salinity


!/////////////////////////////PROGRAM 2/////////////////////////////////////!
subroutine ALGO2_Salinity_to_Conductivity
IMPLICIT NONE

! =====================================================
! Variabel utama
! =====================================================
DOUBLE PRECISION :: S, t, p
DOUBLE PRECISION :: Rt, sqrtRt
DOUBLE PRECISION :: S_calc, dS_dsqrtRt
DOUBLE PRECISION :: R, R0
DOUBLE PRECISION :: rt_t
DOUBLE PRECISION :: A, B, C
DOUBLE PRECISION :: dt
INTEGER :: iter

! =====================================================
! Parameter iterasi
! =====================================================
DOUBLE PRECISION, PARAMETER :: tol = 1.0D-10
INTEGER, PARAMETER :: maxiter = 20

! =====================================================
! Konstanta rt(t)
! =====================================================
DOUBLE PRECISION, PARAMETER :: c0 = 0.6766097D0
DOUBLE PRECISION, PARAMETER :: c1 = 2.00564D-2
DOUBLE PRECISION, PARAMETER :: c2 = 1.104259D-4
DOUBLE PRECISION, PARAMETER :: c3 = -6.9698D-7
DOUBLE PRECISION, PARAMETER :: c4 = 1.0031D-9

! =====================================================
! Konstanta tekanan
! =====================================================
DOUBLE PRECISION, PARAMETER :: e1 = 2.070D-5
DOUBLE PRECISION, PARAMETER :: e2 = -6.370D-10
DOUBLE PRECISION, PARAMETER :: e3 = 3.989D-15

DOUBLE PRECISION, PARAMETER :: d1 = 3.426D-2
DOUBLE PRECISION, PARAMETER :: d2 = 4.464D-4
DOUBLE PRECISION, PARAMETER :: d3 = 4.215D-1
DOUBLE PRECISION, PARAMETER :: d4 = -3.107D-3

! =====================================================
! Konstanta salinitas (PSS-78)
! =====================================================
DOUBLE PRECISION, PARAMETER :: a0 = 0.0080D0
DOUBLE PRECISION, PARAMETER :: a1 = -0.1692D0
DOUBLE PRECISION, PARAMETER :: a2 = 25.3851D0
DOUBLE PRECISION, PARAMETER :: a3 = 14.0941D0
DOUBLE PRECISION, PARAMETER :: a4 = -7.0261D0
DOUBLE PRECISION, PARAMETER :: a5 = 2.7081D0

DOUBLE PRECISION, PARAMETER :: b0 = 0.0005D0
DOUBLE PRECISION, PARAMETER :: b1 = -0.0056D0
DOUBLE PRECISION, PARAMETER :: b2 = -0.0066D0
DOUBLE PRECISION, PARAMETER :: b3 = -0.0375D0
DOUBLE PRECISION, PARAMETER :: b4 = 0.0636D0
DOUBLE PRECISION, PARAMETER :: b5 = -0.0144D0

DOUBLE PRECISION, PARAMETER :: k = 0.0162D0

! =====================================================
! Input
! =====================================================
PRINT *, 'Masukkan salinitas S (PSS-78):'
READ *, S
PRINT *, 'Masukkan suhu t (C):'
READ *, t
PRINT *, 'Masukkan tekanan p (dbar):'
READ *, p

dt = t - 15.0D0

! =====================================================
! Tebakan awal vRt (UNESCO)
! =====================================================
sqrtRt = DSQRT(S / 35.0D0)

! =====================================================
! Newton-Raphson pada vRt
! =====================================================
DO iter = 1, maxiter

   Rt = sqrtRt * sqrtRt

   S_calc = a0 + a1*sqrtRt + a2*Rt + a3*Rt*sqrtRt &
          + a4*Rt*Rt + a5*Rt*Rt*sqrtRt &
          + (dt/(1.0D0 + k*dt)) * &
          ( b0 + b1*sqrtRt + b2*Rt + b3*Rt*sqrtRt &
          + b4*Rt*Rt + b5*Rt*Rt*sqrtRt )

   dS_dsqrtRt = a1 + 2.0D0*a2*sqrtRt + 3.0D0*a3*Rt &
              + 4.0D0*a4*Rt*sqrtRt + 5.0D0*a5*Rt*Rt &
              + (dt/(1.0D0 + k*dt)) * &
              ( b1 + 2.0D0*b2*sqrtRt + 3.0D0*b3*Rt &
              + 4.0D0*b4*Rt*sqrtRt + 5.0D0*b5*Rt*Rt )

   sqrtRt = sqrtRt + (S - S_calc) / dS_dsqrtRt

   IF (DABS(S - S_calc) < tol) EXIT
END DO

Rt = sqrtRt * sqrtRt

! =====================================================
! R pada tekanan nol
! =====================================================
rt_t = c0 + c1*t + c2*t*t + c3*t*t*t + c4*t*t*t*t
R0   = Rt * rt_t

! =====================================================
! Koreksi tekanan
! =====================================================
IF (p .EQ. 0.0D0) THEN
   R = R0
ELSE
   A = e1*p + e2*p*p + e3*p*p*p
   B = 1.0D0 + d1*t + d2*t*t
   C = (d3 + d4*t) * Rt

   R = ( DSQRT((B - A*Rt)**2 + 4.0D0*A*(B + C)) &
       - (B - A*Rt) ) / (2.0D0*A)
END IF

! =====================================================
! Output
! =====================================================
PRINT *
IF (R .LT. 1.0D0) THEN
   WRITE(*,'(A,"0",F7.6)') 'Conductivity ratio R = ', R
ELSE
   WRITE(*,'(A,F9.6)') 'Conductivity ratio R = ', R
END IF

END subroutine ALGO2_Salinity_to_Conductivity


!/////////////////////////////PROGRAM 3/////////////////////////////////////!
subroutine ALGO3_EOS80
IMPLICIT NONE

! ==========================
! INPUT VARIABLES
! ==========================
REAL(8) :: S, t, p

! ==========================
! DENSITY & VOLUME
! ==========================
REAL(8) :: rho, rho0, rho_w
REAL(8) :: V, V350p
REAL(8) :: delta, alpha

! ==========================
! PURE WATER DENSITY COEFF
! ==========================
REAL(8), PARAMETER :: a0 = 999.842594
REAL(8), PARAMETER :: a1 = 6.793952E-2
REAL(8), PARAMETER :: a2 = -9.095290E-3
REAL(8), PARAMETER :: a3 = 1.001685E-4
REAL(8), PARAMETER :: a4 = -1.120083E-6
REAL(8), PARAMETER :: a5 = 6.536332E-9

! ==========================
! SEAWATER DENSITY COEFF
! ==========================
REAL(8), PARAMETER :: b0 = 8.24493E-1
REAL(8), PARAMETER :: b1 = -4.0899E-3
REAL(8), PARAMETER :: b2 = 7.6438E-5
REAL(8), PARAMETER :: b3 = -8.2467E-7
REAL(8), PARAMETER :: b4 = 5.3875E-9

REAL(8), PARAMETER :: c0 = -5.72466E-3
REAL(8), PARAMETER :: c1 = 1.0227E-4
REAL(8), PARAMETER :: c2 = -1.6546E-6

REAL(8), PARAMETER :: d0 = 4.8314E-4

! ==========================
! BULK MODULUS COEFF
! ==========================
REAL(8), PARAMETER :: e0 = 19652.21
REAL(8), PARAMETER :: e1 = 148.4206
REAL(8), PARAMETER :: e2 = -2.327105
REAL(8), PARAMETER :: e3 = 1.360477E-2
REAL(8), PARAMETER :: e4 = -5.155288E-5

REAL(8), PARAMETER :: f0 = 54.6746
REAL(8), PARAMETER :: f1 = -0.603459
REAL(8), PARAMETER :: f2 = 1.09987E-2
REAL(8), PARAMETER :: f3 = -6.1670E-5

REAL(8), PARAMETER :: g0 = 7.944E-2
REAL(8), PARAMETER :: g1 = 1.6483E-2
REAL(8), PARAMETER :: g2 = -5.3009E-4

REAL(8), PARAMETER :: h0 = 3.239908
REAL(8), PARAMETER :: h1 = 1.43713E-3
REAL(8), PARAMETER :: h2 = 1.16092E-4
REAL(8), PARAMETER :: h3 = -5.77905E-7

REAL(8), PARAMETER :: k0 = 8.50935E-5
REAL(8), PARAMETER :: k1 = -6.12293E-6
REAL(8), PARAMETER :: k2 = 5.2787E-8

REAL(8) :: Kw, Aw, Bw, K

! ==========================
! INPUT
! ==========================
PRINT *, 'Enter Salinity (PSU): '
READ *, S
PRINT *, 'Enter Temperature (C): '
READ *, t
PRINT *, 'Enter Pressure (dbar): '
READ *, p

! ==========================
! PURE WATER DENSITY
! ==========================
rho_w = a0 + a1*t + a2*t**2 + a3*t**3 + a4*t**4 + a5*t**5

! ==========================
! SEAWATER DENSITY AT p=0
! ==========================
rho0 = rho_w + &
       (b0 + b1*t + b2*t**2 + b3*t**3 + b4*t**4)*S + &
       (c0 + c1*t + c2*t**2)*S**1.5 + d0*S**2

! ==========================
! BULK MODULUS
! ==========================
Kw = e0 + e1*t + e2*t**2 + e3*t**3 + e4*t**4
Aw = h0 + h1*t + h2*t**2 + h3*t**3
Bw = k0 + k1*t + k2*t**2

K = Kw + &
    (f0 + f1*t + f2*t**2 + f3*t**3)*S + &
    (g0 + g1*t + g2*t**2)*S**1.5 + &
    (Aw + 2.2838E-3*S)*p + &
    Bw*p**2

! ==========================
! FINAL DENSITY
! ==========================
rho = rho0 / (1.0D0 - p/K)

! ==========================
! SPECIFIC VOLUME
! ==========================
V = 1.0D0 / rho

! ==========================
! REFERENCE V(35,0,p)
! ==========================
V350p = 9.7266204E-4 * (1.0D0 - p/(21582.27 + 3.35940552*p + 5.03217E-5*p**2))

! ==========================
! ANOMALIES
! ==========================
delta = (V - V350p) * 1.0D8
alpha = rho - 1000.0D0

! ==========================
! OUTPUT
! ==========================
PRINT *
PRINT *, '============================================='
PRINT *, ' EOS-80 RESULTS'
PRINT *, '============================================='
PRINT '(A,F10.4)', ' Density (kg/m^3)        : ', rho
PRINT '(A,F10.4)', ' Density anomaly (kg/m3): ', alpha
PRINT '(A,E14.6)', ' Specific volume (m3/kg): ', V
PRINT '(A,F12.6)', ' Spec. vol anomaly (x1e-8): ', delta
PRINT *, '============================================='

END subroutine ALGO3_EOS80


!/////////////////////////////PROGRAM 4/////////////////////////////////////!
subroutine ALGO4_Pressure_to_Depth
  implicit none

  ! =========================
  ! Variable declaration
  ! =========================
  real(8) :: p, lat, phi
  real(8) :: g0, gamma_p
  real(8) :: z
  real(8) :: C1, C2, C3, C4
  real(8) :: integral_V

  ! =========================
  ! Constants (EOS-80)
  ! =========================
  C1 = 9.72659D0
  C2 = -2.2512D-5
  C3 = 2.279D-10
  C4 = -1.82D-15

  gamma_p = 2.184D-6   ! m/s^2 per decibar

  ! =========================
  ! Input
  ! =========================
  print *, 'Enter pressure (decibar): '
  read *, p

  print *, 'Enter latitude (degrees): '
  read *, lat

  ! =========================
  ! Gravity calculation
  ! =========================
  phi = lat * 3.141592653589793D0 / 180.0D0

  g0 = 9.780318D0 * &
       (1.0D0 + 5.2788D-3 * sin(phi)**2 &
       + 2.36D-5 * sin(phi)**4)

  ! =========================
  ! Integral of specific volume
  ! =========================
  integral_V = C1*p + C2*p**2 + C3*p**3 + C4*p**4

  ! =========================
  ! Depth calculation
  ! =========================
  z = integral_V / (g0 + gamma_p*p)

  ! =========================
  ! Output
  ! =========================
  print *, '--------------------------------'
  print *, 'Depth (m) = ', z
  print *, '--------------------------------'

end subroutine ALGO4_Pressure_to_Depth


!/////////////////////////////PROGRAM 5/////////////////////////////////////!
subroutine ALGO5_Freezing_Point_Seawater
IMPLICIT NONE

! =====================================================
! Deklarasi variabel
! =====================================================
REAL :: S, p
REAL :: tf
REAL :: S32

! =====================================================
! Konstanta persamaan freezing point
! =====================================================
REAL, PARAMETER :: a0 = -0.0575
REAL, PARAMETER :: a1 =  1.710523E-3
REAL, PARAMETER :: a2 = -2.154996E-4
REAL, PARAMETER :: b  = -7.53E-4

! =====================================================
! Input dari user
! =====================================================
PRINT *, 'Masukkan salinitas S (PSS-78):'
READ *, S

PRINT *, 'Masukkan tekanan p (dbar):'
READ *, p

! =====================================================
! Peringatan batas valid salinitas
! =====================================================
IF (S .LT. 4.0 .OR. S .GT. 40.0) THEN
    PRINT *, 'PERINGATAN: Salinitas di luar range valid (4 - 40)'
END IF

! =====================================================
! Hitung S^(3/2)
! =====================================================
S32 = S * SQRT(S)

! =====================================================
! Hitung freezing point
! =====================================================
tf = a0*S + a1*S32 + a2*S*S + b*p

! =====================================================
! Output dalam bentuk tabel
! =====================================================
PRINT *
WRITE(*,'(A)') '   S        p(dbar)     Freezing Point (C)'
WRITE(*,'(A)') '---------------------------------------------'
WRITE(*,'(F6.1, F12.1, F16.6)') S, p, tf

END subroutine ALGO5_Freezing_Point_Seawater


!/////////////////////////////PROGRAM 6/////////////////////////////////////!
subroutine ALGO6_SPECIFIC_HEAT
IMPLICIT NONE

! ===============================
! DECLARATION
! ===============================
REAL :: S, t, p
REAL :: t2, t3, t4
REAL :: S32
REAL :: Cp0, dCp1, Cp

! ===============================
! CONSTANTS Cp(S,t,0)
! ===============================
REAL, PARAMETER :: cp0_0 = 4217.4
REAL, PARAMETER :: cp0_1 = -3.720283
REAL, PARAMETER :: cp0_2 = 0.1412855
REAL, PARAMETER :: cp0_3 = -2.654387E-3
REAL, PARAMETER :: cp0_4 = 2.093236E-5

REAL, PARAMETER :: As0 = -7.643575
REAL, PARAMETER :: As1 = 0.1072763
REAL, PARAMETER :: As2 = -1.38385E-3

REAL, PARAMETER :: Bs0 = 0.1770383
REAL, PARAMETER :: Bs1 = -4.07718E-3
REAL, PARAMETER :: Bs2 = 5.148E-5

! ===============================
! PRESSURE DEPENDENCE ?1Cp
! ===============================
REAL, PARAMETER :: Ap0 = -4.9592
REAL, PARAMETER :: Ap1 = 1.45747E-2
REAL, PARAMETER :: Ap2 = -3.13885E-4
REAL, PARAMETER :: Ap3 = 2.0357E-6
REAL, PARAMETER :: Ap4 = 1.7168E-8

REAL, PARAMETER :: Bp0 = 2.4931E-4
REAL, PARAMETER :: Bp1 = -1.08645E-5
REAL, PARAMETER :: Bp2 = 2.87533E-7
REAL, PARAMETER :: Bp3 = -4.0027E-9
REAL, PARAMETER :: Bp4 = 2.2956E-11

REAL, PARAMETER :: Cp0p = -5.422E-8
REAL, PARAMETER :: Cp1p = 2.6380E-9
REAL, PARAMETER :: Cp2p = -6.5637E-11
REAL, PARAMETER :: Cp3p = 6.136E-13

! ===============================
! INPUT
! ===============================
PRINT *, '----------------------------------------'
PRINT *, ' SPECIFIC HEAT OF SEAWATER (ALG 6)'
PRINT *, '----------------------------------------'
PRINT *, 'Input Salinity S (psu): '
READ *, S
PRINT *, 'Input Temperature t (C): '
READ *, t
PRINT *, 'Input Pressure p (decibar): '
READ *, p

! ===============================
! PRE-CALCULATION
! ===============================
t2 = t*t
t3 = t2*t
t4 = t3*t
S32 = S*SQRT(S)

! ===============================
! Cp(S,t,0)
! ===============================
Cp0 = (cp0_0 + cp0_1*t + cp0_2*t2 + cp0_3*t3 + cp0_4*t4) &
    + (As0 + As1*t + As2*t2)*S &
    + (Bs0 + Bs1*t + Bs2*t2)*S32

! ===============================
! ?1Cp(0,t,p)
! ===============================
dCp1 = (Ap0 + Ap1*t + Ap2*t2 + Ap3*t3 + Ap4*t4)*p &
     + (Bp0 + Bp1*t + Bp2*t2 + Bp3*t3 + Bp4*t4)*p*p &
     + (Cp0p + Cp1p*t + Cp2p*t2 + Cp3p*t3)*p*p*p

! ===============================
! TOTAL Cp
! ===============================
Cp = Cp0 + dCp1

! ===============================
! OUTPUT TABLE
! ===============================
PRINT *
PRINT *, '---------------------------------------------'
PRINT *, '        SPECIFIC HEAT RESULT TABLE'
PRINT *, '---------------------------------------------'
PRINT '(A,F8.2)', ' Salinity (psu)        : ', S
PRINT '(A,F8.2)', ' Temperature (C)     : ', t
PRINT '(A,F8.2)', ' Pressure (dbar)      : ', p
PRINT *, '---------------------------------------------'
PRINT '(A,F12.3)', ' Cp(S,t,0) [J/kgC]   : ', Cp0
PRINT '(A,F12.3)', ' Cp_pressure         : ', dCp1
PRINT *, '---------------------------------------------'
PRINT '(A,F12.3)', ' Cp(S,t,p) [J/kgC]   : ', Cp
PRINT *, '---------------------------------------------'

END subroutine ALGO6_SPECIFIC_HEAT


!/////////////////////////////PROGRAM 7/////////////////////////////////////!
subroutine ALGO7_ADIABATIC_LAPSE
IMPLICIT NONE

! ==========================
! INPUT VARIABLES
! ==========================
REAL(8) :: S, t, p
REAL(8) :: r

! ==========================
! COEFFICIENTS (Bryden, 1973)
! ==========================
REAL(8), PARAMETER :: a0 = 3.5803E-5
REAL(8), PARAMETER :: a1 = 8.5258E-6
REAL(8), PARAMETER :: a2 = -6.8360E-8
REAL(8), PARAMETER :: a3 = 6.6228E-10

REAL(8), PARAMETER :: b0 = 1.8932E-6
REAL(8), PARAMETER :: b1 = -4.2393E-8

REAL(8), PARAMETER :: c0 = 1.8741E-8
REAL(8), PARAMETER :: c1 = -6.7795E-10
REAL(8), PARAMETER :: c2 = 8.7330E-12
REAL(8), PARAMETER :: c3 = -5.4481E-14

REAL(8), PARAMETER :: d0 = -1.1351E-10
REAL(8), PARAMETER :: d1 = 2.7759E-12

REAL(8), PARAMETER :: e0 = -4.6206E-13
REAL(8), PARAMETER :: e1 = 1.8676E-14
REAL(8), PARAMETER :: e2 = -2.1687E-16

! ==========================
! INPUT
! ==========================
PRINT *, 'Enter Salinity (PSU): '
READ *, S
PRINT *, 'Enter Temperature (C): '
READ *, t
PRINT *, 'Enter Pressure (decibar): '
READ *, p

! ==========================
! ADIABATIC LAPSE RATE
! ==========================
r = a0 + a1*t + a2*t**2 + a3*t**3 &
  + (b0 + b1*t)*(S - 35.0D0) &
  + (c0 + c1*t + c2*t**2 + c3*t**3 &
  + (d0 + d1*t)*(S - 35.0D0))*p &
  + (e0 + e1*t + e2*t**2)*p**2

! ==========================
! OUTPUT
! ==========================
PRINT *
PRINT *, '=========================================='
PRINT *, ' ADIABATIC LAPSE RATE'
PRINT *, '=========================================='
PRINT '(A,ES14.6)', 'r (C/decibar) : ', r
PRINT *, '=========================================='

END subroutine ALGO7_ADIABATIC_LAPSE


!/////////////////////////////PROGRAM 8/////////////////////////////////////!
subroutine ALGO8_POTENTIAL_TEMPERATURE
  IMPLICIT NONE

  REAL :: S, t, p, pr, dp
  REAL :: theta
  REAL :: r1, r2, r3, r4
  REAL :: t1, t2, t3

  ! Koefisien lapse rate
  REAL :: a0,a1,a2,a3,b0,b1
  REAL :: c0,c1,c2,c3,d0,d1
  REAL :: e0,e1,e2

  !---------------------------
  ! Input
  !---------------------------
  PRINT *, 'Enter Salinity (PSU):'
  READ *, S
  PRINT *, 'Enter Temperature (C):'
  READ *, t
  PRINT *, 'Enter Pressure p (decibar):'
  READ *, p
  PRINT *, 'Enter Reference Pressure pr (decibar):'
  READ *, pr

  dp = pr - p
  theta = t

  !---------------------------
  ! Koefisien (Bryden, 1973)
  !---------------------------
  a0 =  3.5803E-5
  a1 =  8.5258E-6
  a2 = -6.8360E-8
  a3 =  6.6228E-10

  b0 =  1.8932E-6
  b1 = -4.2393E-8

  c0 =  1.8741E-8
  c1 = -6.7795E-10
  c2 =  8.7330E-12
  c3 = -5.4481E-14

  d0 = -1.1351E-10
  d1 =  2.7759E-12

  e0 = -4.6206E-13
  e1 =  1.8676E-14
  e2 = -2.1687E-16

  !---------------------------
  ! Runge-Kutta Step 1
  !---------------------------
  r1 = a0 + a1*theta + a2*theta**2 + a3*theta**3 &
     + (b0 + b1*theta)*(S - 35.0) &
     + (c0 + c1*theta + c2*theta**2 + c3*theta**3 &
     + (d0 + d1*theta)*(S - 35.0))*p &
     + (e0 + e1*theta + e2*theta**2)*p**2

  t1 = theta + 0.5*dp*r1

  !---------------------------
  ! Runge-Kutta Step 2
  !---------------------------
  r2 = a0 + a1*t1 + a2*t1**2 + a3*t1**3 &
     + (b0 + b1*t1)*(S - 35.0) &
     + (c0 + c1*t1 + c2*t1**2 + c3*t1**3 &
     + (d0 + d1*t1)*(S - 35.0))*(p + 0.5*dp) &
     + (e0 + e1*t1 + e2*t1**2)*(p + 0.5*dp)**2

  t2 = theta + 0.5*dp*r2

  !---------------------------
  ! Runge-Kutta Step 3
  !---------------------------
  r3 = a0 + a1*t2 + a2*t2**2 + a3*t2**3 &
     + (b0 + b1*t2)*(S - 35.0) &
     + (c0 + c1*t2 + c2*t2**2 + c3*t2**3 &
     + (d0 + d1*t2)*(S - 35.0))*(p + 0.5*dp) &
     + (e0 + e1*t2 + e2*t2**2)*(p + 0.5*dp)**2

  t3 = theta + dp*r3

  !---------------------------
  ! Runge-Kutta Step 4
  !---------------------------
  r4 = a0 + a1*t3 + a2*t3**2 + a3*t3**3 &
     + (b0 + b1*t3)*(S - 35.0) &
     + (c0 + c1*t3 + c2*t3**2 + c3*t3**3 &
     + (d0 + d1*t3)*(S - 35.0))*(p + dp) &
     + (e0 + e1*t3 + e2*t3**2)*(p + dp)**2

  !---------------------------
  ! Potential Temperature
  !---------------------------
  theta = theta + dp*(r1 + 2.0*r2 + 2.0*r3 + r4)/6.0

  PRINT *
  PRINT *, 'Potential Temperature (C):'
  PRINT '(F10.5)', theta

END subroutine ALGO8_POTENTIAL_TEMPERATURE


!/////////////////////////////PROGRAM 9/////////////////////////////////////!
subroutine ALGO9_SOUND_SPEED
IMPLICIT NONE

REAL :: S, t, p, U
REAL :: Cw, A, B, D
REAL :: f0, f1, f2, f3
REAL :: g0, g1, g2, g3

! ==== Constants Cw(t,p)
REAL, PARAMETER :: C00=1402.388, C01=5.03711, C02=-5.80852E-2
REAL, PARAMETER :: C03=3.3420E-4, C04=-1.47800E-6, C05=3.1464E-9
REAL, PARAMETER :: C10=0.153563, C11=6.8982E-4, C12=-8.1788E-6
REAL, PARAMETER :: C13=1.3621E-7, C14=-6.1185E-10
REAL, PARAMETER :: C20=3.1260E-5, C21=-1.7107E-6, C22=2.5974E-8
REAL, PARAMETER :: C23=-2.5335E-10
REAL, PARAMETER :: C30=-9.7729E-8, C31=3.8504E-10, C32=-2.3643E-12

! ==== Constants A(t,p)
REAL, PARAMETER :: A00=1.389, A01=-1.262E-2, A02=7.164E-5
REAL, PARAMETER :: A03=2.006E-6, A04=-3.21E-8
REAL, PARAMETER :: A10=9.4742E-5, A11=-1.2580E-5
REAL, PARAMETER :: A12=-6.4885E-8, A13=1.0507E-8, A14=-2.0122E-10
REAL, PARAMETER :: A20=-3.9064E-7, A21=9.1041E-9
REAL, PARAMETER :: A22=-1.6002E-10, A23=7.988E-12
REAL, PARAMETER :: A30=1.100E-10, A31=6.649E-12, A32=-3.389E-13

! ==== Constants B(t,p)
REAL, PARAMETER :: B00=-1.922E-2, B01=-4.42E-5
REAL, PARAMETER :: B10=7.3637E-5, B11=1.7945E-7

! ==== Constants D(t,p)
REAL, PARAMETER :: D10=-7.9836E-6

! ==== Input
PRINT *, '---------------------------------------'
PRINT *, ' ALGORITHM 9: SOUND SPEED IN SEAWATER'
PRINT *, ' Chen & Millero (1977)'
PRINT *, '---------------------------------------'
PRINT *, 'Input Salinity S (0-40): '
READ *, S
PRINT *, 'Input Temperature t (C): '
READ *, t
PRINT *, 'Input Pressure p (decibars): '
READ *, p

! ==== Cw(t,p)
f0 = C00 + C01*t + C02*t**2 + C03*t**3 + C04*t**4 + C05*t**5
f1 = C10 + C11*t + C12*t**2 + C13*t**3 + C14*t**4
f2 = C20 + C21*t + C22*t**2 + C23*t**3
f3 = C30 + C31*t + C32*t**2
Cw = f0 + f1*p + f2*p**2 + f3*p**3

! ==== A(t,p)
g0 = A00 + A01*t + A02*t**2 + A03*t**3 + A04*t**4
g1 = A10 + A11*t + A12*t**2 + A13*t**3 + A14*t**4
g2 = A20 + A21*t + A22*t**2 + A23*t**3
g3 = A30 + A31*t + A32*t**2
A = g0 + g1*p + g2*p**2 + g3*p**3

! ==== B(t,p)
B = (B00 + B01*t) + (B10 + B11*t)*p

! ==== D(t,p)
D = D10*p

! ==== Sound Speed
U = Cw + A*S + B*S*SQRT(S) + D*S*S

! ==== Output
PRINT *
PRINT *, '------------- RESULT ------------------'
PRINT '(A,F8.2)', ' Salinity (PSU):        ', S
PRINT '(A,F8.2)', ' Temperature (C):      ', t
PRINT '(A,F8.2)', ' Pressure (dbar):       ', p
PRINT *
PRINT '(A,F10.3)', ' Sound Speed (m/s):     ', U
PRINT *, '---------------------------------------'

END subroutine ALGO9_SOUND_SPEED

end program Projek_PDOK_MANTAPJIWA_BANYAKBANGET