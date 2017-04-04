!******************************************************************************

MODULE WTCH
INTEGER                 wnumb,          & ! Number of EXT_WTCH elements
                        cnumb,          & ! Current element number (1:wnumb)
                        wstat,          & ! Working status
                        wleng,          & ! Length of arrays (for turns)
                        wturn,          & ! Number of turns gathered
                        pturn             ! Previous NTURN (to check changes)
REAL(8)                 wshft             ! Coordinate shift (for matrix)
INTEGER, ALLOCATABLE :: whist(:,:,:,:), & ! Histograms (for beam sizes)
                        wpart(:,:,:),   & ! Number of particles (with weights)
                        wtrac(:)          ! Number of turns to track (MEAN)
REAL(8), ALLOCATABLE :: wsize(:,:),     & ! Cell sizes for histograms
                        wstep(:),       & ! Step (accuracy) for MEAN
                        wmean(:,:,:,:), & ! Buffer for center of mass
                        wmatr(:,:,:)      ! Transport matrix (for beam sizes)
END MODULE WTCH

!******************************************************************************

SUBROUTINE WTCH_INIT
USE WTCH
INCLUDE 'commons.f90'

wnumb = COUNT (.NOT.equilibr.AND.elem_type(:elem_numb) == 'EXT_WTCH')
cnumb = wnumb
wleng = NINT (MAXVAL (param(bound_setts+5:bound_setts+6)))
wstat = 0
IF (ALLOCATED(wsize)) THEN
  IF (wnumb /= SIZE(wsize,2)) DEALLOCATE (whist, wpart, wtrac, &
                                          wsize, wstep, wmean)
END IF
IF (wnumb == 0) RETURN

IF (.NOT.ALLOCATED(whist)) ALLOCATE (whist(2,-1001:1001,wnumb,2), &
                              wpart(wleng,wnumb,2), wtrac(wnumb), &
            wsize(2,wnumb), wmean(3,wleng,wnumb,2), wstep(wnumb))
k = 0
DO i = 1, elem_numb
  IF (elem_type(i) /= 'EXT_WTCH') CYCLE
  k = k + 1; j = elem_addr(i)
  wtrac(k) = NINT (param(j))
  wstep(k) = param(j+1)
END DO
WHERE (wstep == 0) wstep = 1.0D-5   ! 0.1 mkm
END

!******************************************************************************

SUBROUTINE WTCH_MATR (shft, wsta)
USE WTCH
REAL(8) shft
INTEGER wsta

IF (wnumb == 0) RETURN
wshft = shft
wstat = wsta
IF (wstat == 7) ALLOCATE (wmatr(6,7,wnumb))
END

!******************************************************************************

SUBROUTINE WTCH_SIZE (emit, sigm, lbeg)
USE WTCH
REAL(8) emit(2), sigm(2), lbeg(2,8), lend(2,9)
INTEGER i

wstat = 0
DO i = 1, wnumb
  CALL LATT_PROPAGATE (lbeg, 1, wmatr(1,1,i), lend)
  wsize(:,i) = SQRT (lend(:,1) * emit + lend(:,5) * emit(2:1:-1) + &
                    (lend(1,3:4) * sigm(2))**2) / 20
END DO
IF (ALLOCATED(wmatr)) DEALLOCATE (wmatr)
END

!******************************************************************************

SUBROUTINE WTCH_START
USE WTCH

IF (wnumb == 0) RETURN
wstat =-1
wturn = 0
pturn = 0
whist = 0
wpart = 0
wmean = 0
END

!******************************************************************************

SUBROUTINE WTCH_GATHER (coord, pwght, nturn)
USE WTCH
REAL(8) coord(6) ! Coordinates of the test particle
INTEGER pwght, & ! Particle's weight
        nturn, & ! Number of turns passed
        i, j

IF (wstat == 0) RETURN
cnumb = MOD (cnumb, wnumb) + 1
IF (wstat > 0) THEN
  wmatr(:,wstat,cnumb) = coord
  IF (wstat <= 6) wmatr(:,wstat,cnumb) = (coord - wmatr(:,7,cnumb)) / wshft
  RETURN
END IF

DO i = 1, 2
  j = MIN (1001, INT(ABS(coord(i*2-1))/wsize(i,cnumb))+1)
  IF (coord(i*2-1) < 0) j = -j
  whist(i,j,cnumb,1) = whist(i,j,cnumb,1) + pwght
END DO
IF (nturn /= pturn) THEN
  pturn = nturn
  wturn = wturn + 1
END IF
IF (wtrac(cnumb) /= 0.AND.wtrac(cnumb) < nturn) RETURN
wpart(  wturn,cnumb,1) = wpart(  wturn,cnumb,1) + pwght
wmean(:,wturn,cnumb,1) = wmean(:,wturn,cnumb,1) + pwght * coord(1:5:2)
END

!******************************************************************************

SUBROUTINE REWR_WTCH (tag, what, IC)
USE WTCH
CHARACTER(1) what   ! What to do: Read or Write
INTEGER      tag, & ! Tag to mark record (for parallel processing)
             IC,  & ! Returned status (0 - OK)
             i, j

IF (tag == 0) THEN
  IF (what == 'R') CALL WTCH_INIT
  IF (wnumb > 0) CALL REWR_PAR (0, what, wsize, 4*wnumb, IC)
END IF
IF (tag == 0.OR.wnumb == 0) RETURN
i = INDEX ('WR', what)
CALL REWR_PAR (tag,   what, whist(1,-1001,1,i), SIZE(whist)/2, IC)
CALL REWR_PAR (tag+1, what, wpart(1,1,i),       SIZE(wpart)/2, IC)
CALL REWR_PAR (tag+2, what, wmean(1,1,1,i),     SIZE(wmean),   IC)
IF (what == 'W') RETURN
whist(:,:,:,1) = whist(:,:,:,1) + whist(:,:,:,2)
wpart(:,:,  1) = wpart(:,:,  1) + wpart(:,:,  2)
wmean(:,:,:,1) = wmean(:,:,:,1) + wmean(:,:,:,2)
END

!******************************************************************************

SUBROUTINE WRITE_WTCH_DATA
USE WTCH
INTEGER       pass, i, j, k, m, IC
CHARACTER(80) header
REAL(8)       hist(-1001:1001), mean, norm, sigm, sig2, pos2, expp, &
              dval, dold, dmax
INTEGER, ALLOCATABLE :: cmean(:,:)

IF (wnumb == 0) RETURN
IF (ANY(wpart(1,:,1) > 0)) ALLOCATE (cmean(3,2*wleng))
DO m = 1, wnumb
  header = 'Sizes'
  IF (wnumb > 1) THEN
    header(6:6) = '_'
    WRITE (header(7:7), '(I1)', IOSTAT=i) m
  END IF
  header = TRIM(header)//' (x,y):'
  DO i = 1, 2
    hist = whist(i,:,m,1) * 40 * SQRT (ACOS(0.0_8)) / SUM (whist(i,:,m,1))
    mean = 0; k = 0
    DO j = 1000, 1, -1
      IF (ALL(hist(-j:j:2*j) == 0)) CYCLE
      IF (k == 0) k = j
      mean = mean + j * (hist(j) - hist(-j))
    END DO
    mean = mean / SUM (hist(-k:k)) / 20
    sigm = 1; dmax = 0.2_8
    DO pass = 1, 20
      dold = dval; dval = 0
      sig2 = sigm ** 2
      norm = -0.05_8
      DO j = -k, k
        IF (j == 0) THEN
          norm = -norm
        ELSE
          pos2 = ((ABS(j) - 0.5_8) * norm - mean) ** 2
          expp = EXP (-pos2/sig2/2)
          dval = dval + expp * (sig2 - pos2) * (expp - sigm * hist(j))
        END IF
      END DO
      IF (pass > 1.AND.dval*dold < 0) dmax = dmax / 2
      sigm = sigm + SIGN (MIN(dmax,ABS(dval)/20), dval)
      IF (pass > 3.AND.MIN(ABS(dval),ABS(dval-dold)) < 0.001) EXIT
    END DO
    sigm = sigm * wsize(i,m) * 20
    CALL WREFF (sigm, header(i*15+2:i*15+9), 'D.E')
  END DO

  IF (wpart(1,m,1) == 0) THEN
    WRITE (1, '(A)') TRIM(header)
    CYCLE
  END IF
  cmean = 0
  DO i = 1, wleng
    IF (wpart(i,m,1) == 0) EXIT
    cmean(:,i) = NINT (wmean(:,i,m,1) / wpart(i,m,1) / wstep(m))
  END DO
  CALL ZIP_D (cmean, 3*wleng, j, 4, IC)
  IF (IC /= 0) CALL FINI ('Error in WRITE_WTCH_DATA: ZIP')
  i = INDEX (header, ':')
  IF (j == 12*wleng) header(i+1:i+1) = ':'
  WRITE (1, '(A)') TRIM(header)
  CALL ZIP_7 (cmean, j, k)
  CALL WRITE_ZTXT (cmean, k)
END DO
IF (ALLOCATED(cmean)) DEALLOCATE (cmean)
END

!******************************************************************************

SUBROUTINE EXPORT_WTCH (turns, nturn, stat)
USE WTCH
INTEGER       turns(3), & ! 1 - flag, 2:3 - turns to output
              nturn,    & ! Total number of tracked turns
              stat,     & ! Returned status, 0 - OK
              ntrc, flag, nall, nold, nwtc, tout(2), i, j, k, IC
CHARACTER(80) str
INTEGER, ALLOCATABLE :: cmean(:,:)

IF (.NOT.ALLOCATED(whist)) CALL WTCH_INIT
stat =-1; IF (wnumb == 0) RETURN
ntrc = COUNT (wtrac >= 0)
flag = 0; nall = 0
REWIND (1)
ALLOCATE (cmean(3,2*wleng))
DO
  READ (1, '(A)', IOSTAT=IC) str
  IF (IC /= 0) EXIT
  IF (flag <= 1.AND.INDEX(str,':') == 0) CYCLE
  IF (flag == 0.AND.str(:12) == 'Turn_count:') THEN
    i = INDEX (str, '(all)=')
    str = ADJUSTL (str(i+6:))
    i = INDEX (str, ' ') - 1
    nold = nall
    READ (str(:i), *) nall
    IF (nall-nold > wleng) CALL FINI ('Error in EXPORT_WTCH: data corrupted')
    flag = 1; nwtc = 0
  ELSE IF (flag == 1.AND.str(:5) == 'Sizes') THEN
    nwtc = nwtc + 1
    IF (nwtc == wnumb) flag = 0
    tout = turns(2:3)
    k = nturn
    IF (wtrac(nwtc) > 0) k = MIN (k, wtrac(nwtc))
    WHERE (tout <= 0) tout = tout + k
    IF (wtrac(nwtc) < 0.OR.tout(1) > nall.OR.tout(2) <= nold) CYCLE
    WHERE (tout >  k) tout = k
    CALL READ_ZTXT (cmean, i)
    IF (nwtc < wnumb) BACKSPACE (1)
    IF (i <= 0) CYCLE
    CALL UNZIP_7 (cmean, j, i, stat)
    IF (stat /= 0) RETURN
    CALL UNZIP_D (cmean, k, j, i, stat)
    IF (k /= 3*wleng.OR.i /= 4) stat = 10
    IF (stat /= 0) RETURN
    DO i = MAX(1,tout(1)-nold), MIN(tout(2),nall)-nold
      str = ' '
      DO j = 1, 3
        CALL WREFF (ABS(cmean(j,i))*wstep(nwtc), str(j*12-10:j*12-3), 'D.Ee0')
        IF (cmean(j,i) < 0) str (j*12-11:j*12-11) = '-'
      END DO
      WRITE (str(37:47), '(I11)') nold+i
      IF (ntrc > 1) WRITE (str(50:52),  '(A,I2.2)') '#', nwtc
      IF (turns(1) > 0) WRITE (16, '(A)') TRIM(str)
      IF (turns(1) < 0) WRITE (*,  '(A)') TRIM(str)
    END DO
  END IF
END DO
DEALLOCATE (cmean)
END

!******************************************************************************

SUBROUTINE PRE_TRACK_EXT (numb, eparm)
INCLUDE 'commons.f90'
INTEGER numb, &       ! Element's number
        aptag         ! Check the aperture after tracking ?
REAL(8) eparm(21), cc ! Parameters of the element

aptag = 0
SELECT CASE (elem_type(numb))
CASE ('EXT_SHFT'); elem_itp(numb) = 101
! Shift of 6D coordinate vector.
  IF (elem_aux(numb) < 0) eparm(1:6)   = -eparm(1:6)
  IF (elem_tag(numb) < 0) eparm(2:6:2) = -eparm(2:6:2)
  aptag = 1
CASE ('EXT_RFQ'); elem_itp(numb) = 102
! RF Quadrupole.
  IF (elem_aux(numb) < 0) eparm(1) = -eparm(1)
CASE ('EXT_TL04'); elem_itp(numb) = 103
! Electron lens for beam-beam compensation in the Tevatron. The space
! charge density is: Rho(R)=Rho_0/(1+(R/r0)**4)
  eparm(2) = 1 / eparm(2)**2                      ! 1/r0**2
  eparm(1) = eparm(1) / eparm(2) * ACOS(0.0_8)*2  ! pi*Rho_0*r0**2
CASE ('EXT_TL08'); elem_itp(numb) = 104
! Electron lens for beam-beam compensation in the Tevatron. The space
! charge density is: Rho(R)=Rho_0/(1+(R/r0)**8)
  eparm(2) = SQRT(2.0_8) / eparm(2)**2            ! sqrt(2)/r0**2
  eparm(1) = eparm(1) / eparm(2) * ACOS(0.0_8)/2  ! pi*Rho_0*r0**2/(4*sqrt(2))
  IF (eparm(7) > 1) eparm(1) = eparm(1) / INT(eparm(7))  ! For curved e-beam
CASE ('EXT_TLA8', 'EXT_TLB8'); elem_itp(numb) = 105
! Electron lens for beam-beam compensation in the Tevatron. The space
! charge density is: Rho(R)=Rho_0*(a*(R/r0)**2+b)/(1+(R/r0)**8)
  eparm(11) = 1 / eparm(2)**4                     ! 1/r0**4
  eparm(5)  = eparm(5) / eparm(6) * SQRT(8.0_8)   ! a/b*2*sqrt(2)
  eparm(2)  = SQRT(2*eparm(11))                   ! sqrt(2)/r0**2
  eparm(1)  = eparm(1) * eparm(6) / eparm(2) * ACOS(0.0_8)/2 ! pi*Rho_0*r0**2*
                                                             ! b/(4*sqrt(2))
  IF (elem_type(numb) == 'EXT_TLB8') THEN
    elem_itp(numb) = 106
    elem_mrk(numb) = 1  ! Mark for time dependence
  END IF
CASE ('EXT_TLGS'); elem_itp(numb) = 107
! Electron lens for beam-beam compensation in the Tevatron. The space
! charge density is: Rho(R)=Rho_0*EXP(-R**2/(2*r0**2))
  eparm(2) = -1 / (2*eparm(2)**2)                 ! -1/(2*r0**2)
  eparm(1) = -eparm(1) / eparm(2) * ACOS(0.0_8)*2 ! 2pi*Pho_0*r0**2
CASE ('EXT_TLHB'); elem_itp(numb) = 108
! Electron lens with Hollow Beam. AV 8/4/2008
! Space charge density is N*e/2/pi/sigma^2*(1-exp(-r^2/2/sigma^2))
! eparm(1) is equal to 4*pi*xi/beta. Positive is focussing
  eparm(2) = -1 / (2*eparm(2)**2)                 ! -1/(2*sigma**2)
CASE ('EXT_TLHC'); elem_itp(numb) = 109
! Electron lens with hollow beam for collimation AV 10/27/2009
! Space charge density is:
! 0     if r < R1
! Const if R1 < r < R2
! 0     if r > R2
! eparm(1) is the maximum kick in radians
! eparm(2) is R2 in cm
! eparm(3) - x offset
! eparm(4) - y offset
! eparm(5) is for pulsing each Nth turn
! eparm(6) = R2 / R1 (by default, 1.5)
  IF (eparm(5) == 0) eparm(5) = 1
  IF (eparm(6) <= 1) eparm(6) = 1.5_8
! eparm(7) - if non-zero, account for entrance/exit e-beam bends (AV 10/7/2013)
! normalization of the kick is from FERMILAB-FN-0972-APC
! updated 3/31/2014 - made the kicks symplectic
  eparm(7) = eparm(7) * eparm(1) / eparm(2) * 7.233D-5 * 2.386_8
! eparm(8) is the random lens current variation, must be [0:1] AV 2/6/2014
! eparm(9) if non-zero, the two entrance/exit bends are treated as S-shape
  eparm(10) = 1                          ! Can be changed if eparm(8) /= 0
  IF (eparm(8) /= 0) elem_mrk(numb) = 1  ! Mark for time dependence
CASE ('EXT_SQL1'); elem_itp(numb) = 110
! SEFT TEL with imperfections by A.Romanov 6/19/09
CASE ('EXT_SQL2'); elem_itp(numb) = 111
! Gaussian TEL with imperfections by A.Romanov 6/25/09
CASE ('EXT_COOL'); elem_itp(numb) = 112
! Electron cooling system.
  eparm( 1) = eparm( 1) * eparm( 2)**3
  eparm( 2) = eparm( 2)**2                        ! (V_eff / c)^2
  eparm(10) = param(gamma_weak)                   ! G(amma)
  eparm(11) = eparm(10)**2 - 1                    ! G^2 - 1
  eparm(12) = eparm(11) / eparm(10)**2            ! B(eta)^2
CASE ('EXT_LONG'); elem_itp(numb) = 113
! Longitudinal shape asymmetry.
  eparm(3) = ACOS(0.0_8)*4 / eparm(3)             ! 2pi/L
CASE ('EXT_RFCV'); elem_itp(numb) = 114
! RF Cavity (nonlinear).
  eparm(2) = ACOS(0.0_8)*4 / eparm(2)             ! 2pi/L
  eparm(5) = ASIN(eparm(4))/ eparm(2)
CASE ('EXT_RFDM'); elem_itp(numb) = 115
! RF Damping.
! No preliminary calculations are required.
CASE ('EXT_DRIF'); elem_itp(numb) = 116
! Drift space, to produce chromaticity of beta-functions up to 3rd order.
! No preliminary calculations are required.
  aptag = 1
CASE ('EXT_WTCH'); elem_itp(numb) = 117
! Watching point for Shottky spectra and physical sizes (non-equilibr. only).
  IF (eparm(2) <= 0) eparm(2) = 1.0D-5            ! 0.1 mkm by default
CASE ('EXT_CRAB'); elem_itp(numb) = 118
! Crab cavity (can be nonlinear) after pi/2 shift (e.g. in the IP).
  IF (elem_aux(numb) < 0) eparm(1) = -eparm(1)
  eparm(1:2) = eparm(1) * [COS(eparm(2)), SIN(eparm(2))]
  eparm(3)   = eparm(3) / param(sigma_weak)
  aptag = 1
CASE ('EXT_CCAV'); elem_itp(numb) = 119
! Different type of crab cavity AV 3/13/2013
  IF (eparm(2) /= 0.OR.eparm(4) == 0) THEN
    eparm(1:2) = eparm(1) * [COS(eparm(2)), SIN(eparm(2))]
  ELSE
    eparm(2) = eparm(4)
  END IF
CASE ('EXT_NLN1'); elem_itp(numb) = 120
! Some nonlinear terms in Hamiltonian:
! 1 - X^2*Px, 2 - X*Y^2, 3 - X*Py^2
! No preliminary calculations are required.
  aptag = 1
CASE ('EXT_DODE'); elem_itp(numb) = 121
! Dodecapole lens.
  eparm(1) = eparm(1) / 120
CASE ('EXT_MULT'); elem_itp(numb) = 122
! Thin multipole.
! First 10 parameters are normal multipole coefficients starting from 1,
! parameters 11-20 are skew-multipole coefficients.
  DO i = 10, 1, -1
    IF (ANY(eparm(i:i+10:10) /= 0)) EXIT
  END DO
  eparm(21) = i + 0.1
CASE ('EXT_NLL1'); elem_itp(numb) = 123
! Non-Linear Lens type 1
! The potential is f(z)=-d*z/(a*z^2+1)
! eparm(1)=d, eparm(2)=a
CASE ('EXT_NLL2'); elem_itp(numb) = 124
! Just octupolar component of NLL1
CASE ('EXT_NLL3'); elem_itp(numb) = 125
! Non-linear Lens type 3
! U(u,v)=(u*sqrt(u**2-1)+v*sqrt(1-v**2)*(-pi/2+acos(v))/(u**2-v**2)
! eparm(1) = lens strength
! eparm(2) = c1 = c * SQRT (beta)
! eparm(3:4) - misalignments (x, y)
! eparm(5) = tilt angle in (x, y) plane
  eparm(1) = eparm(1) / 2
  eparm(6:7) = SIN (eparm(5)) / eparm(2) * [1,-1]
  eparm(5)   = COS (eparm(5)) / eparm(2)
CASE ('EXT_NLLR'); elem_itp(numb) = 126
! Non-linear lens for round beams (integrable, 1st case).
! F(r) = -eparm(1) * r / (eparm(2) * r**2 + 1)
! No preliminary calculations are required.
CASE ('EXT_NLR2'); elem_itp(numb) = 127
! Non-linear lens for round beams (integrable, 2nd case).
! F(r) = -eparm(1) * LOG (SUM (eparm(2:3) * EXP([-r,r])) / &
!                         SUM (eparm(2:3) * EXP([r,-r])))
! No preliminary calculations are required.
CASE ('EXT_NLPO'); elem_itp(numb) = 128
! Polar non-linear Lens.
! U(x,y)= eparm(1)*(x**2-y**2)/(x**2+y**2)**2
  eparm(1) = eparm(1) * 2
CASE ('EXT_CLMR'); elem_itp(numb) = 129
! Collimator (plate shifts from the center, in cm).
! eparm(1:4) = X-, X+, Y-, Y+
  WHERE (eparm(:4) == 0) eparm(:4) = 10000 * [-1, 1, -1, 1] ! by default 100 m
CASE ('EXT_THNL'); elem_itp(numb) = 130
! Thin lens (can be focussing or defocussing in both directions).
! No preliminary calculations are required.
CASE ('EXT_SOLE'); elem_itp(numb) = 131
! Thin solenoid with formulae from madx. AV 2/13/2013
! first parameter is KS, second KSI as in madx
! eparm(11:12) - preceding drift.
  eparm(1:2) = eparm(1:2) / 2
  eparm(3:4) = eparm(1) * [eparm(2), 0.5_8]
  IF (eparm(12) /= 0) eparm(12) = eparm(11) * bet0i
  aptag = 1
CASE ('EXT_DPDG'); elem_itp(numb) = 132
! Thin dipole edge with formulae from madx. AV 2/13/2013
! eparm(1) is H, 2-E1, 3-FINT, 4-HGAP, 5-Tilt
! eparm(11:12) - preceding drift.
  cc = COS (eparm(2))
  eparm(3:4) = [0, 2] * eparm(1) * eparm(3) * eparm(4)
  eparm(6:7) = [1,-1] * eparm(1) * TAN (eparm(2) - eparm(3:4) * (2 / cc - cc))
  eparm(2:3) = (eparm(6) - eparm(7)) * COS (eparm(5)) * SIN (eparm(5))
  eparm([1,4]) = eparm([6,7]) * COS (eparm(5))**2 + &
                 eparm([7,6]) * SIN (eparm(5))**2
  IF (eparm(12) /= 0) eparm(12) = eparm(11) * bet0i
CASE ('EXT_CORR'); elem_itp(numb) = 133
! Thin dipole corrector.
  IF (elem_aux(numb) < 0) eparm(1:2) = -eparm(1:2)
CASE ('EXT_TILT'); elem_itp(numb) = 134
! Tilt in X-Y plane.
  IF (elem_aux(numb) < 0) eparm(1) = -eparm(1)
  eparm([2,5]) = COS (eparm(1))
  eparm([3,4]) = SIN (eparm(1)) * [-1, 1]
  aptag = 1
CASE ('EXT_LEMO'); elem_itp(numb) = 135
! Levichev-Morozov so-called beam-beam compensation.
  eparm(2) = eparm(2) / 2
  eparm(3) = SQRT (eparm(3))
  eparm(4) = eparm(2) * eparm(3)
  eparm(5) = COS  (eparm(4))
  eparm(6) = SIN  (eparm(4))
  eparm(7) = COSH (eparm(4))
  eparm(8) = SINH (eparm(4))
  eparm(11:13:2) = eparm(6:8:2) * eparm(3)
  eparm(12:14:2) = eparm(6:8:2) / eparm(3)
CASE ('EXT_DNSE'); elem_itp(numb) = 136
! Dipole noise - nth turn pulsing AV 4/27/2016
! eparm(1) - horizontal kick value in rad
! eparm(2) - vertical kick value in rad
! eparm(3) - integer, horizontal kick every nth turn
! eparm(4) - integer, vertical kick every nth turn
  WHERE (eparm(3:4) == 0) eparm(3:4) = 1
CASE ('EXT_DNSF'); elem_itp(numb) = 137
! Dipole noise - sinusoidal AV 4/27/2016
! eparm(1) - horizontal kick value in rad
! eparm(2) - vertical kick value in rad
! eparm(3) - frequency of horizontal excitation (in units of rev. freq.)
! eparm(4) - frequency of vertical excitation
  eparm(3:4) = eparm(3:4) * ACOS(0.0_8) * 4 ! 2*pi*f/f0
CASE ('EXT_DNRU'); elem_itp(numb) = 138
! Dipole noise - random uniform MF 04/03/2017
! eparm(1) - max. horizontal kick value in rad
! eparm(2) - max. vertical kick value in rad
! eparm(3) - float [0,1], modulation depth of horizontal kick
! eparm(4) - float [0,1], modulation depth of vertical kick
  elem_mrk(numb) = 1 ! Mark for time dependence
CASE DEFAULT
  CALL FINI ('PRE_TRACK_EXT, unrecognized type: '//elem_type(numb))
END SELECT
elem_tag(numb) = aptag
END

!******************************************************************************

SUBROUTINE TIME_DEP_EXT (type, eparm)
CHARACTER(*) type         ! Element's type
REAL(8)      eparm(21), & ! Parameters of the element
             rnumb        ! Random number
INTEGER      i

SELECT CASE (type)
CASE ('EXT_TLB8')
  DO i = 7, 10
! Noise of bbcomp e-beam's current and separation (horizontal and vertical).
    IF (i == 8.OR.eparm(i) == 0) CYCLE
    CALL RNORR (rnumb, 2)
    eparm(i+5) = eparm(i) * rnumb
  END DO
CASE ('EXT_TLHC')
  IF (eparm(8) == 0) RETURN
  CALL RANCOMB (rnumb, 2)
  eparm(10) = 1 - eparm(8) * rnumb
CASE ('EXT_DNRU')
  IF (eparm(1) == 0 .AND. eparm(2) == 0) RETURN
  CALL RANCOMB (rnumb, 2)
  eparm(5:6) = 1-eparm(3:4)*rnumb
CASE DEFAULT
  CALL FINI ('TIME_DEP_EXT, unrecognized type: '//type)
END SELECT
END

!******************************************************************************

SUBROUTINE TRACK_EXT (numb, type, eparm)
INCLUDE 'commons.f90'
INTEGER    numb,      & ! Element's number
           type,      & ! Element's type
           nsmp
REAL(8) :: eparm(21), & ! Parameters of the element
           xy(2), rr, ra, rb, r2, dd, vv(3), cc(4), yx(2), zz(2), uv(2), &
           pp(2), dp(2), ac(2), rnrm, step(2), kick(2), pi = 2 * ACOS (0.0_8)

SELECT CASE (type)
CASE (101) ! EXT_SHFT
  coord = coord + eparm(1:6)
CASE (102) ! EXT_RFQ
  xy = coord(1:3:2) + eparm(3:4)
  xy(1) = -xy(1)
  coord(2:4:2) = coord(2:4:2) - eparm(1) * xy * coord(5)
CASE (103) ! EXT_TL04
  CALL TL_LENS_XY
  rr = SUM(xy**2); IF (rr == 0) RETURN
  dd = ATAN(rr*eparm(2)) / rr
  coord(2:4:2) = coord(2:4:2) - eparm(1) * dd * xy
CASE (104) ! EXT_TL08
  IF (eparm(7) > 1) THEN
    xy = coord(1:3:2) + eparm(3:4) * COS (pi / eparm(7) - pi)
  ELSE
    CALL TL_LENS_XY
  END IF
  i = 1
  DO WHILE (i == 1.OR.eparm(7) >= i)
    i = i + 1
    rr = SUM(xy**2); dd = 0
    IF (rr > 0) THEN
      ra = rr * eparm(2) + 1    ! sqrt(2)*(r/r0)**2+1
      rb = ra - 2               ! sqrt(2)*(r/r0)**2-1
      dd = (LOG((ra**2+1) / (rb**2+1)) + 2 * (ATAN(ra) + ATAN(rb))) / rr
    END IF
    IF (eparm(7) >= i) THEN
      coord(2:4:2) = coord(2:4:2) - eparm(1) * dd * xy
      xy = coord(1:3:2) + eparm(3:4) * COS(pi*((2*i-1)/eparm(7)-1))
    END IF
  END DO
  coord(2:4:2) = coord(2:4:2) - eparm(1) * dd * xy
CASE (105, 106) ! EXT_TLA8, EXT_TLB8
  CALL TL_LENS_XY
  IF (type == 106) xy = xy + eparm(14:15)
  rr = SUM(xy**2); IF (rr == 0) RETURN
  ra = rr * eparm(2) + 1    ! sqrt(2)*(r/r0)**2+1
  rb = ra - 2               ! sqrt(2)*(r/r0)**2-1
  dd = (LOG((ra**2+1) / (rb**2+1)) + 2 * (ATAN(ra) + ATAN(rb)) + &
        eparm(5) * ATAN(eparm(11)*rr**2)) / rr
  IF (type == 106) dd = dd * (1 + eparm(12))
  coord(2:4:2) = coord(2:4:2) - eparm(1) * dd * xy
CASE (107) ! EXT_TLGS
  CALL TL_LENS_XY
  rr = SUM(xy**2); IF (rr == 0) RETURN
  dd = (1 - EXP(rr*eparm(2))) / rr
  coord(2:4:2) = coord(2:4:2) - eparm(1) * dd * xy
CASE (108) ! EXT_TLHB
  CALL TL_LENS_XY
  rr = SUM(xy**2); IF (rr == 0) RETURN
  dd = 1 + (1 - EXP(rr*eparm(2))) / rr / eparm(2)
  coord(2:4:2) = coord(2:4:2) - eparm(1) * dd * xy
CASE (109) ! EXT_TLHC
  i = eparm(5)
  IF (MOD(nturn,i) /= 0) RETURN
  xy = coord(1:3:2) + eparm(3:4)
  r2 = SUM  (xy**2)
  rr = SQRT (r2)
  ra = eparm(2) / eparm(6)
  IF (rr > ra) THEN ! Main part of the lens - kick from doughnut
    dd = eparm(2) / r2
    IF (rr <= eparm(2)) dd = dd * (r2 / ra**2 - 1) / (eparm(6)**2 - 1)
    coord(2:4:2) = coord(2:4:2) - eparm(1) * eparm(10) * dd * xy
  END IF
  rb = 1.33_8 * eparm(2)
  IF (eparm(7) == 0.OR.rr > rb) RETURN  ! no mapping outside unit circle
  uv = xy / rb                ! normalize to 10sigma = 1.33 outer radius
  CALL TL_LENS_HC
  IF (eparm(9) /= 0) THEN
    pp = dp
    uv(1) = -uv(1)
    CALL TL_LENS_HC
    dp(1) = -dp(1)
    dp = dp + pp
  END IF
  coord(2:4:2) = coord(2:4:2) - eparm(7) * eparm(10) * dp
!PRINT '(I4,6E16.4)', nturn, coord(1), -eparm(7) * eparm(10) * dp(1),&
!                            coord(3), -eparm(7) * eparm(10) * dp(2),&
!                            -0.1/rb*dp(1),-0.1/rb*dp(2)
CASE (110, 111) ! EXT_SQL1, EXT_SQL2: inclined guns
! Inclanation is simulated by slices.
! SQL1 - SEFT, SQL2 - Gaussian.
! eparm(1)   - cm
! eparm(2)   - 1/cm, (4 pi ksi)/beta = G*L*e/(p*c)
! eparm(3:6) - cm
  rr = eparm(1)**2
  IF (type == 111) rr = rr * 2 ! Half of the normalized radius
  nsmp = eparm(7)
  step =(eparm(4:6:2) - eparm(3:5:2)) / nsmp
  xy   = coord(1:3:2) - eparm(3:5:2)  - step / 2
  rnrm = SUM(xy**2) / rr
  kick = 0
  IF (rnrm > 0.0000001) THEN
    IF (type == 110) kick = xy / rnrm * ATAN(rnrm)
    IF (type == 111) kick = xy / rnrm * (1 - EXP(-rnrm))
  END IF
  DO i = 2, nsmp
    xy = xy - step
    rnrm = SUM(xy**2) / rr
    IF (rnrm <= 0.0000001) CYCLE
    IF (type == 110) kick = kick + xy / rnrm * ATAN(rnrm)
    IF (type == 111) kick = kick + xy / rnrm * (1 - EXP(-rnrm))
  END DO
  coord(2:4:2) = coord(2:4:2) + kick * eparm(2) / nsmp
! PRINT *,rnrm
CASE (112) ! EXT_COOL
  dd = coord(6) - SUM(coord(2:4:2)**2) * eparm(11) / 2
  vv(1:2) = coord(2:4:2) * eparm(10)
  vv(3) = dd / eparm(12)
  vv = vv / (1 - dd)
  vv = vv * (1 - eparm(1) / (eparm(2) + eparm(12) * SUM(vv**2))**1.5_8)
  dd = 1 + eparm(12) * vv(3)
  coord(2:4:2) = vv(1:2) / dd / eparm(10)
  coord(6) = (vv(3) + SUM(vv(1:2)**2) / dd / 2) * eparm(12) / dd
CASE (113) ! EXT_LONG
  dd = eparm(1)
  IF (coord(5) > 0) dd = dd * EXP (-coord(5)/eparm(2))
  coord(6) = coord(6) - dd * SIN (coord(5)*eparm(3))
CASE (114) ! EXT_RFCV
  coord(6) = coord(6) + eparm(1) * (SIN ((coord(5) + eparm(5)) * eparm(2)) - &
                                                                 eparm(4))
CASE (115) ! EXT_RFDM
  coord(2:4:2) = coord(2:4:2) * (1 - eparm(1) * coord(5))
CASE (116) ! EXT_DRIF
  coord(1:3:2) = coord(1:3:2) - coord(2:4:2) * coord(6) * &
        (eparm(1:2) + coord(6) * (eparm(3:4) + coord(6) * eparm(5:6)))
  coord(5) = coord(5) - SUM ((eparm(1:2) / 2 + coord(6) * (eparm(3:4) + &
             coord(6) * eparm(5:6) * 1.5_8)) * coord(2:4:2)**2)
CASE (117) ! EXT_WTCH
  CALL WTCH_GATHER (coord, p_wght, nturn)
CASE (118) ! EXT_CRAB
  xy = eparm(1:2) * (1 + eparm(3) * coord(5))
  coord(1:3:2) = coord(1:3:2) - coord(5) * xy
  coord(6) = coord(6) + SUM (coord(2:4:2) * xy)
CASE (119) ! EXT_CCAV
  coord(2:4:2) = coord(2:4:2) - eparm(1:2) * SIN (eparm(3) * coord(5)) * fp0
  coord(6) = coord(6) - eparm(3) * COS (eparm(3) * coord(5)) * &
                                 SUM (eparm(1:2) * coord(1:3:2))
CASE (120) ! EXT_NLN1
  cc(1) = eparm(1) * coord(1)**2
  cc(2) =-eparm(1) * coord(1) * coord(2) * 2 - &
          eparm(2) * coord(3)**2 -             &
          eparm(3) * coord(4)**2
  cc(3) = eparm(3) * coord(1) * coord(4) * 2
  cc(4) =-eparm(2) * coord(1) * coord(3) * 2
  coord(1:4) = coord(1:4) + cc
CASE (121) ! EXT_DODE
  rr = 10 * (coord(1)*coord(3))**2
  xy = coord(1:3:2)**4
  coord(2) = coord(2) - eparm(1) * coord(1) * (xy(1) + 5 * xy(2) - rr)
  coord(4) = coord(4) + eparm(1) * coord(3) * (xy(2) + 5 * xy(1) - rr)
CASE (122) ! EXT_MULT
  xy = coord(1:3:2)
  yx = [-xy(2), xy(1)]
  dp = 0
  DO i = INT(eparm(21)), 1, -1
    pp = dp / (i + 1) + eparm(i:i+10:10)
    dp = pp(1) * xy + pp(2) * yx
  END DO
  coord(2) = coord(2) - dp(1)
  coord(4) = coord(4) + dp(2)
CASE (123) ! EXT_NLL1
  xy = coord(1:3:2)
  dd = (eparm(2) * (xy(1)**2 + xy(2)**2))**2 + &
        eparm(2) * (xy(1)**2 - xy(2)**2) * 2 + 1
  dp = xy * ([1,-1] + eparm(2) * (xy(1)**2 + xy(2)**2))
  coord(2:4:2) = coord(2:4:2) - eparm(1) / dd * dp
CASE (124) ! EXT_NLL2
  xy = coord(1:3:2)
  dp = eparm(2) * (xy**3 - 3 * xy * xy(2:1:-1)**2)
  coord(2:4:2) = coord(2:4:2) + eparm(1) * (dp + xy * [-1,1])
CASE (125) ! EXT_NLL3
  xy = coord(1:3:2) + eparm(3:4)
  xy = xy* eparm(5) + xy(2:1:-1) * eparm(6:7)
  yx = SQRT ((xy(1) + [1,-1])**2 + xy(2)**2)
  zz = 1 / yx
  zz =  zz(1) + [1,-1] * zz(2)
  uv = (yx(1) + [1,-1] * yx(2)) / 2
  IF (ABS(uv(2)) > 1) uv(2) = SIGN (1.0_8, uv(2))
  ra = uv(1)**2 - uv(2)**2
  IF (ABS(ra) < 1.0D-20) ra = SIGN (1.0D-20, ra)
  ac = [ACOSH (uv(1)), ACOS (uv(2)) - pi / 2]
  pp = SQRT ((uv**2 - 1) * [1,-1])
  yx = 0; WHERE (ABS(pp) > 1.0D-12) yx = uv**2 * ac / pp
  dp = ac * pp + [1,-1] * (yx + uv * (1 - 2 * SUM (uv * pp * ac) / ra))
  pp = xy * SUM (dp * zz) + [SUM (dp * zz(2:1:-1)), 0.0_8]
  pp = pp * eparm(5) - pp(2:1:-1) * eparm(6:7)
  coord(2:4:2) = coord(2:4:2) + (eparm(1) / (1 + coord(6)) / ra) * pp
CASE (126) ! EXT_NLLR
  xy = coord(1:3:2)
  coord(2:4:2) = coord(2:4:2) - eparm(1) * xy / (eparm(2) * SUM (xy**2) + 1)
CASE (127) ! EXT_NLR2
  xy = coord(1:3:2)
  rr = SQRT (SUM (xy**2))
  ra = EXP (2 * rr)
  coord(2:4:2) = coord(2:4:2) - eparm(1) * xy / rr * &
        LOG ((eparm(2) + ra) / (eparm(2) * ra + 1))
CASE (128) ! EXT_NLPO
  xy = coord(1:3:2)
  rr = SUM (xy**2)
  coord(2:4:2) = coord(2:4:2) - eparm(1) * xy * ([1,-1] / rr**2 - 2 * &
                                                (xy(1)**2 - xy(2)**2) / rr**3)
CASE (129) ! EXT_CLMR
  IF (nfma(1) > 0.OR.track_out > 0.AND.aptrac == 0) RETURN
  m = 0
  IF (coord(1) <= eparm(1)) m = m + 1
  IF (coord(1) >= eparm(2)) m = m + 2
  IF (coord(3) <= eparm(3)) m = m + 4
  IF (coord(3) >= eparm(4)) m = m + 8
  IF (m == 0) RETURN
  plost = numb + m * 1000000
  IF (watch > 0) CALL PART_DIED
CASE (130) ! EXT_THNL
  coord(2:4:2) = coord(2:4:2) + coord(1:3:2) * eparm(1:2)
CASE (131) ! EXT_SOLE
  zz(1) = COS (eparm(2) * fp0)
  zz(2) = SIN (eparm(2) * fp0)
  xy = coord(1:3:2) + coord(2:4:2) * eparm(11) ! Preceding drift
  rr = SQRT (1 + coord(2)**2 + coord(4)**2)
  pp = coord(2:4:2) / rr
  uv = pp - xy * eparm(3) / fpa
  coord(1) =  xy(1) * zz(1) + xy(2) * zz(2)
     pp(1) =  uv(1) * zz(1) + uv(2) * zz(2)
  coord(3) = -xy(1) * zz(2) + xy(2) * zz(1)
     pp(2) = -uv(1) * zz(2) + uv(2) * zz(1)
  coord(5) = coord(5) + fpp * eparm(2) * fp0 * &
  (xy(1) * uv(2) - xy(2) * uv(1) - fp0 * eparm(4) * (xy(1)**2 + xy(2)**2)) + &
                                                    eparm(12) * (1 - fpp * rr)
  dd = 1 - pp(1)**2 - pp(2)**2
  IF (dd > 0.5_8) THEN
    coord(2:4:2) = pp / SQRT (dd)
  ELSE
    coord(2:4:2) = pp * 2
  END IF
CASE (132) ! EXT_DPDG
  rr = SQRT (1 + coord(2)**2  + coord(4)**2)
  coord(1:3:2) = coord(1:3:2) + coord(2:4:2) * eparm(11) ! Preceding drift
  coord(5) = coord(5) + eparm(12) * (1 - fpp * rr)       ! dz from the drift
  pp = coord(2:4:2) / rr + fp0 * (coord(1) * eparm(:2) + coord(3) * eparm(3:4))
  dd = 1 - pp(1)**2 - pp(2)**2
  IF (dd > 0.5_8) THEN
    coord(2:4:2) = pp / SQRT (dd)
  ELSE
    coord(2:4:2) = pp * 2
  END IF
CASE (133) ! EXT_CORR
  pp = coord(2:4:2) / SQRT (1 + coord(2)**2 + coord(4)**2) + eparm(1:2) * fp0
  coord(2:4:2) = pp / SQRT (1 - pp(1)**2 - pp(2)**2)
CASE (134) ! EXT_TILT
  coord(1:3:2) = coord(1) * eparm(2:3) + coord(3) * eparm(4:5)
  coord(2:4:2) = coord(2) * eparm(2:3) + coord(4) * eparm(4:5)
CASE (135) ! EXT_LEMO
  k = NINT (eparm(1))
  IF (k < 1.OR.k > 4) RETURN
  IF (k /= 2) coord(1:3:2) = coord(1:3:2) + coord(2:4:2) * eparm(2)
  DO i = 1, 2-k/2
    cc(1:3:2) = coord(2:4:2) * eparm(12)
    cc(2:4:2) =-coord(1:3:2) * eparm(11)
    coord(:4) = coord(:4) * eparm(5) - cc
  END DO
  IF (k <= 2) coord(1:3:2) = coord(1:3:2) + coord(2:4:2) * eparm(2) * (3 - k)
  DO i = 1, 2-MIN(1,k/2)
    IF (k == 2) EXIT
    cc(1:3:2) = coord(2:4:2) * eparm(14)
    cc(2:4:2) = coord(1:3:2) * eparm(13)
    coord(:4) = coord(:4) * eparm(7) - cc
  END DO
  IF (MOD(k,2) > 0) coord(1:3:2) = coord(1:3:2) + coord(2:4:2) * eparm(2)
CASE (136) ! EXT_DNSE
  pp = 0
  IF (MOD(nturn,INT(eparm(3))) == 0) pp(1) = 1
  IF (MOD(nturn,INT(eparm(4))) == 0) pp(2) = 1
  coord(2:4:2) = coord(2:4:2) + eparm(1:2) * pp
CASE (137) ! EXT_DNSF
  pp = COS (eparm(3:4) * nturn)
  coord(2:4:2) = coord(2:4:2) + eparm(1:2) * pp
CASE (138) ! EXT_DNRU
  IF (eparm(1) == 0 .AND. eparm(2) == 0) RETURN
  coord(2:4:2) = coord(2:4:2) + eparm(1:2)*eparm(5:6)
CASE DEFAULT
  PRINT *, 'TRACK_EXT, unrecognized type: ', type
  CALL FINI (' ')
END SELECT
CONTAINS
!------------------------------------------------------------------------------

SUBROUTINE TL_LENS_XY
IF (eparm(8) > 0) THEN
  xy = coord(1:3:2) + eparm(3:4) * EXP(-nturn*eparm(8))
ELSE
  xy = coord(1:3:2) + eparm(3:4)
END IF
END SUBROUTINE TL_LENS_XY

!------------------------------------------------------------------------------

SUBROUTINE TL_LENS_HC
REAL(8) :: cv(0:18,0:18) = RESHAPE (                                        &
     [-159.0903610650_8,  0.0000000000_8,  9.5050101765_8,  0.0000000000_8, &
        -0.3436871207_8,  0.0000000000_8, -0.2985366152_8,  0.0000000000_8, &
         0.0985405384_8,  0.0000000000_8,  0.0161206072_8,  0.0000000000_8, &
        -0.0184178080_8,  0.0000000000_8,  0.0043768921_8,  0.0000000000_8, &
        -0.0024667053_8,  0.0000000000_8,  0.0017656902_8, 24.5689120812_8, &
         0.0000000000_8, -4.1629489758_8,  0.0000000000_8,  0.0277619519_8, &
         0.0000000000_8,  0.2922158154_8,  0.0000000000_8, -0.0808736919_8, &
         0.0000000000_8, -0.0136489125_8,  0.0000000000_8,  0.0202978585_8, &
         0.0000000000_8, -0.0023676164_8,  0.0000000000_8,  0.0048290310_8, &
         0.0000000000_8,  0.0000000000_8,  4.4150647143_8,  0.0000000000_8, &
        -1.2048699771_8,  0.0000000000_8, -0.1758690129_8,  0.0000000000_8, &
         0.2122846889_8,  0.0000000000_8, -0.0501043723_8,  0.0000000000_8, &
        -0.0130114351_8,  0.0000000000_8,  0.0125509080_8,  0.0000000000_8, &
        -0.0012017700_8,  0.0000000000_8,  0.0030287292_8,  0.0000000000_8, &
         0.0000000000_8, -0.9559024842_8,  0.0000000000_8,  0.2239272982_8, &
         0.0000000000_8,  0.2243901016_8,  0.0000000000_8, -0.1470003238_8, &
         0.0000000000_8,  0.0270458605_8,  0.0000000000_8,  0.0172007266_8, &
         0.0000000000_8, -0.0086797958_8,  0.0000000000_8,  0.0025126889_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
        -0.1912987845_8,  0.0000000000_8, -0.1086802711_8,  0.0000000000_8, &
         0.2010699929_8,  0.0000000000_8, -0.0967379398_8,  0.0000000000_8, &
         0.0043451753_8,  0.0000000000_8,  0.0168857318_8,  0.0000000000_8, &
        -0.0074661180_8,  0.0000000000_8,  0.0006301034_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8, -0.0109729929_8, &
         0.0000000000_8,  0.1557726732_8,  0.0000000000_8, -0.1386856472_8, &
         0.0000000000_8,  0.0517458756_8,  0.0000000000_8,  0.0101544431_8, &
         0.0000000000_8, -0.0154559657_8,  0.0000000000_8,  0.0064325378_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8, -0.0486994534_8,  0.0000000000_8, &
         0.1235285184_8,  0.0000000000_8, -0.0835135104_8,  0.0000000000_8, &
         0.0162625365_8,  0.0000000000_8,  0.0149054017_8,  0.0000000000_8, &
        -0.0126868460_8,  0.0000000000_8,  0.0017314716_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000000000_8,  0.0373752378_8,  0.0000000000_8, -0.0708156744_8, &
         0.0000000000_8,  0.0395478713_8,  0.0000000000_8,  0.0040759606_8, &
         0.0000000000_8, -0.0151183617_8,  0.0000000000_8,  0.0101165866_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0176475310_8,  0.0000000000_8, -0.0301935457_8,  0.0000000000_8, &
         0.0063198080_8,  0.0000000000_8,  0.0131048984_8,  0.0000000000_8, &
        -0.0140511895_8,  0.0000000000_8,  0.0046796650_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8, -0.0048097188_8, &
         0.0000000000_8,  0.0085630646_8,  0.0000000000_8,  0.0069502206_8, &
         0.0000000000_8, -0.0111886606_8,  0.0000000000_8,  0.0096646598_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0021805188_8,  0.0000000000_8, &
        -0.0055073897_8,  0.0000000000_8,  0.0098023668_8,  0.0000000000_8, &
        -0.0102008844_8,  0.0000000000_8,  0.0053885132_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000000000_8, -0.0022572828_8,  0.0000000000_8,  0.0079247762_8, &
         0.0000000000_8, -0.0067578233_8,  0.0000000000_8,  0.0067532430_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
        -0.0030063404_8,  0.0000000000_8,  0.0050169581_8,  0.0000000000_8, &
        -0.0032598376_8,  0.0000000000_8,  0.0020270143_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0022702579_8, &
         0.0000000000_8, -0.0027121107_8,  0.0000000000_8,  0.0023131595_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0003359492_8,  0.0000000000_8, &
        -0.0006748572_8,  0.0000000000_8, -0.0003956004_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000000000_8, -0.0006581103_8,  0.0000000000_8,  0.0010377451_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000541920_8,  0.0000000000_8,  0.0000470481_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0006895980_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8, -0.0000199036_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000000000_8,  0.0000000000_8,  0.0000000000_8,  0.0000000000_8, &
         0.0000000000_8], SHAPE(cv)),                                       &
          t1(2,0:18) = 1, t2(2,0:18) = 0

t1(:,1) = uv
t2(:,1) = 1
DO j = 2, 18
  t1(:,j) = 2 * uv * t1(:,j-1) - t1(:,j-2)
  t2(:,j) = j * (t1(:,j-1) - uv * t1(:,j)) / (1 - uv**2)
END DO
dp(1) = DOT_PRODUCT (t1(2,:), MATMUL (cv, t2(1,:)))
dp(2) = DOT_PRODUCT (t2(2,:), MATMUL (cv, t1(1,:)))
END SUBROUTINE TL_LENS_HC

!------------------------------------------------------------------------------

SUBROUTINE TL_LENS_HC_OLD
REAL(8) :: cx(6,0:11) = RESHAPE (                                       &
             [1.744584870_8, 3.915534248_8, 1.817264786_8,              & ! 0
              0.515779967_8, 0.080654597_8, 0.005542067_8,              &
              4.905552783_8, 8.822928305_8, 4.452501570_8,              & ! 1
              1.402084759_8, 0.263014562_8, 0.024177168_8,              &
              4.048926045_8, 5.936355689_8, 2.737940944_8,              & ! 2
              0.726798799_8, 0.082721594_8, 0.000000000_8,              &
              3.662363909_8, 5.703729370_8, 2.743715052_8,              & ! 3
              0.790444767_8, 0.102232909_8, 0.000000000_8,              &
              1.925674828_8, 2.847366396_8, 1.164635213_8,              & ! 4
              0.186240796_8, 0.000000000_8, 0.000000000_8,              &
              1.591282990_8, 2.374614686_8, 1.039847407_8,              & ! 5
              0.188566576_8, 0.000000000_8, 0.000000000_8,              &
              0.579240381_8, 0.823497186_8, 0.196477046_8,              & ! 6
              0.000000000_8, 0.000000000_8, 0.000000000_8,              &
              0.426913517_8, 0.625204994_8, 0.177370332_8,              & ! 7
              0.000000000_8, 0.000000000_8, 0.000000000_8,              &
              0.110212063_8, 0.106232618_8, 0.000000000_8,              & ! 8
              0.000000000_8, 0.000000000_8, 0.000000000_8,              &
              0.069472335_8, 0.080552096_8, 0.000000000_8,              & ! 9
              0.000000000_8, 0.000000000_8, 0.000000000_8,              &
              0.011641683_8, 0.000000000_8, 0.000000000_8,              & ! 10
              0.000000000_8, 0.000000000_8, 0.000000000_8,              &
              0.007602317_8, 0.000000000_8, 0.000000000_8,              & ! 11
              0.000000000_8, 0.000000000_8, 0.000000000_8], SHAPE(cx)), &
           cy(6,0:11) = RESHAPE (                                       &
             [4.59969655_8,  3.98293509_8,  1.82523175_8,               & ! 0
              0.48509170_8,  0.10367230_8,  0.02452503_8,               &
              9.90605487_8,  6.11749715_8,  2.54179739_8,               & ! 1
              0.72684607_8,  0.10725870_8,  0.00000000_8,               &
              9.11532492_8,  6.00095724_8,  2.58916168_8,               & ! 2
              0.75848078_8,  0.11901313_8,  0.00000000_8,               &
              5.45091703_8,  3.38838928_8,  1.29298848_8,               & ! 3
              0.21786730_8,  0.00000000_8,  0.00000000_8,               &
              4.46769476_8,  2.80546465_8,  1.10404242_8,               & ! 4
              0.20323442_8,  0.00000000_8,  0.00000000_8,               &
              1.81231385_8,  0.99889104_8,  0.22282669_8,               & ! 5
              0.00000000_8,  0.00000000_8,  0.00000000_8,               &
              1.34555316_8,  0.76304853_8,  0.18833908_8,               & ! 6
              0.00000000_8,  0.00000000_8,  0.00000000_8,               &
              0.32030382_8,  0.11323196_8,  0.00000000_8,               & ! 7
              0.00000000_8,  0.00000000_8,  0.00000000_8,               &
              0.22531207_8,  0.08895966_8,  0.00000000_8,               & ! 8
              0.00000000_8,  0.00000000_8,  0.00000000_8,               &
              0.02408278_8,  0.00000000_8,  0.00000000_8,               & ! 9
              0.00000000_8,  0.00000000_8,  0.00000000_8,               &
              0.01624442_8,  0.00000000_8,  0.00000000_8,               & ! 10
              0.00000000_8,  0.00000000_8,  0.00000000_8,               &
              0.00000000_8,  0.00000000_8,  0.00000000_8,               & ! 11
              0.00000000_8,  0.00000000_8,  0.00000000_8], SHAPE(cy)),  &
           t1(2,0:11) = 1

! CALL t_polynomial (2, 11, uv, t1)
t1(:,1) = uv
DO j = 2, 11
  t1(:,j) = 2 * uv * t1(:,j-1) - t1(:,j-2)
END DO
dp(1) = DOT_PRODUCT (t1(2,0:10:2), MATMUL (cx, t1(1,:)))
dp(2) = DOT_PRODUCT (t1(2,1:11:2), MATMUL (cy, t1(1,:)))
END SUBROUTINE TL_LENS_HC_OLD
END

!******************************************************************************
