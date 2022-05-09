MODULE kinds
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = KIND(0d0)
  PUBLIC :: DP
END MODULE kinds

MODULE constants
  USE kinds
  IMPLICIT NONE
  REAL(KIND=DP) ,PARAMETER :: pi = ATAN(1.0_DP)*4.0_DP
  COMPLEX(KIND=DP), PARAMETER :: ci = (0.0_DP,1.0_DP)
  COMPLEX(KIND=DP), PARAMETER :: czero = (0.0_DP,0.0_DP)
  PUBLIC :: pi,ci,czero
END MODULE constants

MODULE input_variables
  USE kinds
  USE constants
  IMPLICIT NONE
  SAVE
  INTEGER nk_lin
  REAL(KIND=DP) lambda
  REAL(KIND=DP) beta
  REAL(KIND=DP) eps
  REAL(KIND=DP) hubbardu
  LOGICAL use_preset
  !
  PUBLIC :: nk_lin, lambda, beta, eps, hubbardu, use_preset
  !
CONTAINS
  !
  SUBROUTINE read_input
    USE kinds
    IMPLICIT NONE
    INTEGER ios
    REAL(KIND=DP) rdum

    NAMELIST /input/nk_lin, lambda, beta, eps, hubbardu, use_preset
    !
    ! default values 
    lambda = 1.0E5_DP
    beta = 1.0E3_DP
    eps = 1.0E-7_DP
    nk_lin = 256
    hubbardu = 2.0E0_DP
    use_preset = .true.
    !
    READ(5,NML=input,IOSTAT=ios)
    IF(ios.NE.0) THEN 
      PRINT*, 'ERROR in read_input: reading namelist input'
      STOP
    ENDIF
    !
    rdum = LOG(REAL(nk_lin, KIND=DP)) / LOG(2.0E0_DP)
    IF ( (rdum - NINT(rdum)) >  1.0E-8_DP ) THEN
      PRINT*, 'ERROR in read_input: nk_lin is not a integer power of 2.'
      STOP
    ENDIF
    !
    rdum = LOG(lambda) / LOG(1.0E1_DP)
    IF ( (rdum - NINT(rdum)) >  1.0E-8_DP ) THEN
      PRINT*, 'ERROR in read_input: lambda is not a integer power of 10.'
      STOP
    ENDIF
    !
    rdum = LOG(eps) / LOG(1.0E-1_DP)
    IF ( (rdum - NINT(rdum)) >  1.0E-8_DP ) THEN
      PRINT*, 'ERROR in read_input: eps is not a integer power of 0.1.'
      STOP
    ENDIF
    !
  END SUBROUTINE read_input
  !
END MODULE input_variables

MODULE fftw3_run
  USE kinds
  USE constants, ONLY : pi
  IMPLICIT NONE
  INCLUDE "fftw3.f"
  !
CONTAINS
  !
  SUBROUTINE do_fft_2d(n1, n2, in, out, num, isign)
    INTEGER, INTENT(IN) :: n1, n2, num
    COMPLEX(KIND=DP), INTENT(IN) :: in(n1, n2, num)
    COMPLEX(KIND=DP), INTENT(OUT) :: out(n1, n2, num)
    INTEGER, INTENT(IN) ::  isign
    INTEGER, PARAMETER:: rank = 2
    COMPLEX(KIND=DP) in_(n1, n2), out_(n1, n2)
    REAL(KIND=DP) rdum
    INTEGER i
    INTEGER n(rank) ! array of dimension
    INTEGER plan, dist
    !
    rdum = LOG(REAL(n1, KIND=DP)) / LOG(2.0E0_DP)
    IF ( (rdum - NINT(rdum)) >  1.0E-8_DP ) THEN
      PRINT*, 'ERROR in do_fft_2d: n1 is not a integer power of 2.'
      STOP
    ENDIF
    !
    rdum = LOG(REAL(n2, KIND=DP)) / LOG(2.0E0_DP)
    IF ( (rdum - NINT(rdum)) >  1.0E-8_DP ) THEN
      PRINT*, 'ERROR in do_fft_2d: n2 is not a integer power of 2.'
      STOP
    ENDIF
    !
    !IF (isign < 0) THEN
    !  CALL dfftw_plan_dft_2d(plan, n1, n2, in_, out_, FFTW_FORWARD, FFTW_ESTIMATE)
    !ELSEIF (isign > 0) THEN
    !  CALL dfftw_plan_dft_2d(plan, n1, n2, in_, out_, FFTW_BACKWARD, FFTW_ESTIMATE)
    !ELSE
    !  PRINT*, 'ERROR in do_fft_2d: invalid isign.'
    !  STOP
    !ENDIF
    !!
    !DO i = 1, num
    !  in_ = in(:,:,i)
    !  CALL dfftw_execute_dft(plan, in_, out_)
    !  out(:,:,i) = out_(:,:)
    !ENDDO
    !!
    !CALL dfftw_destroy_plan(plan)
    !
    n(1) = n1
    n(2) = n2
    !
    !IF (isign < 0) THEN
    !  CALL dfftw_plan_dft(plan, rank, n, in_, out_, FFTW_FORWARD, FFTW_ESTIMATE)
    !ELSEIF (isign > 0) THEN
    !  CALL dfftw_plan_dft(plan, rank, n, in_, out_, FFTW_BACKWARD, FFTW_ESTIMATE)
    !ELSE
    !  PRINT*, 'ERROR in do_fft_2d: invalid isign.'
    !  STOP
    !ENDIF
    !!
    !DO i = 1, num
    !  in_ = in(:,:,i)
    !  CALL dfftw_execute_dft(plan, in_, out_)
    !  out(:,:,i) = out_(:,:)
    !ENDDO
    !!
    !CALL dfftw_destroy_plan(plan)
    !
    dist = n(1) * n(2)
    !inembed = n
    !onembed = n
    !istride = 1
    !ostride = 1
    !
    IF (isign < 0) THEN
      CALL dfftw_plan_many_dft(plan, rank, n, num, in, n, 1, dist, out, n, 1, dist, FFTW_FORWARD, FFTW_ESTIMATE)
    ELSEIF (isign > 0) THEN
      CALL dfftw_plan_many_dft(plan, rank, n, num, in, n, 1, dist, out, n, 1, dist, FFTW_BACKWARD, FFTW_ESTIMATE)
    ELSE
      PRINT*, 'ERROR in do_fft_2d: invalid isign.'
      STOP
    ENDIF
    !
    CALL dfftw_execute_dft(plan, in, out)
    CALL dfftw_destroy_plan(plan)
    !
  END SUBROUTINE do_fft_2d
  !
END MODULE fftw3_run

PROGRAM main
  USE sparse_ir
  USE sparse_ir_io
  USE sparse_ir_preset
  USE kinds
  USE input_variables
  USE fftw3_run
  IMPLICIT NONE
  !
  TYPE(IR) :: ir_obj
  !
  !! ALLOCATABLE ARRAYS
  INTEGER, ALLOCATABLE :: smpl_matsu(:)
  INTEGER, ALLOCATABLE :: my_freqs(:)
  REAL(KIND=DP), ALLOCATABLE :: smpl_tau(:)
  REAL(KIND=DP), ALLOCATABLE :: k1(:,:), k2(:,:), ek_(:,:)
  REAL(KIND=DP), ALLOCATABLE :: kp(:,:), ek(:)
  COMPLEX(KIND=DP), ALLOCATABLE :: gkf(:,:), gkl(:,:), gkt(:,:), grt(:,:), grt_(:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: srt(:,:), srl(:,:), skl(:,:), uhat(:,:), res(:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: cdummy_k(:,:,:), cdummy_r(:,:,:)
  !
  !! FUNCTIONS
  COMPLEX(KIND=DP), EXTERNAL :: DZNRM2
  !
  INTEGER ntau, nw, size_l, kps1, kps2, nk !, nr
  INTEGER nlambda, ndigit
  REAL(KIND=DP) wmax
  REAL(KIND=DP) smpl_tau_cond, smpl_matsu_cond
  REAL(KIND=DP) inv_beta, inv_kps1, inv_kps2, inv_nk
  !
  INTEGER ios
  INTEGER ix, jx, ii
  INTEGER ik, idum
  REAL(KIND=DP) rdum
  COMPLEX(KIND=DP) iv
  CHARACTER(LEN=100) filename
  !
  CALL read_input
  !
  wmax = lambda / beta
  inv_beta = 1.E0_DP / beta
  !
  rdum = LOG(lambda) / LOG(1.0E1_DP)
  IF ( (rdum - NINT(rdum)) >  1.0E-8_DP ) THEN
    PRINT*, 'ERROR in main: lambda is not a integer power of 10.'
    STOP
  ENDIF
  nlambda = NINT(rdum)
  !
  rdum = LOG(eps) / LOG(1.0E-1_DP)
  IF ( (rdum - NINT(rdum)) >  1.0E-8_DP ) THEN
    PRINT*, 'ERROR in main: lambda is not a integer power of 0.1.'
    STOP
  ENDIF
  ndigit = NINT(rdum)
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!   STEP 1: Store the data of IR basis output by the python modules   !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  IF (use_preset) THEN
    ir_obj = mk_ir_preset(nlambda, ndigit, beta)
  ELSE
    WRITE(filename, *) "ir_nlambda", nlambda, "_ndigit", ndigit, ".dat"
    OPEN(99, FILE=filename, STATUS='old',IOSTAT=ios)
    IF(ios.NE.0) PRINT*, 'ERROR in main while opening file: ', filename
    ir_obj = read_ir(99, beta)
    CLOSE(99)
  ENDIF
  !
  IF (ABS(ir_obj%beta - beta) > 1E-10_DP) THEN
    PRINT*, 'ERROR in main: beta of ir_obj does not match input beta.'
    STOP
  ENDIF
  IF (ABS(ir_obj%wmax - wmax) > 1E-10_DP) THEN
    PRINT*, 'ERROR in main: wmax of ir_obj does not match input wmax.'
    STOP
  ENDIF
  !
  size_l = ir_obj%size
  !
  ntau = ir_obj%ntau
  ALLOCATE(smpl_tau(ntau))
  smpl_tau(:) = ir_obj%tau(:)
  smpl_tau_cond = ir_obj%u%inv_s(ir_obj%u%ns) / ir_obj%u%inv_s(1)
  !
  !PRINT*, 1.0_DP / ir_obj%u%inv_s(:)
  PRINT*, "cond (tau): ", smpl_tau_cond
  !
  nw = ir_obj%nfreq_f
  ALLOCATE(smpl_matsu(nw))
  smpl_matsu(:) = ir_obj%freq_f(:)
  smpl_matsu_cond = ir_obj%uhat_f%inv_s(ir_obj%uhat_f%ns) / ir_obj%uhat_f%inv_s(1)
  !
  PRINT*, "cond (matsu): ", smpl_matsu_cond
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!   STEP 2: Compute the non-interacting band dispersion on a mesh   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  kps1 = nk_lin
  kps2 = nk_lin
  !
  nk = kps1 * kps2
  ! nr = nk
  !
  inv_nk = 1.0E0_DP / REAL(nk, KIND=DP)
  inv_kps1 = 1.0E0_DP / REAL(kps1, KIND=DP)
  inv_kps2 = 1.0E0_DP / REAL(kps2, KIND=DP)
  !
  ALLOCATE(k1(kps1, kps2), k2(kps1, kps2), ek_(kps1, kps2))
  !
  DO jx = 1, kps2
    k2(:, jx) = 2.0E0_DP * pi * REAL((jx - 1), KIND=DP) * inv_kps2
    DO ix = 1, kps1
      k1(ix, jx) = 2.0E0_DP * pi * REAL((ix - 1), KIND=DP) * inv_kps1
    ENDDO
  ENDDO
  !
  DO jx = 1, kps2
    DO ix = 1, kps1
      ek_(ix, jx) = -2.0E0_DP * ( COS(k1(ix, jx)) + COS(k2(ix, jx)) )
    ENDDO
  ENDDO
  !
  PRINT*, "Shape of k1: (", SHAPE(k1), ")"
  PRINT*, "Shape of k1: (", SHAPE(k2), ")"
  PRINT*, "Shape of ek_: (", SHAPE(ek_), ")"
  !
  ALLOCATE(kp(2, nk), ek(nk))
  !
  ik = 0
  DO jx = 1, kps2
    DO ix = 1, kps1
      ik = ik + 1
      kp(1, ik) = k1(ix, jx)
      kp(2, ik) = k2(ix, jx)
      ek(ik) = ek_(ix, jx)
    ENDDO
  ENDDO
  !
  DEALLOCATE(k1, k2, ek_)
  !
  PRINT*, "Shape of k1: (", SHAPE(kp), ")"
  PRINT*, "Shape of ek: (", SHAPE(ek), ")"
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!   STEP 2: Compute non-interacting Green's function on sampling frequencies   !!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !! G(k, iv): (nk, nw)
  ALLOCATE(gkf(nk, nw))
  !
  DO ix = 1, nw
    DO ik = 1, nk
      iv = ci * pi * REAL(ir_obj%freq_f(ix), KIND=DP) * inv_beta
      gkf(ik, ix) = 1.0E0_DP / (iv - ek(ik))
    ENDDO
  ENDDO
  !
  filename = "second_order_gkf_gamma.txt"
  OPEN(100,FILE=filename,STATUS="replace",IOSTAT=ios)
  IF(ios.NE.0) PRINT*, 'ERROR in main while opening file: ', filename
  WRITE(100,'(A)') "#       v      iv              Re(G(k,iv))      Im(G(k,iv))      abs(G(k,iv))"
  DO ix = 1, nw
    iv = pi * REAL(ir_obj%freq_f(ix), KIND=DP) * inv_beta
    WRITE(100,"(I10, 1x, 3(ES16.8, 1x), ES16.8)") ir_obj%freq_f(ix), REAL(iv, KIND=DP), gkf(1, ix), ABS(gkf(1, ix))
  ENDDO
  CLOSE(100)
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!   STEP 3: Evaluate the Green's function on sampling points in the imaginary-time domain   !!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !! G(k, l): (nk, size)
  ALLOCATE(gkl(nk, size_l))
  !
  CALL fit_matsubara_f(ir_obj, gkf, gkl)
  !
  filename = "second_order_gkl_gamma.txt"
  OPEN(101,FILE=filename,STATUS="replace",IOSTAT=ios)
  IF(ios.NE.0) PRINT*, 'ERROR in main while opening file: ', filename
  WRITE(101,'(A)') "#    l      s_l            Re(G(k,l))       Im(G(k,l))      abs(G(k,l))"
  DO ix = 1, size_l
    WRITE(101,"(I6, 3(ES16.8, 1x), ES16.8)") ix, ir_obj%s(ix), gkl(1, ix), ABS(gkl(1, ix))
  ENDDO
  CLOSE(101)
  !
  ALLOCATE(gkt(nk, ntau))
  !
  CALL evaluate_tau(ir_obj, gkl, gkt)
  !
  filename = "second_order_gkt_gamma.txt"
  OPEN(102,FILE=filename,STATUS="replace",IOSTAT=ios)
  IF(ios.NE.0) PRINT*, 'ERROR in main while opening file: ', filename
  WRITE(102,'(A)') "#     tau           Re(G(k,tau))     Im(G(k,tau))    abs(G(k,tau))"
  DO ix = 1, ntau
    WRITE(102,"(3(ES16.8, 1x), ES16.8)") ir_obj%tau(ix), gkt(1, ix), ABS(gkt(1, ix))
  ENDDO
  CLOSE(102)
  !
  filename = "second_order_gkt_m.txt"
  OPEN(103,FILE=filename,STATUS="replace",IOSTAT=ios)
  IF(ios.NE.0) PRINT*, 'ERROR in main while opening file: ', filename
  WRITE(103,'(A)') "#     tau           Re(G(k,tau))     Im(G(k,tau))    abs(G(k,tau))"
  idum = kps1 * (kps2 + 1) / 2 + 1
  DO ix = 1, ntau
    WRITE(103,"(3(ES16.8, 1x), ES16.8)") ir_obj%tau(ix), gkt(idum, ix), ABS(gkt(idum, ix))
  ENDDO
  CLOSE(103)
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!   STEP 4: Transform the Green's function to the real space and evaluate the self-energy   !!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !! Compute G(r, tau): (nk, ntau)
  !!  (1) Reshape gkt into shape of (nk_lin, nk_lin, ntau).
  !!  (2) Apply FFT to the 1st and 2nd axes.
  !!  (3) Reshape the result to (ntau, nk)
  !! G(k, tau): (nk, ntau)
  ALLOCATE(cdummy_k(kps1, kps2, ntau))
  ALLOCATE(cdummy_r(kps1, kps2, ntau))
  !
  DO ii = 1, ntau
    ik = 0
    DO jx = 1, kps2
      DO ix = 1, kps1
        ik = ik + 1
        cdummy_k(ix, jx, ii) = gkt(ik, ii)
      ENDDO
    ENDDO
  ENDDO
  !
  CALL do_fft_2d(kps1, kps2, cdummy_k, cdummy_r, ntau, -1) ! -1: forward FFT transform
  !
  DEALLOCATE(cdummy_k)
  ALLOCATE(grt(nk, ntau))
  !
  DO ii = 1, ntau
    ik = 0
    DO jx = 1, kps2
      DO ix = 1, kps1
        ik = ik + 1
        grt(ik, ii) = cdummy_r(ix, jx, ii) * inv_nk
      ENDDO
    ENDDO
  ENDDO
  !
  DEALLOCATE(cdummy_r)
  !
  filename = "second_order_grt_r_0.txt"
  OPEN(104,FILE=filename,STATUS="replace",IOSTAT=ios)
  IF(ios.NE.0) PRINT*, 'ERROR in main while opening file: ', filename
  WRITE(104,'(A)') "#     tau           Re(G(r,tau))     Im(G(r,tau))    abs(G(r,tau))"
  DO ix = 1, ntau
    WRITE(104,"(3(ES16.8, 1x), ES16.8)") ir_obj%tau(ix), grt(1, ix), ABS(grt(1, ix))
  ENDDO
  CLOSE(104)
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!   STEP 5: Compute the second-order term of the self-energy on the sampling points   !!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !! Sigma(r, tau): (nr, ntau)
  ALLOCATE(srt(nk, ntau))
  !
  !! G(k, beta - tau): (nr, ntau)
  ALLOCATE(grt_(nk, ntau))
  !
  DO ix = 1, ntau
    grt_(:, ntau - ix + 1) = grt(:, ix)
  ENDDO
  !
  srt = grt * grt * grt_ * hubbardu * hubbardu
  !
  DEALLOCATE(grt_)
  !
  filename = "second_order_srt_r_0.txt"
  OPEN(105,FILE=filename,STATUS="replace",IOSTAT=ios)
  IF(ios.NE.0) PRINT*, 'ERROR in main while opening file: ', filename
  WRITE(105,'(A)') "#     tau          Re(Sigma(r,tau)) Im(Sigma(r,tau)) abs(Sigma(r,tau))"
  DO ix = 1, ntau
    WRITE(105,"(3(ES16.8, 1x), ES16.8)") ir_obj%tau(ix), srt(1, ix), ABS(srt(1, ix))
  ENDDO
  CLOSE(105)
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!   STEP 6: Transform the self-energy to the IR basis and then transform it to the k space   !!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !! Sigma(r, l): (nr, size)
  ALLOCATE(srl(nk, size_l))
  !
  CALL fit_tau(ir_obj, srt, srl)
  !
  filename = "second_order_srl_r_0.txt"
  OPEN(106,FILE=filename,STATUS="replace",IOSTAT=ios)
  IF(ios.NE.0) PRINT*, 'ERROR in main while opening file: ', filename
  WRITE(106,'(A)') "#    l   Re(Sigma(r,l))   Im(Sigma(r,l))  abs(Sigma(r,l))"
  DO ix = 1, size_l
    WRITE(106,"(I6, 2(ES16.8, 1x), ES16.8)") ix, srl(1, ix), ABS(srl(1, ix))
  ENDDO
  CLOSE(106)
  !
  ALLOCATE(cdummy_r(kps1, kps2, size_l))
  ALLOCATE(cdummy_k(kps1, kps2, size_l))
  !
  DO ii = 1, size_l
    ik = 0
    DO jx = 1, kps2
      DO ix = 1, kps1
        ik = ik + 1
        cdummy_r(ix, jx, ii) = srl(ik, ii)
      ENDDO
    ENDDO
  ENDDO
  !
  CALL do_fft_2d(kps1, kps2, cdummy_r, cdummy_k, size_l, 1) ! +1: backward FFT transform
  !
  DEALLOCATE(cdummy_r)
  !
  !! Sigma(k, l): (nk, size)
  ALLOCATE(skl(nk, size_l))
  !
  DO ii = 1, size_l
    ik = 0
    DO jx = 1, kps2
      DO ix = 1, kps1
        ik = ik + 1
        skl(ik, ii) = cdummy_k(ix, jx, ii)
      ENDDO
    ENDDO
  ENDDO
  !
  filename = "second_order_skl_gamma.txt"
  OPEN(107,FILE=filename,STATUS="replace",IOSTAT=ios)
  IF(ios.NE.0) PRINT*, 'ERROR in main while opening file: ', filename
  WRITE(107,'(A)') "#    l      s_l          Re(Sigma(k,l))   Im(Sigma(k,l))   abs(Sigma(k,l))"
  DO ix = 1, size_l
    WRITE(107,"(I6, 3(ES16.8, 1x), ES16.8)") ix, ir_obj%s(ix), skl(1, ix), ABS(skl(1, ix))
  ENDDO
  CLOSE(107)
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!   STEP 7: Evaluate the self-energy on arbitrary frequencies/times   !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ALLOCATE(my_freqs(nw))
  my_freqs = ir_obj%freq_f
  !
  ALLOCATE(uhat(nw, size_l))
  uhat = ir_obj%uhat_f%a
  !
  ! NOTE: nw     = ir_obj%nfreq_f = ir_obj%uhat_f%m
  !       size_l = ir_obj%size    = ir_obj%uhat_f%n
  !
  ALLOCATE(res(nk, nw))
  res(:,:) = czero
  !
  DO ii = 1, size_l
    DO ix = 1, nw
      DO ik = 1, nk
        res(ik, ix) = res(ik, ix) + uhat(ix, ii) * skl(ik, ii)
      ENDDO
    ENDDO
  ENDDO
  !
  filename = "second_order_skf_gamma.txt"
  OPEN(108,FILE=filename,STATUS="replace",IOSTAT=ios)
  IF(ios.NE.0) PRINT*, 'ERROR in main while opening file: ', filename
  WRITE(108,'(A)') "#       v      iv            Re(Sigma(k,iv))  Im(Sigma(k,iv))   abs(Sigma(k,iv))"
  DO ix = 1, nw
    iv = pi * REAL(ir_obj%freq_f(ix), KIND=DP) * inv_beta
    WRITE(108,"(I10, 1x, 3(ES16.8, 1x), ES16.8)") ir_obj%freq_f(ix), REAL(iv, KIND=DP), res(1, ix), ABS(res(1, ix))
  ENDDO
  CLOSE(108)
  !
END PROGRAM ! main