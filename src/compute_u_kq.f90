!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE compute_u_kq(ik, q)
  !----------------------------------------------------------------------------
  !
  ! ... diagonalize at k+q
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : RytoeV, tpi
  USE io_global,            ONLY : stdout
  USE io_files,             ONLY : nwordwfcU, iunhub, iunwfc, nwordwfc
  USE mp,                   ONLY : mp_sum
  USE mp_pools,             ONLY : inter_pool_comm, intra_pool_comm, me_pool
  USE mp_bands,             ONLY : intra_bgrp_comm
#ifdef __BANDS
  USE mp_bands,             ONLY : me_bgrp, inter_bgrp_comm
#endif
  USE klist,                ONLY : nkstot, nks, xk, ngk, igk_k, igk_k_d
  USE uspp,                 ONLY : vkb, vkb_d, nkb, okvan
  USE wvfct,                ONLY : et, nbnd, g2kin, g2kin_d, &
                                   current_k, nbndx, btype, npwx
  USE gvecw,                ONLY : gcutw
  USE control_flags,        ONLY : ethr, lscf, istep, max_cg_iter, david
  USE control_flags,        ONLY : cntrl_isolve => isolve
  USE ldaU,                 ONLY : lda_plus_u, wfcU
  USE lsda_mod,             ONLY : current_spin, lsda, isk, nspin
  USE wavefunctions_module, ONLY : evc, evc_d
  USE gvect,                ONLY : g, ngm, gstart
  USE gvecs,                ONLY : doublegrid
  USE dfunct,               ONLY : newd
  USE cell_base,            ONLY : tpiba2
  USE random_numbers,       ONLY : randy
  USE scf,                  ONLY : v, vrs, vltot, kedtau, rho
  USE fft_base,             ONLY : dfftp
  USE noncollin_module,     ONLY : npol
  USE becmod,               ONLY : becp, allocate_bec_type, deallocate_bec_type
  USE buffers
  USE gipaw_module
  USE cpu_gpu_interface
  USE nvtx
  IMPLICIT NONE
  INTEGER :: ik, iter       ! k-point, current iterations
  REAL(DP) :: q(3)          ! q-vector
  REAL(DP) :: avg_iter
  INTEGER :: ig, i
  REAL(DP) :: xkold(3)
  REAL(DP), allocatable :: et_old(:,:)
  REAL(DP) :: rr, arg
  INTEGER :: npw

  CALL start_clock( 'c_bands' )

  ! Initialize the diagonalization
  if (isolve == 1 .or. isolve == 2) then
    nbndx = nbnd ! CG or PPCG
  elseif (isolve == 0) then
    nbndx = david*nbnd ! Davidson TODO: check if 4 times!!!!
  else
    call errore('compute_u_kq', 'wrong isolve', 1)
  endif

  cntrl_isolve = isolve
  max_cg_iter = 200
  iter = 1
  istep = 0
  ethr = conv_threshold
  lscf = .false.
  if (allocated(btype)) deallocate(btype)
  allocate(btype(nbndx,nkstot))
  btype(1:nbnd,:) = 1

  ! save eigenvalues
  allocate( et_old(nbnd,nkstot) )
  et_old = et

  !! debug
  if (iverbosity > 10) &
    write(stdout, '(5X,"compute_u_kq: q = (",F10.4,",",F10.4,",",F10.4,")")') q

  avg_iter = 0.D0

  current_k = ik
  IF ( lsda ) current_spin = isk(ik)
  npw = ngk(ik)

  ! same sorting of G-vector at k+q
  call gk_sort(xk(1,ik),ngm,g,gcutw,npw,igk_k(1,ik),g2kin)

  ! set the k-point
  xkold(:) = xk(:,ik)
  xk(:,ik) = xk(:,ik) + q(:)
  g2kin(1:npw) = ( ( xk(1,ik) + g(1,igk_k(1:npw,ik)) )**2 + &
                   ( xk(2,ik) + g(2,igk_k(1:npw,ik)) )**2 + &
                   ( xk(3,ik) + g(3,igk_k(1:npw,ik)) )**2 ) * tpiba2
#ifdef USE_CUDA
  g2kin_d = g2kin
#endif

  ! various initializations
#ifdef USE_CUDA
  if (nkb > 0) call init_us_2_gpu( npw, igk_k_d(1,ik), xk(1,ik), vkb_d )
#else
  if (nkb > 0) call init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb )
#endif
  if (lda_plus_U) call orthoatwfc1(ik)

  ! read in wavefunctions from the previous iteration
  CALL get_buffer( evc, nwordwfc, iunwfc, ik)
#ifdef USE_CUDA
  evc_d = evc
#endif

#ifdef __BANDS
  ! not needed anymore??
  !!call mp_sum(evc, inter_bgrp_comm)
#endif

#if 0
  ! randomize a little bit in case of CG diagonalization
  if ( isolve == 1 ) then
#ifdef __BANDS
    rr = randy(ik+nks*me_bgrp) ! starting from a well defined k-dependent seed
#else
    rr = randy(ik+nks*me_pool) ! starting from a well defined k-dependent seed
#endif
    do i = 1, nbnd
      do ig = 1, npw
        rr = r_rand*(2.d0*randy() - 1.d0)
        arg = tpi * randy()
        evc(ig,i) = evc(ig,i)*CMPLX(1.d0+rr*cos(arg),rr*sin(arg),kind=DP)
      enddo
    enddo
  endif
#endif

  ! re-update potential
  !call setlocal
  !call plugin_scf_potential(rho, .false., -1d18)
  !call set_vrs (vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid)

  ! compute the D for the pseudopotentials
  !call newd

  ! diagonalization of bands for k-point ik
#if 0
  CALL allocate_bec_type ( nkb, nbnd, becp, intra_bgrp_comm )
  write(stdout,'(5X,''Rotating WFCS'')')
  CALL rotate_wfc ( npwx, npw, nbnd, gstart, nbnd, evc, npol, okvan, evc, et(1,ik) )
  CALL deallocate_bec_type ( becp )
#endif
  avg_iter = 1.d0
  call diag_bands ( iter, ik, avg_iter )

  !! debug
  if (iverbosity > 20) &
    write(stdout,'(5X,"ethr = ",1PE9.2,",  avg # of iterations =",0PF5.1)') &
         ethr, avg_iter

  ! check if diagonalization was ok
  if (iverbosity > 20) then
    write(stdout,'(5X,''eigenvalues at k:'')')
    write(stdout,'(8F9.4)') et_old(1:nbnd,ik)*RytoeV
    write(stdout,'(5X,''eigenvalues at k+q:'')')
    write(stdout,'(8F9.4)') et(1:nbnd,ik)*RytoeV
  endif

  do i = 1, nbnd
    if (abs(et(i,ik) - et_old(i,ik))*RytoeV > 0.2d0) then
      write(stdout,'(5X,''ATTENTION: ik='',I4,''  ibnd='',I3,$)') ik, i
      write(stdout,'(2X,''eigenvalues differ too much!'')')
      write(stdout,'(5X,2(F10.4,2X))') et_old(i,ik)*RytoeV, et(i,ik)*RytoeV
    endif
  enddo

  ! restore the k-point and eigenvalues
  xk(:,ik) = xkold(:)
  etq(:,ik) = et_old(:,ik)
  et = et_old
  deallocate(et_old)

  ! restore wavefunctions
  evq = evc
  CALL get_buffer(evc, nwordwfc, iunwfc, ik)
#ifdef __BANDS
  ! not needed anymore??
  !!call mp_sum(evc, inter_bgrp_comm)
#endif

  CALL stop_clock( 'c_bands' )  
  RETURN

END SUBROUTINE compute_u_kq
