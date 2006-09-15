! File: spanlib.f90
!
! This file is part of the SpanLib library.
! Copyright (C) 2006  Stephane Raynaud
! Contact: stephane dot raynaud at gmail dot com
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.

! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

module spanlib

contains

	! ############################################################
	! ############################################################
	! ## PCA PART ################################################
	! ############################################################
	! ############################################################

	subroutine sl_pca(ff, nkeep, xeof, pc, ev, weights, useteof, bLargeMatrix)

	!	Principal Component Analysis
	!
	! Description:
	!	Perform a decomposition of space-time field in a set of
	!	Empirical Orthogonal Functions (EOFs) and Principal components (PCs).
	!	The input data set can be optionally weighted in space.
	!	By default, the analysis computes  "temporal" (T) or classical
	!	spatial (S) EOFs depending on if the space dimension is greater
	!	than the time dimension. This default behavior can be overridden.
	!
	! Necessary arguments:
	!	- ff:			Space-time array
	!	- nkeep:		Maximum number of modes to keep in outputs
	!
	! Optional arguments:
	!	- xeof:		Space-mode array of EOFs
	!	- pc:			Time-mode array of PCs
	!	- ev:			Mode array of eigen values (variances)
	!	- weights:	Space array of weights
	!	- useteof:	To force the use of T or S EOFs [0 = T, 1 = S, -1 = default]
	!	- bLargeMatrix: Use la_syevd instead of la_syev (faster for large matrices, but uses more workspace) [default:.true.]
	!
	! Dependencies:
	!	sgemm(BLAS) ssyrk(BLAS) la_syev(LAPACK95) la_syevd(LAPACK95)


	! Declarations
	! ============

	use la_precision, only: wp => sp
	use f95_lapack, only: la_syevd, la_syev

	implicit none

	! External
	! --------
	real,    intent(in)           :: ff(:,:)
	integer, intent(in)	         :: nkeep
	real,    intent(out),optional :: pc(size(ff,2),nkeep), &
	&                                xeof(size(ff,1),nkeep), ev(nkeep)
	real,    intent(in), optional :: weights(:)
	integer, intent(in), optional :: useteof
	logical, intent(in), optional :: bLargeMatrix

	! Internal
	! --------
	integer           :: ii,ij,nn,ns,nt
	real, allocatable :: cov(:,:)
	real, allocatable :: wff(:,:), ww(:), zeof(:,:), zff(:,:)
	real, allocatable :: zev(:)
	integer           :: zuseteof, znkeepmax, i,j
	logical           :: zbLargeMatrix

	! Setups
	! ======

	! Sizes
	! -----
	ns = size(ff,1)
	nt = size(ff,2)
	znkeepmax = 100
	if(nkeep>znkeepmax)then
		print*,'[pca] You want to keep a number of PCs '//&
		 & 'greater than ',znkeepmax
		return
	end if

	! What does the user want?
	! ------------------------
	if(.not.present(xeof).and..not.present(pc)&
	  &.and..not.present(ev))then
		print*,'[pca] Nothing to do. Quit.'
		return
	end if

	! By default, T-EOF decompostion if ns > nt
	! -----------------------------------------
	zuseteof = -1
	if(present(useteof))zuseteof = useteof
	if(zuseteof<=0)then
		if(ns>nt)then
			zuseteof=1
		else
			zuseteof=0
		endif
	endif
	znkeepmax=100
	if(zuseteof)then
		if(nkeep>nt)then
			print*,'[pca] You want to keep a number of PCs '//&
				&'greater than the number of EOF:',nt
			return
		end if
	else
		if(nkeep>ns)then
			print*,'[pca] You want to keep a number of PCs '//&
				&'greater than the number of EOF:',ns
			return
		end if
	end if

	! Use ssyevd?
	! -----------
	if(.not.present(bLargeMatrix))then
		zbLargeMatrix = .true.
	else
		zbLargeMatrix = bLargeMatrix
	end if

	! Remove the mean
	! ---------------
	allocate(zff(ns,nt))
	zff = ff - spread(sum(ff,dim=2)/real(nt), ncopies=nt, dim=2)

	! Default weights = 1.
	! --------------------
	allocate(ww(ns))
	allocate(wff(ns,nt))
	ww = 1.
	if(present(weights))then
		ww(:) = weights * real(ns) / sum(weights)
		where(ww==0.)
			ww = 1.
		end where
		do i = 1, nt
			wff(:,i) = zff(:,i) * sqrt(ww)
		end do
	else
		wff = zff
	end if


	! EOF decomposition
	! =================

	if(zuseteof==1)then


		! T-EOF case
		! ----------

		! Covariance
		allocate(cov(nt,nt))
		allocate(zev(nt))
		cov=0.
		call ssyrk('U','T',nt,ns,1.,wff,ns, 0.,cov,nt)
		cov = cov / float(ns)
		deallocate(wff)

		! Diagonalising (cov: input=cov, output=eof)
		if(zbLargeMatrix)then
			call la_syevd(cov,zev,jobz='V')
		else
			call la_syev(cov,zev,jobz='V')
		end if

		! Back to S-EOFs
		if(present(pc).or.present(xeof))then
			allocate(zeof(ns,nkeep))
			call sgemm('N','N',ns,nkeep,nt,1.,zff,ns, &
				& cov(:,nt:nt-nkeep+1:-1),nt,0.,zeof,ns)
			deallocate(cov)
			do i = 1, nkeep
				zeof(:,i) = zeof(:,i) / &
				 &          sqrt(dot_product(ww(:), zeof(:,i)**2))
			end do
			if(.not.present(pc)) deallocate(ww)
		else
			deallocate(cov)
		end if

		! Eigenvalues
		! -----------
		if(present(ev))ev = zev(nt:nt-nkeep+1:-1)

	else

		! S-EOF case (classical)
		! ----------------------

		! Covariance
		allocate(cov(ns,ns))
		allocate(zev(ns))
		cov = 0.
		call ssyrk('U','N',ns,nt,1.,wff,ns, 0.,cov,ns)
		cov = cov / float(nt)
		deallocate(wff)

		! Diagonalisation (cov: input=cov, output=eof)
		if(zbLargeMatrix)then
			call la_syevd(cov,zev,jobz='V')
		else
			call la_syev(cov,zev,jobz='V')
		end if

		! Formatting S-EOFs
		if(present(xeof).or.present(pc))then
			allocate(zeof(ns,nkeep))
			do i = 1, nkeep
				zeof(:,i) = cov(:,ns-i+1) / sqrt(ww(:))
			end do
		end if
	   deallocate(cov)

		! Eigenvalues
		! -----------
		if(present(ev)) ev = zev(ns:ns-nkeep+1:-1)

	end if

	! Free eof array
	! --------------
	if(present(xeof))then
		xeof = zeof
		if(.not.present(pc)) deallocate(zeof)
	end if

	! Finally get PCs
	! ===============
	if(present(pc))then
		if(present(weights))then
			do i=1, nt
				zff(:,i) = zff(:,i) * ww(:)
			end do
		end if
		! pc = matmul( transpose(zff), zeof)
		call sgemm('T','N',nt,nkeep,ns,1.,zff,ns, &
			& zeof,ns,0.,pc,nt)
		do i = 1, nkeep
			pc(:,i) = pc(:,i) / dot_product(zeof(:,i)**2, ww(:))
		end do
	end if

	end subroutine sl_pca


	subroutine sl_pcarec(xeof, pc, ffrec, istart, iend)

	! Title:
	!	Reconstruction of a set of PCA components
	!
	! Description:
	!	Perform a reconstruction using a set of components previously
	!	computed with a PCA. All the reconstructed components are summed.
	!	A reconstructed component is simply the "product" of an EOF
	!	by its PC. The sum of all reconstructed component is the original field.
	!
	! Necessary arguments:
	!	- xeof:		Space-mode array of EOFs
	!	- pc:		Time-mode array of PCs
	!	- ffrec:	Space-time array of the reconstructed field
	!
	! Optional arguments:
	!	- istart:	Index of the first component to use
	!	- iend:		Index of the last component to use

	implicit none


	! Declarations
	! ============

	! External
	! --------
	real,	intent(in)              :: xeof(:,:), pc(:,:)
	real,	intent(out)             :: ffrec(size(xeof,1),size(pc,1))
	integer,intent(in),	optional	:: istart, iend

	! Internal
	! --------
	integer           :: nkept, itmp, zistart, ziend, nt, ns, i
	real, allocatable	:: zpc(:,:)


	! Setup
	! =====
	nkept = size(xeof,2)
	zistart=1
	if(present(istart))zistart=istart
	if(present(iend))then
		ziend=iend
	else
		ziend=nkept
	end if
	if(zistart.lt.1.or.zistart.gt.nkept)then
		zistart=1
		print*,'[pcarec] istart lower than 1 => set to 1'
	end if
	if(ziend.lt.1.or.ziend.gt.nkept)then
		ziend=nkept
		print*,'[pcarec] iend greater than the number '//&
			&'of avalaible modes => reduced to ',ziend
	end if
	if(zistart>ziend)then
		itmp=ziend
		ziend=zistart
		zistart=itmp
		print*,'[pcarec] istart > iend => inversion'
	end if
	ns = size(xeof,1)
	nt = size(pc,1)


	! Computation
	! ===========
	ffrec = 0.
	if(nt<ns) then
		do i = 1, nt
			ffrec(:, i) = ffrec(:, i) + &
				&	matmul(xeof(:, zistart:ziend), pc(i, zistart:ziend))
		end do
	else
		allocate(zpc(ziend-zistart+1, nt))
		zpc = transpose(pc(:, zistart:ziend))
		do i = 1, ns
			ffrec(i, :) = ffrec(i, :) + &
				&	matmul(xeof(i, zistart:ziend), zpc)
		end do
	end if

	end subroutine sl_pcarec




  !############################################################
  !############################################################
  !## MSSA PART ###############################################
  !############################################################
  !############################################################

	subroutine sl_mssa(ff, nwindow, nkeep, steof, stpc, ev, bLargeMatrix)

	! Title:
	!	Multi-channel Singular Spectrum Analysis
	!
	! Description:
	!	Perform a decomposition of space-time field in a set of
	!	space-time Empirical Orthogonal Functions (EOFs) and
	!	time Principal components (PCs), according to a window
	!	parameter.
	!
	! Necessary arguments:
	!	- ff:      Space-time array
	!	- nwindow: Window size
	!	- nkeep:   Maximum number of modes to keep in outputs
	!
	! Optional arguments:
	!	- steof: SpaceXwindow-mode array of EOFs
	!	- stpc:  Time-mode array of PCs
	!	- ev:    Mode array of eigen values (variances)
	!	- bLargeMatrix: Use ssyevd instead of ssyev (faster for large matrices, but uses more workspace) [default:.true.]
	!
	! Dependencies:
	!	sl_diasym


	! Declarations
	! ============

	use la_precision, only: wp => sp
	use f95_lapack, only: la_syevd, la_syev

	implicit none

	! External
	! --------
	real,   intent(in)            :: ff(:,:)
	integer,intent(in)            :: nwindow, nkeep
	real,   intent(out), optional :: steof(size(ff,1)*nwindow, nkeep), &
		& stpc(size(ff,2)-nwindow+1, nkeep), ev(nkeep)
	logical, intent(in), optional :: bLargeMatrix

	! Internal
	! --------
	real, allocatable :: cov(:,:), zev(:), zff(:,:), zsteof(:,:), wpc(:)
	real :: wsteof
	integer :: nchan, nsteof, nt, znkeepmax
	integer :: iw, iw1, iw2, i1, i2, im, ic1, ic2
 	logical :: zbLargeMatrix


	! Setup
	! =====

	! Sizes
	! -----
	nchan = size(ff,1)
	nsteof = nchan * nwindow
	nt = size(ff,2)
	znkeepmax = 100
	if(nkeep>znkeepmax)then
		print*,'[mssa] You want to keep a number of PCs '//&
		 & 'greater than ',znkeepmax
		return
	else if(nkeep>nsteof) then
		print*,'[mssa] You want to keep a number of PCs greater '// &
			& 'than the number of ST-EOFs:',nsteof
		return
	end if

	! Use ssyevd?
	! -----------
	if(.not.present(bLargeMatrix))then
		zbLargeMatrix = .true.
	else
		zbLargeMatrix = bLargeMatrix
	end if

	! Remove the mean
	! ---------------
	allocate(zff(nchan, nt))
	zff = ff - spread(sum(ff,dim=2)/real(nt), ncopies=nt, dim=2)

	! Set the block-Toeplitz covariance matrix
	! ========================================
	allocate(cov(nsteof, nsteof))
	do ic1 = 1, nchan
		do ic2 = 1, nchan
			do iw2 = 1, nwindow
				do iw1 = 1, iw2
					i1 = (ic1-1) * nwindow + iw1
					i2 = (ic2-1) * nwindow + iw2
					iw = iw2 - iw1 + 1
					cov(i1,i2) = &
						& dot_product(zff(ic1, 1  : nt-iw+1),  &
						&             zff(ic2, iw : nt	 )) / &
						& real(nt-iw+1)
					cov(i2,i1) = cov(i1,i2)
				end do
			end do
		end do
	end do

	! Diagonalisation
	! ===============
	allocate(zev(nsteof))
	if(zbLargeMatrix)then
		call la_syevd(cov,zev,jobz='V')
	else
		call la_syev(cov,zev,jobz='V')
	end if



	! Get ST-EOFs and eigenvalues
	! ===========================
	if(present(steof).or.present(stpc))then
		allocate(zsteof(nsteof, nkeep))
		zsteof = cov(:, nsteof : nsteof-nkeep+1 : -1)
		deallocate(cov)
		if(present(steof))then
			steof = zsteof
			deallocate(zsteof)
		end if
	end if

	! Eigen values
	! ------------
	if(present(ev))then
		ev = zev(nsteof : nsteof-nkeep+1 : -1)
	end if
	deallocate(zev)


	! Get ST-PCs
	! ==========
	if(present(stpc))then
		allocate(wpc(nt-nwindow+1))
		stpc = 0.
		do im = 1, nkeep
			do iw = 1, nwindow
				call sgemm('T','N', nt-nwindow+1, 1, nchan, 1.,&
					& zff(:,iw:iw+nt-nwindow), nchan, &
					& steof(iw:iw+(nchan-1)*nwindow:nwindow, im),nchan,&
					& 0., wpc, nt-nwindow+1)
					stpc(:, im)  =  stpc(:, im) + wpc
			end do
			stpc(:, im) = stpc(:, im) / sum(steof(:,im)**2)
		end do
	end if

  end subroutine sl_mssa



	subroutine sl_mssarec(steof, stpc, nwindow, ffrec, istart, iend)

	! Title:
	!	Reconstruction of a set of MSSA components
	!
	! Description:
	!	Same as for the reconstruction of PCA components, but for MSSA.
	!
	! Necessary arguments:
	!	- steof:   SpaceXwindow-mode array of EOFs
	!	- stpc:    Time-mode array of PCs
	!	- nwindow: Window size
	!	- ffrec:   Space-time array of the reconstructed field
	!
	! Optional arguments:
	!	- istart: Index of the first component to use
	!	- iend:   Index of the last component to use

	implicit none


	! Declarations
	! ============

	! External
	! --------
	real,   intent(in)           :: steof(:,:), stpc(:,:)
	real,   intent(out)          :: ffrec(size(steof, 1)/nwindow,&
	 &                                    size(stpc, 1)+nwindow-1)
	integer,intent(in)           :: nwindow
	integer,intent(in), optional :: istart, iend

	! Internal
	! --------
	integer :: ntpc, nchan, nt, ic, im, iw, nkept, &
	 &         itmp, zistart, ziend
	real, allocatable :: reof(:), epc(:,:)


	! Setup
	! =====

	! Sizes
	! -----
	ntpc  = size(stpc, 1)
	nt    = ntpc+nwindow-1
	nchan = size(steof, 1)/nwindow
	nkept = size(steof, 2)
	allocate(reof(nwindow))
	allocate(epc(nwindow, ntpc-nwindow+1))
	ffrec = 0.

	! Range
	! -----
	if(present(istart))zistart=istart
	if(present(iend))then
		ziend=iend
	else
		ziend=nkept
	end if
	zistart = 1
	if(zistart.lt.1.or.zistart.gt.nkept)then
		zistart = 1
		print*,'[mssarec] istart lower than 1 => set to 1'
	end if
	if(ziend.lt.1.or.ziend.gt.nkept)then
		ziend = nkept
		print*,'[mssarec] iend greater than the number of '// &
		 &     'avalaible modes => reduced to',iend
	end if
	if(zistart>ziend)then
		itmp    = ziend
		ziend   = zistart
		zistart = itmp
		print*,'[mssarec] istart > iend => inversion'
	end if

	! Computation
	! ===========
	ffrec = 0.
	do im = zistart, ziend ! sum over the selection of modes

		! (ntpc-nwindow+1) length slices
		do iw = 1, nwindow
			epc(iw,:) = stpc(iw : iw+ntpc-nwindow, im)
		end do

		do ic = 1, nchan ! sum over the channels (= space or PCs from simple PCA)

			! reversed eof
			reof = steof(nwindow+(ic-1)*nwindow : 1+(ic-1)*nwindow : -1, im)

			! * middle * [nwindow length projections]
			ffrec(ic, nwindow : ntpc) =  ffrec(ic, nwindow : ntpc) + &
				& matmul(reof, epc) / real(nwindow)

		  do iw = 1, nwindow-1

			 ! * beginning * [iw length projections]
			 ffrec(ic, iw) = ffrec(ic, iw) + &
				  & dot_product(reof(nwindow-iw+1:nwindow), &
				  &	stpc(1:iw, im)           ) / real(iw)
!
			 ! * end * [iw length projections]
			 ffrec(ic, nt-iw+1) = ffrec(ic, nt-iw+1) + &
				  & dot_product(reof(1:iw), &
				  &	stpc(ntpc-iw+1:ntpc, im) ) / real(iw)

			end do

		end do

	end do

	end subroutine sl_mssarec



  !############################################################
  !############################################################
  !## SVD PART ################################################
  !############################################################
  !############################################################

	subroutine sl_svd(ll,rr,nkeep,leof,reof,lpc,rpc,ev,bLargeMatrix)

	! Title:
	!	Singular Value Decomposition
	!
	! Description:
	!	Singular value decomposition between two datasets having
	!	the same length in time.
	!
	! Necessary arguments:
	!	- ll:    Left space-time array
	!	- rr:    Right space-time array
	!	- nkeep: Maximum number of modes to keep in outputs
	!
	! Optional arguments:
	!	- leof:  Left EOFs
	!	- reof:  Right EOFs
	!	- lpc:   Left PCs
	!	- rpc:   Right PCs
	!	- bLargeMatrix: Use la_sgesdd instead of la_sgesvd (faster for large matrices, but uses more workspace) [default:.true.]


	! Declarations
	! ============

	use la_precision, only: wp => sp
	use f95_lapack, only: la_gesdd, la_gesvd
	implicit none

	! External
	! --------
	real,    intent(in)           :: ll(:,:),rr(:,:)
	integer, intent(in)	         :: nkeep
	real,    intent(out),optional :: lpc(size(ll,2),nkeep), &
	&                                leof(size(ll,1),nkeep), &
	&                                rpc(size(ll,2),nkeep), &
	&                                reof(size(rr,1),nkeep), &
	&                                ev(nkeep)
	logical, intent(in), optional :: bLargeMatrix

	! Internal
	! --------
	integer           :: ii,ij,nn,ns,nsl,nsr,nt
	real, allocatable :: zll(:,:), zrr(:,:), cov(:,:)
	real, allocatable :: zev(:), zleof(:,:)
	integer           :: znkeepmax, i,j
	logical           :: zbLargeMatrix


	! Sizes
	! -----
	nsl = size(ll,1)
	nsr = size(rr,1)
	nt = size(ll,2)
	if(nsl/=nsr.or.nt/=size(rr,2))then
		print*,'[svd] Left and right array have incompatible sizes'
		return
	end if
	ns = min(nsr,nsl)
	znkeepmax = 100
	if(nkeep>znkeepmax)then
		print*,'[svd] You want to keep a number of PCs '//&
		 & 'greater than ',znkeepmax
		return
	end if
	if(nkeep>ns)then
		print*,'[svd] You want to keep a number of PCs '//&
			&'greater than the number of EOF:',ns
		return
	end if

	! What does the user want?
	! ------------------------
	if(.not.present(leof).and..not.present(lpc).and.&
	  &.not.present(reof).and..not.present(rpc).and.&
	  &.not.present(ev))then
		print*,'[svd] Nothing to do. Quit.'
		return
	end if

	! Use ssyevd?
	! -----------
	if(.not.present(bLargeMatrix))then
		zbLargeMatrix = .true.
	else
		zbLargeMatrix = bLargeMatrix
	end if


	! Computations
	! ============

	! Remove the mean
	! ---------------
	allocate(zll(nsl,nt))
	zll = ll - spread(sum(ll,dim=2)/real(nt), ncopies=nt, dim=2)
	allocate(zrr(nsr,nt))
	zrr = rr - spread(sum(rr,dim=2)/real(nt), ncopies=nt, dim=2)

	! Cross-covariances
	! -----------------
	allocate(cov(nsl,nsr))
	call sgemm('N','T',nsl,nsr,nt,1.,zll,nsl, &
				& zrr,nsr,0.,cov,nsl)
	cov = cov / float(nt)
	if(.not.present(leof).and..not.present(lpc)) deallocate(zll)
	if(.not.present(reof).and..not.present(rpc)) deallocate(zrr)

	! SVD
	! ---
	allocate(zleof(nsr,ns), zev(ns))
	if(bLargeMatrix)then
		call la_gesdd(cov, zev, u=zleof, job='V')
	else
		call la_gesvd(cov, zev, u=zleof, job='V')
	end if
	if(.not.present(reof).and..not.present(rpc)) deallocate(cov)


	! Get output arrays
	! =================

	! Eigen values
	! ------------
	if(present(ev)) ev = zev(1:nkeep)
	print*,'svd ev:',zev(1:nkeep)
	deallocate(zev)

	! EOFs
	! ----
	if(present(leof).or.present(lpc))then
		leof = zleof(:,1:nkeep)
		deallocate(zleof)
	end if
	if(present(reof).or.present(rpc))then
		do i = 1, nkeep
			reof(:, i) = cov(i,:)
		end do
		deallocate(cov)
	end if

	! PCs
	! ---
	if(present(lpc))then
		lpc = matmul(transpose(zll), leof)
		do i = 1, nkeep
			lpc(:,i) = lpc(:,i) / sum(leof(:,i)**2)
		end do
	end if
	if(present(rpc))then
		rpc = matmul(transpose(zrr), reof)
		do i = 1, nkeep
			rpc(:,i) = rpc(:,i) / sum(reof(:,i)**2)
		end do
	end if

	end subroutine sl_svd




  ! ############################################################
  ! ############################################################
  ! ## TOOLS PART ##############################################
  ! ############################################################
  ! ############################################################

	subroutine sl_phasecomp(ffrec,np,phases,weights,offset,firstphase)

	! Title:
	!	Phase composites
	!
	! Description:
	!	Performs phase composites of S-T oscillatory field.
	!	This field is typically a reconstructed pair of MSSA modes.
	!	Composites are evaluated according to an index defined by the
	!	first PC of the input field and its derivative.
	!	Space weights can be optionally used to compute the PC.
	!	A minimal normalized amplitude can be also used: when the
	!	index is under value, data are not used to compute phases.
	!	It is also possible so specify the angle of the first phase
	!	in the 360 degrees phase diagram circle: zero means the
	!	the first phase conincides with the maximmum.
	!
	!
	! Necessary arguments:
	!	- ffrec: Space-time array
	!	- np:    Number of requested phases over the 360 degrees cycle [default:8]
	!
	! Optional arguments:
	!	- weights:    Space array of weights
	!	- offset:     Minimal normalized amplitude of the index [default:0.]
	!	- firstphase: Value in degrees of the first phase [default:0.]
	!
	! Dependencies:
	!	sl_pca


	implicit none

	! Declarations
	! ============

	! External
	! --------
	integer, intent(in)           :: np
	real,    intent(in)           :: ffrec(:,:)
	real,    intent(in), optional :: weights(:)
	real,    intent(in), optional :: offset, firstphase
	real,    intent(out)          :: phases(size(ffrec, 1),np)

	! Internal
	! --------
	real, allocatable :: xeof(:,:), pc(:,:)
	real :: dpc(size(ffrec,2)), amp(size(ffrec,2))
	integer :: nt, iphase
	real :: angles(np), projection(size(ffrec,2))
	real :: pi, deltarad, pcos, psin, zoffset, zfirstphase
	logical :: select_amplitude(size(ffrec,2)), &
	 &         select_phase(size(ffrec,2))
	integer :: itime(size(ffrec,2)), nsel, i, ns
	integer, allocatable :: isel(:)


	! Setup
	! =====
	nt = size(ffrec,2)
	pi = acos(-1.)
	itime = (/ (i, i=1, nt) /)
	ns = size(ffrec, 1)
	if(present(offset))then
		zoffset=offset
	else
		zoffset=0.
	end if

	! Find the first PC and its derivative
	! ====================================
	allocate(pc(nt,1))
	call sl_pca(ffrec, 1, pc=pc, weights=weights)
	pc = pc * sqrt(real(nt)/sum(pc**2))
	dpc = 0.5 * (eoshift(pc(:,1),  1, pc(nt,1)) - &
	 &           eoshift(pc(:,1), -1, pc(1,1)))
	dpc((/1,nt/)) = dpc((/1,nt/)) * 2.
	dpc = dpc * sqrt(real(nt)/sum(dpc**2))
	amp = sqrt(pc(:,1)**2 + dpc**2)


	! Compute the maps
	! ================

	! Define the marks
	! ----------------
	deltarad = 2 * pi / real(np)
	if(present(firstphase))then
		zfirstphase = modulo(firstphase * 2 * pi / 360., 2 * pi)
	else
	   zfirstphase = 0.
	end if
	angles = (/ (real(iphase), iphase=0,np-1) /) * deltarad + &
	 &       zfirstphase

	! Compute the phase maps
	! ----------------------
	phases = 0.
	select_amplitude = amp >= zoffset
	do iphase = 1, np
		pcos = cos(angles(iphase))
		psin = sin(angles(iphase))
		projection =  (pc(:,1)*pcos+dpc*psin) / amp
		select_phase = ( projection >= cos(0.5*deltarad) ) &
		 &             .and. select_amplitude
		if(any(select_phase))then
			nsel = count(select_phase)
			allocate(isel(nsel))
			isel = pack(itime, select_phase)
			phases(:,iphase) = sum(ffrec(:,isel), dim=2) / real(nsel)
			deallocate(isel)
		end if
	end do

	end subroutine sl_phasecomp


end module spanlib
