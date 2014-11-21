

subroutine optimec_pca(data, eof, ec, mv, ns, nt, nm)

    implicit none
    integer, intent(in) :: ns, nt, nm
    real(kind=8), intent(in) :: eof(ns,nm), data(ns,nt), mv
    real(kind=8), intent(inout) :: ec(nt,nm)

    integer :: m, IPRINT(2), IFLAG,icall, it
    parameter(m=7)
    real(kind=8) :: diff(ns), F, G(nm), W(nm*(2*m+1)+2*m), mvtol, EPS, XTOL, norm
    logical :: diagco, beof(ns), bec(nt), bdata(ns,nt)

    integer :: MP,LP
    real(kind=8) :: GTOL,STPMIN,STPMAX,DIAG(nm)

      EXTERNAL LB2
      COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX

    mvtol = epsilon(1d0)
    bdata = abs((data-mv)/mv)>mvtol
    beof = abs((eof(:,1)-mv)/mv)>mvtol
    bec = abs((ec(:,1)-mv)/mv)>mvtol

    norm = sum(data**2, mask=bdata)/dble(size(data,1))

    DIAGCO= .FALSE.
    EPS= 1.0D-5*norm
    XTOL= 1.0D-16
    IPRINT(1)= -1
    IPRINT(2)= 0
    IFLAG=0

    !$OMP PARALLEL &
    !$OMP SHARED(eof,ec,data,bdata,beof,bec,diagco,eps,xtol,nt,nm) &
    !$OMP PRIVATE(it,icall,diff,F,G,W,IFLAG)
    !$OMP DO

    do it=1,nt

        if(bec(it))then

            IFLAG = 0
            do icall=1,2000

                diff = matmul(eof, ec(it,:))
                diff = merge(diff-data(:,it), 0d0, bdata(:,it).and.beof)

                F = sum(diff**2)
                G = 2*matmul(transpose(eof), diff)

                call LBFGS(nm,M,ec(it,:),F,G,DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IFLAG)

                if(IFLAG==0)exit

            end do

        end if

    end do
    !$OMP END DO
    !$OMP END PARALLEL

end subroutine optimec_pca
