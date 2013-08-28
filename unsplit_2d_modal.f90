! ===================================== 
! Unsplit 2D modal algorithm for tracer transport
! Uses tensor product Legendre polynomial basis to simulate 2D tracer transport equations with variable windspeeds
!
! Dependencies:
! 	netCDF
!	LAPACK
!	mDGmod.f90 ; mDG_elem_update.f90
! By: Devin Light ; Aug. 2013
! =====================================

PROGRAM execute
	IMPLICIT NONE
	USE mDGmod
!	USE netcdf

	INTEGER :: ntest,start_res

	write(*,*) '======'
	write(*,*) 'TEST 1: Uniform advection (u=v=1)'
	write(*,*) '======'

	start_res = 8 ! Number of elements in each direction
	call test2d_modal(1,start_res,start_res,2,3,20,0.2D0)

	CONTAINS

	SUBROUTINE test2d_modal(ntest,nex0,ney0,nscale,nlevel,noutput,maxcfl)
		IMPLICIT NONE
		! Inputs
		INTEGER, INTENT(IN) :: ntest,nex0,ney0,nscale,nlevel,noutput
		REAL(KIND=8), INTENT(IN) :: maxcfl

		! Local variables
		INTEGER, DIMENSION(10) :: tmp_method
		INTEGER :: nmethod, nmethod_final, imethod,ierr
		INTEGER :: dgorder,p

		CHARACTER(len=40) :: cdf_out

		INTEGER :: nex,ney,nxplot,nyplot
		REAL(KIND=8) :: dxel,dyel,dxplot,dyplot
		REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: q, L, dL
		REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: qnodes,qweights, x_elcent, y_elcent, xplot,yplot
		REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: A, &   ! Coefficent array
													   u,v, & ! Velocites within each element
													  uedge,vedge ! Velocities along each edge of the element
		INTEGER :: i,j,k

		CHARACTER(len=9) :: outdir
		REAL(KIND=4), DIMENSION(2) :: tstart,tend
		REAL(KIND=4) :: t0,tf

		if(nlevel.lt.1) STOP 'nlev should be at least 1 in test2d_modal'

		nmethod_final = 1
		tmp_method = 0
		tmp_method(1) = 1 ! 2D Modal, no limiting

		DO nmethod=1,nmethod_final
			imethod = tmp_method(nmethod)

			SELECT CASE(imethod)
				CASE(1)
				  WRITE(*,*) '2D Modal, Unsplit, No limiting'
				  outdir = 'mdgunlim/'
				  dgorder = 3
				  WRITE(*,*) 'Max order of Legendre basis=',dgorder
			END SELECT

			! Initialize quadrature weights and nodes (only done once per call)
			ALLOCATE(qnodes(0:dgorder),qweights(0:dgorder), L(0:dgorder,0:dgorder), dL(0:dgorder,0:dgorder), STAT=ierr)
	        CALL quad_nodes(dgorder+1,qnodes)
			CALL quad_weights(dgorder+1,qnodes,qweights)

			DO i=0,dgorder
				DO j=0,dgorder
					L(i,j) = legendre(qnodes(i),j)
					dL(i,j) = dlegendre(qnodes(i),j)
				ENDDO
			ENDDO

			DO p=1,nlevel
				nex = nex0*nscale**(p-1)
				ney = ney0*nscale**(p-1)
				dxel = 1D0/DBLE(nex)
				dyel = 1D0/DBLE(ney)

				nxplot = (dgorder+1)*nex
				nyplot = (dgorder+1)*ney
				dxplot = 1D0/DBLE(nxplot)
				dyplot = 1D0/DBLE(nyplot)


				ALLOCATE(q(1:nxplot,1:nyplot),x_elcent(1:nex),y_elcent(1:ney), xplot(1:nxplot), yplot(1:nyplot), &
						 A(1:nex,1:ney,0:dgorder,0:dgorder), u(1:nex,1:ney,0:dgorder,0:dgorder), v(1:nex,1:ney,0:dgorder,0:dgorder),&
						 uedge(1:nex,1:ney,1:2,0:dgorder), vedge(1:nex,1:ney,1:2,0:dgorder), STAT=ierr)

				! Note that elements are ordered row-wise within the domain

				! q(i,j) : Approx. solution on plotting grid at point (xi,yj) based on series expansion
				! A(i,j,k) : Coefficent array ; i = element, (j,k) = which Legendre polynomial (A(i,j,k)*P_j(xi)*P_k(eta))
				! u(i,j,k), v(i,j,k) : Horizontal/vertical vel. array ; i = element, (j,k) = horiz,vertical location within element
				! uedge(i,j,k), vedge(i,j,k) : Edge horiz./vert. vel. array ; i = element, j = which edge (lower index closer to origin), k = which node
				
				! Initialize x- and y- grids
				x_elcent(1) = dxel/2D0
				DO i=2,nex
					x_elcent(i) = x_elcent(i-1)+dxel
				ENDDO

				y_elcent(1) = dyel/2D0
				DO i=2,ney
					y_elcent(i) = y_elcent(i-1)+dyel
				ENDDO

				xplot(1) = dxplot/2D0
				DO i=2,nxplot
					xplot(i) = xplot(i-1)+dxplot
				ENDDO

				yplot(1) = dyplot/2D0
				DO i=2,nyplot
					yplot(i) = yplot(i-1)+dyplot
				ENDDO

				! Initialize q, A, u, and v
				CALL init2d(ntest,nex,ney,dgorder,A,u,v,uedge,vedge,x_elcent,y_elcent,qnodes,qweights,L)

				DEALLOCATE(q,x_elcent,y_elcent,xplot,yplot,A,u,v,uedge,vedge, STAT=ierr)
			ENDDO
		ENDDO

	END SUBROUTINE test2d_modal

	SUBROUTINE init2d(ntest,nex,ney,dgorder,A,u,v,uedge,vedge,x_elcent,y_elcent,qnodes,qweights,L)
		IMPLICIT NONE
		! Inputs
		INTEGER, INTENT(IN) :: ntest,nex,ney,dgorder
		REAL(KIND=8), DIMENSION(1:nex), INTENT(IN) :: x_elcent
		REAL(KIND=8), DIMENSION(1:ney), INTENT(IN) :: y_elcent
		REAL(KIND=8), DIMENSION(0:dgorder), INTENT(IN) :: qnodes,qweights
		REAL(KIND=8), DIMENSION(0:dgorder,0:dgorder), INTENT(IN) :: L

		! Ouputs
		REAL(KIND=8), DIMENSION(1:nex,1:ney,0:dgorder,0:dgorder), INTENT(OUT) ::  A
		REAL(KIND=8), DIMENSION(1:(dgorder+1)*nex,1:(dgorder+1)*ney) , INTENT(OUT):: u,v
		REAL(KIND=8), DIMENSION(0:nex,1:(dg+1)*ney), INTENT(OUT) :: uedge
		REAL(KIND=8), DIMENSION(1:(dg+1)*nex,0:ney), INTENT(OUT) :: vedge
		CHARACTER(len=40), INTENT(OUT) :: cdf_out

		! Local variables
		REAL(KIND=8) :: dxmin,dymin,dxel,dyel
		REAL(KIND=8), DIMENSION(1:(dgorder+1)*nex) :: x
		REAL(KIND=8), DIMENSION(1:(dgorder+1)*ney) :: y

		REAL(KIND=8), DIMENSION(1:(dgorder+1)*nex,1:(dgorder+1)*ney,2) :: psiu
		REAL(KIND=8), DIMENSION(1:(dgorder+1)*nex,1:(dgorder+1)*ney,2) :: psiv
		REAL(KIND=8), DIMENSION(0:nex,1:(dg+1)*ney,2) :: psiu_edge
		REAL(KIND=8), DIMENSION(1:(dg+1)*nex,1:0:ney,2) :: psiv_edge

		REAL(KIND=8), DIMENSION(0:dgorder,0:dgorder) :: foo

		REAL(KIND=8) :: PI
		INTEGER :: i,j,l,m,s,t

		PI = DACOS(-1D0)
		dxel = x_elcent(2)-x_elcent(1)
		dyel = y_elcent(2)-y_elcent(1)

		! Minimum internode spacing, mapped to physical domain
		dxmin = (dxel/2D0)*MINVAL(qnodes(1:dgorder)-qnodes(0:dgorder-1))
		dymin = (dyel/2D0)*MINVAL(qnodes(1:dgorder)-qnodes(0:dgorder-1))

		! Compute x- and y- values where we want velocities within element
		DO i=1,nex
			x(1+(i-1)*(dgorder+1):i*(dgorder+1)) = (dxel/2D0)*qnodes(0:dgorder)+x_elcent(i)
		ENDDO

		DO j=1,ney
			y(1+(i-1)*(dgorder+1):i*(dgorder+1)) = (dyel/2D0)*qnodes(0:dgorder)+y_elcent(i)
		ENDDO

		! Fill stream function arrays
		SELECT CASE(ntest)
			CASE(1) ! Uniform u=v=1 ; no time dependence

				DO j=1,(dgorder+1)*ney
					psiu(:,j,1) = (y(j)-dymin/2D0) - x(:)
					psiu(:,j,2) = (y(j)+dymin/2D0) - x(:)
			
					psiu_edge(0,j,1) = (y(j)-dymin/2D0) - (x_elcent(:) - dxel/2D0) ! Vel thru left face 1st element col
					psiu_edge(0,j,2) = (y(j)+dymin/2D0) - (x_elcent(:) - dxel/2D0)
					psiu_edge(1:nex,j,1) = (y(j)-dymin/2D0) - (x_elcent(:) + dxel/2D0) ! Vel thru right face all others
					psiu_edge(1:nex,j,2) = (y(j)+dymin/2D0) - (x_elcent(:) + dxel/2D0)
				ENDDO


				DO i=1,(dgorder+1)*nex
					psiv(i,:,1) = y(:) - (x(i)-dxmin/2D0)
					psiv(i,:,2) = y(:) - (x(i)+dxmin/2D0)

					psiv_edge(i,0,1) = (y_elcent(:) - dyel/2D0) - (x(i)-dxmin/2D0) ! Vel thru bot face 1st element row
					psiv_edge(i,0,2) = (y_elcent(:) - dyel/2D0) - (x(i)+dxmin/2D0)
					psiv_edge(i,1:ney,1) = (y_elcent(:) + dyel/2D0) - (x(i)-dxmin/2D0) ! Vel thru top face all others
					psiv_edge(i,1:ney,2) = (y_elcent(:) + dyel/2D0) - (x(i)+dxmin/2D0)					
				ENDDO

		ENDSELECT

		! Compute velocity from stream functions
		DO j=1,(dgorder+1)*ney
			u(:,j) = (psiu(:,j,2)-psiu(:,j,1))/dymin
		ENDDO
		DO i=1,(dgorder+1)*nex
			v(i,:) = -(psiv(i,:,2)-psiv(i,:,1))/dxmin
		ENDDO
		uedge(:,:) = (psiu_edge(:,:,2)-psiu_edge(:,:,1))/dymin
		vedge(:,:) = -(psiv_edge(:,:,2)-psiv_edge(:,:,1))/dxmin

		! Initialize A array
		SELECT CASE(ntest)
			CASE(1) ! sine wave advection
				cdf_out = 'dg2d_sine_adv.nc'
!          q(:,j) = sin(2.d0*pi*x)*sin(2.d0*pi*y(j))
				DO i=1:nex
					DO j=1:ney
						DO l=0:dgorder
						DO m=0:dgorder
							DO s=0:dgorder
							DO t=0:dgorder
		FOO(s,t) = qweights(s)*qweights(t)*DSIN(2D0*PI*x(1+s+(i-1)*(dgorder+1)))*DSIN(2D0*PI*y(1+t+(j-1)*(dgorder+1)))* &
					L(s,l)*L(t,m)
							ENDDO
							ENDDO
		A(i,j,l,m) = ((2D0*l+1)*(2D0*m+1)/4D0)*SUM(FOO)
						ENDDO
						ENDDO
					ENDDO
				ENDDO
		ENDSELECT

	END SUBROUTINE init2d
END PROGRAM execute
