! ===================================== 
! Unsplit 2D modal algorithm for tracer transport
! Uses tensor product Legendre polynomial basis to simulate 2D tracer transport equations with variable windspeeds
!
! Dependencies:
! 	netCDF
!	LAPACK
!	mDGmod.f90 ; coeff_update.f90
! By: Devin Light ; Aug. 2013
! =====================================

PROGRAM execute
	USE mDGmod
	USE netcdf

	IMPLICIT NONE

	INTEGER :: ntest,start_res
	LOGICAL :: transient

	write(*,*) '======'
	write(*,*) 'TEST 1: Uniform advection (u=v=1)'
	write(*,*) '======'

	transient = .FALSE.
!	transient = .TRUE.
	start_res = 20 ! Number of elements in each direction
!	CALL test2d_modal(1,start_res,start_res,2,3,20,0.1D0) !1D0/(2D0*4D0-1D0) !0.3D0/sqrt(2d0)

	write(*,*) '======'
	write(*,*) 'TEST 2: Smooth cosbell deformation'
	write(*,*) '======'

	transient = .TRUE.
!	transient = .FALSE.
	CALL test2d_modal(6,start_res,start_res,2,3,20,0.1D0) !1D0/(2D0*4D0-1D0)

	write(*,*) '======'
	write(*,*) 'TEST 3: Standard cosbell deformation'
	write(*,*) '======'

	transient = .TRUE.
!	CALL test2d_modal(5,start_res,start_res,2,3,20,1D0/(2D0*4D0-1D0))

	write(*,*) '======'
	write(*,*) 'TEST 4: Square wave deformation'
	write(*,*) '======'
	
	transient = .TRUE.
!	CALL test2d_modal(7,start_res,start_res,2,3,20,0.01D0) !1D0/(2D0*4D0-1D0)

CONTAINS

	SUBROUTINE test2d_modal(ntest,nex0,ney0,nscale,nlevel,noutput,maxcfl)
		IMPLICIT NONE
		! Inputs
		INTEGER, INTENT(IN) :: ntest,nex0,ney0,nscale,nlevel,noutput
		REAL(KIND=8), INTENT(IN) :: maxcfl

		! Local variables
		INTEGER, DIMENSION(10) :: tmp_method
	    REAL(KIND=8), DIMENSION(nlevel) :: e1, e2, ei
		REAL(KIND=8) :: cnvg1, cnvg2, cnvgi
		INTEGER :: nmethod, nmethod_final, imethod,ierr,nstep, nout
		INTEGER :: dgorder, norder,p

		CHARACTER(len=40) :: cdf_out

		INTEGER :: nex,ney,nxplot,nyplot,nxiplot,netaplot
		REAL(KIND=8) :: dxel,dyel,dxplot,dyplot,tfinal, tmp_umax, tmp_vmax, dxm, dym,dt, time, dxiplot,detaplot
		REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: q, q0, Leg, dLeg, L_xi_plot, L_eta_plot, &
													 u,v,u_tmp,v_tmp, & ! Velocites within each element
													 uedge,vedge,uedge_tmp, vedge_tmp, & ! Velocities along each edge of the element
													 foo,C0,C
		REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: qnodes, qweights, x_elcent, y_elcent, xplot,yplot,xiplot,etaplot
		REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: A,&  ! Coefficent array
														LL,dLL
		INTEGER :: i,j,l,m,k,s,t,n

		CHARACTER(len=9) :: outdir
		REAL(KIND=4), DIMENSION(2) :: tstart,tend,t1,t2
		REAL(KIND=4) :: t0,tf,tick,tock

		REAL(KIND=8) :: tmp_qmax,tmp_qmin

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
				  norder = 3
				  dgorder = 2*(norder+1)
				  WRITE(*,*) 'N=',norder,'Uses a total of',(norder+1)**2,'Legendre basis polynomials'
			END SELECT

			! Initialize quadrature weights and nodes (only done once per call)
			ALLOCATE(qnodes(0:dgorder),qweights(0:dgorder), Leg(0:dgorder,0:norder), dLeg(0:dgorder,0:norder), &
					LL(0:norder,0:norder,0:dgorder,0:dgorder), dLL(0:norder,0:norder,0:dgorder,0:dgorder), &
					STAT=ierr)
!			ALLOCATE(foo(0:dgorder,0:dgorder), STAT=ierr)
	        CALL quad_nodes(dgorder+1,qnodes)
			CALL quad_weights(dgorder+1,qnodes,qweights)

			nxiplot = norder+1
			netaplot = norder+1
			dxiplot = 2D0/DBLE(nxiplot) ! 2D0 because each element gets mapped to [-1,1]
			detaplot = 2D0/DBLE(netaplot)

			ALLOCATE(xiplot(1:nxiplot), etaplot(1:netaplot), L_xi_plot(1:nxiplot,0:norder),L_eta_plot(1:netaplot,0:norder), STAT=ierr)

			xiplot(1) = dxiplot/2D0 - 1D0
			DO i=2,nxiplot
				xiplot(i) = xiplot(i-1)+dxiplot
			ENDDO

			etaplot(1) = detaplot/2D0 - 1D0
			DO i=2,netaplot
				etaplot(i) = etaplot(i-1)+detaplot
			ENDDO

!			xiplot = qnodes(:)
!			etaplot = qnodes(:)

			! Precompute Legendre basis evaluated at quadrature points and plotting points
			DO i=1,nxiplot
				DO j=0,norder
					L_xi_plot(i,j) = legendre(xiplot(i),j)
				ENDDO
			ENDDO
			
			DO i=1,netaplot
				DO j=0,norder
					L_eta_plot(i,j) = legendre(etaplot(i),j)
				ENDDO
			ENDDO

			DO i=0,dgorder
				DO j=0,norder
					Leg(i,j) = legendre(qnodes(i),j)
					dLeg(i,j) = dlegendre(qnodes(i),j)
				ENDDO
			ENDDO

			DO l=0,norder
				DO m=0,norder
					DO s=0,dgorder
						DO t=0,dgorder
							LL(l,m,s,t) = Leg(s,l)*Leg(t,m)
							dLL(l,m,s,t) = dLeg(s,l)*Leg(t,m)
						ENDDO
					ENDDO
				ENDDO
			ENDDO


			DO p=1,nlevel
				
				t0 = etime(tstart)

				nex = nex0*nscale**(p-1)
				ney = ney0*nscale**(p-1)
				dxel = 1D0/DBLE(nex)
				dyel = 1D0/DBLE(ney)

				nxplot = (nxiplot)*nex
				nyplot = (netaplot)*ney
				dxplot = 1D0/DBLE(nxplot)
				dyplot = 1D0/DBLE(nyplot)


				ALLOCATE(q(1:nxplot,1:nyplot),q0(1:nxplot,1:nyplot),x_elcent(1:nex),y_elcent(1:ney), xplot(1:nxplot), yplot(1:nyplot), &
						 A(1:nex,1:ney,0:norder,0:norder), u(1:(dgorder+1)*nex,1:(dgorder+1)*ney), & 
						 u_tmp(1:(dgorder+1)*nex,1:(dgorder+1)*ney), v(1:(dgorder+1)*nex,1:(dgorder+1)*ney),	&
						 v_tmp(1:(dgorder+1)*nex,1:(dgorder+1)*ney), uedge(1:nex,1:(dgorder+1)*ney), &
						 uedge_tmp(1:nex,1:(dgorder+1)*ney), vedge(1:(dgorder+1)*nex,1:ney), vedge_tmp(1:(dgorder+1)*nex,1:ney),&
						 C0(1:nex,1:ney), C(1:nex,1:ney), STAT=ierr)

				! Note that elements are ordered row-wise within the domain

				! q(i,j) : Approx. solution on output grid at point (x_i,y_j) based on element expansion
				! A(i,j,l,m) : Coefficent array ; (i,j) = element, (l,m) = which Legendre polynomial (q_ij = SUM(SUM(A(i,j,l,m)*P_l(xi)*P_m(eta))))
				! u(i,j), v(i,j) : Horizontal/vertical vel. array ; (i,j) = horiz,vertical location
				! uedge(i,k), vedge(i,j,k) : Edge horiz./vert. vel. array ; i = element, k = which node
				
				! Initialize x- and y- grids and xi- and eta- plotting grids
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

				CALL init2d(ntest,nex,ney,nxplot,nyplot,dgorder,norder,A,q0,u,v,uedge,vedge,x_elcent,y_elcent,xplot,yplot, &
							qnodes,qweights,Leg,cdf_out,tfinal)

				! Store element averages for conservation estimation
				DO i=1,nex
					DO j=1,ney
						C0(i,j) = A(i,j,0,0) 
					ENDDO
				ENDDO

				u_tmp = u
				v_tmp = v
				uedge_tmp = uedge
				vedge_tmp = vedge
!				u_tmp = 1D0
!				v_tmp = 1D0
!				uedge_tmp = 1D0
!				vedge_tmp = 1D0

				cdf_out = outdir // cdf_out

				! Set up timestep
!				dxm = dxel/DBLE(dgorder+1)
!				dym = dyel/DBLE(dgorder+1)
				dxm = dxel
				dym = dyel

				tmp_umax = MAX(MAXVAL(DABS(u)),MAXVAL(DABS(uedge)))
				tmp_vmax = MAX(MAXVAL(DABS(v)),MAXVAL(DABS(vedge)))

				IF(noutput .eq. -1) THEN
!					nstep = CEILING( (tfinal/maxcfl)*(tmp_umax/dxm + tmp_vmax/dym) )
					nstep = CEILING( (tfinal/maxcfl)*(tmp_umax + tmp_vmax)/(dxm) )
					nout = nstep
				ELSE
!					nstep = noutput*CEILING( (tfinal/maxcfl)*(tmp_umax/dxm + tmp_vmax/dym)/DBLE(noutput) )
					nstep = noutput*CEILING( (tfinal/maxcfl)*(tmp_umax + tmp_vmax)/(dxm)/DBLE(noutput) )
					nout = noutput
				ENDIF

				dt = tfinal/DBLE(nstep)

				IF(p .eq. 1) THEN ! Set up netCDF file
					CALL output2d(q0,xplot,yplot,nxplot,nyplot,tfinal,cdf_out,nout,-1)
				ENDIF

				CALL output2d(q0,xplot,yplot,nxplot,nyplot,0D0,cdf_out,p,0) ! Set up variables for this value of p ; Write x, y, and initial conditions

				! Time integration
				tmp_qmax = MAXVAL(q0)
				tmp_qmin = MINVAL(q0)

				time = 0D0
				DO n=1,nstep

					CALL coeff_update(q,A,u_tmp,v_tmp,uedge_tmp,vedge_tmp,qnodes,qweights,Leg,dLL,LL,L_xi_plot,L_eta_plot,dxel,dyel,& 
									  dt,dgorder,norder,nxplot,nyplot,nex,ney,nxiplot,netaplot,transient,time)

					! Store element averages for conservation estimation (for modal DG these are just the 0th order coeffs)
					DO i=1,nex
						DO j=1,ney
							C(i,j) = A(i,j,0,0) 
						ENDDO
					ENDDO


					time = time + dt
	
					IF((MOD(n,nstep/nout).eq.0).OR.(n.eq.nstep)) THEN
						! Write output
						write(*,*) 'Outputting at time t=',time
						CALL output2d(q,xplot,yplot,nxplot,nyplot,time,cdf_out,p,2)
					ENDIF
					
					tmp_qmax = MAX(tmp_qmax,MAXVAL(q))
					tmp_qmin = MIN(tmp_qmin,MINVAL(q))

				ENDDO

				e1(p) = SUM(ABS(q(1:nxplot,1:nyplot)-q0))/DBLE(nxplot)/DBLE(nyplot)
				e2(p) = SQRT(SUM((q(1:nxplot,1:nyplot)-q0)**2)/DBLE(nxplot)/DBLE(nyplot))
!e2(p)=0D0
!do i = 1,nex
!do j = 1,ney
!	do l=0,dgorder
!	do m=0,dgorder
!foo(l,m) = qweights(l)*qweights(m)*(q(1+l+(i-1)*nxiplot,1+m+(j-1)*netaplot)-q0(1+l+(i-1)*nxiplot,1+m+(j-1)*netaplot))**2
!	enddo
!	enddo
!e2(p) = e2(p)+sum(foo)
!enddo
!enddo
!e2(p) = sqrt(e2(p))

				ei(p) = MAXVAL(ABS(q(1:nxplot,1:nyplot)-q0))
				tf = etime(tend) - t0
				if (p.eq.1) then
					write(UNIT=6,FMT='(A117)') &
'nex    ney      E1        E2       Einf      convergence rate  overshoot  undershoot    cons   cputime   time step'
					cnvg1 = 0.d0
					cnvg2 = 0.d0
					cnvgi = 0.d0
		       else
        		  cnvg1 = -log(e1(p)/e1(p-1))/log(dble(nscale))
        		  cnvg2 = -log(e2(p)/e2(p-1))/log(dble(nscale))
        		  cnvgi = -log(ei(p)/ei(p-1))/log(dble(nscale))
		       end if
       write(*,990) nex, ney, e1(p), e2(p), ei(p), &
            cnvg1, cnvg2, cnvgi, &
            tmp_qmax-MAXVAL(q0), &
            MINVAL(q0)-tmp_qmin, &
            SUM(C-C0)/DBLE(nex*ney), tf, nstep


				IF(p .eq. nlevel) THEN
					CALL output2d(q,xplot,yplot,nxplot,nyplot,tfinal,cdf_out,p,1) ! Close netCDF files
				ENDIF
				DEALLOCATE(q,q0,x_elcent,y_elcent,xplot,yplot,A,u,v,uedge,vedge,u_tmp,v_tmp,uedge_tmp,vedge_tmp, STAT=ierr)

			ENDDO
		ENDDO
		DEALLOCATE(qnodes,qweights,Leg,dLeg,LL,dLL,xiplot,etaplot,L_xi_plot,L_eta_plot, STAT=ierr)

990    format(2i6,3e12.4,3f5.2,3e12.4,f8.2,i8)

	END SUBROUTINE test2d_modal

	SUBROUTINE init2d(ntest,nex,ney,nxplot,nyplot,dgorder,norder,A,q,u,v,uedge,vedge,x_elcent,y_elcent,xplot,yplot, &
					  qnodes,qweights,Leg,cdf_out,tfinal)
		IMPLICIT NONE
		! Inputs
		INTEGER, INTENT(IN) :: ntest,nex,ney,nxplot,nyplot,dgorder,norder
		REAL(KIND=8), DIMENSION(1:nex), INTENT(IN) :: x_elcent
		REAL(KIND=8), DIMENSION(1:ney), INTENT(IN) :: y_elcent
		REAL(KIND=8), DIMENSION(1:nxplot), INTENT(IN) :: xplot
		REAL(KIND=8), DIMENSION(1:nyplot), INTENT(IN) :: yplot
		REAL(KIND=8), DIMENSION(0:dgorder), INTENT(IN) :: qnodes,qweights
		REAL(KIND=8), DIMENSION(0:dgorder,0:norder), INTENT(IN) :: Leg

		! Ouputs
		REAL(KIND=8), INTENT(OUT) :: tfinal
		REAL(KIND=8), DIMENSION(1:nex,1:ney,0:norder,0:norder), INTENT(OUT) ::  A
		REAL(KIND=8), DIMENSION(1:nxplot,1:nyplot), INTENT(OUT) :: q
		REAL(KIND=8), DIMENSION(1:(dgorder+1)*nex,1:(dgorder+1)*ney) , INTENT(OUT):: u,v
		REAL(KIND=8), DIMENSION(1:nex,1:(dgorder+1)*ney), INTENT(OUT) :: uedge
		REAL(KIND=8), DIMENSION(1:(dgorder+1)*nex,1:ney), INTENT(OUT) :: vedge
		CHARACTER(len=40), INTENT(OUT) :: cdf_out

		! Local variables
		REAL(KIND=8) :: dxmin,dymin,dxel,dyel,r_s,ph
		REAL(KIND=8), DIMENSION(1:(dgorder+1)*nex) :: x
		REAL(KIND=8), DIMENSION(1:(dgorder+1)*ney) :: y

		REAL(KIND=8), DIMENSION(1:(dgorder+1)*nex,1:(dgorder+1)*ney,2) :: psiu
		REAL(KIND=8), DIMENSION(1:(dgorder+1)*nex,1:(dgorder+1)*ney,2) :: psiv
		REAL(KIND=8), DIMENSION(1:nex,1:(dgorder+1)*ney,2) :: psiu_edge
		REAL(KIND=8), DIMENSION(1:(dgorder+1)*nex,1:ney,2) :: psiv_edge
	    REAL(KIND=8), DIMENSION(1:nxplot,1:nyplot) :: r
		REAL(KIND=8), DIMENSION(0:dgorder,0:dgorder) :: foo,g,r_el
	
		REAL(KIND=8), DIMENSION(1:(dgorder+1)*nex) :: erru
		REAL(KIND=8), DIMENSION(1:(dgorder+1)*ney) :: errv

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
			x(1+(i-1)*(dgorder+1):i*(dgorder+1)) = dxel*qnodes(0:dgorder)/2D0+x_elcent(i)
		ENDDO

		DO j=1,ney
			y(1+(j-1)*(dgorder+1):j*(dgorder+1)) = dyel*qnodes(0:dgorder)/2D0+y_elcent(j)
		ENDDO

		! Fill stream function arrays
		SELECT CASE(ntest)
			CASE(1,10) ! Uniform u=v=1 ; no time dependence
				DO j=1,(dgorder+1)*ney
					psiu(:,j,1) = (y(j)-dymin/2D0) - x(:)
					psiu(:,j,2) = (y(j)+dymin/2D0) - x(:)
			
!					psiu_edge(0,j,1) = (y(j)-dymin/2D0) - (x_elcent(1) - dxel/2D0) ! Vel thru left face 1st element col
!					psiu_edge(0,j,2) = (y(j)+dymin/2D0) - (x_elcent(1) - dxel/2D0)
					psiu_edge(1:nex,j,1) = (y(j)-dymin/2D0) - (x_elcent(:) + dxel/2D0) ! Vel thru right face all others
					psiu_edge(1:nex,j,2) = (y(j)+dymin/2D0) - (x_elcent(:) + dxel/2D0)
				ENDDO


				DO i=1,(dgorder+1)*nex
					psiv(i,:,1) = y(:) - (x(i)-dxmin/2D0)
					psiv(i,:,2) = y(:) - (x(i)+dxmin/2D0)

!					psiv_edge(i,0,1) = (y_elcent(1) - dyel/2D0) - (x(i)-dxmin/2D0) ! Vel thru bot face 1st element row
!					psiv_edge(i,0,2) = (y_elcent(1) - dyel/2D0) - (x(i)+dxmin/2D0)
					psiv_edge(i,1:ney,1) = (y_elcent(:) + dyel/2D0) - (x(i)-dxmin/2D0) ! Vel thru top face all others
					psiv_edge(i,1:ney,2) = (y_elcent(:) + dyel/2D0) - (x(i)+dxmin/2D0)					
				ENDDO

			CASE(5:7) ! Leveque deformation flow
				!(1/pi)*sin(pi*xf(i))**2 * sin(pi*yf(j))**2
				DO j=1,(dgorder+1)*ney
					psiu(:,j,1) = (1/PI)*DSIN(PI*x(:))**2 * DSIN(PI*(y(j)-dymin/2D0))**2
					psiu(:,j,2) = (1/PI)*DSIN(PI*x(:))**2 * DSIN(PI*(y(j)+dymin/2D0))**2

					psiu_edge(1:nex,j,1) = (1/PI)*DSIN(PI*(x_elcent(:)+dxel/2D0))**2 * DSIN(PI*(y(j)-dymin/2D0))**2
					psiu_edge(1:nex,j,2) = (1/PI)*DSIN(PI*(x_elcent(:)+dxel/2D0))**2 * DSIN(PI*(y(j)+dymin/2D0))**2
				ENDDO

				DO i=1,(dgorder+1)*nex
					psiv(i,:,1) = (1/PI)*DSIN(PI*(x(i)-dxmin/2D0))**2 * DSIN(PI*y(:))**2
					psiv(i,:,2) = (1/PI)*DSIN(PI*(x(i)+dxmin/2D0))**2 * DSIN(PI*y(:))**2

					psiv_edge(i,1:ney,1) = (1/PI)*DSIN(PI*(x(i)-dxmin/2D0))**2 * DSIN(PI*(y_elcent(:)+dyel/2D0))**2
					psiv_edge(i,1:ney,2) = (1/PI)*DSIN(PI*(x(i)+dxmin/2D0))**2 * DSIN(PI*(y_elcent(:)+dyel/2D0))**2
				ENDDO
		ENDSELECT

		! Compute velocity from stream functions
		DO j=1,(dgorder+1)*ney
			u(:,j) = (psiu(:,j,2)-psiu(:,j,1))/dymin
!			u(:,j) = 2D0*DSIN(PI*x(:))**2 * DCOS(PI*y(j))
		ENDDO
		DO i=1,(dgorder+1)*nex
			v(i,:) = -(psiv(i,:,2)-psiv(i,:,1))/dxmin
!			v(i,:) = -2D0*DSIN(PI*y(:))**2 * DCOS(PI*x(i))
		ENDDO
		uedge(:,:) = (psiu_edge(:,:,2)-psiu_edge(:,:,1))/dymin
		vedge(:,:) = -(psiv_edge(:,:,2)-psiv_edge(:,:,1))/dxmin

		! Initialize A array and q array
		SELECT CASE(ntest)
			CASE(1) ! sine wave advection
				cdf_out = 'dg2d_sine_adv.nc'
				tfinal = 1D0

			DO i=1,nxplot
				DO j=1,nyplot
					q(i,j) = DSIN(2D0*PI*xplot(i))*DSIN(2D0*PI*yplot(j))
				ENDDO
			ENDDO

				DO i=1,nex
					DO j=1,ney
						DO l=0,norder
						DO m=0,norder
							DO s=0,dgorder
 							DO t=0,dgorder
		FOO(s,t) = qweights(s)*qweights(t)*DSIN(2D0*PI*x(1+s+(i-1)*(dgorder+1)))*DSIN(2D0*PI*y(1+t+(j-1)*(dgorder+1)))* &
					Leg(s,l)*Leg(t,m)

							ENDDO
							ENDDO
		A(i,j,l,m) = ((2D0*l+1)*(2D0*m+1)/4D0)*SUM(FOO)
						ENDDO
						ENDDO
					ENDDO
				ENDDO

			CASE(5) ! standard cosbell deformation
				cdf_out = 'dg2d_def_cosbell.nc'
				tfinal = 5D0

				DO j = 1,nyplot
				r(:,j) = 4.d0*sqrt((xplot-0.25d0)**2 + (yplot(j)-0.25d0)**2)
				ENDDO
				q = 0.d0
				WHERE (r.lt.1.d0)
				q = (0.5d0*(1.d0 + cos(pi*r)))
				END WHERE

				DO i=1,nex
					DO j=1,ney
						DO l=0,norder
						DO m=0,norder
							DO s=0,dgorder
 							DO t=0,dgorder

						r_s = 4D0*SQRT( (x(1+s+(i-1)*(dgorder+1))-0.25D0)**2 + (y(1+t+(j-1)*(dgorder+1))-0.25D0)**2 )
						ph = 0D0
						IF(r_s .lt. 1D0) THEN
							ph = (0.5d0*(1.d0 + cos(pi*r_s)))
						ENDIF
		
		FOO(s,t) = qweights(s)*qweights(t)*ph*Leg(s,l)*Leg(t,m)

							ENDDO
							ENDDO
		A(i,j,l,m) = ((2D0*l+1)*(2D0*m+1)/4D0)*SUM(FOO)
						ENDDO
						ENDDO
					ENDDO
				ENDDO

			CASE(6,10)	! smoothed cosbell deformation
				cdf_out = 'dg2d_smth_cosbell_rk4.nc'
				tfinal = 5D0
!				tfinal = (sqrt(3.D0)-1D0)

				DO j = 1,nyplot
				r(:,j) = 3.d0*sqrt((xplot(:)-0.4d0)**2 + (yplot(j)-0.4d0)**2)
				ENDDO
				q = 0.d0
				WHERE (r.lt.1.d0)
				q = (0.5d0*(1.d0 + cos(pi*r)))**3
				END WHERE

				DO i=1,nex
					DO j=1,ney
						DO l=0,norder
						DO m=0,norder
							DO s=0,dgorder
 							DO t=0,dgorder
						r_s = 3D0*SQRT( (x(1+s+(i-1)*(dgorder+1))-0.4D0)**2 + (y(1+t+(j-1)*(dgorder+1))-0.4D0)**2 )
						ph = 0D0
						IF(r_s .lt. 1D0) THEN
							ph = (0.5d0*(1.d0 + cos(pi*r_s)))**3
						ENDIF
		
		FOO(s,t) = qweights(s)*qweights(t)*ph*Leg(s,l)*Leg(t,m)

							ENDDO
							ENDDO
		A(i,j,l,m) = ((2D0*l+1)*(2D0*m+1)/4D0)*SUM(FOO)
						ENDDO
						ENDDO
					ENDDO
				ENDDO

			CASE(7) ! square wave
				cdf_out = 'dg2d_fdef_sqwave.nc'
				tfinal = 5D0

				DO j=1,nyplot
					r(:,j) = MAX(ABS(x-0.3d0),ABS(y(j)-0.5d0))/0.15d0
				ENDDO
			    q = 0.d0
			    WHERE (r.lt.1.d0)
			    		q = 1.d0
			    END WHERE

				DO i=1,nex
					DO j=1,ney
						DO l=0,norder
						DO m=0,norder
							DO s=0,dgorder
 							DO t=0,dgorder
						r_s = MAX(ABS(x(1+s+(i-1)*(dgorder+1))-0.3d0),ABS(y(1+t+(j-1)*(dgorder+1))-0.5d0))/0.15d0
						ph = 0D0
						IF(r_s .lt. 1D0) THEN
							ph = 1D0
						ENDIF
		
		FOO(s,t) = qweights(s)*qweights(t)*ph*Leg(s,l)*Leg(t,m)

							ENDDO
							ENDDO
		A(i,j,l,m) = ((2D0*l+1)*(2D0*m+1)/4D0)*SUM(FOO)
						ENDDO
						ENDDO
					ENDDO
				ENDDO
			
		ENDSELECT
	END SUBROUTINE init2d

	SUBROUTINE output2d(q,x,y,nxplot,nyplot,tval_in,cdf_out,ilvl,stat)
		IMPLICIT NONE

		! Inputs
		INTEGER, INTENT(IN) :: nxplot,nyplot,stat,ilvl
		CHARACTER(len=40), INTENT(IN) :: cdf_out
		REAL(KIND=8), INTENT(IN) :: tval_in
		REAL(KIND=8), DIMENSION(1:nxplot), INTENT(IN) :: x
		REAL(KIND=8), DIMENSION(1:nyplot), INTENT(IN) :: y
		REAL(KIND=8), DIMENSION(1:nxplot,1:nyplot), INTENT(IN) :: q
		
		! Outputs

		! Local variables
		INTEGER :: cdfid ! ID for netCDF file
		INTEGER, PARAMETER :: NDIMS = 3
		INTEGER :: ierr
	    INTEGER :: idq,idt,idx,idy,dimids(NDIMS)
	    INTEGER :: x_dimid, y_dimid, t_dimid
		INTEGER, DIMENSION(1:NDIMS) :: start, count
		CHARACTER(len=8) :: nxname,xname,nyname,yname,qname

		REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: tmp
		REAL(KIND=8), DIMENSION(1:nxplot) :: temp	
		INTEGER :: i,j

	    SAVE cdfid, idq, t_dimid, start, count

		IF(stat .eq. -1) THEN
			! Create netCDF file and time variables
			ierr = NF90_CREATE(TRIM(cdf_out),NF90_CLOBBER,cdfid)

			ierr = NF90_REDEF(cdfid)
			ierr = NF90_DEF_DIM(cdfid, "nt", ilvl+1, t_dimid)
			ierr = NF90_DEF_VAR(cdfid, "time", NF90_FLOAT, t_dimid,idt)
			ierr = NF90_ENDDEF(cdfid)

			! Calculate time at output levels (note ilvl=noutput)
			ALLOCATE(tmp(1:ilvl+1), STAT=ierr)
			DO i=0,ilvl
				tmp(i+1) = DBLE(i)*tval_in/DBLE(ilvl)
			ENDDO

			! Write t values
			ierr = NF90_PUT_VAR(cdfid,idt,tmp)
			DEALLOCATE(tmp, STAT=ierr)

			RETURN

		ELSEIF(stat .eq. 0) THEN
			! Create dimensions and variables for this level of runs (ilvl = p)
			start = 1
			count = 1

			! Define names of variables
			WRITE(nxname,'(a2,i1)') 'nx',ilvl
			WRITE(nyname,'(a2,i1)') 'ny',ilvl
			WRITE(xname, '(a1,i1)') 'x',ilvl
			WRITE(yname, '(a1,i1)') 'y',ilvl
			WRITE(qname, '(a1,i1)') 'Q',ilvl

			ierr = NF90_REDEF(cdfid)
			ierr = NF90_DEF_DIM(cdfid, TRIM(nxname), nxplot, x_dimid)
			ierr = NF90_DEF_DIM(cdfid, TRIM(nyname), nyplot, y_dimid)

			dimids(1) = x_dimid
			dimids(2) = y_dimid
			dimids(3) = t_dimid

			ierr = NF90_DEF_VAR(cdfid, TRIM(qname),NF90_FLOAT,dimids,idq)
			ierr = NF90_DEF_VAR(cdfid, TRIM(xname),NF90_FLOAT,x_dimid,idx)
			ierr = NF90_DEF_VAR(cdfid, TRIM(yname),NF90_FLOAT,y_dimid,idy)

			ierr = NF90_enddef(cdfid)

			! Write x and y values
			ierr = NF90_PUT_VAR(cdfid, idx, x)
			ierr = NF90_PUT_VAR(cdfid, idy, y)

			start(3) = 1

		ELSEIF(stat .eq. 1) THEN
			ierr = NF90_CLOSE(cdfid)
			RETURN
		ENDIF

		! Write out concentration field q
		count(1) = nxplot
		DO j=1,nyplot
			start(2) = j
			temp(:) = q(:,j)
			ierr = NF90_PUT_VAR(cdfid,idq,temp,start,count)
		ENDDO
		
		! Increment t level 
		start(3) = start(3) + 1 

	END SUBROUTINE output2d

END PROGRAM execute
