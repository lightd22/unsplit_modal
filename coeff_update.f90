! ==========================================
!
! Coefficent update for 2D Modal DG 
! Transport equation with cartesian gridding
!
! Updates coefficents A(i,j,l,m) for each element 
! See paper: Nair et al (2005) "A Discontinuous Galerkin Transport Scheme.."
!
! By: Devin Light ; Aug 2013
! ==========================================

SUBROUTINE coeff_update(q,A,u_tmp,v_tmp,uedge_tmp,vedge_tmp,qnodes,qweights,Leg,dLL,LL,L_xi_plot,L_eta_plot,dxel,dyel,&
						dt,dgorder,norder,nxplot,nyplot,nex,ney,nxiplot,netaplot,transient,time)
	IMPLICIT NONE

	! External functions
	REAL(KIND=8), EXTERNAL :: dadt ! RHS function for evolution ODE for expansion coefficent
	REAL(KIND=8), EXTERNAL :: vel_update ! Velocity update function
	
	! Inputs
	INTEGER, INTENT(IN) :: dgorder,norder,nxplot,nyplot,nex,ney,nxiplot,netaplot
	REAL(KIND=8), INTENT(IN) :: dt, dxel, dyel,time
	REAL(KIND=8), DIMENSION(0:dgorder), INTENT(IN) :: qnodes,qweights
	REAL(KIND=8), DIMENSION(1:(dgorder+1)*nex,1:(dgorder+1)*ney), INTENT(IN) :: u_tmp,v_tmp
	REAL(KIND=8), DIMENSION(1:nex,1:(dgorder+1)*ney), INTENT(IN) :: uedge_tmp
	REAL(KIND=8), DIMENSION(1:(dgorder+1)*nex,1:ney), INTENT(IN) :: vedge_tmp
	REAL(KIND=8), DIMENSION(0:dgorder,0:norder), INTENT(IN) :: Leg
	REAL(KIND=8), DIMENSION(0:norder,0:norder,0:dgorder,0:dgorder), INTENT(IN) :: LL,dLL
	REAL(KIND=8), DIMENSION(1:nxiplot,0:norder), INTENT(IN) :: L_xi_plot
	REAL(KIND=8), DIMENSION(1:netaplot,0:norder), INTENT(IN) :: L_eta_plot
	LOGICAL, INTENT(IN) :: transient

	! Outputs
	REAL(KIND=8), DIMENSION(1:nxplot,1:nyplot), INTENT(OUT) :: q
	REAL(KIND=8), DIMENSION(1:nex,1:ney,0:norder,0:norder), INTENT(INOUT) :: A

	! Local variables
	REAL(KIND=8), DIMENSION(1:nex,1:ney,0:norder,0:norder) :: A1,A2,A3,A4,A_tmp
	REAL(KIND=8), DIMENSION(1:(dgorder+1)*nex,1:(dgorder+1)*ney) :: u,v
	REAL(KIND=8), DIMENSION(1:nex,1:(dgorder+1)*ney) :: uedge
	REAL(KIND=8), DIMENSION(1:(dgorder+1)*nex,1:ney) :: vedge
	REAL(KIND=8), DIMENSION(0:dgorder,0:dgorder) :: u_int, v_int, R
	REAL(KIND=8), DIMENSION(0:nex,1:ney,0:dgorder) :: Fhat
	REAL(KIND=8), DIMENSION(1:nex,0:ney,0:dgorder) :: Ghat
	REAL(KIND=8), DIMENSION(0:norder,0:norder) :: FOO2D
	REAL(KIND=8) :: coef1
	INTEGER :: i,j,s,t,l,m

	REAL(KIND=4), DIMENSION(2) :: tstart,tend,t1,t2
	REAL(KIND=4) :: t0,tf,tick,tock

!go to 100
!---------

!	coef1 = dt/dxel/dyel
	coef1 = 1d0/dxel/dyel
	! ######################
	! Time step using Strong-Stability Preserving RK3
	! ######################

	! Update velocities
		u = u_tmp *vel_update(time)
		uedge = uedge_tmp*vel_update(time)
		v = v_tmp*vel_update(time)
		vedge = vedge_tmp*vel_update(time)

	! Update fluxes
	CALL numflx(A,Ghat,Fhat,uedge,vedge,Leg,dgorder,norder,nex,ney)

	! Perform first step of ssprk3
	DO i=1,nex
	DO j=1,ney
		u_int = u(1+(i-1)*(dgorder+1):i*(dgorder+1),1+(j-1)*(dgorder+1):j*(dgorder+1))
		v_int = v(1+(i-1)*(dgorder+1):i*(dgorder+1),1+(j-1)*(dgorder+1):j*(dgorder+1))
		CALL R_fill(R,A(i,j,:,:),LL,dgorder,norder)

		DO l=0,norder
		DO m=0,norder
	!		A1(i,j,l,m) = A(i,j,l,m)+ &
	!				coef1*dadt(i,j,l,m,R,Fhat,Ghat,u_int,v_int,qweights,LL,dLL,Leg,dxel,dyel,dgorder,norder,nex,ney)

			A1(i,j,l,m) = coef1*dadt(i,j,l,m,R,Fhat,Ghat,u_int,v_int,qweights,LL,dLL,Leg,dxel,dyel,dgorder,norder,nex,ney)
		ENDDO
		ENDDO
	ENDDO
	ENDDO

	IF( transient ) THEN
	! Update velocities
		u = u_tmp *vel_update(time+dt/2d0)
		uedge = uedge_tmp*vel_update(time+dt/2d0)
		v = v_tmp*vel_update(time+dt/2d0)
		vedge = vedge_tmp*vel_update(time+dt/2d0)
	ENDIF

	A_tmp = A + dt*A1/2d0
	! Update fluxes
	CALL numflx(A_tmp,Ghat,Fhat,uedge,vedge,Leg,dgorder,norder,nex,ney)
	! Perform second step of ssprk3
	DO i=1,nex
	DO j=1,ney
		u_int = u(1+(i-1)*(dgorder+1):i*(dgorder+1),1+(j-1)*(dgorder+1):j*(dgorder+1))
		v_int = v(1+(i-1)*(dgorder+1):i*(dgorder+1),1+(j-1)*(dgorder+1):j*(dgorder+1))
		CALL R_fill(R,A_tmp(i,j,:,:),LL,dgorder,norder)

		DO l=0,norder
		DO m=0,norder
!			A2(i,j,l,m) = (0.75D0)*A(i,j,l,m) + (0.25D0)*(A1(i,j,l,m) + &
!						coef1*dadt(i,j,l,m,R,Fhat,Ghat,u_int,v_int,qweights,LL,dLL,Leg, &
!						dxel,dyel,dgorder,norder,nex,ney))
			A2(i,j,l,m) = coef1*dadt(i,j,l,m,R,Fhat,Ghat,u_int,v_int,qweights,LL,dLL,Leg, &
						dxel,dyel,dgorder,norder,nex,ney)
		ENDDO
		ENDDO
	ENDDO
	ENDDO

	IF( transient ) THEN
	! Update velocities
		u = u_tmp*vel_update(time+dt/2d0)
		uedge = uedge_tmp*vel_update(time+dt/2d0)
		v = v_tmp*vel_update(time+dt/2d0)
		vedge = vedge_tmp*vel_update(time+dt/2d0)
	ENDIF
	A_tmp = A + dt*A2/2d0

	! Update fluxes, fill R array
	CALL numflx(A_tmp,Ghat,Fhat,uedge,vedge,Leg,dgorder,norder,nex,ney)
	! Perform final step of ssprk3
	DO i=1,nex
	DO j=1,ney
		u_int = u(1+(i-1)*(dgorder+1):i*(dgorder+1),1+(j-1)*(dgorder+1):j*(dgorder+1))
		v_int = v(1+(i-1)*(dgorder+1):i*(dgorder+1),1+(j-1)*(dgorder+1):j*(dgorder+1))
		CALL R_fill(R,A_tmp(i,j,:,:),LL,dgorder,norder)

		DO l=0,norder
		DO m=0,norder
!			A(i,j,l,m) = A(i,j,l,m)/3D0 + (2D0/3D0)*(A2(i,j,l,m) + &
!						coef1*dadt(i,j,l,m,R,Fhat,Ghat,u_int,v_int,qweights,LL,dLL,Leg, &
!						dxel,dyel,dgorder,norder,nex,ney))
			A3(i,j,l,m) = coef1*dadt(i,j,l,m,R,Fhat,Ghat,u_int,v_int,qweights,LL,dLL,Leg, &
						      dxel,dyel,dgorder,norder,nex,ney)
		ENDDO
		ENDDO
	ENDDO
	ENDDO

	!OPTIONAL RK4 STEP
	IF( transient ) THEN
	! Update velocities
		u = u_tmp*vel_update(time+dt)
		uedge = uedge_tmp*vel_update(time+dt)
		v = v_tmp*vel_update(time+dt)
		vedge = vedge_tmp*vel_update(time+dt)
	ENDIF
	A_tmp = A + dt*A3

	! Update fluxes, fill R array
	CALL numflx(A_tmp,Ghat,Fhat,uedge,vedge,Leg,dgorder,norder,nex,ney)
	! Perform final step of ssprk3
	DO i=1,nex
	DO j=1,ney
		u_int = u(1+(i-1)*(dgorder+1):i*(dgorder+1),1+(j-1)*(dgorder+1):j*(dgorder+1))
		v_int = v(1+(i-1)*(dgorder+1):i*(dgorder+1),1+(j-1)*(dgorder+1):j*(dgorder+1))
		CALL R_fill(R,A_tmp(i,j,:,:),LL,dgorder,norder)

		DO l=0,norder
		DO m=0,norder
!			A(i,j,l,m) = A(i,j,l,m)/3D0 + (2D0/3D0)*(A2(i,j,l,m) + &
!						coef1*dadt(i,j,l,m,R,Fhat,Ghat,u_int,v_int,qweights,LL,dLL,Leg, &
!						dxel,dyel,dgorder,norder,nex,ney))
			A4(i,j,l,m) = coef1*dadt(i,j,l,m,R,Fhat,Ghat,u_int,v_int,qweights,LL,dLL,Leg, &
						      dxel,dyel,dgorder,norder,nex,ney)
		ENDDO
		ENDDO
	ENDDO
	ENDDO

	A = A + dt/6d0*(A1+2D0*A2+2D0*A3+A4)

	! #########
	! END SPPRK3 TIMESTEP ; BEGIN FILLING q ARRAY
	! #########

100 continue

	DO i=1,nex
	DO j=1,ney

		DO s=1,nxiplot
		DO t=1,netaplot

			DO l=0,norder
			DO m=0,norder

				FOO2D(l,m) = A(i,j,l,m)*L_xi_plot(s,l)*L_eta_plot(t,m)

			ENDDO
			ENDDO
	
			q(s+(i-1)*nxiplot,t+(j-1)*netaplot) = SUM(FOO2D)

		ENDDO
		ENDDO
	ENDDO
	ENDDO

END SUBROUTINE coeff_update

REAL(KIND=8) FUNCTION dadt(i,j,l,m,Rin,Fhat,Ghat,u_int,v_int,qweights,LL,dLL,Leg,dxel,dyel,dgorder,norder,nex,ney)
! Computes dadt for element (i,j) for use in SSPRK3 update step
! Does this in 3 steps: dadt = [(2s+1)*(2t+1)/(2*dxel*dyel)](A+B+C) where
!	1) A = Int(Int( dyel*dPs/dxi*Pt*F+dxel*dPt/deta*Ps*G )) -- Interior contribution to change
!	2) B = -dxel*Int( Ps*(Ghat(eta=1) - (-1)**t*Ghat(eta=-1)) ) -- Contribution through North/South faces
!	3) C = -dyel*Int( Pt*(Fhat(xi=1) - (-1)**s*Fhat(xi=-1)) ) -- Contribution through East/West faces
!
	IMPLICIT NONE

	! Inputs
	INTEGER, INTENT(IN) :: i,j,l,m,dgorder,norder,nex,ney
	REAL(KIND=8), INTENT(IN) :: dxel,dyel
 	REAL(KIND=8), DIMENSION(0:dgorder,0:dgorder), INTENT(IN) :: u_int,v_int,Rin
	REAL(KIND=8), DIMENSION(0:dgorder,0:norder), INTENT(IN) :: Leg
	REAL(KIND=8), DIMENSION(0:norder,0:norder,0:dgorder,0:dgorder), INTENT(IN) :: LL,dLL
	REAL(KIND=8), DIMENSION(0:dgorder), INTENT(IN) :: qweights
	REAL(KIND=8), DIMENSION(0:nex,1:ney,0:dgorder), INTENT(IN) :: Fhat
	REAL(KIND=8), DIMENSION(1:nex,0:ney,0:dgorder), INTENT(IN) :: Ghat

	! Local variables
	INTEGER :: s,t,k,n
	REAL(KIND=8) :: A,B,C
	REAL(KIND=8), DIMENSION(0:dgorder) :: FOO1D
	REAL(KIND=8), DIMENSION(0:dgorder,0:dgorder) :: FOO2D

	! Step 1: Compute A
!	DO t=0,dgorder
!	DO s=0,dgorder
!		FOO2D(s,t) = qweights(s)*qweights(t)*Rin(s,t)*( dyel*dLL(l,m,s,t)*u_int(s,t) + dxel*dLL(m,l,t,s)*v_int(s,t) )
!	ENDDO
!	ENDDO
!	A = SUM(FOO2D)
!	DO s=0,dgorder
!	DO t=0,dgorder
!		FOO2D(s,t) = dLL(m,l,t,s)
!	ENDDO
!	ENDDO

	FOO2D = TRANSPOSE(dLL(m,l,:,:))

	A = SUM(spread(qweights(0:dgorder),dim=2,ncopies=dgorder+1)*spread(qweights(0:dgorder),dim=1,ncopies=dgorder+1)*Rin &
			*(dyel*dLL(l,m,:,:)*u_int+dxel*FOO2D*v_int))

	! Step 2: Compute B
!	DO s=0,dgorder
!		FOO1D(s) = qweights(s)*Leg(s,l)*( Ghat(i,j,s) - ((-1D0)**m)*Ghat(i,j-1,s) )
!	ENDDO
!	B = -1D0*dxel*SUM(FOO1D)
	B = -1D0*dxel*SUM(qweights*Leg(:,l)*(Ghat(i,j,:)-((-1D0)**m)*Ghat(i,j-1,:)))

	! Step 3: Compute C
!	DO t=0,dgorder
!		FOO1D(t) = qweights(t)*Leg(t,m)*( Fhat(i,j,t) - ((-1D0)**l)*Fhat(i-1,j,t) )
!	ENDDO
!	C = -1D0*dyel*SUM(FOO1D)
	C = -1D0*dyel*SUM(qweights*Leg(:,m)*(Fhat(i,j,:)-((-1D0)**l)*Fhat(i-1,j,:)))

	dadt = ( (2D0*DBLE(l)+1D0)*(2D0*DBLE(m)+1D0)/2D0 )*(A+B+C)

END FUNCTION dadt

SUBROUTINE numflx(A,Ghat,Fhat,uedge,vedge,Leg,dgorder,norder,nex,ney)
! Computes approximate solution at quadrature locations interior to element and fills numerical flux arrays
	IMPLICIT NONE
	! Inputs
	INTEGER, INTENT(IN) :: dgorder, norder,nex,ney
	REAL(KIND=8), DIMENSION(0:dgorder,0:norder), INTENT(IN) :: Leg
	REAL(KIND=8), DIMENSION(1:nex,1:ney,0:norder,0:norder), INTENT(IN) :: A
	REAL(KIND=8), DIMENSION(1:nex,1:(dgorder+1)*ney), INTENT(IN) :: uedge
	REAL(KIND=8), DIMENSION(1:(dgorder+1)*nex,1:ney), INTENT(IN) :: vedge

	! Outputs
	REAL(KIND=8), DIMENSION(0:nex,1:ney,0:dgorder), INTENT(OUT) :: Fhat
	REAL(KIND=8), DIMENSION(1:nex,0:ney,0:dgorder), INTENT(OUT) :: Ghat

	! Local Variables
	INTEGER :: i,j,k,l,m,s,t
	REAL(KIND=8), DIMENSION(0:norder,0:norder) :: FOO
	REAL(KIND=8), DIMENSION(0:dgorder) :: u_tmp, v_tmp
	REAL(KIND=8), DIMENSION(0:norder) :: arr1,arr2
	REAL(KIND=8), DIMENSION(1:nex+1,1:ney+1,0:norder,0:norder) :: A_tmp
	
	INTEGER :: which_el
	REAL(KIND=8), DIMENSION(0:norder,0:norder) :: which_sign

	! Set up arrays for computing P_l(+-1)
	arr1 = 1D0 ! P_l(1) = 1 for all orders l
	arr2 = (/ ((-1D0)**i , i=0,norder) /) ! P_l(-1) = (-1)**l

	! A_tmp is the periodically extended version of A
	A_tmp(:,:,:,:) = 0D0
	A_tmp(1:nex,1:ney,0:norder,0:norder) = A(1:nex,1:ney,0:norder,0:norder)
	A_tmp(nex+1,1:ney,:,:) = A(1,1:ney,:,:)
	A_tmp(1:nex,ney+1,:,:) = A(1:nex,1,:,:)

	DO i=1,nex
	DO j=1,ney

		u_tmp = uedge(i,1+(j-1)*(dgorder+1):j*(dgorder+1))
		v_tmp = vedge(1+(i-1)*(dgorder+1):i*(dgorder+1),j)

		! Compute numerical flux through North edges (Ghat) (There are dgorder+1 nodes per element)
		DO k=0,dgorder
			which_sign(0,:) = arr1 - (arr1-arr2)*(1-INT(SIGN(1D0,v_tmp(k))) )/2
			DO l=1,norder
				which_sign(l,:) = which_sign(0,:)
			ENDDO
			which_el = j + (1-INT(SIGN(1D0,v_tmp(k))) )/2

			Ghat(i,j,k) = v_tmp(k)*SUM(Leg(k,:)*SUM(A_tmp(i,which_el,:,:)*which_sign,2))
!			IF(v_tmp(k) .ge. 0D0) THEN

!				DO m=0,norder
!				DO l=0,norder
!					FOO(l,m) = A_tmp(i,j,l,m)*Leg(k,l)
!				ENDDO
!				ENDDO	
					
!			ELSE

!				DO m=0,norder
!				DO l=0,norder
!					FOO(l,m) = A_tmp(i,j+1,l,m)*Leg(k,l)*(-1D0)**(m)
!				ENDDO
!				ENDDO	
!			ENDIF
!			Ghat(i,j,k) = v_tmp(k)*SUM(FOO)
		ENDDO
		! Extend periodically
		Ghat(1:nex,0,:) = Ghat(1:nex,ney,:)

		! Compute numerical flux through East edges (Fhat) (There are dgorder+1 nodes per element)
		DO k=0,dgorder
			which_sign(:,0) = arr1 - (arr1-arr2)*(1-INT(SIGN(1D0,u_tmp(k))) )/2
			DO m=1,norder
				which_sign(:,m) = which_sign(:,0)
			ENDDO
			which_el = i + (1-INT(SIGN(1D0,u_tmp(k))) )/2
			Fhat(i,j,k) = u_tmp(k)*SUM(Leg(k,:)*SUM(A_tmp(which_el,j,:,:)*which_sign,1))

!			IF(u_tmp(k) .ge. 0D0) THEN

!				DO m=0,norder
!				DO l=0,norder
!					FOO(l,m) = A_tmp(i,j,l,m)*Leg(k,m)
!				ENDDO
!				ENDDO	
					
!			ELSE

!				DO m=0,norder
!				DO l=0,norder
!					FOO(l,m) = A_tmp(i+1,j,l,m)*Leg(k,m)*(-1D0)**(l)
!				ENDDO
!				ENDDO	
!			ENDIF
!			Fhat(i,j,k) = u_tmp(k)*SUM(FOO)
		ENDDO
		! Extend periodically
		Fhat(0,1:ney,:) = Fhat(nex,1:ney,:)

	ENDDO
	ENDDO
	
END SUBROUTINE numflx

SUBROUTINE R_fill(R,Ain,LL,dgorder,norder)
! Compute approx soln at interior quad points (xi_s,eta_t) for element (i,j)
	IMPLICIT NONE
	! Inputs
	INTEGER, INTENT(IN) :: dgorder,norder
	REAL(KIND=8), DIMENSION(0:norder,0:norder,0:dgorder,0:dgorder), INTENT(IN) :: LL
	REAL(KIND=8), DIMENSION(0:norder,0:norder) :: Ain
	! Outputs
	REAL(KIND=8), DIMENSION(0:dgorder,0:dgorder), INTENT(OUT) :: R
	! Local variables
	INTEGER :: s,t

	DO t=0,dgorder
		DO s=0,dgorder
			R(s,t) = SUM(Ain*LL(:,:,s,t))
		ENDDO
	ENDDO

END SUBROUTINE R_fill

REAL(KIND=8) FUNCTION vel_update(t)
		IMPLICIT NONE
		REAL(KIND=8), INTENT(IN) :: t
	    REAL(KIND=8) :: pi
	    REAL(KIND=8), parameter :: t_period = 5.d0

	    pi = DACOS(-1D0)
	    vel_update = DCOS(pi*t/t_period)
END FUNCTION vel_update

