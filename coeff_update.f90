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

SUBROUTINE coeff_update(q,A,u,v,uedge,vedge,qnodes,qweights,Leg,dLeg,L_xi_plot,L_eta_plot,dxel,dyel,dt,dgorder,nxplot,&
						nyplot,nex,ney,nxiplot,netaplot)
	IMPLICIT NONE

	! External functions
	REAL(KIND=8), EXTERNAL :: dadt ! RHS function for evolution ODE for kth expansion coefficent
	
	! Inputs
	INTEGER, INTENT(IN) :: dgorder,nxplot,nyplot,nex,ney,nxiplot,netaplot
	REAL(KIND=8), INTENT(IN) :: dt, dxel, dyel
	REAL(KIND=8), DIMENSION(0:dgorder), INTENT(IN) :: qnodes,qweights
	REAL(KIND=8), DIMENSION(1:(dgorder+1)*nex,1:(dgorder+1)*ney), INTENT(IN) :: u,v
	REAL(KIND=8), DIMENSION(1:nex,1:(dgorder+1)*ney), INTENT(IN) :: uedge
	REAL(KIND=8), DIMENSION(1:(dgorder+1)*nex,1:ney), INTENT(IN) :: vedge
	REAL(KIND=8), DIMENSION(0:dgorder,0:dgorder), INTENT(IN) :: Leg,dLeg	
	REAL(KIND=8), DIMENSION(1:nxiplot,0:dgorder), INTENT(IN) :: L_xi_plot
	REAL(KIND=8), DIMENSION(1:netaplot,0:dgorder), INTENT(IN) :: L_eta_plot

	! Outputs
	REAL(KIND=8), DIMENSION(1:nxplot,1:nyplot), INTENT(OUT) :: q
	REAL(KIND=8), DIMENSION(1:nex,1:ney,0:dgorder,0:dgorder), INTENT(INOUT) :: A

	! Local variables
	REAL(KIND=8), DIMENSION(1:nex,1:ney,0:dgorder,0:dgorder) :: A1,A2, R
	REAL(KIND=8), DIMENSION(0:dgorder,0:dgorder) :: u_tmp, v_tmp, R_tmp
	REAL(KIND=8), DIMENSION(0:nex,1:ney,0:dgorder) :: Fhat
	REAL(KIND=8), DIMENSION(1:nex,0:ney,0:dgorder) :: Ghat
	REAL(KIND=8), DIMENSION(0:dgorder,0:dgorder) :: FOO2D
	REAL(KIND=8) :: coef1
	INTEGER :: i,j,s,t,l,m

	REAL(KIND=4), DIMENSION(2) :: tstart,tend,t1,t2
	REAL(KIND=4) :: t0,tf,tick,tock

	coef1 = dt/dxel/dyel
	! ######################
	! Time step using Strong-Stability Preserving RK3
	! ######################

	tick = etime(t1)

	! Update fluxes, fill R array
	CALL phitld_numflx_fill(A,R,Ghat,Fhat,uedge,vedge,Leg,dgorder,nex,ney)

	tock = etime(t2) - tick
	write(*,*) 'Numflx:',tock,'s'

	! Perform first step of ssprk3
	DO i=1,nex
	DO j=1,ney
		R_tmp = R(i,j,:,:)
		u_tmp = u(1+(i-1)*(dgorder+1):i*(dgorder+1),1+(j-1)*(dgorder+1):j*(dgorder+1))
		v_tmp = v(1+(i-1)*(dgorder+1):i*(dgorder+1),1+(j-1)*(dgorder+1):j*(dgorder+1))

		DO s=0,dgorder
		DO t=0,dgorder
			A1(i,j,s,t) = A(i,j,s,t)+ &
					coef1*dadt(i,j,s,t,R_tmp,Fhat,Ghat,u_tmp,v_tmp,qweights,Leg,dLeg,dxel,dyel,dgorder,nex,ney)
		ENDDO
		ENDDO
	ENDDO
	ENDDO

	tock = etime(t2) - tick
	write(*,*) 'RK3 1:',tock,'s'


	! Update fluxes, fill R array
	CALL phitld_numflx_fill(A1,R,Ghat,Fhat,uedge,vedge,Leg,dgorder,nex,ney)
	! Perform second step of ssprk3
	DO i=1,nex
	DO j=1,ney
		R_tmp = R(i,j,:,:)
		u_tmp = u(1+(i-1)*(dgorder+1):i*(dgorder+1),1+(j-1)*(dgorder+1):j*(dgorder+1))
		v_tmp = v(1+(i-1)*(dgorder+1):i*(dgorder+1),1+(j-1)*(dgorder+1):j*(dgorder+1))

		DO s=0,dgorder
		DO t=0,dgorder
			A2(i,j,s,t) = (0.75D0)*A(i,j,s,t) + (0.25D0)*(A1(i,j,s,t) + &
						coef1*dadt(i,j,s,t,R_tmp,Fhat,Ghat,u_tmp,v_tmp,qweights,Leg,dLeg, &
						dxel,dyel,dgorder,nex,ney))
		ENDDO
		ENDDO
	ENDDO
	ENDDO

	! Update fluxes, fill R array
	CALL phitld_numflx_fill(A2,R,Ghat,Fhat,uedge,vedge,Leg,dgorder,nex,ney)
	! Perform final step of ssprk3
	DO i=1,nex
	DO j=1,ney
		R_tmp = R(i,j,:,:)
		u_tmp = u(1+(i-1)*(dgorder+1):i*(dgorder+1),1+(j-1)*(dgorder+1):j*(dgorder+1))
		v_tmp = v(1+(i-1)*(dgorder+1):i*(dgorder+1),1+(j-1)*(dgorder+1):j*(dgorder+1))

		DO s=0,dgorder
		DO t=0,dgorder
			A(i,j,s,t) = A(i,j,s,t)/3D0 + (2D0/3D0)*(A2(i,j,s,t) + &
						coef1*dadt(i,j,s,t,R_tmp,Fhat,Ghat,u_tmp,v_tmp,qweights,Leg,dLeg, &
						dxel,dyel,dgorder,nex,ney))
		ENDDO
		ENDDO
	ENDDO
	ENDDO

	! #########
	! END SPPRK3 TIMESTEP ; BEGIN FILLING q ARRAY
	! #########

	DO i=1,nex
	DO j=1,ney

		DO s=1,nxiplot
		DO t=1,netaplot

			DO l=0,dgorder
			DO m=0,dgorder

				FOO2D(l,m) = A(i,j,l,m)*L_xi_plot(s,l)*L_eta_plot(t,m)

			ENDDO
			ENDDO
	
			q(s+(i-1)*nxiplot,t+(j-1)*netaplot) = SUM(FOO2D)

		ENDDO
		ENDDO
	ENDDO
	ENDDO

END SUBROUTINE coeff_update

REAL(KIND=8) FUNCTION dadt(i,j,s,t,Rin,Fhat,Ghat,u_int,v_int,qweights,Leg,dLeg,dxel,dyel,dgorder,nex,ney)
! Computes dadt for element (i,j) for use in SSPRK3 update step
! Does this in 3 steps: dadt = [(2s+1)*(2t+1)/(2*dxel*dyel)](A+B+C) where
!	1) A = Int(Int( dyel*dPs/dxi*Pt*F+dxel*dPt/deta*Ps*G )) -- Interior contribution to change
!	2) B = -dxel*Int( Ps*(Ghat(eta=1) - (-1)**t*Ghat(eta=-1)) ) -- Contribution through North/South face
!	3) C = -dyel*Int( Pt*(Fhat(xi=1) - (-1)**s*Fhat(xi=-1)) ) -- Contribution through East/West face
!
	IMPLICIT NONE

	! Inputs
	INTEGER, INTENT(IN) :: i,j,s,t,dgorder,nex,ney
	REAL(KIND=8), INTENT(IN) :: dxel,dyel
 	REAL(KIND=8), DIMENSION(0:dgorder,0:dgorder), INTENT(IN) :: Leg,dLeg,u_int,v_int,Rin
	REAL(KIND=8), DIMENSION(0:dgorder), INTENT(IN) :: qweights
	REAL(KIND=8), DIMENSION(0:nex,1:ney,0:dgorder), INTENT(IN) :: Fhat
	REAL(KIND=8), DIMENSION(1:nex,0:ney,0:dgorder), INTENT(IN) :: Ghat


	! Local variables
	INTEGER :: l,m
	REAL(KIND=8) :: A,B,C
	REAL(KIND=8), DIMENSION(0:dgorder) :: FOO1D
	REAL(KIND=8), DIMENSION(0:dgorder,0:dgorder) :: FOO2D

	! Step 1: Compute A
	DO l=0,dgorder
	DO m=0,dgorder
		FOO2D(l,m) = qweights(l)*qweights(m)*Rin(l,m)*( dyel*dLeg(l,s)*Leg(m,t)*u_int(l,m) + dxel*dLeg(m,t)*Leg(l,s)*v_int(l,m) )
	ENDDO
	ENDDO
	A = SUM(FOO2D)

	! Step 2: Compute B
	DO l=0,dgorder
		FOO1D(l) = qweights(l)*Leg(l,s)*( Ghat(i,j,l) - ((-1D0)**t)*Ghat(i,j-1,l) )
	ENDDO
	B = -1D0*dxel*SUM(FOO1D)

	! Step 3: Compute C
	DO m=0,dgorder
		FOO1D(m) = qweights(m)*Leg(m,t)*( Fhat(i,j,m) - ((-1D0)**s)*Fhat(i-1,j,m) )
	ENDDO
	C = -1D0*dyel*SUM(FOO1D)

	dadt = ( (2D0*DBLE(s)+1D0)*(2D0*DBLE(t)+1D0)/2D0 )*(A+B+C)

END FUNCTION dadt

SUBROUTINE phitld_numflx_fill(A,R,Ghat,Fhat,uedge,vedge,Leg,dgorder,nex,ney)
! Computes approximate solution at quadrature locations interior to element and fills numerical flux arrays
	IMPLICIT NONE
	! Inputs
	INTEGER, INTENT(IN) :: dgorder, nex,ney
	REAL(KIND=8), DIMENSION(0:dgorder,0:dgorder), INTENT(IN) :: Leg
	REAL(KIND=8), DIMENSION(1:nex,1:ney,0:dgorder,0:dgorder), INTENT(IN) :: A
	REAL(KIND=8), DIMENSION(1:nex,1:(dgorder+1)*ney), INTENT(IN) :: uedge
	REAL(KIND=8), DIMENSION(1:(dgorder+1)*nex,1:ney), INTENT(IN) :: vedge

	! Outputs
	REAL(KIND=8), DIMENSION(1:nex,1:ney,0:dgorder,0:dgorder), INTENT(OUT) :: R
	REAL(KIND=8), DIMENSION(0:nex,1:ney,0:dgorder), INTENT(OUT) :: Fhat
	REAL(KIND=8), DIMENSION(1:nex,0:ney,0:dgorder), INTENT(OUT) :: Ghat

	! Local Variables
	INTEGER :: i,j,k,l,m,s,t
	REAL(KIND=8), DIMENSION(0:dgorder,0:dgorder) :: FOO
	REAL(KIND=8), DIMENSION(0:dgorder) :: u_tmp, v_tmp
	REAL(KIND=8), DIMENSION(1:nex+1,1:ney+1,0:dgorder,0:dgorder) :: A_tmp

	REAL(KIND=4), DIMENSION(2) :: tstart,tend,t1,t2
	REAL(KIND=4) :: t0,tf,tick,tock

	! A_tmp is the periodically extended version of A
	A_tmp(:,:,:,:) = 0D0
	A_tmp(1:nex,1:ney,0:dgorder,0:dgorder) = A(1:nex,1:ney,0:dgorder,0:dgorder)
	A_tmp(nex+1,1:ney,:,:) = A(1,1:ney,:,:)
	A_tmp(1:nex,ney+1,:,:) = A(1:nex,1,:,:)

	DO i=1,nex
	DO j=1,ney

		u_tmp = uedge(i,1+(j-1)*(dgorder+1):j*(dgorder+1))
		v_tmp = vedge(1+(i-1)*(dgorder+1):i*(dgorder+1),j)

		! Compute approx soln at interior quad points (xi_s,eta_t)
		DO s=0,dgorder
		DO t=0,dgorder

			DO m=0,dgorder
			DO l=0,dgorder
				FOO(l,m) = A(i,j,l,m)*Leg(s,l)*Leg(t,m)
			ENDDO
			ENDDO

		R(i,j,s,t) = SUM(FOO)

		ENDDO
		ENDDO		
		tick = etime(t1)
		! Compute numerical flux through North edges (Ghat) (There are dgorder+1 nodes per element)
		DO k=0,dgorder
			IF(v_tmp(k) .ge. 0D0) THEN

				DO m=0,dgorder
				DO l=0,dgorder
					FOO(l,m) = A_tmp(i,j,l,m)*Leg(k,l)
				ENDDO
				ENDDO	
					
			ELSE

				DO m=0,dgorder
				DO l=0,dgorder
					FOO(l,m) = A_tmp(i,j+1,l,m)*Leg(k,l)*(-1D0)**(m)
				ENDDO
				ENDDO	
			ENDIF
			Ghat(i,j,k) = v_tmp(k)*SUM(FOO)
		ENDDO
		! Extend periodically
		Ghat(1:nex,0,:) = Ghat(1:nex,ney,:)

		tock = etime(t2) - tick
!		write(*,*) 'Ghat fill:',tock,'s'

		tick = etime(t1)

		! Compute numerical flux through East edges (Fhat) (There are dgorder+1 nodes per element)
		DO k=0,dgorder
			IF(u_tmp(k) .ge. 0D0) THEN

				DO m=0,dgorder
				DO l=0,dgorder
					FOO(l,m) = A_tmp(i,j,l,m)*Leg(k,m)
				ENDDO
				ENDDO	
					
			ELSE

				DO m=0,dgorder
				DO l=0,dgorder
					FOO(l,m) = A_tmp(i+1,j,l,m)*Leg(k,m)*(-1D0)**(l)
				ENDDO
				ENDDO	
			ENDIF
			Fhat(i,j,k) = u_tmp(k)*SUM(FOO)
		ENDDO
		! Extend periodically
		Fhat(0,1:ney,:) = Fhat(nex,1:ney,:)

		tock = etime(t2)-tick
		tf = etime(t2)-t0
!		write(*,*) 'Fhat fill:',tock,'s'
!		write(*,*) 'Total:',tf,'s'
		
	ENDDO
	ENDDO
	
END SUBROUTINE phitld_numflx_fill
