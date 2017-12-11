!Solving 1D-Sod's problem by finite volume approach
!Weno scheme is used for spatial evolution and 3rd order RK is used for temporal evolution 
program weno5
	implicit none
	INTERFACE flux
		Function flux (q)
		double precision,intent(in):: q(:,:)
		double precision flux (size(q,1),size(q,2))
		End Function
		
		Function WENO5LF1d (q,lambda,dx)
		double precision,intent(in):: q(:,:)
		double precision,intent(in):: lambda
		double precision,intent(in):: dx
		double precision WENO5LF1d (size(q,1),size(q,2))
		End Function 
	End INTERFACE
!
!          Interface triple
!             Function triple (m)
!                Real, intent(in) :: m(:,:)
!                Real triple (size(m,1),size(m,2))
!             End Function
!          End Interface
!          Real, Dimension(3,3) :: z,x
!          x(:,:) = 0.0
!          x(1,1) = 1.0
!          x(2,2) = 1.0
!          x(3,3) = 1.0
!          z = triple(x)
!          Print *, shape(z)
!          Print *, z
!    End Program
!    Function triple(m)
!        Real, intent(in) :: m(:,:)
!        Real triple(size(m,1),size(m,2))
!        triple =  m(:,:) * 3.0
!    End Function
	
	!Declaration of variables---------------------------------------------------------------------
	double precision                  :: lambda,dt,dt0,lambda0,PRL,alpha,CFL,tFinal,t=0,nE,gama,dx,x_middle,xstart,xend
	integer                           :: it=0,i,j
	integer, parameter                :: N=201
	double precision, dimension(N)    :: x,rho0,p0,u0,E0,a0,rho,u,E,a,p
	double precision, dimension(3,N)  :: q0,q,F,dF
	double precision, dimension(2)    :: pi,ui,rhoi
	real                              :: cfl1,tEnd
        character(len=1024) :: filename,ff,ff2
	!real, dimension(3,3)             :: temp, temp_shift
	!Assigning values-----------------------------------------------------------------------------
	CFL=0.6
	tFinal=2.0
	nE=200
	gama=1.4
	xstart=-1
	xend=1
	!Testing cshift function
		!temp=reshape((/1,2,3,4,5,6,7,8,9/),(/ 3,3 /))
		!temp_shift=cshift(temp,-1,dim=2)
		!print *,'(3i3)', temp(1,:)
		!print *,'(3i3)', temp(2,:)
		!print *,'(3i3)', temp(3,:)
		!
		!print *,'(3i3)', temp_shift(1,:)
		!print *,'(3i3)', temp_shift(2,:)
		!print *,'(3i3)', temp_shift(3,:)
		
	dx=(xend-xstart)/nE
	x(1)=xstart
	do i=2,N
	x(i)=x(i-1)+dx
	end do
	!Defining initial conditions-------------------------------------------------------------------
	x_middle=(xend+xstart)/2
	!print *," x_middle=",x_middle
	!Sod's Problem
	rhoi=(/1.0,0.125/)
	pi=(/1.0,0.1/)
	ui=(/0.0,0.0/)
	
	do i=1,N
	if (x(i)>x_middle) then
	p0(i)=pi(2)
	rho0(i)=rhoi(2)
	u0(i)=ui(2)
	else
	p0(i)=pi(1)
	rho0(i)=rhoi(1)
	u0(i)=ui(1)
	end if
	end do
	do i=1,N
	E0(i)=p0(i)/(rho0(i)*(gama-1))+0.5*u0(i)*u0(i)
	!print *,i,E0(i) 
	a0(i)=SQRT(gama*p0(i)/rho0(i))
	
	q0(1,i)=rho0(i)
	q0(2,i)=rho0(i)*u0(i)
	q0(3,i)=rho0(i)*E0(i)
	end do

	!F=flux(q0)
	!print *, F
	!print *,"a0(1)=", a0(1),"E0(1)=",E0(1)
	!print *,"u0(1)=",u0(1)
	!print *,"F0(1,1)=",F(1,1)
	!print *,"F0(2,1)=",F(2,1)
	!print *,"F0(3,1)=",F(3,1)
	!Exact solution----------------------------------------------------------------------------------
	alpha=(gama+1)/(gama-1)
	!PRL=p(2)/p(1)
	!cright=SQRT(gama*p(2)/rho(2))
	!cleft=SQRT(gama*p(1)/rho(1))
	!CRL=cright/cleft
	!MACHLEFT=(u(2)-u(1))/cleft
	lambda0=MAXVAL(ABS(u0)+a0)
	dt0=CFL*dx/lambda0
	!Solver loop-------------------------------------------------------------------------------------
	!Load initial conditions
	q=q0
	dt=dt0
	lambda=lambda0
        dF=WENO5LF1d(q,lambda,dx) 
        !print *,q
	!print *,dF
	do while(t<tFinal)
		!RK Initial step
		!1st stage
                       q0=q
                       dF=WENO5LF1d(q,lambda,dx) 
			q = q0-dt*dF
			do j=1,3
				q(j,1)=q0(j,1)
				q(j,N)=q0(j,N)
			end do
		!RK second stage
			dF=WENO5LF1d(q,lambda,dx)
			q =0.75*q0+0.25*(q-dt*dF)
			do j=1,3
				q(j,1)=q0(j,1)
				q(j,N)=q0(j,N)
			end do
		!Third stage
			dF=WENO5LF1d(q,lambda,dx)
			q =(q0+2*(q-dt*dF))/3.0
			do j=1,3
				q(j,1)=q0(j,1)
				q(j,N)=q0(j,N)
			end do
		!Computing primary properties
			do i=1,N
				rho(i)=q(1,i)
				u(i)=q(2,i)/rho(i)
				E(i)=q(3,i)/rho(i)
				p(i)=(gama-1)*rho(i)*(E(i)-0.5*u(i)*u(i))
				a(i)=SQRT(gama*p(i)/rho(i))
			end do
		!Update dt and lambda
			lambda=maxval(abs(u)+a)
			dt=CFL*dx/lambda
                	if(t+dt>tFinal)then
				dt=tFinal-t
                	end if



		!Update time and iteration count
			t=t+dt
			it=it+1
		!Plot------------------------------------------------------------------------------------------

       write (filename, "(A,I4)") "datarho",it+1000 
	OPEN(UNIT=1000,FILE=filename)
       DO i=1,n
       WRITE(1000,*) x(i),rho(i)
       ENDDO
       CLOSE(1000)
	!write(ff,*) "filename=",filename
	!write(ff2,*) "gnuplot -e ",ff, " -p drawrho_plot.plt"
	end do

	print *,"it=",it,"t=",t 

!print *,u


!        CALL SYSTEM("gnuplot -p drawrho_plotmulti.plt")        
!       OPEN(UNIT=1000,FILE='datau')
!       DO i=1,n
!       WRITE(1000,*) x(i),u(i)
!       ENDDO
!       CLOSE(1000)
!       CALL SYSTEM('gnuplot -p datau_plot.plt')       
!       OPEN(UNIT=1000,FILE='datap')
!       DO i=1,n
!       WRITE(1000,*) x(i),p(i)
!       ENDDO
!       CLOSE(1000)
!       CALL SYSTEM('gnuplot -p datap_plot.plt')       

 

end program weno5
!Function triple(m)
!    Real, intent(in) :: m(:,:)
!    Real triple(size(m,1),size(m,2))
!    triple =  m(:,:) * 3.0
!End Function
Function flux(q)
       integer, parameter   :: N=201
       double precision, intent(in)     :: q(:,:)
       double precision flux(size(q,1),size(q,2))
       double precision, dimension(N)   :: rho,E,p,u
       real                 :: gama=1.4
       integer              :: i
       do i=1,N
       	rho(i)=q(1,i)
       	u(i)=q(2,i)/rho(i)
       	E(i)=q(3,i)/rho(i)
       	p(i)=(gama-1)*rho(i)*(E(i)-0.5*u(i)*u(i))
       	flux(1,i)=rho(i)*u(i)
       	flux(2,i)=rho(i)*u(i)*u(i)+p(i)
       	flux(3,i)=u(i)*(rho(i)*E(i)+p(i))
       end do
       !print *,"u(1)=",u(1)
       !print *,"t=",t
       !print *,"rho(1)=",rho(1)
       !print *,"p(1)=",p(1)
       !print *,"E(1)=",E(1)
       !print *,"F0(1,1)=",flux(1,1)
       !print *,"F0(2,1)=",flux(2,1)
       !print *,"F0(3,1)=",flux(3,1)
end Function

Function WENO5LF1d(q,lambda,dx)
        INTERFACE flux
                Function flux (q)
                double precision,intent(in):: q(:,:)
                double precision flux (size(q,1),size(q,2))
                End Function
	end INTERFACE
        integer, parameter   :: N=201
	double precision, intent(in)     :: q(:,:)
	double precision, intent(in)     :: lambda,dx
	double precision                 :: alpha,gama=1.4,eps=1E-006,kp1=0.1,kp2=0.6,kp3=0.3,kn1=0.3,kn2=0.6,kn3=0.1
	double precision, dimension(3,N) :: F,Fp,Fn,Fcin1,Fcin2,Fcin3,Fcin,Fcip1,Fcip2,Fcip3,Fcip,Fpf,Fpff,Fpb,Fpbb,Fnf,Fnff,Fnbb,Fnb
        double precision, dimension(3,N) :: wn1,wn2,wn3,wwn1,wwn2,wwn3,betan1,betan2,betan3 
	double precision, dimension(3,N) :: wp1,wp2,wp3,wwp1,wwp2,wwp3,betap1,betap2,betap3
	real,dimension(3,3)		 :: temp
	double precision WENO5LF1d(size(q,1),size(q,2))
        !print *,"I entered WENO5LF1d"
        !print *,"epsilon=",eps
	F=flux(q)
        !print *,F
	Fp=0.5*(F+lambda*q)
	Fn=0.5*(F-lambda*q)
	!print *,"Fp and Fn generated"
	Fn=cshift(Fn,1,dim=2)
 	!print *,Fn 
	!Positive flux
	Fpf=cshift(Fp,1,dim=2)
	Fpff=cshift(Fpf,1,dim=2)
	Fpb=cshift(Fp,-1,dim=2)
	Fpbb=cshift(Fpb,-1,dim=2)
	
	Fcip1=(2*Fpbb-7*Fpb+11*Fp)/6.0
	Fcip2=(-1*Fpb+5*Fp+2*Fpf)/6.0
	Fcip3=(2*Fp+5*Fpf-1*Fpff)/6.0
	!print *,Fcip1	
	betap1=(13.0/12.0)*(Fpbb-2*Fpb+Fp)*(Fpbb-2*Fpb+Fp)+(1.0/4.0)*(Fpbb-4*Fpb+3*Fp)*(Fpbb-4*Fpb+3*Fp)
	betap2=(13.0/12.0)*(Fpb-2*Fp+Fpf)*(Fpb-2*Fp+Fpf)+(1.0/4.0)*(Fpb-Fpf)*(Fpb-Fpf)
	betap3=(13.0/12.0)*(Fp-2*Fpf+Fpff)*(Fp-2*Fpf+Fpff)+(1.0/4.0)*(3*Fp-4*Fpf+Fpff)*(3*Fp-4*Fpf+Fpff)
	!print *,betap1
	!temp=reshape((/1,2,3,4,5,6,7,8,9/),(/3,3/))
	!print *,temp*temp
	!print *,temp+5
	!print *,1/(temp*temp)
	!print *,1/temp
	!print *,k1,k2,k3
	wwp1=kp1/((eps+betap1)*(eps+betap1))
	wwp2=kp2/((eps+betap2)*(eps+betap2))
	wwp3=kp3/((eps+betap3)*(eps+betap3))
        !print *,wwp1
	wp1=wwp1/(wwp1+wwp2+wwp3)
	wp2=wwp2/(wwp1+wwp2+wwp3)
	wp3=wwp3/(wwp1+wwp2+wwp3)
	!print *,wp1
	Fcip=wp1*Fcip1+wp2*Fcip2+wp3*Fcip3
	!print *,Fcip
	
	!!Negative flux
	Fnf=cshift(Fn,1,dim=2)
	Fnff=cshift(Fnf,1,dim=2)
	Fnb=cshift(Fn,-1,dim=2)
	Fnbb=cshift(Fnb,-1,dim=2)
	
	Fcin1=(-1*Fnbb+5*Fnb+2*Fn)/6.0
	Fcin2=(2*Fnb+5*Fn-1*Fnf)/6.0
	Fcin3=(11*Fn-7*Fnf+2*Fnff)/6.0
	
	betan1=(13.0/12.0)*(Fnbb-2*Fnb+Fn)*(Fnbb-2*Fnb+Fn)+(1.0/4.0)*(Fnbb-4*Fnb+3*Fn)*(Fnb-4*Fnb+3*Fn)
	betan2=(13.0/12.0)*(Fnb- 2*Fn+ Fnf)*(Fnb-2*Fn+Fnf)+(1.0/4.0)*(Fnb-Fnf)*(Fnb-Fnf)
	betan3=(13.0/12.0)*(Fn-  2*Fnf+Fnff)*(Fn-2*Fnf+Fnff)+(1.0/4.0)*(3*Fn-4*Fnf+Fnff)*(3*Fn-4*Fnf+Fnff)
	wwn1=kn1/((eps+betan1)*(eps+betan1))
	wwn2=kn2/((eps+betan2)*(eps+betan2))
	wwn3=kn3/((eps+betan3)*(eps+betan3))
	wn1=wwn1/(wwn1+wwn2+wwn3)
	wn2=wwn2/(wwn1+wwn2+wwn3)
	wn3=wwn3/(wwn1+wwn2+wwn3)
	Fcin=wn1*Fcin1+wn2*Fcin2+wn3*Fcin3
	WENO5LF1d=((Fcip-cshift(Fcip,-1,dim=2))+(Fcin-cshift(Fcin,-1,dim=2)))/dx
	!print *,WENO5LF1d
end Function WENO5LF1d
