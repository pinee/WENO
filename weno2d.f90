!Solving 1D-Sod's problem by finite volume approach
!Weno scheme is used for spatial evolution and 3rd order RK is used for temporal evolution 
program weno5
	INTERFACE flux
		Function flux (q,normal)
		double precision,intent(in) :: q(:,:,:)
		integer,intent(in)          :: normal(:)
		double precision flux (size(q,1),size(q,2),size(q,3))
		End Function
		
		Function WENO5LF2d (q,a,dx,dy)
		double precision,intent(in) :: q(:,:,:)
		double precision,intent(in) :: a
		double precision,intent(in) :: dx,dy
		double precision WENO5LF2d (size(q,1),size(q,2),size(q,3))
		End Function 
	End INTERFACE
	
	!Declaration of variables---------------------------------------------------------------------
	double precision                      :: dt,dt0,a0,a01,a02,a,a1,a2,CFL,tFinal,t=0,gama
	double precision                      :: dx,dy,x_middle,y_middle,xstart,xend,ystart,yend
	integer                               :: it=0,i,j
	integer, parameter                    :: nx=51,ny=51
	double precision, dimension(nx,ny)    :: rho0,p0,u0,v0,vn0,E0,c0,rho,u,E,c,p,v,vn,lambda01,lambda02,lambda1,lambda2
	double precision, dimension(nx)       :: x
	double precision, dimension(ny)       :: y
	double precision, dimension(nx,ny)    :: reg1=0.0,reg2=0.0,reg3=0.0,reg4=0.0 
	double precision, dimension(nx,ny,4)  :: q0,q,F,dF
	double precision, dimension(4)        :: pi,ui,rhoi,vi
	real                                  :: cfl1,tEnd
        character(len=1024)                   :: filename,ff,ff2
	!Assigning values-----------------------------------------------------------------------------
	CFL=0.475
	tFinal=0.2
	gama=1.4
	xstart=0
	xend=1
	ystart=0
        yend=1
		
	dx=(xend-xstart)/nx
	dy=(yend-ystart)/ny
	x(1)=xstart
	do i=2,nx
	x(i)=x(i-1)+dx
	end do
	y(1)=ystart
	do i=2,ny
	y(i)=y(i-1)+dy
	end do
	!Defining initial conditions-------------------------------------------------------------------
	x_middle=(xend+xstart)/2
	y_middle=(yend+ystart)/2
	!print *," x_middle=",x_middle
	!Configuration 5
	rhoi=(/1.0,0.5197,0.1072,0.2579/)
	pi=(/1.0,0.4,0.0439,0.15/)
	ui=(/0.0,-0.7259,-0.7259,0.0/)
	vi=(/0.0,0.0,-1.4045,-1.4045/)
		
	do i=1,nx
	do j=1,ny	
	if (x(i)>=x_middle) then
		if(y(j)>=y_middle) then			
                        reg1(i,j)=1.0
	        else 		
			reg4(i,j)=1.0
		end if
	else 
		if(y(j)>=y_middle) then
			reg2(i,j)=1.0
	        else
			reg3(i,j)=1.0
		end if
	end if
	rho0(i,j)=reg1(i,j)*rhoi(1)+reg2(i,j)*rhoi(2)+reg3(i,j)*rhoi(3)+reg4(i,j)*rhoi(4)
	p0(i,j)=reg1(i,j)*pi(1)+reg2(i,j)*pi(2)+reg3(i,j)*pi(3)+reg4(i,j)*pi(4)
	u0(i,j)=reg1(i,j)*ui(1)+reg2(i,j)*ui(2)+reg3(i,j)*ui(3)+reg4(i,j)*ui(4)
	v0(i,j)=reg1(i,j)*vi(1)+reg2(i,j)*vi(2)+reg3(i,j)*vi(3)+reg4(i,j)*vi(4)
	end do
	end do


       write (filename, "(A,I4)") "datarho",it+1000 
	OPEN(UNIT=1000,FILE=filename)
        DO i=1,nx
	  DO j=1,ny
          WRITE(1000,*) x(i),y(j),rho0(i,j)
          ENDDO
        ENDDO
       CLOSE(1000)


	do i=1,nx
		do j=1,ny
			E0(i,j)=p0(i,j)/(rho0(i,j)*(gama-1))+0.5*((u0(i,j)*u0(i,j))+(v0(i,j)*v0(i,j)))
			!print *,i,E0(i,j) 
			c0(i,j)=SQRT(gama*p0(i,j)/rho0(i,j))
			q0(i,j,1)=rho0(i,j)
			q0(i,j,2)=rho0(i,j)*u0(i,j)
			q0(i,j,3)=rho0(i,j)*v0(i,j)
			q0(i,j,4)=rho0(i,j)*E0(i,j)
		end do
	end do
	!F=flux(q0)
	!print *, F
	!print *,"a0(1)=", a0(1),"E0(1)=",E0(1)
	!print *,"u0(1)=",u0(1)
	!print *,"F0(1,1)=",F(1,1)
	!print *,"F0(2,1)=",F(2,1)
	!print *,"F0(3,1)=",F(3,1)
	vn0=sqrt((u0*u0)+(v0*v0))
	lambda01=vn0+c0
	lambda02=vn0-c0
	
	a01=MAXVAL(ABS(lambda01))
	a02=MAXVAL(ABS(lambda02))
	a0=max(a01,a02)
	
!	print *,'Value of a0 is',a0
	dt0=CFL*min(dx/a0,dy/a0)
	!Solver loop-------------------------------------------------------------------------------------
	!Load initial conditions
	q=q0
	dt=dt0
	a=a0
        !dF=WENO5LF2d(q,a,dx,dy) 
        !print *,rho0
	!print *,dF
	do while(t<tFinal)
!		!RK Initial step
!		!1st stage
                       q0=q
!		print*,"q before is"
!		do i=1,nx
!		print*,"i=",i
!		print*,"-----------------------------------------------"
!		do j=1,ny
!		write(*,*),q(i,j,:)
!		end do
!		end do
!		print*,"-----------------------------------------------"
                        dF=WENO5LF2d(q,a,dx,dy) 
			q = q0-dt*dF
!		print*,"The value for dt is -------------------"
!		write(*,*),dt
!		print*,"q after is"
!		do i=1,nx
!		print*,"i=",i
!		print*," -----------------------------------------------"
!		do j=1,ny
!		write(*,*),q(i,j,:)
!		end do
!		end do
!		print*,"-----------------------------------------------"
			do i=1,nx
			do k=1,4
				q(i,1,k)=q0(i,1,k)
				q(i,ny,k)=q0(i,ny,k)
			end do
			end do
			do j=1,ny
			do k=1,4
				q(1,j,k)=q0(1,j,k)
				q(nx,j,k)=q0(nx,j,k)
			end do
			end do

		!RK second stage
                       dF=WENO5LF2d(q,a,dx,dy) 
			q =0.75*q0+0.25*(q-dt*dF)
			do i=1,nx
			do k=1,4
				q(i,1,k)=q0(i,1,k)
				q(i,ny,k)=q0(i,ny,k)
			end do
			end do
			do j=1,ny
			do k=1,4
				q(1,j,k)=q0(1,j,k)
				q(nx,j,k)=q0(nx,j,k)
			end do
			end do
		!Third stage
                       dF=WENO5LF2d(q,a,dx,dy) 
			q =(q0+2*(q-dt*dF))/3.0
			do i=1,nx
			do k=1,4
				q(i,1,k)=q0(i,1,k)
				q(i,ny,k)=q0(i,ny,k)
			end do
			end do
			do j=1,ny
			do k=1,4
				q(1,j,k)=q0(1,j,k)
				q(nx,j,k)=q0(nx,j,k)
			end do
			end do
		!Computing primary properties
			do i=1,nx
			do j=1,ny
				rho(i,j)=q(i,j,1)
				u(i,j)=q(i,j,2)/rho(i,j)
				v(i,j)=q(i,j,3)/rho(i,j)
				E(i,j)=q(i,j,4)/rho(i,j)
				p(i,j)=(gama-1)*rho(i,j)*(E(i,j)-0.5*(u(i,j)*u(i,j))+(v(i,j)*v(i,j)))
				c(i,j)=SQRT(gama*p(i,j)/rho(i,j))
			end do
			end do
		!Update dt and lambda
                	vn=sqrt((u*u)+(v*v))
			lambda1=vn+c
			lambda2=vn-c
			a1=maxval(abs(lambda1))
		        a2=maxval(abs(lambda2))
			a=max(a1,a2)
			dt=CFL*min(dx/a,dy/a)
                	if(t+dt>tFinal)then
				dt=tFinal-t
                	end if


	!Update time and iteration count
			t=t+dt
			it=it+1
		!Plot------------------------------------------------------------------------------------------
!
!       write (filename, "(A,I4)") "datarho",it+1000 
!	OPEN(UNIT=1000,FILE=filename)
!        DO i=1,nx
!	  DO j=1,ny
!          WRITE(1000,*) x(i),y(j),rho(i,j)
!          ENDDO
!        ENDDO
!       CLOSE(1000)
!	!write(ff,*) "filename=",filename
!	!write(ff2,*) "gnuplot -e ",ff, " -p drawrho_plot.plt"
!	end do

!	print *,"it=",it,"t=",t 
!
!	!print *,u
	end do

!	do i=1,nx
!	write(*,*),rho(i,:)
!	end do
!        CALL SYSTEM("gnuplot -p drawrho_plotmulti.plt")        
!        OPEN(UNIT=1000,FILE='datarho')
!        DO i=1,nx
! 	DO j=1,ny      
!	WRITE(1000,*) x(i),y(j),rho(i,j)
!	ENDDO
!        ENDDO
!        CLOSE(1000)
!       CALL SYSTEM('gnuplot -p drawrho_plot.plt')       
!       OPEN(UNIT=1000,FILE='datap')
!       DO i=1,n
!       WRITE(1000,*) x(i),p(i)
!       ENDDO
!       CLOSE(1000)
!       CALL SYSTEM('gnuplot -p datap_plot.plt')       

 

end program weno5
Function flux(q,normal)
        integer, parameter                        :: nx=51,ny=51
	integer, intent(in)                       :: normal(:)
        double precision, intent(in)              :: q(:,:,:)
        double precision flux(size(q,1),size(q,2),size(q,3))
        double precision, dimension(nx,ny)        :: rho,E,p,u,v,vm
        real                                      :: gama=1.4
        integer                                   :: i
        nmx=normal(1)
        nmy=normal(2)
        do i=1,nx
	do j=1,ny
       	rho(i,j)=q(i,j,1)
       	u(i,j)=q(i,j,2)/rho(i,j)
       	v(i,j)=q(i,j,3)/rho(i,j)
       	E(i,j)=q(i,j,4)/rho(i,j)
	vm(i,j)=(u(i,j)*nmx)+(v(i,j)*nmy)
       	p(i,j)=(gama-1)*rho(i,j)*(E(i,j)-0.5*((u(i,j)*u(i,j))+(v(i,j)*v(i,j))))
       	flux(i,j,1)=rho(i,j)*vm(i,j)
       	flux(i,j,2)=rho(i,j)*vm(i,j)*u(i,j)+p(i,j)*nmx
       	flux(i,j,3)=rho(i,j)*vm(i,j)*v(i,j)+p(i,j)*nmy
       	flux(i,j,4)=vm(i,j)*(rho(i,j)*E(i,j)+p(i,j))
        end do
	end do
end Function

Function WENO5LF2d(q,a,dx,dy)
        INTERFACE flux
                Function flux (q,normal)
                double precision,intent(in)  :: q(:,:,:)
                integer,intent(in)           :: normal(:)
		double precision flux (size(q,1),size(q,2),size(q,3))
                End Function
	end INTERFACE
        integer, parameter                   :: nx=51,ny=51
	double precision, intent(in)         :: q(:,:,:)
	double precision, intent(in)         :: a,dx,dy
	integer, dimension(2)                :: normal
	double precision                     :: alpha,gama=1.4,eps=1E-006,kp1=0.1,kp2=0.6,kp3=0.3,kn1=0.3,kn2=0.6,kn3=0.1
	double precision, dimension(nx,ny,4) :: F,Fp,Fn,Fcin1,Fcin2,Fcin3,Fcin,Fcip1,Fcip2,Fcip3,Fcip,Fpf,Fpff,Fpb,Fpbb,Fnf,Fnff,Fnbb,Fnb
	double precision, dimension(nx,ny,4) :: G,Gp,Gn,Gcin1,Gcin2,Gcin3,Gcin,Gcip1,Gcip2,Gcip3,Gcip,Gpf,Gpff,Gpb,Gpbb,Gnf,Gnff,Gnbb,Gnb
        double precision, dimension(nx,ny,4) :: wn1,wn2,wn3,wwn1,wwn2,wwn3,betan1,betan2,betan3 
	double precision, dimension(nx,ny,4) :: wp1,wp2,wp3,wwp1,wwp2,wwp3,betap1,betap2,betap3
	double precision WENO5LF2d(size(q,1),size(q,2),size(q,3))
        !print *,"I entered WENO5LF2d--------------------------------------------------------------"
        !print *,"epsilon=",eps
	normal=(/ 1,0 /)
	F=flux(q,normal)
	!print *,F
        !print *,F
	Fp=0.5*(F+a*q)
	Fn=0.5*(F-a*q)
	!print *,"Fp and Fn generated"
	Fn=cshift(Fn,1,dim=2)
 	!print *,Fn 
	!Positive flux
	Fpf=cshift(Fp,1,dim=2)
	Fpff=cshift(Fpf,1,dim=2)
	Fpb=cshift(Fp,-1,dim=2)
	Fpbb=cshift(Fpb,-1,dim=2)
	
	Fcip1=(2*Fpbb-7*Fpb+11*Fp)/6.0
	!print *,"--------------------------------------------------"
	!print *,"Fcip1 is",Fcip1
	Fcip2=(-1*Fpb+5*Fp+2*Fpf)/6.0
	Fcip3=(2*Fp+5*Fpf-1*Fpff)/6.0
	!print *,Fcip1	
	betap1=(13.0/12.0)*(Fpbb-2*Fpb+Fp)*(Fpbb-2*Fpb+Fp)+(1.0/4.0)*(Fpbb-4*Fpb+3*Fp)*(Fpbb-4*Fpb+3*Fp)
	betap2=(13.0/12.0)*(Fpb-2*Fp+Fpf)*(Fpb-2*Fp+Fpf)+(1.0/4.0)*(Fpb-Fpf)*(Fpb-Fpf)
	betap3=(13.0/12.0)*(Fp-2*Fpf+Fpff)*(Fp-2*Fpf+Fpff)+(1.0/4.0)*(3*Fp-4*Fpf+Fpff)*(3*Fp-4*Fpf+Fpff)
	!print *,"--------------------------------------------------"
	!print *,betap1
	wwp1=kp1/((eps+betap1)*(eps+betap1))
	wwp2=kp2/((eps+betap2)*(eps+betap2))
	wwp3=kp3/((eps+betap3)*(eps+betap3))
	!print *,"--------------------------------------------------"
        !print *,wwp1
	wp1=wwp1/(wwp1+wwp2+wwp3)
	wp2=wwp2/(wwp1+wwp2+wwp3)
	wp3=wwp3/(wwp1+wwp2+wwp3)
	!print *,"--------------------------------------------------"
	!print *,wp1
	Fcip=wp1*Fcip1+wp2*Fcip2+wp3*Fcip3
	!print *,"--------------------------------------------------"
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
	WENO5LF2d=((Fcip-cshift(Fcip,-1,dim=2))+(Fcin-cshift(Fcin,-1,dim=2)))/dx
	!print *,"-------------------------------------------------------"
	!print *,WENO5LF2d

!Right extrapolation
	normal=(/ 0,1 /)
	G=flux(q,normal)
        !print *,G
	Gp=0.5*(G+a*q)
	Gn=0.5*(G-a*q)
	!print *,"Gp and Gn generated"
	Gn=cshift(Gn,1,dim=1)
 	!print *,Gn 
	!Positive flux
	Gpf=cshift(Gp,1,dim=1)
	Gpff=cshift(Gpf,1,dim=1)
	Gpb=cshift(Gp,-1,dim=1)
	Gpbb=cshift(Gpb,-1,dim=1)
	
	Gcip1=(2*Gpbb-7*Gpb+11*Gp)/6.0
	Gcip2=(-1*Gpb+5*Gp+2*Gpf)/6.0
	Gcip3=(2*Gp+5*Gpf-1*Gpff)/6.0
	!print *,Gcip1	
	betap1=(13.0/12.0)*(Gpbb-2*Gpb+Gp)*(Gpbb-2*Gpb+Gp)+(1.0/4.0)*(Gpbb-4*Gpb+3*Gp)*(Gpbb-4*Gpb+3*Gp)
	betap2=(13.0/12.0)*(Gpb-2*Gp+Gpf)*(Gpb-2*Gp+Gpf)+(1.0/4.0)*(Gpb-Gpf)*(Gpb-Gpf)
	betap3=(13.0/12.0)*(Gp-2*Gpf+Gpff)*(Gp-2*Gpf+Gpff)+(1.0/4.0)*(3*Gp-4*Gpf+Gpff)*(3*Gp-4*Gpf+Gpff)
	!print *,betap1
	wwp1=kp1/((eps+betap1)*(eps+betap1))
	wwp2=kp2/((eps+betap2)*(eps+betap2))
	wwp3=kp3/((eps+betap3)*(eps+betap3))
        !print *,wwp1
	wp1=wwp1/(wwp1+wwp2+wwp3)
	wp2=wwp2/(wwp1+wwp2+wwp3)
	wp3=wwp3/(wwp1+wwp2+wwp3)
	!print *,wp1
	Gcip=wp1*Gcip1+wp2*Gcip2+wp3*Gcip3
	!print *,Gcip
	
	!!Negative flux
	Gnf=cshift(Gn,1,dim=1)
	Gnff=cshift(Gnf,1,dim=1)
	Gnb=cshift(Gn,-1,dim=1)
	Gnbb=cshift(Gnb,-1,dim=1)
	
	Gcin1=(-1*Gnbb+5*Gnb+2*Gn)/6.0
	Gcin2=(2*Gnb+5*Gn-1*Gnf)/6.0
	Gcin3=(11*Gn-7*Gnf+2*Gnff)/6.0
	
	betan1=(13.0/12.0)*(Gnbb-2*Gnb+Gn)*(Gnbb-2*Gnb+Gn)+(1.0/4.0)*(Gnbb-4*Gnb+3*Gn)*(Gnb-4*Gnb+3*Gn)
	betan2=(13.0/12.0)*(Gnb- 2*Gn+ Gnf)*(Gnb-2*Gn+Gnf)+(1.0/4.0)*(Gnb-Gnf)*(Gnb-Gnf)
	betan3=(13.0/12.0)*(Gn-  2*Gnf+Gnff)*(Gn-2*Gnf+Gnff)+(1.0/4.0)*(3*Gn-4*Gnf+Gnff)*(3*Gn-4*Gnf+Gnff)
	wwn1=kn1/((eps+betan1)*(eps+betan1))
	wwn2=kn2/((eps+betan2)*(eps+betan2))
	wwn3=kn3/((eps+betan3)*(eps+betan3))
	wn1=wwn1/(wwn1+wwn2+wwn3)
	wn2=wwn2/(wwn1+wwn2+wwn3)
	wn3=wwn3/(wwn1+wwn2+wwn3)
	Gcin=wn1*Gcin1+wn2*Gcin2+wn3*Gcin3
	WENO5LF2d=WENO5LF2d+((Gcip-cshift(Gcip,-1,dim=1))+(Gcin-cshift(Gcin,-1,dim=1)))/dy
!		print*,"dF is"
!		do i=1,nx
!		print*,"i=",i
!		print*,"-----------------------------------------------"
!		do j=1,ny
!		write(*,*),WENO5LF2d(i,j,:)
!		end do
!		end do
!		print*,"-----------------------------------------------"
!	print *,WENO5LF2d
end Function WENO5LF2d
