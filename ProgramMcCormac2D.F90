!**********************************************************
!                                                         !
!	ПРОГРАММА ДЛЯ РАСЧЕТА ПЛОСКИХ ВЯЗКИХ ТЕЧЕНИЙ ПО СХЕМЕ !
!	    МАК-КОРМАКА 2-ГО ПОРЯДКА ТОЧНОСТИ                 !
!            Ver. 1.0.1 (28.01.2012)                        !
!**********************************************************

program McCormac2D
implicit none

real,parameter :: N=25, M=25 ! Расчетная сетка
real,parameter :: tmax = 300 ! количество шагов по времени

real rho(0:N,0:M), u(0:N,0:M), v(0:N,0:M), E(0:N,0:M), p(0:N,0:M), T(0:N,0:M)
real txx(0:N,0:M), txy(0:N,0:M), tyx(0:N,0:M), tyy(0:N,0:M)
real Qdx(0:N,0:M), Qdy(0:N,0:M)
real q(4,0:N,0:M), F(4,0:N,0:M), G(4,0:N,0:M), qp(4,0:N,0:M)
real hx(0:N), x(0:N), hy(0:M), y(0:M), hxmin, hymin
real Nstep, time, timeshow, taup, tauc, tstep
real Lx, Ly, Lh, h0
real gamma, R, pij

real mu, muT, cp, kte, Pr, PrT, k  ! параметры связанные с вязкостью и теплопроводностью

integer i,j,l ! целочисленные индексы
real rho00, p00, M00, u00, v00, c00, E00, T00  ! параметры набегающего потока
real rho0, p0, M0, u0, v0, c0, E0, T0, nr, tpulse ! параметры вдуваемой струи
real Ku, DBL, cv, rcv

!--------------------------------------!
! ПАРАМЕТРЫ НАБЕГАЮЩЕГО ПОТОКА И СТРУИ !
!--------------------------------------!
! Общие параметры газа
R = 287.0; gamma = 1.4; mu = 18.27*0.001; k = 0.02

! 1. Параметры набегающего потока
p00=0.0663*10**6; M00=2.9; v00=0; T00 = 108; 
rho00 = p00/(T00*R); c00=sqrt(gamma*p00/rho00); u00= M00*c00; E00=p00/(gamma-1) + 0.5*rho00*(u00*u00 + v00*v00)
rcv = gamma*(gamma - 1)*M00*M00; cv = 1/rcv 

! 2. Параметры вдуваемой струи
nr=18.7; M0=1.0; u0 = 0.0; T0 = 217; 
rho0 = p0/(T0*R); p0=p00*nr; c0=sqrt(gamma*p0/rho0); v0=M0*c0; E0=p0/(gamma-1) + 0.5*rho0*(u0*u0 + v0*v0)
tpulse=0.01

! 3. Параметры турбулентнсти
muT=1.0e-10; cp = 1.0e-10; kte = 1.0e-10; PrT = 0.9
 
! 4. Расчет шагов по сетке
Lx=15*0.01; Ly=7.62*0.01; Lh=2*Lx/3; h0=0.1*Lx; DBL = 0.1*Ly

x(0) = 0
do i=1, N-1
    hx(i) = Lx/(N)
    x(i+1) = x(i) + hx(i)
end do
hxmin = Lx/N

y(0) = 0
do j=1, M-1
    hy(j)=Ly/(M)
    y(j+1) = y(j) + hy(j)
end do
hymin = Ly/M

! 5. Параметры схемы
Ku = 0.01

! Очистка файлов для записи (см. ниже)
open(1,FILE='p.txt');   close(1, status='DELETE')
open(2,FILE='v.txt');   close(2, status='DELETE')
open(3,FILE='M.txt');   close(3, status='DELETE')

call startwrite()

!============================================================================================================
!	1. НАЧАЛЬНЫЕ УСЛОВИЯ (задаются через массивы РЕАЛЬНЫХ переменных - давления, плотности и т.п.)
!============================================================================================================

! Вектор q описывает поток:
! q1 = rho, q2 = rho*u, q3 = rho*v, q4 = E
! где (u,v,w) - компоненты вектора скорости.

do i=1,N-1
	do j=1,M-1
!		if (y(j).le.DBL) then
!			q(1,i,j) = rho00
!!			q(2,i,j) = rho00*u00*(y(j))**2/(DBL**2)
!			q(2,i,j) = rho00*u00
!			q(3,i,j) = rho00*v00
!			q(4,i,j) = p00/(gamma-1) + 0.5*rho00*((q(2,i,j)/rho00)**2 + v00*v00)		
!		else
		    q(1,i,j) = rho00
			q(2,i,j) = rho00*u00
			q(3,i,j) = rho00*v00
			q(4,i,j) = E00
!		end if
	end do
end do

call boundary(q) ! вызов процедуры, применяющей граничные условия

!============================================================================================================
!	2. Основной вычислительный цикл прогрммы (метод Мак-Кормака)
!============================================================================================================
time = 0        ! обнуление реального времени моделирования
timeshow = 0    ! обнулние счетчика кадров
Nstep = 0       ! обнуление счетчика шагов по времени

do while (Nstep.lt.tmax)

!==========================!
!   2.1. ШАГ "ПРЕДИКТОР"   !
!==========================!
    
    taup=tau(q) ! вычисление шага по времени
    
    ! Вычисление исходных функций по вектору q на всей сетке
    do i=1,N-1
        do j=1,M-1
            rho(i,j) = q(1,i,j)
            u(i,j) = q(2,i,j)/rho(i,j)
            v(i,j) = q(3,i,j)/rho(i,j)
            E(i,j) = q(4,i,j)
            p(i,j) = (gamma-1)*(E(i,j)-0.5*rho(i,j)*(u(i,j)*u(i,j)+v(i,j)*v(i,j)))
            T(i,j) = p(i,j)/(R*rho(i,j))
        enddo
    enddo

    ! вычисление напряжений, тепловых потоков и функций F и G на шаге "предиктор"
    do i=1,N-1
        do j=1,M-1
            ! Вычисление напряжений на шаге "предиктор"        
            txx(i,j) = (4/3)*(mu+muT)*(u(i+1,j)-u(i,j))/hx(i) - (2/3)*(mu+muT)*0.5*(v(i,j+1)-v(i,j-1))/hy(j) - (2/3)*rho(i,j)*kte
            tyy(i,j) = (4/3)*(mu+muT)*(v(i,j+1)-v(i,j))/hy(j) - (2/3)*(mu+muT)*0.5*(u(i+1,j)-u(i-1,j))/hx(i) - (2/3)*rho(i,j)*kte  ! ЧТО ТАКОЕ kte?
            txy(i,j) = (mu+muT)*((u(i+1,j)-u(i,j))/hx(i) + 0.5*(v(i,j+1)-v(i,j-1))/hy(j))
            tyx(i,j) = (mu+muT)*(0.5*(u(i+1,j)-u(i-1,j))/hx(i) + (v(i,j+1)-v(i,j))/hy(j))
            
            ! Вычисление тепловых потоков на шаге "предиктор"
            Qdx(i,j) = -(k + cp*muT/PrT)*(T(i+1,j) - T(i,j))/hx(i)
            Qdy(i,j) = -(k + cp*muT/PrT)*(T(i,j+1) - T(i,j))/hy(j)
            
            ! Вычисление функции F на шаге "предиктор"
            F(1,i,j) = rho(i,j)*u(i,j)
            F(2,i,j) = rho(i,j)*u(i,j)*u(i,j) + p(i,j) - txx(i,j)
            F(3,i,j) = rho(i,j)*u(i,j)*v(i,j) - txy(i,j)
            F(4,i,j) = (E(i,j)+p(i,j)-txx(i,j))*u(i,j) - txy(i,j)*v(i,j) + Qdx(i,j)
            
            ! Вычисление функции G на шаге "предиктор"
            G(1,i,j) = rho(i,j)*v(i,j)
            G(2,i,j) = rho(i,j)*u(i,j)*v(i,j) - tyx(i,j)
            G(3,i,j) = rho(i,j)*v(i,j)*v(i,j) + p(i,j) - tyy(i,j)
            G(4,i,j) = (E(i,j) + p(i,j) - tyy(i,j))*v(i,j) - tyx(i,j)*u(i,j) + Qdy(i,j)
        enddo
    enddo

    ! Вычисление вектора q для перехода от предиктора к корректору
    do i=1,N-1
        do j=1,M-1
            do l=1,4
                qp(l,i,j) = q(l,i,j) - (taup/hx(i))*(F(l,i,j) - F(l,i-1,j))-(taup/hy(j))*(G(l,i,j) - G(l,i,j-1))
            end do
        end do
    end do
    
    call boundary(qp)
    
    !====================!
    !   ШАГ "КОРРЕКТОР   !
    !====================!
    
    tauc=tau(qp)
    
    ! Вычисление исходных функций по вектору q на всей сетке
    do i=1,N-1
        do j=1,M-1
            rho(i,j) = qp(1,i,j)
            u(i,j) = qp(2,i,j)/rho(i,j)
            v(i,j) = qp(3,i,j)/rho(i,j)
            E(i,j) = qp(4,i,j)
            p(i,j) = (gamma-1)*(E(i,j)-0.5*rho(i,j)*(u(i,j)*u(i,j)+v(i,j)*v(i,j)))
            T(i,j) = p(i,j)/(R*rho(i,j))
        enddo
    enddo

    ! Вычисление напряжений, тепловых потоков и функций F и G на шаге "корректор"
    do i=1,N-1
        do j=1,M-1
            ! Вычисление напряжений на шаге "корректор"        
            txx(i,j) = (4/3)*(mu+muT)*(u(i,j)-u(i-1,j))/hx(i) - (2/3)*(mu+muT)*0.5*(v(i,j+1)-v(i,j-1))/hy(j) - (2/3)*rho(i,j)*kte 
            tyy(i,j) = (4/3)*(mu+muT)*(v(i,j)-v(i,j-1))/hy(j) - (2/3)*(mu+muT)*0.5*(u(i+1,j)-u(i-1,j))/hx(i) - (2/3)*rho(i,j)*kte 
            txy(i,j) = (mu+muT)*((u(i,j)-u(i-1,j))/hx(i) + 0.5*(v(i,j+1)-v(i,j-1))/hy(j))
            tyx(i,j) = (mu+muT)*(0.5*(u(i+1,j)-u(i-1,j))/hx(i) + (v(i,j)-v(i,j-1))/hy(j))
            
            ! Вычисление тепловых потоков на шаге "корректор"
            Qdx(i,j) = -(k + cp*muT/PrT)*(T(i,j) - T(i-1,j))/hx(i)
            Qdy(i,j) = -(k + cp*muT/PrT)*(T(i,j) - T(i,j-1))/hy(j)
            
            ! Вычисление функции F на шаге "корректор"
            F(1,i,j) = rho(i,j)*u(i,j)
            F(2,i,j) = rho(i,j)*u(i,j)*u(i,j) + p(i,j) - txx(i,j)
            F(3,i,j) = rho(i,j)*u(i,j)*v(i,j) - txy(i,j)
            F(4,i,j) = (E(i,j)+p(i,j)-txx(i,j))*u(i,j) - txy(i,j)*v(i,j) + Qdx(i,j)
            
            ! Вычисление функции G на шаге "корректор"
            G(1,i,j) = rho(i,j)*v(i,j)
            G(2,i,j) = rho(i,j)*u(i,j)*v(i,j) - tyx(i,j)
            G(3,i,j) = rho(i,j)*v(i,j)*v(i,j) + p(i,j) - tyy(i,j)
            G(4,i,j) = (E(i,j) + p(i,j) - tyy(i,j))*v(i,j) - tyx(i,j)*u(i,j) + Qdy(i,j)
        enddo
    enddo

    ! Вычисление вектора q для перехода на следующий временной слой
    do i=1,N-1
        do j=1,M-1
            do l=1,4
                q(l,i,j) = 0.5*(q(l,i,j) + qp(l,i,j)) - 0.5*(tauc/hx(i))*(F(l,i+1,j) - F(l,i,j))-0.5*(tauc/hy(j))*(G(l,i,j+1) - G(l,i,j))
            end do
        end do
    end do
    
    call boundary(q)
    
    time=time+taup+tauc
    
    ! Процедура вывода результатов в случае нужного шага
	if((timeshow-Nstep).eq.0) then 
		if (Nstep.lt.100) then ; timeshow = timeshow + 1
		else if (tstep.lt.1000) then ; timeshow = timeshow + 25
		else ; timeshow = timeshow + 50
		end if
		open(4,FILE='M.txt',access='APPEND')
		open(5,FILE='p.txt',access='APPEND')      ! открытие файла для записи
		open(6,FILE='v.txt',access='APPEND')
		do i = 1,N-1
			do j = 1,M-1
			    pij = (gamma-1)*(q(4,i,j)-0.5*q(1,i,j)*((q(2,i,j)/q(1,i,j))**2 + (q(3,i,j)/q(1,i,j))**2))
		!	   write(4,'(1x,f26.6, f26.3,f26.3,f36.3)') time, x(i),y(j), sqrt(((q(2,n,m)/q(1,n,m))**2+(q(3,n,m)/q(1,n,m)))/(gamma*pij))
			!   write(5,'(1x,f26.6, f26.3,f26.3,f36.3)') time, x(i),y(j), pij
			write(5,'(1x, f7.0, i3, i3,f36.3)') Nstep, i, j, pij
		!	   write(6,'(1x,f26.6, f26.3,f26.3,f36.3,f36.3)') time, x(i),y(j), q(2,n,m)/q(1,n,m), q(3,n,m)/q(1,n,m)
			end do
		end do
		close(4); close(5); close(6)
	end if    


    Nstep = Nstep + 1
end do

!------------------------------------------------------
!------------------------------------------------------
!--------------- ПРОЦЕДУРЫ ----------------------------
!------------------------------------------------------
!------------------------------------------------------

contains

!=================================!
!	6. Вычисление шага по времени !
!=================================!

real function tau(q)
	real vel,vmax, pij, q(4,0:N,0:M)
	vmax = 0
	do i = 1,N-1
		do j = 1,M-1		    
		    pij = (gamma-1)*(q(4,i,j)-0.5*q(1,i,j)*((q(2,i,j)/q(1,i,j))**2 + (q(3,i,j)/q(1,i,j))**2))
			vel=sqrt(gamma*pij/q(1,i,j))+sqrt((q(2,i,j)/q(1,i,j))**2+(q(3,i,j)/q(1,i,j))**2)
			if (vel.gt.vmax) then ; vmax = vel ; end if
		end do
	end do

	tau=sqrt(0.5-mu)*Ku*min(hxmin,hymin)/vmax
	tau = 0.0000001
end function	
!========================!
!	7.Граничные условия  !
!========================!

subroutine boundary(q)
    real q(4,0:N,0:M)
    do i=0,N    ! по координате x
        if ((x(i).lt.Lh).or.(x(i).gt.(Lh+h0))) then
            q(1,i,0)= q(1,i,1)
            q(2,i,0)=-q(2,i,1)
            q(3,i,0)=-q(3,i,1)
            q(4,i,0)= q(4,i,1)
        else !if (time.ge.tpulse) then
            q(1,i,0) = rho0
            q(2,i,0) = rho0*u0
            q(3,i,0) = rho0*v0
            q(4,i,0) = E0
!		else
!		     q(1,i,0)= q(1,i,1)
!            q(2,i,0)=-q(2,i,1)
!            q(3,i,0)=-q(3,i,1)
!            q(4,i,0)= q(4,i,1)
        end if

		! условие свободной границы сверху
		do l=1,4
			q(l,i,M)=q(l,i,M-1)
		end do        
    end do
    
    do j=0,M    ! по координате y
!		if (y(j).le.DBL) then
		! степенной закон в пограничном слое
!			q(1,0,j) = rho00
	!		q(2,0,j) = rho00*u00*(y(j))**2/(DBL**2)
!			q(2,0,j) = rho00*u00
!			q(3,0,j) = rho00*v00
!			q(4,0,j) = p00/(gamma-1) + 0.5*rho00*((q(2,0,j)/rho00)**2 + v00*v00)		
!		else
			! параметры набегающего потока
!			q(1,0,j) = rho00
	!		q(2,0,j) = rho00*u00*(y(j))**2/(DBL**2)
!	        q(2,0,j) = rho00*u00
!			q(3,0,j) = rho00*v00
!			q(4,0,j) = p00/(gamma-1) + 0.5*rho00*((q(2,0,j)/rho00)**2 + v00*v00)		
!		end if
		! условие свободной границы слева
		do l=1,4
			q(l,0,j) = q(l,1,j)
		end do

		! условие свободной границы справа
		do l=1,4
			q(l,N,j) = q(l,N-1,j)
		end do
    end do
end subroutine

!==============================================!
!	8. Процедуры для вывода на экран сообщений !
!==============================================!

subroutine startwrite()
	write(*,*) "========================================="
	write(*,*) "McCormac scheme for 2D gasdynamics"
	write(*,*) "========================================="
	write(*,*)
	write(*,*) "Calculation is started..."
	write(*,*)
end subroutine

subroutine endwrite()
	write(*,*)
	write(*,*) "Calculation if finished."
	write(*,*)
end subroutine

end program