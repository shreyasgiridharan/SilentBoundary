    module modulefem
    !*********************************************************************
    !    Function:
    !**********************************************************************


    implicit none

    logical, dimension(:), allocatable :: IsFixDof
    integer :: ndtn = 0, ITime = 0
    integer :: NDivX = 0, NDivY = 0
    integer :: NEl = 0, NNod = 0
    integer, dimension(:,:), allocatable :: ICon
    double precision :: gx = 0.D0, gy = 0.D0
    double precision :: Meshdx = 0.D0, Meshdy = 0.D0
    double precision :: Maxdx = 0.D0, Maxdy = 0.D0
    double precision :: eledx = 0.D0, eledy = 0.D0
    double precision, dimension(:), allocatable :: GrvF, InrF, Mas, Area, Areai, ExtF, v0, v, dis, DamF
    double precision, dimension(:,:), allocatable :: NodCo
    double precision, dimension(:,:,:), allocatable :: SigG, EpsG, F, HS, B_trial, EpsP, EpsE, edP, ETA, Sig0, &
                                                       Statevar, Swp, Dswp, Plastind, FBARNEW
    double precision, dimension(:,:,:,:), allocatable :: B
    double precision :: MatProp(32), dt = 0.D0, FBAR, epsv, delta = 0.d0
    double precision :: BulkW = 2200000.d0 ! Bulk Modulus of Water 2200000 kN/m2
    
    integer :: nplast = 0.d0
    double precision :: PHI, PSI, COH
        
    integer, parameter :: DATUnit = 2
    integer, parameter :: LOGUnit = 2
    integer, parameter :: MSHUnit = 3
    integer, parameter :: RESUnit = 4

    contains

    subroutine solve()
    !**********************************************************************
    !    Function:
    !**********************************************************************

    implicit none
    integer :: IIter, IStep = 0, IPTime, dPsi
    integer :: iprint, imeshwrite
    real :: start, finish
    iprint = 1000.0
    imeshwrite = 100.0
    open(LOGUnit, file = 'output.log')
    call initial()
    call WtMSH(MSHUnit)

    call CPU_TIME(start)
    do ITime = 1, ndtn ! physical time scale
        if (0.eq.1) then
            call Map2nod()
            call update()
        else
            call Map2nod()
            call update_small()
        endif

        !if (itime.eq.1 .or. int(ITime/5000.0) .eq. (ITime/5000.0)) then
        if(itime.eq.1 .or. (itime/imeshwrite)*imeshwrite == itime .or. itime == ndtn) then
            IStep = IStep + 1
            call WtRES(IStep)
        endif
        if(itime.eq.1 .or. (itime/iprint)*iprint == itime .or. itime == ndtn) then

        ! Switch on for bi-axial element test - writes vertical strain and q
        !write (LOGUnit, *) -epsg(2,1,1),',',-sigg(2,1,1)+sigg(1,1,1),',',epsg(1,1,1)+epsg(2,1,1)
        !write (*, *) -epsg(2,1,1),',', sigg(1,1,1)-sigg(2,1,1),',', itime*dt
        
        ! Switch on for cyclic tests - writes shear stress vs. shear strain
        !write (LOGUnit, *) epsg(3,1,1)*2.d0,',',sigg(3,1,1),',', itime*dt ! Multiplied by 2 to get gamma_xy
        !write (*, *) epsg(3,1,1)*2.d0,',',sigg(3,1,1),',', itime*dt
        
        ! Looking at vertical stress and strain on 1st element
        !write (LOGUnit, *) epsg(2,1,1),',',sigg(2,1,1),',',  itime*dt ! Multiplied by 2 to get gamma_xy
        !write (LOGUnit, *) sigg(2,1,1),',',dis(2),',',  itime*dt ! Multiplied by 2 to get gamma_xy
        write (LOGUnit, *) sigg(2,1,1),',',dis(2),',',  itime*dt ! Multiplied by 2 to get gamma_xy
        write (*, *) epsg(2,1,1),',',sigg(2,1,1),',', itime*dt
                    
        endif
    enddo
    call CPU_TIME(finish)
    write(*,*) 'Time = ',finish-start,' seconds'
    write(*,*) 'Time = ',(finish-start)/60.d0,' minutes'
    write(*,*) 'Time = ',(finish-start)/3600.d0,' hours'
    read(*,*)
    close(2)

    end subroutine solve

    subroutine readdata()
    !**********************************************************************
    !    Function: Make a file
    !**********************************************************************

    implicit none
    double precision :: r1(2), r2(2), LPos(2), temp, Gpy, PI
    integer :: I, J, IEl, INod(4), ig

    open(DATUnit, file = 'input.dat')
    read(DATUnit, *) NDivX, NDivY
    read(DATUnit, *) Meshdx, Meshdy
    NNod = (NDivX + 1)*(NDivY + 1) 
    NEl = NDivX * NDivY
    allocate (NodCo(2, NNod)); nodco = 0.d0
    allocate (ICon(4, NEl)); icon = -1
    allocate (IsFixDof(2 * NNod))
    allocate(Mas(2 * NNod)); Mas = 0.d0
    allocate(GrvF(2 * NNod)); GrvF = 0.d0
    allocate(InrF(2 * NNod)); InrF = 0.d0
    allocate(ExtF(2 * NNod)); ExtF = 0.d0
    allocate(DamF(2 * NNod)); DamF = 0.d0
    allocate(v(2 * NNod)); V = 0.d0
    allocate(v0(2 * NNod)); V0 = 0.d0
    allocate(dis(2 * NNod)); dis = 0.d0
    allocate(B(2, 4, NEl, 4)); B = 0.d0
    allocate(Area(NEl)); Area = 0.d0
    allocate(Areai(NEl)); Areai = 0.d0
    allocate(Sigg(4, 4, NEl)); Sigg = 0.d0
    allocate(Sig0(4, 4, NEl));  Sig0 = 0.d0
    !Sigg(1,:,:) = -100.d0   ! Stress value xx initialisation
    !Sigg(2,:,:) = -100.d0   ! Stress value yy initialisation
    !Sigg(4,:,:) = -100.d0   ! Stress value zz initialisation
    allocate(Swp(4, 4, NEl));  Swp = 0.d0
    allocate(FBARNEW(1, 4, NEl));  FBARNEW = 0.d0
    allocate(DSwp(4, 4, NEl)); DSwp = 0.d0
    allocate(EpsP(4, 4, NEl)); EpsP = 0.d0
    allocate(edP(1 ,4, NEl)); edP = 0.d0
    allocate(EpsE(3, 4, NEl)); EpsE = 0.d0
    allocate(eta(1,4,nel)); eta =0.d0
    allocate(EpsG(3, 4, NEl)); EpsG = 0.d0
    allocate(F(4, 4, nel)); F(1,:,:) = 1.d0; F(2,:,:) = 1.d0; F(3,:,:) = 0.d0; F(4,:,:) = 0.d0
    allocate(B_trial(4, 4, nel)); B_trial(1,:,:) = 1.d0; B_trial(2,:,:) = 1.d0; B_trial(3,:,:) = 0.d0; B_trial(4,:,:) = 0.d0

    allocate(Statevar(32,4,Nel)); Statevar = 0.d0 ! State variable
    allocate(PlastInd(1,4,Nel)); PlastInd = 0.d0 ! State variable

    allocate(HS(4, NEL, 4)); HS = 0.d0

    maxdx=NodCo(1, 1); maxdy=NodCo(2, 1)
    do I = 1, NNod
        read(DATUnit, *) temp, NodCo(1, I), NodCo(2, I), temp
        if(NodCo(1, I).gt.maxdx) then; maxdx=NodCo(1, I); endif; if(NodCo(2, I).gt.maxdy) then; maxdy=NodCo(2, I); endif;
    end do
        
    eledx=maxdx/NDivX; eledy=maxdy/NDivY;

    do I = 1, NEl
        read(DATUnit, *) temp, ICon(1, I), ICon(2, I), ICon(3, I), ICon(4, I)
    end do
    
    read(DATUnit, *) MatProp(1), MatProp(2), MatProp(3), MatProp(4), MatProp(5), MatProp(6), MatProp(7) ! Phi and Psi added
    read(DATUnit, *) MatProp(8), MatProp(9), MatProp(10) !DUMMY VALUES
    read(DATUnit, *) MatProp(11), MatProp(12), MatProp(13), MatProp(14), MatProp(15), MatProp(16)
    read(DATUnit, *) MatProp(17), MatProp(18), MatProp(19), MatProp(20), MatProp(21), MatProp(22), MatProp(23), MatProp(24), MatProp(25)
    !read(DATUnit, *) MatProp(26), MatProp(27), MatProp(28), MatProp(29), MatProp(30), MatProp(31), MatProp(32) !for clay only
    read(DATUnit, *) MatProp(26)!for sand constitutive model only, bulk mod of water
    read(DATUnit, *) gx, gy
    read(DATUnit, *) ndtn, dt
    close(1)
    !Bulkw = MatProp(4)
    do IEl = 1, NEl
        INod(:) = ICon(:, IEl)
        !Area(IEl) = abs(NodCo(1,INod(3))-NodCo(1,INod(1)))* abs(NodCo(2,INod(3))-NodCo(2,INod(1)))
    end do

    do I=1,nnod
        NodCo(1,I) = NodCo(1,I) * Meshdx
        NodCo(2,I) = NodCo(2,I) * Meshdy
    enddo
    PI = 4.d0 * atan(1.d0) ! Approximation for PI
    PHI  = SIN(PI*MatProp(5)/180.D0)
    PSI  = SIN(PI*MatProp(6)/180.D0)
    COH = MatProp(7)*COS(PI*PHI/180.D0)
   

    !factorg = dt / (ndtn*dt)

!    write(*,*) 'Gravity Load Factor = ', factorg
!    read(*,*) 
    !Initial stresses, gravity stress
    !do IEl = 1, nel
    !    INod(:) = ICon(:, IEl)
    !    do ig = 1, 4
    !        GPy=(NodCo(2,INod(1))+NodCo(2,INod(4)))/2.d0
    !        Sig0(2,ig,IEl)=(0.5d0-GPy)*MatProp(3)*gy
    !        Sig0(1,ig,IEl)=Sig0(2,ig,IEl)/2.0
    !        Sig0(3,ig,IEl)=(Sig0(1,ig,IEl)+Sig0(2,ig,IEl))/2.0
    !    enddo !gauss
    !end do !nel
    
    ! Boundary Conditions
    IsFixDof = .false.
    !do I = 1, NNod ! vertical bar problem
    !    if ((NodCo(2, I) .eq. 0.d0)) then
    !        IsFixDof((I - 1) * 2 + 2) = .true.
    !    end if
    !    if ((NodCo(2, I) .eq. 0.d0)) then
    !        ! IsFixDof((I - 1) * 2 + 1) = .true.
    !    end if
    !    if ((NodCo(2, I) .eq. 1.d0)) then
    !        !  IsFixDof((I - 1) * 2 + 2) = .true.
    !    end if
    !end do
    !IsFixDof(1) = .true.   ! switch off for biaxial test
    !IsFixDof(2) = .true.    
    !IsFixDof(3) = .true.   ! switch off for biaxial test
    !IsFixDof(4) = .true. 
    
    do I = 1,Nnod
        if(NodCo(2,i).eq.0.d0) then
            IsFixDof((I - 1) * 2 + 2) = .true. 
            IsFixDof((I - 1) * 2 + 1) = .true.  
        end if !Base fixed
        if(NodCo(1,i).eq.0.d0) then
            IsFixDof((I - 1) * 2 + 1) = .true. 
        end if !Left side roller support
        if(NodCo(1,i).eq.0.2d0) then
            IsFixDof((I - 1) * 2 + 1) = .true. 
        end if !Right side roller support
    enddo 
    
    !WRITE(*,*) 'Youngs Modulus = ', MatProp(1)
    !WRITE(*,*) 'Poissons ration = ', MatProp(2)
    WRITE(*,*) 'Mass = ', MatProp(3)
    WRITE(*,*) 'Bulk Modulus of Water = ', MatProp(4)
    !WRITE(*,*) 'PHI = ', PHI
    !WRITE(*,*) 'PSI = ', PSI
    !WRITE(*,*) 'Cohesion = ', COH
    WRITE(*,*) 'Gravity in Y Direction = ', gy
    !READ(*,*)    
    
    end subroutine readdata

    subroutine Initial()

    !**********************************************************************
    !    Function: Calculates B-Strain Displacement matrix (1 and 4 gauss points)
    !**********************************************************************

    implicit none

    integer :: IEl, I, J, K, INod(4), Id, ig
    double precision :: LPos(2, 1), dNxi(2, 4), Ja(2, 2), JaI(2, 2), A
    double precision :: xi, eta, rp, rm, sp, sm, temp 

    temp = 1.d0/sqrt(3.d0)  ! Gauﬂ points - 4 
    !temp=0.d0              ! Gauﬂ points - 1

    B = 0.0d0

    do IEl = 1, nel
        INod(:) = ICon(:, IEl)

        do ig = 1, 4

            select case (ig)
            case(1)
                xi = -temp
                eta = -temp

            case (2)
                xi = temp
                eta = -temp

            case(3)
                xi = temp
                eta = temp

            case(4)
                xi = -temp
                eta = temp
            end select

            rp = 1.0 + xi
            rm = 1.0 - xi;
            sp = 1.0 + eta;
            sm = 1.0 - eta;

            dNxi(1, 1) = -0.25D0 * sm; dNxi(1, 2) = +0.25D0 * sm; dNxi(1, 3) = +0.25D0 * sp
            dNxi(1, 4) = -0.25D0 * sp; dNxi(2, 1) = -0.25D0 * rm; dNxi(2, 2) = -0.25D0 * rp
            dNxi(2, 3) = +0.25D0 * rp; dNxi(2, 4) = +0.25D0 * rm

            HS(1, iel, ig) = (1.D0 - xi)*(1.D0 - eta)/4.D0
            HS(2, iel, ig) = (1.D0 + xi)*(1.D0 - eta)/4.D0
            HS(3, iel, ig) = (1.D0 + xi)*(1.D0 + eta)/4.D0
            HS(4, iel, ig) = (1.D0 - xi)*(1.D0 + eta)/4.D0

            Area(IEl) = abs(NodCo(1, INod(3)) - NodCo(1, INod(1))) * abs(NodCo(2, INod(3)) - NodCo(2, INod(1)))

            Ja = 0.0D0

            do I = 1, 2
                do J = 1, 2
                    do K = 1, 4
                        Ja(I, J) = Ja(I, J) + dNxi(I, K) * NodCo(J, INod(K))
                    end do
                end do
            end do
            A = Ja(1, 1) * Ja(2, 2) - Ja(1, 2) * Ja(2, 1)

            if (A .gt. 0.D0) then
                JaI(1, 1) = +Ja(2, 2)/A; JaI(1, 2) = -Ja(1, 2)/A
                JaI(2, 1) = -Ja(2, 1)/A; JaI(2, 2) = +Ja(1, 1)/A
            else
                write(LOGUnit, *) 'negative or zero Jacobian !!'; stop
            end if

            do J = 1, 4
                B(1, J, IEl, ig) = dNxi(1, J) * JaI(1, 1) + dNxi(2, J) * JaI(1, 2)
                B(2, J, IEl, ig) = dNxi(1, J) * JaI(2, 1) + dNxi(2, J) * JaI(2, 2)
            end do
        enddo
    enddo

    ! call Map2Nod()

    end subroutine Initial

    subroutine Map2Nod()

    !**********************************************************************
    !    Function: Internal and External Forces, Mass at each nodes
    !**********************************************************************

    implicit none
    integer :: I, IEl, Id, INod(4), ig, J, count = 0
    double precision :: factor, factor1, kp, ks, c, G_mod, K_mod, cp, consta,constb, E_mod, cs
    GrvF = 0.D0; InrF = 0.D0; ExtF = 0.D0; Mas = 0.D0; Damf=0.d0
       
    !!!!!!!!!!!!!!!!! ******     
    !k = 100.d0
    !c = 5000.d0 
    !k = 0.d0; c = 0.d0
    consta=1.d0; constb=1.d0;
    delta = min(eledx/(2.d0*consta),eledy/(2.d0*constb)) ! distance of the supposed damping layer
    G_mod = MatProp(1) / (2.d0 * (1 + MatProp(2)))
    K_mod = MatProp(1) / (3.d0 * (1 - 2 * MatProp(2)))
    E_mod = MatProp(1)*(1.d0-MatProp(2))/((1+MatProp(2))*(1.d0-2.d0*MatProp(2)))
    cp=sqrt(E_mod/MatProp(3))
    cs=sqrt(G_mod/MatProp(3))
    kp = MatProp(3) * cp**2.d0 / delta
    ks = MatProp(3) * cs**2.d0 / delta
    !k = 1000
    
    !!!!!!!!!!*************
     if(itime*dt.ge.0.d0) then 
        factor1 = itime * dt / 1.d0 
    endif 
    
     do I = 1,Nnod
        if(NodCo(2,i).eq.1.d0) then
            if(NodCo(1,i).ge.0.0d0.and.NodCo(1,i).le.0.2d0) then
                count = count + 1
            end if 
        end if
     end do
     
    !if(itime.eq.300001.or.itime.eq.350001.or.itime.eq.400001.or.itime.eq.450001.or.itime.eq.500001.or.itime.eq.550001.or.itime.eq.600001.or.itime.eq.650001.or.itime.eq.700001.or.itime.eq.750001) then 
    !if(itime.eq.300001.or.itime.eq.400001.or.itime.eq.500001.or.itime.eq.600001.or.itime.eq.700001) then 
    !if(itime.eq.30001.or.itime.eq.40001.or.itime.eq.50001.or.itime.eq.60001.or.itime.eq.70001.or.itime.eq.80001.or.itime.eq.90001) then 
    !if(itime.eq.1) then 
        do I = 1,Nnod
            if(NodCo(2,i).eq.1.d0) then
                if(NodCo(1,i).eq.0.0d0.or.NodCo(1,i).le.0.2d0) then
                    !ExtF((I - 1) * 2 + 2) = -100.d0/count * factor1
                    !ExtF((I - 1) * 2 + 2) = -0.50d0  !* 1000.d0 !/ count * 2.d0
                    ExtF((I - 1) * 2 + 2) = -100.d0 / count !* 1000.d0
                else 
                    !ExtF((I - 1) * 2 + 2) = -0.25d0 !* 1000.d0 !/ count * 2.d0
                end if
            end if
        end do     
    !end if         
    
   if(itime*dt.gt.1.d0) then 
        factor = 1.d0
    else
        factor = itime * dt / 1.d0    
    endif 
    
    if(dt.gt.min((2.d0*consta*delta/cp),(2.d0*constb*delta/cs))) write(*,*) 'reduce'
   !write(*,*) itime*dt,factor
    
    do IEl = 1, nel
        INod(:) = ICon(:, IEl)
        do ig = 1, 4 ! gauss
            do I = 1, 4
                Id = (INod(I) - 1) * 2
                !InrF(Id + 1) = InrF(Id + 1)+ ((Sigg(1,ig,IEl)+Swp(1,ig,Iel)) * B(1,I,iel,ig) +&
                !    (Sigg(3,ig,IEl))               * B(2, I, iel, ig)) * Area(iel)/4.d0  ! For undrained conditions
                !InrF(Id + 2) = InrF(Id + 2)+ ((Sigg(3,ig,IEl))               * B(1,I,iel,ig) +& 
                !    (Sigg(2,ig,IEl)+Swp(2,ig,iel)) * B(2, I, iel, ig)) * Area(iel)/4.d0  ! For undrained conditions
                InrF(Id + 1) = InrF(Id + 1)+ ((Sigg(1,ig,IEl)) * B(1,I,iel,ig) +&
                   (Sigg(3,ig,IEl)) * B(2, I, iel, ig)) * Area(iel)/4.d0 ! For drained conditions
                InrF(Id + 2) = InrF(Id + 2)+ ((Sigg(3,ig,IEl)) * B(1,I,iel,ig) +& 
                   (Sigg(2,ig,IEl)) * B(2, I, iel, ig)) * Area(iel)/4.d0 ! For drained conditions
                GrvF(Id + 1) = GrvF(Id + 1) + Area(iel) * MatProp(3) * hs(i, iel, ig) * gx /4.d0 
                GrvF(Id + 2) = GrvF(Id + 2) + Area(iel) * MatProp(3) * hs(i, iel, ig) * gy /4.d0 !* factor ! factor for gravity load application
                Mas(Id + 1) = Mas(Id + 1) + Area(iel) * MatProp(3) * hs(i, iel, ig) /4.d0
                Mas(Id + 2) = Mas(Id + 2) + Area(iel) * MatProp(3) * hs(i, iel, ig) /4.d0
                
                !!!!!***********
                
                if(IsFixDof(Id+1).and.IsFixDof(Id+2)) then
                DamF(Id + 1) = DamF(Id + 1) + constb*MatProp(3)*cs*Area(Iel)*v0(Id+1) + ks*Area(Iel)*dis(Id+1)
                DamF(Id + 2) = DamF(Id + 2) + consta*MatProp(3)*cp*Area(Iel)*v0(Id+2) + kp*Area(Iel)*dis(Id+2)
                end if
                
                if(IsFixDof(Id+1).and..not.IsFixDof(Id+2)) then
                DamF(Id + 1) = DamF(Id + 1) + consta*MatProp(3)*cp*Area(Iel)*v0(Id+1) + kp*Area(Iel)*dis(Id+1)
                DamF(Id + 2) = DamF(Id + 2) + constb*MatProp(3)*cs*Area(Iel)*v0(Id+2) + ks*Area(Iel)*dis(Id+2)
                end if
                
                !!!!!***********
                
            end do
        enddo !gauss
    end do !nel
    end subroutine Map2Nod

    subroutine Update()

    !**********************************************************************/************************
    !    Function: Large deformation formulation with the Deformation gradient, Left Cauchy Tensors
    !***********************************************************************************************

    implicit none
    integer :: IEl, INod(4), Id, I, J, ig, n = 2
    double precision :: delV(4), deps(3), tem, dampf = 0.01d0, df(4), B_temp(4)
    double precision :: B_eigen(2, 2), B_ev(2, 2), abserr = 1.0e-09


    V = 0.D0
    do I = 1, NNod
        Id = (I - 1) * 2
        if (.not.IsFixDof(Id + 1).and.(Mas(Id + 1) .gt. 0.D0)) then

            tem = V0(Id + 1)+(GrvF(Id + 1) + ExtF(Id + 1) - InrF(Id + 1))/Mas(Id + 1) * dt
            V(Id + 1) = tem - sign(1.d0, tem) * dampf * abs(tem)
        end if

        if (.not.IsFixDof(Id + 2).and.(Mas(Id + 2) .gt. 0.D0)) then

            tem = V0(Id + 2)+(GrvF(Id + 2) + ExtF(Id + 2) - InrF(Id + 2))/Mas(Id + 2) * dt
            V(Id + 2) = tem - sign(1.d0, tem) * dampf * abs(tem)
        end if
    end do

    !V = 0.d0
    !V(7) = 0.1d0
    !V(5) = 0.1d0

    !                  Prescribed velocity
    !        do I = 1, NNod
    !            Id = (I - 1) * 2
    !            if (NodCo(1, I) .eq. 0.d0) then
    !                V(Id + 1) = -0.05d0
    !            end if
    !            if (NodCo(1, I) .eq. 2.d0) then
    !                V(Id + 1) = 0.05d0
    !            end if
    !        enddo

    do IEl = 1, NEl
        INod(:) = ICon(:, IEl)
        do ig = 1, 4
            delV = 0.0
            do I = 1, 4
                Id = (INod(I) - 1) * 2
                delV(1) = delV(1) + B(1, i, iel, ig) * V(Id + 1)
                delV(2) = delV(2) + B(2, i, iel, ig) * V(Id + 2)
                delV(3) = delV(3) + B(2, i, iel, ig) * V(Id + 1)
                delV(4) = delV(4) + B(1, i, iel, ig) * V(Id + 2)
            enddo

            !-----------nonlinear elastic---------------------!
            !--------------deformation gradient------------------!

            dF(1) = 1 + delV(1) * dt
            dF(2) = 1 + delV(2) * dt
            dF(3) = delV(3) * dt
            dF(4) = delV(4) * dt

            B_temp(1) = B_trial(1, ig, iel) * dF(1) + B_trial(3, ig, iel) * dF(3)
            B_temp(2) = B_trial(3, ig, iel) * dF(4) + B_trial(2, ig, iel) * dF(2)
            B_temp(3) = B_trial(1, ig, iel) * dF(4) + B_trial(3, ig, iel) * dF(2)
            B_temp(4) = B_trial(4, ig, iel) * dF(1) + B_trial(2, ig, iel) * dF(3)

            B_trial(1, ig, iel) = df(1) * B_temp(1) + df(3) * B_temp(3)
            B_trial(2, ig, iel) = df(4) * B_temp(3) + df(2) * B_temp(2)
            B_trial(3, ig, iel) = df(1) * B_temp(3) + df(3) * B_temp(2)
            B_trial(4, ig, iel) = df(4) * B_temp(1) + df(2) * B_temp(4)

            B_eigen(1, 1) = B_trial(1, ig, iel)
            B_eigen(2, 2) = B_trial(2, ig, iel)
            B_eigen(1, 2) = B_trial(3, ig, iel)
            B_eigen(2, 1) = B_trial(4, ig, iel)

            call Eigen(B_eigen, B_ev, abserr, n)

            EpsG(1, ig, iel) = (0.5d0 * log(b_eigen(1, 1)) * b_ev(1, 1)**2) + (0.5d0 * log(b_eigen(2, 2)) * b_ev(1, 2)**2)
            EpsG(2, ig, iel) = (0.5d0 * log(b_eigen(1, 1)) * b_ev(2, 1)**2) + (0.5d0 * log(b_eigen(2, 2)) * b_ev(2, 2)**2)
            EpsG(3, ig, iel) = (0.5d0 * log(b_eigen(1, 1)) * b_ev(1, 1) * b_ev(2, 1))+(0.5d0 * log(b_eigen(2, 2)) * b_ev(1, 2) * b_ev(2, 2))

            !                      ! to check
            !                      F(1, ig,iel) = dF(1) * F(1, ig,iel) + dF(3) * F(4, ig,iel)
            !                      F(2, ig,iel) = dF(4) * F(3, ig,iel) + dF(2) * F(2, ig,iel)
            !                      F(3, ig,iel) = dF(1) * F(3, ig,iel) + dF(3) * F(2, ig,iel)
            !                      F(4, ig,iel) = dF(4) * F(1, ig,iel) + dF(2) * F(4, ig,iel)

            !------------neo-Hookean  model (old)  -------------------------!
            !  call Finitedef(MatProp(1), MatProp(2), SigG(:, ig, iel), EpsG(:, ig, iel))

        enddo !gauss
    end do ! elements

    dis = dis + v * dt
    v0 = v

    ! Node-co updated
    !              do J = 1, Nnod
    !                  Id = (J - 1) * 2
    !                  NodCo(1, J) = NodCo(1, J) + v(Id + 1) * dt
    !                  NodCo(2, J) = NodCo(2, J) + v(Id + 2) * dt
    !              enddo
    !
    !              call initial()

    end subroutine Update


    subroutine  Update_small()
    !**********************************************************************
    !    Function: Small deformation formulation
    !*************************** *******************************************

    implicit none
    integer :: IEl, INod(4), Id, I, ig, J
    double precision :: delV(4), deps(4), tem, dampf = 0.d0, check, epsilon, rad, WaterP
    !Initial Velocity
    !V = 0.d0
    
        !IsFixDof = .false.
        !do I = 1,Nnod
        !if(NodCo(2,i).eq.0.d0) then
        !    IsFixDof((I - 1) * 2 + 2) = .true. 
        !end if !Base fixed in y direction
        !enddo
     
    !ExtF(5) =  50.d0
    !ExtF(7) = -50.d0
    !ExtF(1) =  50.d0
    !ExtF(3) = -50.d0
    
    do I = 1, NNod
        Id = (I - 1) * 2
        !if (.not.IsFixDof(Id + 1).and.(Mas(Id + 1) .gt. 0.D0)) then
        if ((Mas(Id + 1) .gt. 0.D0)) then
            tem = V0(Id + 1)+(GrvF(Id + 1) + ExtF(Id + 1) - InrF(Id + 1) - Damf(Id+1))/Mas(Id + 1) * dt 
            V(Id + 1) = tem - sign(1.d0, tem) * dampf * abs(tem) 
        end if
    
        !!if (.not.IsFixDof(Id + 2).and.(Mas(Id + 2) .gt. 0.D0)) then
        if ((Mas(Id + 2) .gt. 0.D0)) then
            tem = V0(Id + 2)+(GrvF(Id + 2) + ExtF(Id + 2) - InrF(Id + 2) - Damf(Id + 2))/Mas(Id + 2) * dt
            V(Id + 2) = tem - sign(1.d0, tem) * dampf * abs(tem) 
        end if
        !dampf = 0.d0
    end do
    
    do IEl = 1, NEl
        INod(:) = ICon(:, IEl)
        do ig = 1, 4 ! gauss
            delV = 0.0
            do I = 1, 4
                Id = (INod(I) - 1) * 2
                delV(1) = delV(1) + B(1, i, iel, ig) * V(Id + 1)
                delV(2) = delV(2) + B(2, i, iel, ig) * V(Id + 2)
                delV(3) = delV(3) + B(2, i, iel, ig) * V(Id + 1)
                delV(4) = delV(4) + B(1, i, iel, ig) * V(Id + 2)
            end do
            dEps(1:2) = delV(1:2) * dt
            dEps(3) = (delV(3) + delV(4)) * dt
        
            EpsG(1:2, ig, IEl) = EpsG(1:2, ig, IEl) +  dEps(1:2)
            EpsG(3, ig, IEl) = EpsG(3, ig, IEl) +  dEps(3)

            EpsV = deps(1) + deps(2)                              ! Volumetric strain, for Undrained analysis
            Dswp(1:4,ig,iel) = BulkW * EpsV                       ! Incremental water pressure, for Undrained analysis
            Swp(1:4,ig,iel) = Swp(1:4,ig,iel) + Dswp(1:4,ig,iel)  ! Total water pressure, for Undrained analysis
        !enddo !gauss
            
        !do ig = 1, 4
           call Elastic(MatProp(1),MatProp(2), deps(:), Sigg(:,ig,IEl)) ! Elastic material
           !call Intf_DLL_Plaxis(dt, MatProp(:), dEps(:), Sigg(:,ig,iel), Statevar(:,ig,iel))
        enddo !gauss
    end do ! particles
    !sig0 = sigg
    !plastind(1,:,:) = statevar(31,:,:)
    dis = dis + v * dt
    v0 = v

    end subroutine Update_small

    !*********************************************************
    !Elastic Hookes model
    !**********************************************************
    subroutine Elastic(E, nu, eps, Sig)
    
    implicit none
    
    double precision, intent(in) :: E, nu
    double precision, intent(inout) :: Sig(8)
    double precision, intent(in) :: eps(3)
    ! local variables
    double precision :: G_mod, K_mod, Eps_tr
    
    G_mod = E / (2.d0 * (1 + nu))
    K_mod = E / (3.d0 * (1 - 2 * nu))
    
    Eps_tr = eps(1) + eps(2)
    Sig(1) = sig(1) +  ((K_mod * Eps_tr) + 2 * G_mod * (eps(1) - (Eps_tr/3.d0)))
    Sig(2) = sig(2) +  ((K_mod * Eps_tr) + 2 * G_mod * (eps(2) - (Eps_tr/3.d0)))
    Sig(3) = sig(3) +  (2 * G_mod * eps(3))
    Sig(4) = sig(4) +  ((K_mod * Eps_tr) + 2 * G_mod * (0.d0 - (Eps_tr/3.d0)))
    
    !Sig(3) = sig(3) +  ((K_mod * Eps_tr) + 2 * G_mod * (0.d0 - (Eps_tr/3.d0)))
    !Sig(4) = sig(4) +  (2 * G_mod * eps(3))
    
    endsubroutine elastic



    subroutine VonMises(E, enu, Yield, Eps3, EpsP, s, EpsE)
    !*********************************************************************
    ! Von Mises Elasto-Plastic Model
    !*******************************************************************************

    implicit double precision (a - h, o - z)

    double precision, intent(in) :: E, enu, Yield, Eps3(3)
    double precision, intent(inout) :: S(7)
    double precision, intent(inout) :: EpsP(6), EpsE(3)
    integer :: i
    double precision :: Eps(6), devt(6), dir(6), dev(6), eps_v(3), EpsPT(3)

    B_K = E / (3.d0 * (1 - 2.d0 * enu)) ! Bulk-modulus K
    G_mod = E /(2.d0 * (1.d0 + enu)) ! 2nd lame - mu (G) shear mod

    gamma=0.d0
    eps = 0.d0
    eps(1:2) = Eps3(1:2)
    eps(4) = Eps3(3)

    ! write (LOgUnit, * ) Eps3
    treps = eps(1) + eps(2) + eps(3)

    ! dev of strain
    do i=1,3
        dev(i) = eps(i) - treps/3.d0
        dev(i+3) = eps(i+3) / 2.d0
    enddo

    !deviatoric trial force
    do i=1,6
        devt(i) = 2.d0 * G_mod * (dev(i) - EpsP(i))
    enddo

    !norm
    devtnorm = dsqrt(devt(1)**2.d0+ devt(2)**2.d0 + devt(3)**2.d0 + 2.d0* devt(4)**2.d0  &
        + 2.d0* devt(5)**2.d0 + 2.d0* devt(6)**2.d0)
    !write (LOGUnit, * ) devt(1), devt(2), devt(3)
    if (devtnorm.eq.0.d0)  devtnorm =0.0000001d0

    !direction of plastic flow
    do i = 1, 6
        dir(i) = devt(i)/(devtnorm)
    enddo

    !determine yield criterion
    Yn = (dsqrt(2.d0/3.d0) * (Yield ))
    phi = devtnorm - Yn

    !to compute stresses
    if ((phi .lt. 0.d0) ) then !elastic
        s(1) = B_K * treps + devt(1) !  stresses
        s(2) = B_K * treps + devt(2)
        s(3) = B_K * treps + devt(3)
        s(4) =  devt(4) !   plastic strains
        s(5) = EpsP(2)
        s(6) = EpsP(4)
        s(7) = gamma
        EpsE(1) = eps(1)
        EpsE(2) = eps(2)
        EpsE(3) = eps(3)

    else !plastic
        flag = 1.d0
        Yn = (dsqrt(2.d0/3.d0) * (Yield ))

        gamma = phi / (2.d0  * G_mod)

        !update plastic strain
        do i = 1, 6
            EpsP(i) = EpsP(i)+ (dir(i)*(gamma))
        enddo


        EpsPv = (Eps(1) - EpsE(1) + Eps(2) - EpsE(2) + Eps(3) - EpsE(3)) / 3.d0

        do i=1,3
            EpsPT(i) = EpsP(i) + EpsPv
        enddo

        treps1 =Eps(1) - EpsPT(1) + Eps(2) - EpsPT(2) + Eps(3) - EpsPT(3)

        s(1) = B_K * treps1 + devt(1) - ( 2.d0 * G_mod * flag * gamma * dir(1)) !  stresses
        s(2) = B_K * treps1 + devt(2) - ( 2.d0 * G_mod *  flag * gamma * dir(2))
        s(3) =  B_K * treps1 + devt(3) - ( 2.d0 * G_mod *  flag * gamma * dir(3))
        s(4) = devt(4) - 2.d0 * G_mod * (flag * gamma * dir(4)) !   plastic strains
        s(5) = EpsP(2)
        s(6) = EpsP(4)
        sigeq = (0.5d0 * ((s(1)-s(2))**2 + (s(2)-s(3))**2  +  (s(3)-s(1))**2 + 6*s(4)**2))
        s(7) = dsqrt (sigeq)

    endif

    end subroutine VonMises



    subroutine Eigen(a, x, abserr, n)
    !===========================================================
    ! Evaluate eigenvalues and eigenvectors
    ! of a real symmetric matrix a(n,n): a*x = lambda*x
    ! method: Jacoby method for symmetric matrices
    ! Alex G. (December 2009)
    !-----------------------------------------------------------
    ! input ...
    ! a(n,n) - array of coefficients for matrix A
    ! n      - number of equations
    ! abserr - abs tolerance [sum of (off-diagonal elements)^2]
    ! output ...
    ! a(i,i) - eigenvalues
    ! x(i,j) - eigenvectors
    ! comments ...
    !===========================================================
    implicit none
    integer i, j, k, n
    double precision a(n, n), x(n, n)
    double precision abserr, b2, bar
    double precision beta, coeff, c, s, cs, sc

    ! initialize x(i,j)=0, x(i,i)=1
    ! *** the array operation x=0.0 is specific for Fortran 90/95
    x = 0.0
    do i = 1, n
        x(i, i) = 1.0
    end do

    ! find the sum of all off-diagonal elements (squared)
    b2 = 0.0
    do i = 1, n
        do j = 1, n
            if (i .ne. j) b2 = b2 + a(i, j)**2
        end do
    end do

    if (b2 <= abserr) return

    ! average for off-diagonal elements /2
    bar = 0.5 * b2/float(n * n)

    do while (b2 .gt. abserr)
        do i = 1, n - 1
            do j = i + 1, n
                if (a(j, i)**2 <= bar) cycle ! do not touch small elements
                b2 = b2 - 2.0 * a(j, i)**2
                bar = 0.5 * b2/float(n * n)
                ! calculate coefficient c and s for Givens matrix
                beta = (a(j, j) - a(i, i))/(2.0 * a(j, i))
                coeff = 0.5 * beta/sqrt(1.0 + beta**2)
                s = sqrt(max(0.5 + coeff, 0.0))
                c = sqrt(max(0.5 - coeff, 0.0))
                ! recalculate rows i and j
                do k = 1, n
                    cs = c * a(i, k) + s * a(j, k)
                    sc = -s * a(i, k) + c * a(j, k)
                    a(i, k) = cs
                    a(j, k) = sc
                end do
                ! new matrix a_{k+1} from a_{k}, and eigenvectors
                do k = 1, n
                    cs = c * a(k, i) + s * a(k, j)
                    sc = -s * a(k, i) + c * a(k, j)
                    a(k, i) = cs
                    a(k, j) = sc
                    cs = c * x(k, i) + s * x(k, j)
                    sc = -s * x(k, i) + c * x(k, j)
                    x(k, i) = cs
                    x(k, j) = sc
                end do
            end do
        end do
    end do
    return
    end subroutine Eigen

    subroutine MkOpFiles()
    !**********************************************************************
    !    Function: Output files
    !**********************************************************************
    implicit none

    call MkFile(MSHUnit, 'Model.post.msh')
    call MkFile(RESUnit, 'Model.post.res')

    end subroutine MkOpFiles


    subroutine MkFile(Unit, flNam)
    !**********************************************************************
    !
    !    Function: Make a file
    !
    !**********************************************************************

    implicit none

    integer Unit
    character flNam * (*)

    if (FlExist(flNam)) then
        open(Unit, file = flNam)
        close(Unit, Status = 'Delete')
    endif
    open(Unit, file = flNam)

    end subroutine MkFile

    logical function FlExist(flNam)
    !**********************************************************************
    !
    !    Function: To check the existence of a file
    !
    !**********************************************************************

    implicit none

    logical lfil
    character flNam * (*)

    lfil = .false.
    inquire(file = flNam, exist = lfil)
    if (lfil) then
        FlExist = .true.
    else
        FlExist = .false.
    endif

    end function FlExist


    subroutine WtMSH(Unit)
    !**********************************************************************
    !    Function: Writing GiD *.msh file
    !!**********************************************************************

    implicit none
    integer Unit
    ! local variables
    integer :: IPart, INod, IEl, J, K
    double precision, dimension(2) :: Pos, r1, r2, VPos

    write(Unit, *) 'MESH dimension 2 ElemType Quadrilateral Nnode 4'
    write(Unit, *) 'Coordinates'
    do INod = 1, NNod
        write(Unit, 111) INod, NodCo(1, INod), NodCo(2, INod)
    end do
    write(Unit, *) 'End Coordinates'
    write(Unit, *) 'Elements'
    do IEl = 1, NEl
        write(Unit, "(5I7)") IEl, ICon(1, IEl), ICon(2, IEl), ICon(3, IEl), ICon(4, IEl)
    end do
    write(Unit, *) 'End Elements'
    close(Unit)
111 format(I7, 2E16.6E3)

    end subroutine WtMSH

    subroutine WtRES(IStep)
    !**********************************************************************
    !
    !    Function:
    !
    !**********************************************************************

    implicit none

    integer IStep
    ! local variables
    integer :: IEl, J, K, INod, Id, Unit, ig

    Unit = ResUnit
    if (IStep .eq. 1) then
        write(Unit, *) 'GiD Post Results File 1.0'
        write(Unit, *) 'GaussPoints "Material_Point" Elemtype Quadrilateral'
        write(Unit, *) 'Number of Gauss Points: 4'
        write(Unit, *) 'Natural Coordinates: Internal'
        write(Unit, *) 'end gausspoints'
        write(Unit, *) 'Result "Boundary" "MPM"', IStep, 'Vector OnNodes'
        write(Unit, *) 'ComponentNames "X-fix", "Y-fix"'
        write(Unit, *) 'values'
        do INod = 1, NNod
            Id = (INod - 1) * 2
            J = 0; K = 0
            if (IsFixDof(Id + 1)) J = 1
            if (IsFixDof(Id + 2)) K = 1
            write(Unit, "(3I7)") INod, J, K
        end do
        write(Unit, *) 'end values'
    end if
    write(Unit, *) 'Result "displacement" "MPM"', IStep, 'Vector OnNodes'
    write(Unit, *) 'ComponentNames "comp. x", "comp. y"'
    write(Unit, *) 'values'
    do INod = 1, NNod
        Id = (INod - 1) * 2
        write(Unit, "(I7, 2E16.6E3)") INod, Dis(Id + 1), Dis(Id + 2)
    end do
    write(Unit, *) 'end values'

    write(Unit, *) 'Result "Stress" "MPM"', IStep, 'Vector OnGaussPoints "Material_Point"'
    write(Unit, *) 'ComponentNames "sigma xx", "sigma yy", "sigma zz", "sigma xy"'
    write(Unit, *) 'values'
    do IEl = 1, NEl
        write(Unit, "(I7, 4E16.6E4)") IEl, Sigg(1, 1, IEl), Sigg(2, 1, IEl), Sigg(4, 1, IEl),Sigg(3, 1, IEl)
        write(Unit, "(4E16.6E4)") Sigg(1, 2, IEl), Sigg(2, 2, IEl), Sigg(4, 1, IEl),Sigg(3, 2, IEl)
        write(Unit, "(4E16.6E4)") Sigg(1, 3, IEl), Sigg(2, 3, IEl),Sigg(4, 1, IEl), Sigg(3, 3, IEl)
        write(Unit, "(4E16.6E4)") Sigg(1, 4, IEl), Sigg(2, 4, IEl), Sigg(4, 1, IEl),Sigg(3, 4, IEl)
    end do
    write(Unit, *) 'end values'

    write(Unit, *) 'Result "Strain" "MPM"', IStep, 'Vector OnGaussPoints "Material_Point"'
    write(Unit, *) 'ComponentNames "eps xx", "eps yy", "eps xy"'
    write(Unit, *) 'values'
    do IEl = 1, NEl
        write(Unit, "(I7, 3E16.6E3)") IEl, EpsG(1, 1, IEl), EpsG(2, 1, IEl), EpsG(3, 1, IEl)
        write(Unit, "( 3E16.6E3)") Epsg(1, 2, IEl), EpsG(2, 2, IEl), EpsG(3, 2, IEl)
        write(Unit, "( 3E16.6E3)") Epsg(1, 3, IEl), Epsg(2, 3, IEl), Epsg(3, 3, IEl)
        write(Unit, "( 3E16.6E3)") Epsg(1, 4, IEl), Epsg(2, 4, IEl), Epsg(3, 4, IEl)
    end do
    write(Unit, *) 'end values'


    write(Unit, *) 'Result "Plastic strain" "MPM"', IStep, 'Vector OnGaussPoints "Material_Point"'
    write(Unit, *) 'ComponentNames "EpsP xx", "EpsP yy", "EpsP zz", "EpsP xy"'
    write(Unit, *) 'values'
    do IEl = 1, NEl
        write(Unit, "(I7, 4E16.6E4)") IEl, Sigg(1, 1, IEl), Sigg(2, 1, IEl), Sigg(3, 1, IEl), Sigg(4, 1, IEl)
        write(Unit, "(4E16.6E4)") Sigg(1, 2, IEl), Sigg(2, 2, IEl), Sigg(3, 2, IEl), Sigg(4, 1, IEl)
        write(Unit, "(4E16.6E4)") Sigg(1, 3, IEl), Sigg(2, 3, IEl), Sigg(3, 3, IEl),Sigg(4, 1, IEl)
        write(Unit, "(4E16.6E4)") Sigg(1, 4, IEl), Sigg(2, 4, IEl), Sigg(3, 4, IEl),Sigg(4, 1, IEl)
    end do
    write(Unit, *) 'end values'
    
    write(Unit, *) 'Result "Failure" "MPM"', IStep, 'Scalar OnGaussPoints "Material_Point"'
    write(Unit, *) 'ComponentNames "Failure Indicator"'
    write(Unit, *) 'values'
    do IEl = 1, NEl
        write(Unit, "(I7, 4E16.6E4)") IEl, Plastind(1,1,IEl)
        write(Unit, "(4E16.6E4)") Plastind(1,2,IEl)
        write(Unit, "(4E16.6E4)") Plastind(1,3,IEl)
        write(Unit, "(4E16.6E4)") Plastind(1,4,IEl)
    end do
    write(Unit, *) 'end values'    
    
    end subroutine WtRES
    
    end module modulefem

    !************program *****************
    program MPM2D

    use modulefem
    implicit none
    
  
    call MkOpFiles()
    call readdata()

    call solve()

    end program
    !*************program*****************