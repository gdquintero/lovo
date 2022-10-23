program lovo

    use sort
    use bmalgencan, only: algencan
    use iso_c_binding, only: c_ptr, c_loc,c_f_pointer
  
    implicit none
  
    ! Re-define this type (pdata_type) anyway you want. Algencan
    ! receives a pointer to a 'structure of this type'. Algencan has no
    ! access to the structure. It simple passes the pointer back to the
    ! user defined subroutines evalf, evalg, evalc, evalj, and
    ! evalhl. So, this is a trade safe way (no common blocks) of passing
    ! to the user-provided routines any information related to the
    ! problem. In this example, it is only being used for the user to
    ! count by itself the number of calls to each routine.
    type :: pdata_type
       integer :: counters(5) = 0
    end type pdata_type
  
    ! LOCAL SCALARS
    logical :: corrin,extallowed,rhoauto,scale
    integer :: allocerr,hlnnzmax,ierr,inform,istop,jnnzmax,m,maxoutit,n,nwcalls,nwtotit,outiter,p,totiter
    real(kind=8) :: bdsvio,csupn,epsfeas,epscompl,epsopt,f,finish,nlpsupn,rhoini,ssupn,start
    type(pdata_type), target :: pdata
    
    ! LOCAL ARRAYS
    logical, allocatable :: lind(:),uind(:)
    real(kind=8), allocatable :: c(:),lbnd(:),ubnd(:),lambda(:),x(:)

    !--> LOVO Algorithm variables <--

    integer, pointer :: samples=>null(),samples_train=>null(),samples_validation=>null()
    integer, pointer :: r_comb=>null(),order_lovo=>null(),rows_train=>null()
    real(kind=8), pointer :: t(:)=>null(),y(:)=>null(),train(:,:)=>null(),validation(:,:)=>null()
    real(kind=8) :: Fmin
    real(kind=8), allocatable :: Fmin_aux(:),indices(:)
    integer, allocatable :: Imin(:),combi(:)
    integer :: i,ind_train,n_Imin,i4_choose
    

    !--> End LOVO Algorithm variables <--

    allocate(order_lovo,rows_train,samples,samples_train,samples_validation,r_comb)

    ! Reading data and storing it in the variables t and y
    Open(Unit = 100, File = "output/data.txt", ACCESS = "SEQUENTIAL")

    read(100,*) samples

    samples_train = 10
    samples_validation = 30
    rows_train = samples - (samples_train + samples_validation) + 1
    order_lovo = samples_train - 2
    r_comb = i4_choose(samples_train,order_lovo)

    allocate(t(samples_train),y(samples),Fmin_aux(samples),indices(samples),Imin(r_comb),&
            combi(order_lovo),train(rows_train,samples_train),validation(rows_train,samples_validation),stat=allocerr)

    if ( allocerr .ne. 0 ) then
        write(*,*) 'Allocation error in main program'
        stop
    end if

    do i = 1, samples
        read(100,*) y(i)
    enddo

    t(:) = (/(i, i = 1, samples_train)/)

    close(100)

    !--> End Reading data and storing it in the variables t and y <--
  
    ! Number of variables
  
    n = 3
    
    allocate(x(n),lind(n),lbnd(n),uind(n),ubnd(n),stat=allocerr)
  
    if ( allocerr .ne. 0 ) then
       write(*,*) 'Allocation error.'
       stop
    end if
  
    ! Initial guess and bound constraints
    
    x(1:n) = 1.0d0
  
    lind(1:n) = .false.
    ! lbnd(1:n) = - 100.0d0
  
    uind(1:n) = .false.
    ! ubnd(1:n) =   100.0d0
  
    ! Number equality (m) and inequality (p) constraints.
    
    m = 0
    p = 0

    ind_train = 1

    call train_test_split()

    call compute_Fmin(x,n,ind_train,Fmin_aux,Fmin)

    call mount_Imin(x,n,Fmin,ind_train,combi,Imin,n_Imin)
    
    allocate(lambda(m+p),c(m+p),stat=allocerr)
  
    if ( allocerr .ne. 0 ) then
       write(*,*) 'Allocation error.'
       stop
    end if
  
    ! Initial guess for the Lagrange multipliers
    
    lambda(1:m+p) = 0.0d0
  
    ! Number of entries in the Jacobian of the constraints
    
    jnnzmax = n
  
    ! This should be the number of entries in the Hessian of the
    ! Lagrangian. But, in fact, some extra space is need (to store the
    ! Hessian of the Augmented Lagrangian, whose size is hard to
    ! predict, and/or to store the Jacobian of the KKT system). Thus,
    ! declare it as large as possible.
    
    hlnnzmax = 100000
  
    ! Feasibility, complementarity, and optimality tolerances
    
    epsfeas  = 1.0d-08
    epscompl = 1.0d-08
    epsopt   = 1.0d-08
  
    ! Maximum number of outer iterations
    
    maxoutit = 50
  
    ! rhoauto means that Algencan will automatically set the initial
    ! value of the penalty parameter. If you set rhoauto = .false. then
    ! you must set rhoini below with a meaningful value.
    rhoauto = .true.
  
    if ( .not. rhoauto ) then
       rhoini = 1.0d-08
    end if
  
    ! scale = .true. means that you allow Algencan to automatically
    ! scale the constraints. In any case, the feasibility tolerance
    ! (epsfeas) will be always satisfied by the UNSCALED original
    ! constraints.
    scale = .false.
  
    ! extallowed = .true. means that you allow Gencan (the active-set
    ! method used by Algencan to solve the bound-constrained
    ! subproblems) to perform extrapolations. This strategy may use
    ! extra evaluations of the objective function and the constraints
    ! per iterations; but it uses to provide overal savings. You should
    ! test both choices for the problem at hand.
    extallowed = .true.
  
    ! extallowed = .true. means that you allow the inertia of the
    ! Jacobian of the KKT system to be corrected during the acceleration
    ! process. You should test both choices for the problem at hand.
    corrin = .false.
  
    ! call cpu_time(start)
  
    ! call algencan(evalf,evalg,evalc,evalj,evalhl,jnnzmax,hlnnzmax, &
    !      n,x,lind,lbnd,uind,ubnd,m,p,lambda,epsfeas,epscompl,epsopt,maxoutit, &
    !      scale,rhoauto,rhoini,extallowed,corrin,f,csupn,ssupn,nlpsupn,bdsvio, &
    !      outiter,totiter,nwcalls,nwtotit,ierr,istop,c_loc(pdata))
  
    ! call cpu_time(finish)
  
    ! write(*,*)
    ! write(*,*) 'Number of variables                                   = ',n
    ! write(*,*) 'Number of equality constraints                        = ',m
    ! write(*,*) 'Number of inequality constraints                      = ',p
    
    ! write(*,*)
    ! write(*,*) '(REPORTED BY SOLVER) istop                            = ',istop
    ! write(*,*) '(REPORTED BY SOLVER) ierr                             = ',ierr
    ! write(*,*) '(REPORTED BY SOLVER) f                                = ',f
    ! write(*,*) '(REPORTED BY SOLVER) csupn                            = ',csupn
    ! write(*,*) '(REPORTED BY SOLVER) ssupn                            = ',ssupn
    ! write(*,*) '(REPORTED BY SOLVER) nlpsupn                          = ',nlpsupn
    ! write(*,*) '(REPORTED BY SOLVER) bounds violation                 = ',bdsvio
    ! write(*,*) '(REPORTED BY SOLVER) Number of outer iterations       = ',outiter
    ! write(*,*) '(REPORTED BY SOLVER) Number of inner iterations       = ',totiter
    ! write(*,*) '(REPORTED BY SOLVER) Number of Newton-KKT trials      = ',nwcalls
    ! write(*,*) '(REPORTED BY SOLVER) Number of Newton-KKT iterations  = ',nwtotit
    
    ! write(*,*)
    ! write(*,*) '(COMPUTED BY CALLER) Number of calls to evalf         = ',pdata%counters(1)
    ! write(*,*) '(COMPUTED BY CALLER) Number of calls to evalg         = ',pdata%counters(2)
    ! write(*,*) '(COMPUTED BY CALLER) Number of calls to evalc         = ',pdata%counters(3)
    ! write(*,*) '(COMPUTED BY CALLER) Number of calls to evalj         = ',pdata%counters(4)
    ! write(*,*) '(COMPUTED BY CALLER) Number of calls to evalhl        = ',pdata%counters(5)
    ! write(*,*) '(COMPUTED BY CALLER) CPU time in seconds              = ',finish - start
  
    ! ! *****************************************************************
    ! ! *****************************************************************
    ! ! Just checking ...
  
    ! inform = 0
    
    ! call evalf(n,x,f,inform,c_loc(pdata))
    
    ! if ( inform .ne. 0 ) then
    !    write(*,*) 'error when calling evalf in the main file. '
    !    stop
    ! end if
  
    ! call evalc(n,x,m,p,c,inform,c_loc(pdata))
  
    ! if ( inform .ne. 0 ) then
    !    write(*,*) 'error when calling evalc in the main file. '
    !    stop
    ! end if
  
    ! csupn = max( 0.0d0, max( maxval( abs( c(1:m) ) ), maxval( c(m+1:m+p) ) ) )
  
    ! bdsvio = max( 0.0d0, max( maxval( lbnd(1:n) - x(1:n), lind(1:n) ), maxval( x(1:n) - ubnd(1:n), uind(1:n) ) ) )
  
    ! write(*,*)
    ! write(*,*) '(COMPUTED BY CALLER) f                                = ',f
    ! write(*,*) '(COMPUTED BY CALLER) csupn                            = ',csupn
    ! write(*,*) '(COMPUTED BY CALLER) bounds violation                 = ',bdsvio
  
    ! write(*,*)
    ! write(*,*) 'When a quantity appears as computed by solver and computed by caller, they must coincide.'
    ! write(*,*) '(In case they do not coincide, please report it as a bug.)'
    ! ! *****************************************************************
    ! ! *****************************************************************
    
    ! deallocate(lind,lbnd,uind,ubnd,x,lambda,c,stat=allocerr)
    ind_train = 100
    ! if ( allocerr .ne. 0 ) then
    !    write(*,*) 'Deallocation error.'
    !    stop
    ! end if


    stop
  
    contains  

    ! *****************************************************************
    ! *****************************************************************

    subroutine train_test_split()

        implicit none
  
        integer :: i,j,k

        ! Mounting train and validation matrices
        do i = 1, rows_train

            k = i
            do j = 1, samples_train
                train(i,j) = y(k)
                k = k + 1
            enddo

            k = i
            do j = 1, samples_validation
                validation(i,j) = y(k + samples_train)
                k = k + 1
            enddo
        enddo

    end subroutine train_test_split

    ! *****************************************************************
    ! *****************************************************************

    subroutine mount_Imin(x,n,Fmin,ind_train,combi,Imin,n_Imin)

        implicit none

        integer,        intent(in) :: ind_train,n
        real(kind=8),   intent(in) :: Fmin,x(n)
        integer,        intent(inout) :: combi(order_lovo)
        integer,        intent(out) :: n_Imin,Imin(r_comb)
        integer :: i,j
        real(kind=8) :: F_i,Fi_aux

        n_Imin = 0

        do i = 1, r_comb
            call comb_unrank(samples_train,order_lovo,i,combi)
            
            F_i = 0.0d0
            do j = 1, order_lovo
                call fi(x,n,combi(j),ind_train,Fi_aux)
                F_i = F_i + Fi_aux
            enddo

            if (F_i .eq. Fmin) then
                Imin(n_Imin + 1) = i
                n_Imin = n_Imin + 1
            endif
        enddo

    end subroutine mount_Imin

    ! *****************************************************************
    ! *****************************************************************

    subroutine compute_grad_Fi(x,n,C_i,combi,grad_Fi)
        
        implicit none

        integer,        intent(in) :: C_i,n
        real(kind=8),   intent(in) :: x(n)
        integer,        intent(inout) :: combi(order_lovo)
        real(kind=8),   intent(out) :: grad_Fi(n)
        real(kind=8) :: zi



    end subroutine compute_grad_Fi

    ! *****************************************************************
    ! *****************************************************************

    subroutine compute_Fmin(x,n,ind_train,Fmin_aux,fun)

        implicit none

        integer,        intent(in) :: n,ind_train
        real(kind=8),   intent(in) :: x(n)
        real(kind=8),   intent(inout) :: Fmin_aux(samples_train)
        real(kind=8),   intent(out) :: fun
        integer :: i,kflag

        Fmin_aux(:) = 0.0d0

        kflag = 2

        indices(:) = (/(i, i = 1, samples_train)/)
        
        do i = 1, samples_train
            call fi(x,n,i,ind_train,Fmin_aux(i))
        end do

        ! Sorting
        call DSORT(Fmin_aux,indices,samples_train,kflag)

        ! Lovo function 
        fun = sum(Fmin_aux(1:order_lovo))

    end subroutine compute_Fmin

    ! *****************************************************************
    ! *****************************************************************

    subroutine fi(x,n,i,ind_train,fun)

        implicit none

        integer,        intent(in) :: n,i,ind_train
        real(kind=8),   intent(in) :: x(n)
        real(kind=8),   intent(out) :: fun   
        
        call model(x,n,i,ind_train,fun)
        fun = 0.5d0 * ((fun - train(ind_train,i))**2)

    end subroutine fi

    ! *****************************************************************
    ! *****************************************************************

    subroutine model(x,n,i,ind_train,fun)

        implicit none 

        integer,        intent(in) :: n,i,ind_train
        real(kind=8),   intent(in) :: x(n)
        real(kind=8),   intent(out) :: fun

        fun = train(ind_train,samples_train) + &
            x(1) * (t(i) - t(samples_train)) + &
            x(2) * (t(i) - t(samples_train))**2 + &
            x(3) * (t(i) - t(samples_train))**3
        
    end subroutine model

    ! *****************************************************************
    ! *****************************************************************

    subroutine evalf(n,x,f,inform,pdataptr)

        implicit none

        ! SCALAR ARGUMENTS
        integer, intent(in) :: n
        integer, intent(inout) :: inform
        real(kind=8), intent(out) :: f
        type(c_ptr), optional, intent(in) :: pdataptr

        ! ARRAY ARGUMENTS
        real(kind=8), intent(in) :: x(n)

        ! This routine must compute the objective function.
        
        ! LOCAL SCALARS
        type(pdata_type), pointer :: pdata
        
        call c_f_pointer(pdataptr,pdata)
        pdata%counters(1) = pdata%counters(1) + 1
        
        f = ( x(1) + 4.0d0 ) ** 4.0d0 + x(2) ** 2.0d0
        
    end subroutine evalf

    ! *****************************************************************
    ! *****************************************************************

    subroutine evalg(n,x,g,inform,pdataptr)

        implicit none

        ! SCALAR ARGUMENTS
        integer, intent(in) :: n
        integer, intent(inout) :: inform
        type(c_ptr), optional, intent(in) :: pdataptr

        ! ARRAY ARGUMENTS
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(out) :: g(n)
        
        ! This routine must compute the gradient of the objective
        ! function.
        
        ! LOCAL SCALARS
        type(pdata_type), pointer :: pdata
        
        call c_f_pointer(pdataptr,pdata)
        pdata%counters(2) = pdata%counters(2) + 1
        
        g(1) = 4.0d0 * ( x(1) + 4.0d0 ) ** 3.0d0
        g(2) = 2.0d0 * x(2)

    end subroutine evalg

    ! *****************************************************************
    ! *****************************************************************

    subroutine evalc(n,x,m,p,c,inform,pdataptr)

        implicit none

        ! SCALAR ARGUMENTS
        integer, intent(in) :: m,n,p
        integer, intent(inout) :: inform
        type(c_ptr), optional, intent(in) :: pdataptr

        ! ARRAY ARGUMENTS
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(out) :: c(m+p)

        ! This routine must compute all the m+p constraints.
        
        ! LOCAL SCALARS
        type(pdata_type), pointer :: pdata
        
        call c_f_pointer(pdataptr,pdata)
        pdata%counters(3) = pdata%counters(3) + 1
        
        c(1) = x(1) ** 3.0d0 - x(2) - 1.0d0
        
    end subroutine evalc

    ! *****************************************************************
    ! *****************************************************************

    subroutine evalj(n,x,m,p,ind,sorted,jsta,jlen,lim,jvar,jval,inform,pdataptr)

        implicit none
        
        ! SCALAR ARGUMENTS
        integer, intent(in) :: lim,m,n,p
        integer, intent(inout) :: inform
        type(c_ptr), optional, intent(in) :: pdataptr

        ! ARRAY ARGUMENTS
        logical, intent(in) :: ind(m+p)
        real(kind=8), intent(in) :: x(n)
        logical, intent(out) :: sorted(m+p)
        integer, intent(out) :: jsta(m+p),jlen(m+p),jvar(lim)
        real(kind=8), intent(out) :: jval(lim)
        
        ! This routine must compute the Jacobian of the constraints. In
        ! fact, only gradients of constraints j such that ind(j) =
        ! .true. need to be computed.
        
        ! LOCAL SCALARS
        integer :: i
        type(pdata_type), pointer :: pdata
        
        call c_f_pointer(pdataptr,pdata)
        pdata%counters(4) = pdata%counters(4) + 1

        ! Only gradients of constraints j such that ind(j) = .true. need
        ! to be computed.
        
        if ( ind(1) ) then
            if ( lim .lt. n ) then
            inform = -94
            return
            end if
            
            jsta(1) = 1
            jlen(1) = n
            
            jvar(1:n) = (/ (i,i=1,n) /)

            jval(1) = 3.0d0 * x(1) ** 2.0d0
            jval(2) = - 1.0d0

            ! Says whether the variables' indices in jvar (related to this
            ! constraint) are in increasing order. In case they are,
            ! Algencan takes advantage of this. Implement sorted gradients
            ! of constraints if you can do this in a natural (cheap)
            ! way. Under no circumnstance use a sorting algorithm. (It is
            ! better to set sorted(1) = .false. in this case.)
            
            sorted(1) = .true.
        end if
        
    end subroutine evalj

    ! *****************************************************************
    ! *****************************************************************

    subroutine evalhl(n,x,m,p,lambda,lim,inclf,hlnnz,hlrow,hlcol,hlval,inform,pdataptr)

        implicit none
        
        ! SCALAR ARGUMENTS
        logical, intent(in) :: inclf
        integer, intent(in) :: m,n,lim,p
        integer, intent(out) :: hlnnz
        integer, intent(inout) :: inform
        type(c_ptr), optional, intent(in) :: pdataptr

        ! ARRAY ARGUMENTS
        real(kind=8), intent(in) :: lambda(m+p),x(n)
        integer, intent(out) :: hlrow(lim),hlcol(lim)
        real(kind=8), intent(out) :: hlval(lim)

        ! This routine must compute the Hessian of the Lagrangian. The
        ! Hessian of the objective function must NOT be included if inclf
        ! = .false.
        
        ! LOCAL SCALARS
        type(pdata_type), pointer :: pdata
        
        call c_f_pointer(pdataptr,pdata)
        pdata%counters(5) = pdata%counters(5) + 1

        hlnnz = 0

        ! If .not. inclf then the Hessian of the objective function must not be included
        
        if ( inclf ) then
            if ( hlnnz + 2 .gt. lim ) then
            inform = -95
            return
            end if
        
            hlnnz = hlnnz + 1
            
            hlrow(hlnnz) = 1
            hlcol(hlnnz) = 1
            hlval(hlnnz) = 12.0d0 * ( x(1) + 4.0d0 ) ** 2.0d0
        
            hlnnz = hlnnz + 1
            
            hlrow(hlnnz) = 2
            hlcol(hlnnz) = 2
            hlval(hlnnz) = 2.0d0
        end if

        ! Note that entries of the Hessian of the Lagrangian can be
        ! repeated. If this is case, them sum of repeated entrances is
        ! considered. This feature simplifies the construction of the
        ! Hessian of the Lagrangian.
        
        if ( hlnnz + 1 .gt. lim ) then
            inform = -95
            return
        end if
        
        hlnnz = hlnnz + 1
        
        hlrow(hlnnz) = 1
        hlcol(hlnnz) = 1
        hlval(hlnnz) = lambda(1) * 6.0d0 * x(1)
        
    end subroutine evalhl

    end program lovo