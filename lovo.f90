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

    integer :: samples,samples_train,samples_validation,r_comb,order_lovo,rows_train,i,ind_train,n_Imin,i4_choose,nuk
    integer, allocatable :: Imin(:),combi(:)
    real(kind=8), allocatable :: xtrial(:),xk(:),t(:),y(:),train(:,:),validation(:,:),Fmin_aux(:),&
                                 indices(:),grad_Fi(:),hess_Fi(:,:)
    real(kind=8) :: Fmin,sigma
    real(kind=8), parameter :: sigmin = 1.d0

    !--> End LOVO Algorithm variables <--

    ! Reading data and storing it in the variables t and y
    Open(Unit = 100, File = "output/data.txt", ACCESS = "SEQUENTIAL")

    read(100,*) samples

    samples_train = 10
    samples_validation = 30
    rows_train = samples - (samples_train + samples_validation) + 1
    order_lovo = samples_train - 2
    r_comb = i4_choose(samples_train,order_lovo)

    allocate(t(samples_train),y(samples),Fmin_aux(samples),indices(samples),&
            Imin(r_comb),combi(order_lovo),train(rows_train,samples_train),&
            validation(rows_train,samples_validation),stat=allocerr)

    if ( allocerr .ne. 0 ) then
        write(*,*) 'Allocation error in main program'
        stop
    end if

    do i = 1, samples
        read(100,*) y(i)
    enddo

    t(:) = (/(i, i = 1, samples_train)/)

    close(100)

    ! Number of variables
  
    n = 3

    allocate(xtrial(n),xk(n),grad_Fi(n),hess_Fi(n,n),stat=allocerr)

    if ( allocerr .ne. 0 ) then
        write(*,*) 'Allocation error.'
        stop
     end if
    
    allocate(x(n),lind(n),lbnd(n),uind(n),ubnd(n),stat=allocerr)
  
    if ( allocerr .ne. 0 ) then
       write(*,*) 'Allocation error.'
       stop
    end if
  
    ! Initial guess and bound constraints
  
    lind(1:n) = .false.
    ! lbnd(1:n) = - 100.0d0
  
    uind(1:n) = .false.
    ! ubnd(1:n) =   100.0d0
  
    ! Number equality (m) and inequality (p) constraints.
    
    m = 0
    p = 0

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

    call train_test_split()
    
    allocate(lambda(m+p),c(m+p),stat=allocerr)
  
    if ( allocerr .ne. 0 ) then
       write(*,*) 'Allocation error.'
       stop
    end if

    ind_train = 1

    call cpu_time(start)
    call lovo_algorithm(ind_train)
    call cpu_time(finish)

    print '("Time = ",f10.4," seconds.")',finish-start

    deallocate(lind,lbnd,uind,ubnd,x,lambda,c,stat=allocerr)

    if ( allocerr .ne. 0 ) then
        write(*,*) 'Deallocation error.'
        stop
    end if

    stop
  
    contains  

    ! *****************************************************************
    ! *****************************************************************

    subroutine lovo_algorithm(ind_train)

        implicit none

        integer, intent(in) :: ind_train
        integer, parameter :: max_iter=1000,max_iter_sub=100
        real(kind=8), parameter :: theta=1.0d0,tol=1.0d-6
        real(kind=8) :: fxk, fxtrial
        integer :: iter,iter_sub
    
        xk(:) = 1.0d0

        iter = 0

        call compute_Fmin(xk,n,ind_train,Fmin)

        call mount_Imin(xk,n,Fmin,ind_train,combi,Imin,n_Imin)

        nuk = Imin(n_Imin)

        x(1:n) = xk(1:n)

        sigma = sigmin

        call algencan(evalf,evalg,evalc,evalj,evalhl,jnnzmax,hlnnzmax, &
            n,x,lind,lbnd,uind,ubnd,m,p,lambda,epsfeas,epscompl,epsopt,maxoutit, &
            scale,rhoauto,rhoini,extallowed,corrin,f,csupn,ssupn,nlpsupn,bdsvio, &
            outiter,totiter,nwcalls,nwtotit,ierr,istop,c_loc(pdata))

        xtrial(1:n) = x(1:n)

        print*, xtrial


    end subroutine lovo_algorithm


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

    subroutine regularized_Taylor(x,n,ind_train,nuk,sigma,res)

        implicit none

        integer,        intent(in) :: n,nuk,ind_train
        real(kind=8),   intent(in) :: x(n),sigma
        real(kind=8),   intent(out) :: res

        call compute_Fmin(x,n,ind_train,res)
        call compute_grad_Fi(x,n,nuk,grad_Fi)
        res = res + dot_product(grad_Fi,x(1:n) - xk(1:n))
        res = res + 0.5d0 * sigma * (norm2(x(1:n) - xk(1:n))**2)

    end subroutine regularized_Taylor

    ! *****************************************************************
    ! *****************************************************************

    subroutine grad_regularized_Taylor(x,n,ind_train,nuk,sigma,res)

        implicit none

        integer,        intent(in) :: n,nuk,ind_train
        real(kind=8),   intent(in) :: x(n),sigma
        real(kind=8),   intent(out) :: res(n)

        call compute_grad_Fi(x,n,nuk,res)
        res = res + sigma * (x(1:n) - xk(1:n))

    end subroutine grad_regularized_Taylor

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

    subroutine compute_grad_Fi(x,n,ind_Ci,res)
        
        implicit none

        integer,        intent(in) :: ind_Ci,n
        real(kind=8),   intent(in) :: x(n)
        real(kind=8),   intent(out) :: res(n)
        real(kind=8) :: zi
        integer :: i,j

        combi(:) = 0

        call comb_unrank(samples_train,order_lovo,ind_Ci,combi)
        
        zi = 0.0d0
        res(:) = 0.0d0

        do i = 1, order_lovo
            call fit_model(x,n,combi(i),ind_train,zi)
            zi = zi - train(ind_train,combi(i))

            do j = 1, n
                res(j) = res(j) + zi * ((t(combi(i)) - t(samples_train))**j)
            enddo

        enddo

    end subroutine compute_grad_Fi

    ! *****************************************************************
    ! *****************************************************************

    subroutine compute_hess_Fi(n,ind_Ci,combi,res)

        implicit none

        integer,        intent(in) :: ind_Ci,n
        integer,        intent(inout) :: combi(order_lovo)
        real(kind=8),   intent(out) :: res(n,n)
        integer :: i,j,k

        call comb_unrank(samples_train,order_lovo,ind_Ci,combi)

        res(:,:) = 0.0d0

        do k = 1, order_lovo
            do i = 1, n
                do j = 1, n
                    res(i,j) = res(i,j) + (t(combi(k)) - t(samples_train))**(i + j)
                enddo
            enddo
        enddo

    end subroutine compute_hess_Fi

    ! *****************************************************************
    ! *****************************************************************

    subroutine compute_Fmin(x,n,ind_train,res)

        implicit none

        integer,        intent(in) :: n,ind_train
        real(kind=8),   intent(in) :: x(n)
        real(kind=8),   intent(out) :: res
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
        res = sum(Fmin_aux(1:order_lovo))

    end subroutine compute_Fmin

    ! *****************************************************************
    ! *****************************************************************

    subroutine fi(x,n,i,ind_train,fun)

        implicit none

        integer,        intent(in) :: n,i,ind_train
        real(kind=8),   intent(in) :: x(n)
        real(kind=8),   intent(out) :: fun   
        
        call fit_model(x,n,i,ind_train,fun)
        fun = 0.5d0 * ((fun - train(ind_train,i))**2)

    end subroutine fi

    ! *****************************************************************
    ! *****************************************************************

    subroutine fit_model(x,n,i,ind_train,res)

        implicit none 

        integer,        intent(in) :: n,i,ind_train
        real(kind=8),   intent(in) :: x(n)
        real(kind=8),   intent(out) :: res

        res = train(ind_train,samples_train) + &
            x(1) * (t(i) - t(samples_train)) + &
            x(2) * (t(i) - t(samples_train))**2 + &
            x(3) * (t(i) - t(samples_train))**3
        
    end subroutine fit_model

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
        
        call regularized_Taylor(x,n,ind_train,nuk,sigma,f)
        
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
        
        call grad_regularized_Taylor(x,n,ind_train,nuk,sigma,g)

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
        
        ! c(1) = x(1) ** 3.0d0 - x(2) - 1.0d0
        
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

        integer :: i

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
            hlnnz = n
            if ( hlnnz .gt. lim ) then
            inform = -95
            return
            end if
            
            hlrow(1:n) = (/(i, i = 1, n)/)
            hlcol(1:n) = (/(i, i = 1, n)/)
            hlval(1:n) = sigma
    
        end if

        ! Note that entries of the Hessian of the Lagrangian can be
        ! repeated. If this is case, them sum of repeated entrances is
        ! considered. This feature simplifies the construction of the
        ! Hessian of the Lagrangian.
        
        if ( hlnnz .gt. lim ) then
            inform = -95
            return
        end if
        
        ! hlnnz = hlnnz + 1
        
        ! hlrow(hlnnz) = 1
        ! hlcol(hlnnz) = 1
        ! hlval(hlnnz) = lambda(1) * 6.0d0 * x(1)
        
    end subroutine evalhl

    end program lovo