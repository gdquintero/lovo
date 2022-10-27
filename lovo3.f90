Program lovo
    use sort

    implicit none 
    
    ! LOCAL SCALARS
    logical :: checkder
    integer :: hnnzmax,inform,jcnnzmax,m,n,nvparam,allocerr
    real(kind=8) :: cnorm,efacc,efstain,eoacc,eostain,epsfeas,epsopt,f,nlpsupn,snorm

    ! LOCAL ARRAYS
    character(len=80) :: specfnm,outputfnm,vparam(10)
    logical :: coded(11)
    logical, pointer :: equatn(:),linear(:)
    real(kind=8), pointer :: l(:),lambda(:),u(:),x(:)

    !--> LOVO Algorithm variables <--

    integer :: samples,samples_train,samples_validation,r_comb,order_lovo,rows_train,i,ind_train,n_Imin,i4_choose,nuk
    integer, allocatable :: Imin(:),combi(:),t(:)
    real(kind=8), allocatable :: xtrial(:),xk(:),y(:),train(:,:),validation(:,:),Fmin_aux(:),&
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

    allocate(t(samples_train),y(samples),Fmin_aux(samples_train),indices(samples_train),&
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

    call train_test_split()

    ! Set parameters
    n = 3
    m = 0

    allocate(x(n),l(n),u(n),equatn(m),linear(m),lambda(m),&
            xtrial(n),xk(n),grad_Fi(n),hess_Fi(n,n),stat=allocerr)

    if ( allocerr .ne. 0 ) then
        write(*,*) 'Allocation error.'
        stop
    end if

    ! stop

    l(:) = - 1.0d-20
    u(:) = 1.0d+20

    ! Coded subroutines

    coded(1:6)  = .true.  ! evalf, evalg, evalh, evalc, evaljac, evalhc
    coded(7:11) = .false. ! evalfc,evalgjac,evalgjacp,evalhl,evalhlp


    ! Upper bounds on the number of sparse-matrices non-null elements
    jcnnzmax = 10000
    hnnzmax  = 10000

    ! Checking derivatives?
    checkder = .true.

    ! Parameters setting
    epsfeas   = 1.0d-08
    epsopt    = 1.0d-08
  
    efstain   = sqrt( epsfeas )
    eostain   = epsopt ** 1.5d0
  
    efacc     = sqrt( epsfeas )
    eoacc     = sqrt( epsopt )

    outputfnm = ''
    specfnm   = ''

    nvparam   = 0
    vparam(1) = 'ITERATIONS-OUTPUT-DETAIL 0' 

    ind_train = 1

    call lovo_algorithm(ind_train)

    CONTAINS

    ! *****************************************************************
    ! MAIN ALGORITHM (WEAKLY CRITICAL POINTS)
    ! *****************************************************************

    subroutine lovo_algorithm(ind_train)

        implicit none

        integer, intent(in) :: ind_train
        real(kind=8) :: fxk,fxtrial,tol,theta
        integer :: iter,iter_sub,max_iter,max_iter_sub,i

        tol = 1.0d-4
        max_iter = 1
        max_iter_sub = 1
    
        xk(:) = 1.0d0

        iter = 0

        ! call compute_Fmin(xk,n,Fmin)

        ! call mount_Imin(xk,n,Fmin,combi,Imin,n_Imin)

        ! fxk = Fmin

        ! nuk = Imin(n_Imin)

        ! call grad_regularized_Taylor(xk,n,nuk,sigma,grad_Fi)

        ! print*, norm2(grad_Fi), ind_train

        stop

        do
            iter = iter + 1
            
            x(1:n) = xk(1:n)

            sigma = sigmin

            ! do 

            call algencan(myevalf,myevalg,myevalh,myevalc,myevaljac,myevalhc, &
            myevalfc,myevalgjac,myevalgjacp,myevalhl,myevalhlp,jcnnzmax, &
            hnnzmax,epsfeas,epsopt,efstain,eostain,efacc,eoacc,outputfnm, &
            specfnm,nvparam,vparam,n,x,l,u,m,lambda,equatn,linear,coded, &
            checkder,f,cnorm,snorm,nlpsupn,inform)

            xtrial(1:n) = x(1:n)

            ! call grad_regularized_Taylor(xtrial,n,nuk,sigma,grad_Fi)

            ! print*, norm2(grad_Fi), ind_train, Imin(1:n_Imin)

            if (iter .ge. max_iter) exit

        enddo

    end subroutine lovo_algorithm

    ! *****************************************************************
    ! DIVISION OF THE DATA SET
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
    ! REGULARIZED MODEL
    ! *****************************************************************

    subroutine regularized_Taylor(x,n,nuk,sigma,res)

        implicit none

        integer,        intent(in) :: n,nuk
        real(kind=8),   intent(in) :: x(n),sigma
        real(kind=8),   intent(out) :: res

        call compute_Fmin(x,n,res)
        call compute_grad_F_i(x,n,nuk,grad_Fi)
        res = res + dot_product(grad_Fi,x(1:n) - xk(1:n))
        res = res + 0.5d0 * sigma * (norm2(x(1:n) - xk(1:n))**2)

    end subroutine regularized_Taylor

    ! *****************************************************************
    ! REGULARIZED MODEL GRADIENT
    ! *****************************************************************

    subroutine grad_regularized_Taylor(x,n,nuk,sigma,res)

        implicit none

        integer,        intent(in) :: n,nuk
        real(kind=8),   intent(in) :: x(n),sigma
        real(kind=8),   intent(out) :: res(n)

        call compute_grad_F_i(x,n,nuk,res)
        res = res + sigma * (x(1:n) - xk(1:n))

    end subroutine grad_regularized_Taylor

    ! *****************************************************************
    ! GRADIENT OF ERROR FUNCTIONS Fi
    ! *****************************************************************

    subroutine compute_grad_F_i(x,n,ind_Ci,res)
        
        implicit none

        integer,        intent(in) :: ind_Ci,n
        real(kind=8),   intent(in) :: x(n)
        real(kind=8),   intent(out) :: res(n)
        real(kind=8) :: zi,yi
        integer :: i,j,ti,tm

        combi(:) = 0

        call comb_unrank(samples_train,order_lovo,ind_Ci,combi)
        
        zi = 0.0d0
        res(:) = 0.0d0

        tm = t(samples_train)

        do i = 1, order_lovo
            ti = t(combi(i))
            yi = train(ind_train,combi(i))

            call fit_model(x,n,combi(i),zi)
            zi = zi - yi

            do j = 1, n
                res(j) = res(j) + zi * ((ti - tm)**j)
            enddo

        enddo

        combi(:) = 0

    end subroutine compute_grad_F_i

    ! *****************************************************************
    ! HESSIAN OF ERROR FUNCTIONS
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
    ! GRADIENT OF ERROR FUNCTIONS
    ! *****************************************************************

    subroutine compute_grad_fi(x,n,ind,res)
        
        implicit none

        integer,        intent(in) :: ind,n
        real(kind=8),   intent(in) :: x(n)
        real(kind=8),   intent(out) :: res(n)
        real(kind=8) :: zi,yi
        integer :: ti,tm

        zi = 0.0d0
        ti = t(ind)
        tm = t(samples_train)
        yi = train(ind_train,ind)
        res(:) = 0.0d0

        call fit_model(x,n,ind,zi)

        zi = zi - yi

        res(1) = ti - tm
        res(2) = (ti - tm)**2
        res(3) = (ti - tm)**3

        res(1:n) = zi * res(1:n)

    end subroutine compute_grad_fi

    ! *****************************************************************
    ! I_min(x)
    ! *****************************************************************

    subroutine mount_Imin(x,n,Fmin,combi,Imin,n_Imin)

        implicit none

        integer,        intent(in) :: n
        real(kind=8),   intent(in) :: Fmin,x(n)
        integer,        intent(inout) :: combi(order_lovo)
        integer,        intent(out) :: n_Imin,Imin(r_comb)
        integer :: i,j
        real(kind=8) :: F_i,Fi_aux

        n_Imin = 0
        combi(:) = 0

        do i = 1, r_comb
            call comb_unrank(samples_train,order_lovo,i,combi)
            
            F_i = 0.0d0
            do j = 1, order_lovo
                Fi_aux = 0.0d0
                call fi(x,n,combi(j),Fi_aux)
                F_i = F_i + Fi_aux
            enddo

            if (abs(F_i - Fmin) .le. 1.0d-6) then
                Imin(n_Imin + 1) = i
                n_Imin = n_Imin + 1
            endif
            combi(:) = 0
        enddo

    end subroutine mount_Imin

    ! *****************************************************************
    ! LOVO FUNCTION
    ! *****************************************************************

    subroutine compute_Fmin(x,n,res)

        implicit none

        integer,        intent(in) :: n
        real(kind=8),   intent(in) :: x(n)
        real(kind=8),   intent(out) :: res
        integer :: i,kflag

        Fmin_aux(:) = 0.0d0

        kflag = 2

        indices(:) = (/(i, i = 1, samples_train)/)
        
        do i = 1, samples_train
            call fi(x,n,i,Fmin_aux(i))
        end do

        ! Sorting
        call DSORT(Fmin_aux,indices,samples_train,kflag)

        ! Lovo function 
        res = sum(Fmin_aux(1:order_lovo))

    end subroutine compute_Fmin

    ! *****************************************************************
    ! QUADRATIC ERROR FUNCTIONS
    ! *****************************************************************

    subroutine fi(x,n,i,fun)

        implicit none

        integer,        intent(in) :: n,i
        real(kind=8),   intent(in) :: x(n)
        real(kind=8),   intent(out) :: fun 
        real(kind=8) :: yi

        yi = train(ind_train,i)
        
        call fit_model(x,n,i,fun)
        fun = 0.5d0 * ((fun - yi)**2)

    end subroutine fi

    ! *****************************************************************
    ! MODEL TO BE ADJUSTED
    ! *****************************************************************

    subroutine fit_model(x,n,i,res)

        implicit none 

        integer,        intent(in) :: n,i
        real(kind=8),   intent(in) :: x(n)
        real(kind=8),   intent(out) :: res
        real(kind=8) :: ti,tm,ym

        ym = train(ind_train,samples_train)
        ti = t(i)
        tm = t(samples_train)

        res = x(1) * (ti - tm)
        res = res + x(2) * ((ti - tm)**2)
        res = res + x(3) * ((ti - tm)**3)
        res = res + ym
        
    end subroutine fit_model

    !==============================================================================
    ! SUBROUTINES FOR ALGENCAN
    !==============================================================================

    !******************************************************************************
    ! OBJECTIVE FUNCTION
    !******************************************************************************
    subroutine myevalf(n,x,f,flag)
        implicit none

        ! SCALAR ARGUMENTS
        integer, intent(in) :: n
        integer, intent(out) :: flag
        real(kind=8), intent(out) :: f

        ! ARRAY ARGUMENTS
        real(kind=8), intent(in) :: x(n)

        ! Compute objective function

        flag = 0

        call regularized_Taylor(x,n,nuk,sigma,f)

    end subroutine myevalf

    !******************************************************************************
    ! GRADIENT OF THE OBJECTIVE FUNCTION
    !******************************************************************************
    subroutine myevalg(n,x,g,flag)
        implicit none

        ! SCALAR ARGUMENTS
        integer, intent(in) :: n
        integer, intent(out) :: flag

        ! ARRAY ARGUMENTS
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(out) :: g(n)

        ! Compute gradient of the objective function

        flag = 0

        call grad_regularized_Taylor(x,n,nuk,sigma,g)

    end subroutine myevalg

    !******************************************************************************
    ! HESSIAN FOR THE OBJECTIVE FUNCTION
    !******************************************************************************
    subroutine myevalh(n,x,hrow,hcol,hval,hnnz,lim,lmem,flag)
        implicit none

        ! SCALAR ARGUMENTS
        logical, intent(out) :: lmem
        integer, intent(in) :: lim,n
        integer, intent(out) :: flag,hnnz

        ! ARRAY ARGUMENTS
        integer, intent(out) :: hcol(lim),hrow(lim)
        real(kind=8), intent(in)  :: x(n)
        real(kind=8), intent(out) :: hval(lim)

        integer :: i

        ! Compute (lower triangle of the) Hessian of the objective function
        flag = 0
        lmem = .false.
        hnnz = n

        hrow(1:n) = (/(i, i = 1, n)/)
        hcol(1:n) = (/(i, i = 1, n)/)
        hval(1:n) = sigma
    end subroutine myevalh

    !******************************************************************************
    ! CONSTRAINTS
    !******************************************************************************
    subroutine myevalc(n,x,ind,c,flag)
        implicit none

        ! SCALAR ARGUMENTS
        integer, intent(in) :: ind,n
        integer, intent(out) :: flag
        real(kind=8), intent(out) :: c

        ! ARRAY ARGUMENTS
        real(kind=8), intent(in) :: x(n)

        ! Compute ind-th constraint
        flag = -1

    end subroutine myevalc

    !******************************************************************************
    ! JACOBIAN OF THE CONSTRAINTS
    !******************************************************************************
    subroutine myevaljac(n,x,ind,jcvar,jcval,jcnnz,lim,lmem,flag)

        implicit none

        ! SCALAR ARGUMENTS
        logical, intent(out) :: lmem
        integer, intent(in) :: ind,lim,n
        integer, intent(out) :: flag,jcnnz

        ! ARRAY ARGUMENTS
        integer, intent(out) :: jcvar(lim)
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(out) :: jcval(lim)

        flag = -1

    end subroutine myevaljac

    !******************************************************************************
    ! HESSIAN OF THE CONSTRAINTS
    !******************************************************************************
    subroutine myevalhc(n,x,ind,hcrow,hccol,hcval,hcnnz,lim,lmem,flag)

        implicit none

        ! SCALAR ARGUMENTS
        logical, intent(out) :: lmem
        integer, intent(in) :: ind,lim,n
        integer, intent(out) :: flag,hcnnz

        ! ARRAY ARGUMENTS
        integer, intent(out) :: hccol(lim),hcrow(lim)
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(out) :: hcval(lim)

        flag = -1

    end subroutine myevalhc

    ! ******************************************************************
    ! ******************************************************************

    subroutine myevalfc(n,x,f,m,c,flag)

        implicit none

        ! SCALAR ARGUMENTS
        integer, intent(in) :: m,n
        integer, intent(out) :: flag
        real(kind=8), intent(out) :: f

        ! ARRAY ARGUMENTS
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(out) :: c(m)

        flag = - 1

    end subroutine myevalfc

    ! ******************************************************************
    ! ******************************************************************

    subroutine myevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,lim,lmem,flag)

        implicit none

        ! SCALAR ARGUMENTS
        logical, intent(out) :: lmem
        integer, intent(in) :: lim,m,n
        integer, intent(out) :: flag,jcnnz

        ! ARRAY ARGUMENTS
        integer, intent(out) :: jcfun(lim),jcvar(lim)
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(out) :: g(n),jcval(lim)

        flag = - 1

    end subroutine myevalgjac

    ! ******************************************************************
    ! ******************************************************************

    subroutine myevalgjacp(n,x,g,m,p,q,work,gotj,flag)

        implicit none

        ! SCALAR ARGUMENTS
        logical, intent(inout) :: gotj
        integer, intent(in) :: m,n
        integer, intent(out) :: flag
        character, intent(in) :: work

        ! ARRAY ARGUMENTS
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(inout) :: p(m),q(n)
        real(kind=8), intent(out) :: g(n)

        flag = - 1

    end subroutine myevalgjacp

    ! ******************************************************************
    ! ******************************************************************

    subroutine myevalhl(n,x,m,lambda,sf,sc,hlrow,hlcol,hlval,hlnnz,lim,lmem,flag)

        implicit none

        ! SCALAR ARGUMENTS
        logical, intent(out) :: lmem
        integer, intent(in) :: lim,m,n
        integer, intent(out) :: flag,hlnnz
        real(kind=8), intent(in) :: sf

        ! ARRAY ARGUMENTS
        integer, intent(out) :: hlcol(lim),hlrow(lim)
        real(kind=8), intent(in) :: lambda(m),sc(m),x(n)
        real(kind=8), intent(out) :: hlval(lim)

        flag = - 1

    end subroutine myevalhl

    ! ******************************************************************
    ! ******************************************************************

    subroutine myevalhlp(n,x,m,lambda,sf,sc,p,hp,goth,flag)

        implicit none

        ! SCALAR ARGUMENTS
        logical, intent(inout) :: goth
        integer, intent(in) :: m,n
        integer, intent(out) :: flag
        real(kind=8), intent(in) :: sf

        ! ARRAY ARGUMENTS
        real(kind=8), intent(in) :: lambda(m),p(n),sc(m),x(n)
        real(kind=8), intent(out) :: hp(n)

        flag = - 1

    end subroutine myevalhlp
end Program lovo
