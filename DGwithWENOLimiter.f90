! DG_2D_New 重新开始的DG
 
    program main

    implicit none

    ! 变量

    real pi,CFL,x,y,t,xa,xb,ya,yb,dt,tend,gamma,hx,hy,hx1,hy1
    integer i,j,k,n,Nx,Ny,RKorder,dim,d,frame,skiptime,dimpol,i1,j1
    integer Qlength,countstep,count,skip,quad
    real stoptime1,stoptime2,stoptime3

    parameter(dim = 4)
    parameter(dimpol = 6)
    parameter(quad = 4)
    parameter(pi = 4*atan(1.0d0))
    parameter(CFL = 0.05d0)
    parameter(Nx = 40)
    parameter(Ny = 40)
    parameter(xa = 0d0)
    parameter(xb = 4d0)
    parameter(ya = 0d0)
    parameter(yb = 4d0)
    parameter(hx = (xb - xa)/Nx)
    parameter(hy = (yb - ya)/Ny)
    parameter(hx1 = 0.5d0*hx)
    parameter(hy1 = 0.5d0*hy)
    parameter(tend = 1.5d0/pi)
    parameter(RKorder = 3)
    parameter(gamma = 1.4d0)
    parameter(frame = 200)
    parameter(Qlength = Nx*Ny*frame)
    parameter(skiptime = 1)
    parameter(stoptime1 = 0.5d0)
    parameter(stoptime2 = 1.0d0)
    parameter(stoptime3 = 1.5d0)
    real XX(Nx + 1),YY(Ny + 1),Xc(Nx),Yc(Ny),Q(Nx,Ny,dim),QR(Nx,Ny,dim),Linfty,Linfty1,L2,L2V(dim)
    real TT(frame),QF(Nx,Ny,dim,frame),Error(Nx,Ny)
    real QQ1(1,Qlength),QQ2(1,Qlength),QQ3(1,Qlength),QQ4(1,Qlength),alphax,alphay,alpha1,alpha2
    real Pol(Nx,Ny,dimpol,dim),PolValue(quad,quad,dimpol,dim),lambda(quad),weight(quad),lambdai(quad),lambdaj(quad)
    real PolDx(quad,quad,dimpol,dim),PolDy(quad,quad,dimpol,dim),Poledge(4,quad,dimpol,dim),diag(dimpol),diagpol(Nx,Ny,dimpol,dim)
    real PolDxx(quad,quad,dimpol,dim),PolDxy(quad,quad,dimpol,dim),PolDyy(quad,quad,dimpol,dim),weight2(quad,quad)
    real QRquad(quad*Nx,quad*Ny),Qquad(quad*Nx,quad*Ny),Xcquad(quad*Nx),Ycquad(quad*Ny),Errorquad(quad*Nx,quad*Ny)
    real,external :: heaviside,g
    integer :: countstop = 0
    
    ! 初值
    real p
    !p(x,y) = (1 - (10*exp(1 - ((x - 5)**2 + (y - 5)**2)))/(11.2*pi**2))**(2.5*1.4) ! Vortex
    !p(x,y) = 2.5d0 ! cloud
    !p(x,y) = 1 ! sin
    p(x,y) = sin(0.5d0*pi*(x + y))
    real rho
    !rho(x,y) = (1 - (10*exp(1 - ((x - 5)**2 + (y - 5)**2)))/(11.2*pi**2))**2.5 ! Vortex
    !rho(x,y) = 1 + heaviside(0.5 - y) ! cloud
    !rho(x,y) = 1 + 0.2d0*sin(pi*(x + y)) ! sin
    rho(x,y) = sin(0.5d0*pi*(x + y))
    real v1
    !v1(x,y) = 1 + 5/(2*pi)*exp(0.5*(1 - ((x - 5)**2 + (y - 5)**2)))*(-(y - 5)) ! Vortex
    !v1(x,y) = -0.5 + heaviside(0.5 - y) ! cloud
    !v1(x,y) = 0.7d0 ! sin
    v1(x,y) = sin(0.5d0*pi*(x + y))
    real v2
    !v2(x,y) = 1 + 5/(2*pi)*exp(0.5*(1 - ((x - 5)**2 + (y - 5)**2)))*(x - 5) ! Vortex
    !v2(x,y) = 0.5d0*x*sin(4d0*pi*x)*exp(-(y - 0.5d0)**800) ! cloud
    !v2(x,y) = 0.3d0 ! sin
    v2(x,y) = sin(0.5d0*pi*(x + y))
    
    ! 守恒量
    real U1
    U1(x,y) = sin(0.5d0*pi*(x + y))
    real U2
    U2(x,y) = sin(0.5d0*pi*(x + y))
    real U3
    U3(x,y) = sin(0.5d0*pi*(x + y))
    real U4
    U4(x,y) = sin(0.5d0*pi*(x + y))
    
    ! 基函数
    real Phi1
    Phi1(x,y) = 1d0
    real Phi2
    Phi2(x,y) = x/hx1
    real Phi3
    Phi3(x,y) = y/hy1
    real Phi4
    Phi4(x,y) = (x/hx1)**2 - 1d0/3d0
    real Phi5
    Phi5(x,y) = (x*y)/(hx1*hy1)
    real Phi6
    Phi6(x,y) = (y/hy1)**2 - 1d0/3d0
    
    ! 导数
    real Phix1
    Phix1(x,y) = 0
    real Phix2
    Phix2(x,y) = 1d0/hx1
    real Phix3
    Phix3(x,y) = 0
    real Phix4
    Phix4(x,y) = 2d0*x/(hx1**2)
    real Phix5
    Phix5(x,y) = y/(hx1*hy1)
    real Phix6
    Phix6(x,y) = 0
    
    real Phiy1
    Phiy1(x,y) = 0
    real Phiy2
    Phiy2(x,y) = 0
    real Phiy3
    Phiy3(x,y) = 1d0/hy1
    real Phiy4
    Phiy4(x,y) = 0
    real Phiy5
    Phiy5(x,y) = x/(hx1*hy1)
    real Phiy6
    Phiy6(x,y) = 2d0*y/(hy1**2)
    
    ! 二阶导
    real Phixx1
    Phixx1(x,y) = 0
    real Phixx2
    Phixx2(x,y) = 0
    real Phixx3
    Phixx3(x,y) = 0
    real Phixx4
    Phixx4(x,y) = 2d0/(hx1**2)
    real Phixx5
    Phixx5(x,y) = 0
    real Phixx6
    Phixx6(x,y) = 0
    
    real Phixy1
    Phixy1(x,y) = 0
    real Phixy2
    Phixy2(x,y) = 0
    real Phixy3
    Phixy3(x,y) = 0
    real Phixy4
    Phixy4(x,y) = 0
    real Phixy5
    Phixy5(x,y) = 1/(hx1*hy1)
    real Phixy6
    Phixy6(x,y) = 0
    
    real Phiyy1
    Phiyy1(x,y) = 0
    real Phiyy2
    Phiyy2(x,y) = 0
    real Phiyy3
    Phiyy3(x,y) = 0
    real Phiyy4
    Phiyy4(x,y) = 0
    real Phiyy5
    Phiyy5(x,y) = 0
    real Phiyy6
    Phiyy6(x,y) = 2d0/(hy1**2)
    
    
    
    ! Gauss积分点和权重
    lambda(1) = 0.3399810435848562648026658     
    weight(1) = 0.6521451548625461426269361
    lambda(2) = 0.8611363115940525752239465     
    weight(2) = 0.3478548451374538573730639
    lambda(3) = -0.3399810435848562648026658     
    weight(3) = 0.6521451548625461426269361
    lambda(4) = -0.8611363115940525752239465     
    weight(4) = 0.3478548451374538573730639
    
    lambdai = hx1*lambda
    lambdaj = hy1*lambda
    
    diag(1) = 1d0
    diag(2) = 1d0/3d0
    diag(3) = 1d0/3d0
    diag(4) = 4d0/45d0
    diag(5) = 1d0/9d0
    diag(6) = 4d0/45d0
    
    diag = hx*hy*diag
    
    do i = 1,Nx
        do j = 1,Ny
            do d = 1,dim
                diagpol(i,j,:,d) = diag
            end do
        end do
    end do
    
    ! 构造多项式值、边界值和梯度的矩阵
    PolDxx = 0
    PolDxy = 0
    PolDyy = 0
    do i = 1,quad
        do j = 1,quad
            do d = 1,dim
                PolValue(i,j,1,d) = Phi1(lambdai(i),lambdaj(j))
                PolValue(i,j,2,d) = Phi2(lambdai(i),lambdaj(j))
                PolValue(i,j,3,d) = Phi3(lambdai(i),lambdaj(j))
                PolValue(i,j,4,d) = Phi4(lambdai(i),lambdaj(j))
                PolValue(i,j,5,d) = Phi5(lambdai(i),lambdaj(j))
                PolValue(i,j,6,d) = Phi6(lambdai(i),lambdaj(j))
                
                PolDx(i,j,1,d) = Phix1(lambdai(i),lambdaj(j))
                PolDx(i,j,2,d) = Phix2(lambdai(i),lambdaj(j))
                PolDx(i,j,3,d) = Phix3(lambdai(i),lambdaj(j))
                PolDx(i,j,4,d) = Phix4(lambdai(i),lambdaj(j))
                PolDx(i,j,5,d) = Phix5(lambdai(i),lambdaj(j))
                PolDx(i,j,6,d) = Phix6(lambdai(i),lambdaj(j))
                
                PolDy(i,j,1,d) = Phiy1(lambdai(i),lambdaj(j))
                PolDy(i,j,2,d) = Phiy2(lambdai(i),lambdaj(j))
                PolDy(i,j,3,d) = Phiy3(lambdai(i),lambdaj(j))
                PolDy(i,j,4,d) = Phiy4(lambdai(i),lambdaj(j))
                PolDy(i,j,5,d) = Phiy5(lambdai(i),lambdaj(j))
                PolDy(i,j,6,d) = Phiy6(lambdai(i),lambdaj(j))
                
                PolDxx(i,j,4,d) = Phixx4(lambdai(i),lambdaj(j))
                
                PolDxy(i,j,5,d) = Phixy5(lambdai(i),lambdaj(j))
                
                PolDyy(i,j,6,d) = Phiyy6(lambdai(i),lambdaj(j))
            end do
        end do
    end do
    
    do i = 1,quad
        do d = 1,dim
            Poledge(1,i,1,d) = Phi1(hx1,lambdaj(i))
            Poledge(1,i,2,d) = Phi2(hx1,lambdaj(i))
            Poledge(1,i,3,d) = Phi3(hx1,lambdaj(i))
            Poledge(1,i,4,d) = Phi4(hx1,lambdaj(i))
            Poledge(1,i,5,d) = Phi5(hx1,lambdaj(i))
            Poledge(1,i,6,d) = Phi6(hx1,lambdaj(i))
            
            Poledge(2,i,1,d) = Phi1(-hx1,lambdaj(i))
            Poledge(2,i,2,d) = Phi2(-hx1,lambdaj(i))
            Poledge(2,i,3,d) = Phi3(-hx1,lambdaj(i))
            Poledge(2,i,4,d) = Phi4(-hx1,lambdaj(i))
            Poledge(2,i,5,d) = Phi5(-hx1,lambdaj(i))
            Poledge(2,i,6,d) = Phi6(-hx1,lambdaj(i))
            
            Poledge(3,i,1,d) = Phi1(lambdai(i),hy1)
            Poledge(3,i,2,d) = Phi2(lambdai(i),hy1)
            Poledge(3,i,3,d) = Phi3(lambdai(i),hy1)
            Poledge(3,i,4,d) = Phi4(lambdai(i),hy1)
            Poledge(3,i,5,d) = Phi5(lambdai(i),hy1)
            Poledge(3,i,6,d) = Phi6(lambdai(i),hy1)
            
            Poledge(4,i,1,d) = Phi1(lambdai(i),-hy1)
            Poledge(4,i,2,d) = Phi2(lambdai(i),-hy1)
            Poledge(4,i,3,d) = Phi3(lambdai(i),-hy1)
            Poledge(4,i,4,d) = Phi4(lambdai(i),-hy1)
            Poledge(4,i,5,d) = Phi5(lambdai(i),-hy1)
            Poledge(4,i,6,d) = Phi6(lambdai(i),-hy1)
        end do
    end do
    
    ! 网格
    do i = 1,Nx + 1
        XX(i) = xa + hx*(i - 1)
    end do

    do j = 1,Ny + 1
        YY(j) = ya + hy*(j - 1)
    end do
    
    ! 整格点
    Xc = 0.5d0*(XX(1:Nx) + XX(2:Nx + 1))
    Yc = 0.5d0*(YY(1:Ny) + YY(2:Ny + 1))

    ! 初始时间层的值
    Pol = 0
    do i = 1,Nx
        do j = 1,Ny
            lambdai = Xc(i) + hx1*lambda
            lambdaj = Yc(j) + hy1*lambda
            do d = 1,dimpol
                do i1 = 1,quad
                    do j1 = 1,quad
                        Pol(i,j,d,1) = Pol(i,j,d,1) + weight(i1)*weight(j1)*U1(lambdai(i1),lambdaj(j1))*PolValue(i1,j1,d,1)
                        Pol(i,j,d,2) = Pol(i,j,d,2) + weight(i1)*weight(j1)*U2(lambdai(i1),lambdaj(j1))*PolValue(i1,j1,d,1)
                        Pol(i,j,d,3) = Pol(i,j,d,3) + weight(i1)*weight(j1)*U3(lambdai(i1),lambdaj(j1))*PolValue(i1,j1,d,1)
                        Pol(i,j,d,4) = Pol(i,j,d,4) + weight(i1)*weight(j1)*U4(lambdai(i1),lambdaj(j1))*PolValue(i1,j1,d,1)
                    end do
                end do
            end do
        end do
    end do
    Pol = (hx1*hy1*Pol)/diagpol
    
    call Pol_to_Q(Pol,Q,Nx,Ny,dimpol,dim)

    QR = Q ! 真解
    QF(:,:,:,1) = Q ! 动画

    t = 0
    TT(1) = t
    count = 2
    skip = 0

    ! 开始迭代
    do while (t < tend)
        
        alphax = 1
        alphay = 1

        dt = CFL*(hx/alphax + hy/alphay)

        if (t + dt <= tend) then
            t = t + dt
        else
            dt = tend - t
            t = tend
        end if
        
        if (dt < 1e-6) then
            dt = tend - t
            t = tend
        end if
        
        ! 记录指定时刻的图像
        if ((t >= stoptime1) .and. (countstop == 0)) then
            dt = stoptime1 + dt - t
            t = stoptime1
            countstop = countstop + 1
            skip = skiptime - 1
        end if
        
        if ((t >= stoptime2) .and. (countstop == 1)) then
            dt = stoptime2 + dt - t
            t = stoptime2
            countstop = countstop + 1
            skip = skiptime - 1
        end if
        
        if ((t >= stoptime3) .and. (countstop == 2)) then
            dt = stoptime3 + dt - t
            t = stoptime3
            countstop = countstop + 1
            skip = skiptime - 1
        end if
        
        if (t == tend) then
            skip = skiptime - 1
        end if
        
        skip = skip + 1
        
        
        
        if (skip == skiptime) then
            skip = 0
        end if
        
        print *,"已经迭代到 t = ",t
        
        if (RKorder == 3) then
            call RK3DG(Pol,PolValue,PolDx,PolDy,PolDxx,PolDxy,PolDyy,Poledge,Q,hx,hy,dt,Nx,Ny,dim,dimpol,diagpol,lambdai,lambdaj,weight)
        end if
        
        if (skip == 0) then
    
            TT(count) = t
            call Calculate_Real_Solution(Q(:,:,2),Xc,Yc,Nx,Ny,t)
            Q(:,:,3) = abs(Q(:,:,2) - Q(:,:,1))
            QF(:,:,:,count) = Q
            print *,"将 t = ",t,"时刻的信息记入第",count,"帧"
            count = count + 1
            
        end if
     
    end do

    ! 计算误差
    Error = abs(Q(:,:,1) - Q(:,:,2))
    Linfty = 0
    L2V = 0

    
    L2 = 0
    do i = 1,Nx
        do j = 1,Ny
            L2 = L2 + Error(i,j)**2
            if (Error(i,j) > Linfty) then
                Linfty = Error(i,j)
            end if
        end do
    end do
    L2 = (L2/(Nx*Ny))**0.5d0
    
    ! 计算更准确的误差
    
    ! 构造积分点的坐标
    do i = 1,Nx
        lambdai = hx1*lambda + Xc(i)
        Xcquad(quad*(i - 1) + 1:quad*i) = lambdai
    end do
    do j = 1,Ny
        lambdaj = hy1*lambda + Yc(j)
        Ycquad(quad*(j - 1) + 1:quad*j) = lambdaj
    end do
    
    ! 计算真解
    call Calculate_Real_Solution(QRquad,Xcquad,Ycquad,quad*Nx,quad*Ny,tend)
    
    ! 计算数值解
    Qquad = 0
    do i = 1,Nx
        do j = 1,Ny
            do d = 1,dimpol
                Qquad(quad*(i - 1) + 1:quad*i,quad*(j - 1) + 1:quad*j) = Qquad(quad*(i - 1) + 1:quad*i,quad*(j - 1) + 1:quad*j) + Pol(i,j,d,1)*PolValue(:,:,d,1)
            end do
        end do
    end do
    
    Errorquad = QRquad - Qquad
    
    do i = 1,quad
        do j = 1,quad
            weight2(i,j) = weight(i)*weight(j)
        end do
    end do
    
    L2 = 0
    do i = 1,Nx
        do j = 1,Ny
            L2 = L2 + sum(weight2*Errorquad(quad*(i - 1) + 1:quad*i,quad*(j - 1) + 1:quad*j)**2)
        end do
    end do
    L2 = (hx1*hy1*L2)**0.5d0
    
    Linfty = 0
    do i = 1,quad*Nx
        do j = 1,quad*Ny
            if (abs(Errorquad(i,j)) > Linfty) then
                Linfty = abs(Errorquad(i,j))
            end if
        end do
    end do
    
    print *,L2
    print *,Linfty

    QQ1 = reshape(QF(:,:,1,:),(/1,Qlength/))
    QQ2 = reshape(QF(:,:,2,:),(/1,Qlength/))
    QQ3 = reshape(QF(:,:,3,:),(/1,Qlength/))
    QQ4 = reshape(QF(:,:,4,:),(/1,Qlength/))
    
    open(unit = 2,file = 'T.txt')
        do i = 1,count - 1
            write(2,*) TT(i)
        end do
        
    print *,"正在写入文件"
    
    open(unit = 3,file = 'Q2.txt')
    open(unit = 4,file = 'Q3.txt')
    open(unit = 5,file = 'Q4.txt')
    open(unit = 1,file = 'Q1.txt')
        do i = 1,Nx*Ny*(count - 1)
            if (i <= Nx*Ny*count) then
                write(1,*) QQ1(1,i)
                write(3,*) QQ2(1,i)
                write(4,*) QQ3(1,i)
                write(5,*) QQ4(1,i)
                if (mod(i,5000) == 0) then
                    print *,"已完成",i,"/",Nx*Ny*(count - 1)
                end if
            end if
        end do
        

    end program main
    
    
    
    
    
    
    

    subroutine f(u,y)


    real u(4),y(4,2),gamma,v1,v2,p,rho

    gamma = 1.4d0

    rho = u(1)
    v1 = u(2)/u(1)
    v2 = u(3)/u(1)
    p = (gamma - 1)*(u(4) - 0.5d0*(u(2)**2 + u(3)**2)/u(1))

    y(1,1) = 0.5d0*u(1)**2
    y(1,2) = 0.5d0*u(1)**2
    y(2,1) = 0.5d0*u(2)**2
    y(2,2) = 0.5d0*u(2)**2
    y(3,1) = 0.5d0*u(3)**2
    y(3,2) = 0.5d0*u(3)**2
    y(4,1) = 0.5d0*u(4)**2
    y(4,2) = 0.5d0*u(4)**2

    end subroutine f
    
    
    
    
    
    subroutine f1(u,y)
    
    real u(4),y(4),y2(4,2)
    
    call f(u,y2)
    
    y = y2(:,1)
    
    end subroutine f1
    
    
    
    
    
    subroutine f2(u,y)
    
    real u(4),y(4),y2(4,2)
    
    call f(u,y2)
    
    y = y2(:,2)
    
    end subroutine f2
    
    
    
    
    
    function heaviside(x)
    
    implicit none
    
    real x,heaviside
    
    if (x > 0) then
        heaviside = 1
    else
        heaviside = 0
    end if
    
    return
    
    end
    
    
    
    
    ! 从多项式的值得到整格点上的值
    subroutine Pol_to_Q(Pol,Q,Nx,Ny,dimpol,dim)
    
    integer Nx,Ny,dimpol,dim,i,j,k,d
    real Pol(Nx,Ny,dimpol,dim),Q(Nx,Ny,dim)
    
    do i = 1,Nx
        do j = 1,Ny
            do d = 1,dim
                Q(i,j,d) = Pol(i,j,1,d) - (1d0/3d0)*Pol(i,j,4,d) - (1d0/3d0)*Pol(i,j,6,d)
            end do
        end do
    end do
    
    end subroutine Pol_to_Q
    
    
    ! 最大波速
    subroutine max_wavespeed(alpha,rho,u,v,p,gamma,direction)
    
    real u,v,gamma,p,rho
    integer direction
    
    if (direction == 1) then
        alpha = (abs(gamma*p/rho))**0.5d0 + abs(u)
    else if (direction == 2) then
        alpha = (abs(gamma*p/rho))**0.5d0 + abs(v)
    end if
    
    end subroutine max_wavespeed
    
    
    function g(z,x,y,t)
    
    real z,g,x,y,t,pi
    parameter(pi = 4d0*atan(1d0))
    
    
    g = z - (z + (2**1d0)*t*sin(0.5d0*pi*z) - x - y)/(1 + (2**1d0)*0.5d0*pi*t*cos(0.5d0*pi*z))
    
    end function g
    
    
    
    subroutine Newtondd(xi,z0,x1,y1,t,Nmax,Tol)
    
    real xi,z,x1,y1,Tol,x,t,pi
    parameter(pi = 4d0*atan(1d0))
    integer Nmax,N
    
    ! 对于xi = x + y, u 满足 1D burgers 方程 : u_t + u*u_xi = 0
    ! 这时 xi + t*sin(pi*xi/2) = x + y
    ! 相当于解 f(z) = 0, 其中f(z) = z + t*sin(pi*z/2) - x - y
    ! f'(z) = 1 + (pi/2)*t*cos(pi*z/2)
    ! 迭代函数 g(z) = z - f(z)/f'(z)
    
    ! g(z) = z - (z + t*sin(0.5d0*pi*z) - x - y)/(1 + 0.5d0*pi*t*cos(0.5d0*pi*z))
    
    z = g(z0,x1,y1,t)
    N = 1
    
    do while ((abs(z - z0) >= Tol) .and. N <= Nmax)
    
        z0 = z
        z = g(z0,x1,y1,t)
        N = N + 1
    
    end do
    
    xi = g(z,x1,y1,t)
    
    end subroutine Newtondd
    
    
    
    subroutine Calculate_Real_Solution(Q,Xc,Yc,Nx,Ny,t)
    
    integer Nx,Ny,Nmax
    real Q(Nx,Ny),Xc(Nx),Yc(Ny),t,x0,xstar,Tol,pi,z0,xi
    parameter(Nmax = 1000)
    parameter(Tol = 1e-15)
    parameter(pi = 4d0*atan(1.0d0))
    
    do i = 1,Nx
        do j = 1,Ny
            
            if (Xc(i) + Yc(j) > 6d0) then
                z0 = 8d0
            else if ((Xc(i) + Yc(j) <= 6d0) .and. (Xc(i) + Yc(j) >= 2d0)) then
                z0 = 4d0
            else
                z0 = 0d0
            end if
            
            call Newtondd(xi,z0,Xc(i),Yc(j),t,Nmax,Tol)
            
            Q(i,j) = sin(0.5d0*pi*xi)
            
        end do
    end do
    
    end subroutine Calculate_Real_Solution
    
    
    
    
    



