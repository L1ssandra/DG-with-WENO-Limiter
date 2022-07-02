    subroutine WENO_Limiter(Pol,PolValue,PolDx,PolDy,PolDxx,PolDxy,PolDyy,Poledge,Q,hx,hy,dt,Nx,Ny,dim,dimpol,diagpol,lambdai,lambdaj,weight)
    
    integer i,j,Nx,Ny,dim,dimpol,quad,k1,i1,j1
    parameter(quad = 4)
    real Pol(Nx,Ny,dimpol,dim),PolValue(quad,quad,dimpol,dim),PolDx(quad,quad,dimpol,dim),PolDy(quad,quad,dimpol,dim),Polnew(Nx,Ny,dimpol,dim)
    real Q(Nx,Ny,dim),diagpol(Nx,Ny,dimpol,dim),Poledge(4,quad,dimpol,dim)
    real LPol(Nx,Ny,dimpol,dim),lambdai(quad),lambdaj(quad),weight(quad)
    real PolDxx(quad,quad,dimpol,dim),PolDxy(quad,quad,dimpol,dim),PolDyy(quad,quad,dimpol,dim)
    real PolR(Nx,Ny,dimpol,dim),PolL(Nx,Ny,dimpol,dim),PolU(Nx,Ny,dimpol,dim),PolD(Nx,Ny,dimpol,dim),Is_Trouble_Cell(Nx,Ny)
    real uintR(quad),uextR(quad),uintL(quad),uextL(quad),uintU(quad),uextU(quad),uintD(quad),uextD(quad)
    real pC(dimpol,dim),pR(dimpol,dim),pL(dimpol,dim),pU(dimpol,dim),pD(dimpol,dim),pnew(dimpol,dim)
    real w0(quad),w1(quad),w2(quad),w3(quad),w4(quad),beta0(quad),beta1(quad),beta2(quad),beta3(quad),beta4(quad),epsilon
    real gamma0(quad),gamma1(quad),gamma2(quad),gamma3(quad),gamma4(quad),S(quad)
    real ux(quad,quad,dim),uy(quad,quad,dim),uxx(quad,quad,dim),uxy(quad,quad,dim),uyy(quad,quad,dim),ubar(dim),Ii
    integer boundconditionX,boundconditionY
    
    epsilon = 1e-12
    
    gamma0 = 0.6d0
    gamma1 = 0.1d0
    gamma2 = 0.1d0
    gamma3 = 0.1d0
    gamma4 = 0.1d0
    
    boundconditionX = 1
    boundconditionY = 1
    
    hx1 = 0.5d0*hx
    hy1 = 0.5d0*hy
    
    Polnew = 0
    
    ! 赋边界条件
    if (boundconditionX == 1) then ! 周期边界
        PolR(1:Nx - 1,:,:,:) = Pol(2:Nx,:,:,:)
        PolR(Nx,:,:,:) = Pol(1,:,:,:)
        
        PolL(2:Nx,:,:,:) = Pol(1:Nx - 1,:,:,:)
        PolL(1,:,:,:) = Pol(Nx,:,:,:)
    end if
    
    if (boundconditionY == 1) then ! 周期边界    
        PolU(:,1:Ny - 1,:,:) = Pol(:,2:Ny,:,:)
        PolU(:,Ny,:,:) = Pol(:,1,:,:)
        
        PolD(:,2:Ny,:,:) = Pol(:,1:Ny - 1,:,:)
        PolD(:,1,:,:) = Pol(:,Ny,:,:)
    end if
    
    Is_Trouble_Cell = 0
    
    ! 获取整格点的值
    call Pol_to_Q(Pol,Q,Nx,Ny,dimpol,dim)
    
    ! KXRCF检测器
    do i = 1,Nx
        do j = 1,Ny
            
            Ii = 1e-9
            Pmax = 1e-9
            
            ! 判断流入
            if (Q(i,j,1) < 0) then ! 右边和上边流入
                uintR = 0
                uextR = 0
                uintU = 0
                uextU = 0
                do d = 1,dimpol
                    uintR = uintR + Pol(i,j,d,1)*Poledge(1,:,d,1)
                    uextR = uextR + PolR(i,j,d,1)*Poledge(2,:,d,1)
                    uintU = uintU + Pol(i,j,d,1)*Poledge(3,:,d,1)
                    uextU = uextU + PolU(i,j,d,1)*Poledge(4,:,d,1)
                end do
                do k1 = 1,quad
                    
                    Ii = Ii + weight(k1)*((uintR(k1) - uextR(k1))*hy1 + (uintU(k1) - uextU(k1))*hx1)
                    if (abs(uintR(k1)) > Pmax) then
                        Pmax = abs(uintR(k1))
                    end if
                    if (abs(uintU(k1)) > Pmax) then
                        Pmax = abs(uintU(k1))
                    end if
                end do
            else ! 左边和下边流入
                uintL = 0
                uextL = 0
                uintD = 0
                uextD = 0
                do d = 1,dimpol
                    uintL = uintL + Pol(i,j,d,1)*Poledge(2,:,d,1)
                    uextL = uextL + PolL(i,j,d,1)*Poledge(1,:,d,1)
                    uintD = uintD + Pol(i,j,d,1)*Poledge(4,:,d,1)
                    uextD = uextD + PolD(i,j,d,1)*Poledge(3,:,d,1)
                end do
                do k1 = 1,quad
                    Ii = Ii + weight(k1)*((uintL(k1) - uextL(k1))*hy1 + (uintD(k1) - uextD(k1))*hx1)
                    if (abs(uintL(k1)) > Pmax) then
                        Pmax = abs(uintL(k1))
                    end if
                    if (abs(uintD(k1)) > Pmax) then
                        Pmax = abs(uintD(k1))
                    end if
                end do
            end if
            
            Ii = abs(Ii)/(max(hx1,hy1)**1.5d0*(hx + hy)*Pmax)
            
            if (Ii >= 1) then
                Is_Trouble_Cell(i,j) = 1
            end if
            
        end do
    end do
    
    
    ! WENO Limiter
    do i = 1,Nx
        do j = 1,Ny
            if (Is_Trouble_Cell(i,j) == 1) then
                
                pC = Pol(i,j,:,:)
                pR = PolR(i,j,:,:)
                pL = PolL(i,j,:,:)
                pU = PolU(i,j,:,:)
                pD = PolD(i,j,:,:)
                
                ! 修正邻居多项式的单元平均值
                ubar = pC(1,:)
                pR(1,:) = ubar
                pL(1,:) = ubar
                pU(1,:) = ubar
                pD(1,:) = ubar
                
                ! 计算光滑指示器
                beta0 = 0
                beta1 = 0
                beta2 = 0
                beta3 = 0
                beta4 = 0
                
                ! 中心
                ux = 0
                uy = 0
                uxx = 0
                uxy = 0
                uyy = 0
                do m = 1,dim
                    do d = 1,dimpol
                        ux(:,:,m) = ux(:,:,m) + pC(d,m)*PolDx(:,:,d,m)
                        uy(:,:,m) = uy(:,:,m) + pC(d,m)*PolDy(:,:,d,m)
                        uxx(:,:,m) = uxx(:,:,m) + pC(d,m)*PolDxx(:,:,d,m)
                        uxy(:,:,m) = uxy(:,:,m) + pC(d,m)*PolDxy(:,:,d,m)
                        uyy(:,:,m) = uyy(:,:,m) + pC(d,m)*PolDyy(:,:,d,m)
                    end do
                end do
                
                do i1 = 1,quad
                    do j1 = 1,quad
                        beta0 = beta0 + weight(i1)*weight(j1)*(ux(i1,j1,:)**2 + uy(i1,j1,:)**2 + hx*hy*(uxx(i1,j1,:)**2 + uxy(i1,j1,:)**2 + uyy(i1,j1,:)**2))
                    end do
                end do
                
                beta0 = beta0*hx1*hy1
                
                ! 右边
                ux = 0
                uy = 0
                uxx = 0
                uxy = 0
                uyy = 0
                do m = 1,dim
                    do d = 1,dimpol
                        ux(:,:,m) = ux(:,:,m) + pR(d,m)*PolDx(:,:,d,m)
                        uy(:,:,m) = uy(:,:,m) + pR(d,m)*PolDy(:,:,d,m)
                        uxx(:,:,m) = uxx(:,:,m) + pR(d,m)*PolDxx(:,:,d,m)
                        uxy(:,:,m) = uxy(:,:,m) + pR(d,m)*PolDxy(:,:,d,m)
                        uyy(:,:,m) = uyy(:,:,m) + pR(d,m)*PolDyy(:,:,d,m)
                    end do
                end do
                
                do i1 = 1,quad
                    do j1 = 1,quad
                        beta1 = beta1 + weight(i1)*weight(j1)*(ux(i1,j1,:)**2 + uy(i1,j1,:)**2 + hx*hy*(uxx(i1,j1,:)**2 + uxy(i1,j1,:)**2 + uyy(i1,j1,:)**2))
                    end do
                end do
                
                beta1 = beta1*hx1*hy1
                
                ! 左边
                ux = 0
                uy = 0
                uxx = 0
                uxy = 0
                uyy = 0
                do m = 1,dim
                    do d = 1,dimpol
                        ux(:,:,m) = ux(:,:,m) + pL(d,m)*PolDx(:,:,d,m)
                        uy(:,:,m) = uy(:,:,m) + pL(d,m)*PolDy(:,:,d,m)
                        uxx(:,:,m) = uxx(:,:,m) + pL(d,m)*PolDxx(:,:,d,m)
                        uxy(:,:,m) = uxy(:,:,m) + pL(d,m)*PolDxy(:,:,d,m)
                        uyy(:,:,m) = uyy(:,:,m) + pL(d,m)*PolDyy(:,:,d,m)
                    end do
                end do
                
                do i1 = 1,quad
                    do j1 = 1,quad
                        beta2 = beta2 + weight(i1)*weight(j1)*(ux(i1,j1,:)**2 + uy(i1,j1,:)**2 + hx*hy*(uxx(i1,j1,:)**2 + uxy(i1,j1,:)**2 + uyy(i1,j1,:)**2))
                    end do
                end do
                
                beta2 = beta2*hx1*hy1
                
                ! 上边
                ux = 0
                uy = 0
                uxx = 0
                uxy = 0
                uyy = 0
                do m = 1,dim
                    do d = 1,dimpol
                        ux(:,:,m) = ux(:,:,m) + pU(d,m)*PolDx(:,:,d,m)
                        uy(:,:,m) = uy(:,:,m) + pU(d,m)*PolDy(:,:,d,m)
                        uxx(:,:,m) = uxx(:,:,m) + pU(d,m)*PolDxx(:,:,d,m)
                        uxy(:,:,m) = uxy(:,:,m) + pU(d,m)*PolDxy(:,:,d,m)
                        uyy(:,:,m) = uyy(:,:,m) + pU(d,m)*PolDyy(:,:,d,m)
                    end do
                end do
                
                do i1 = 1,quad
                    do j1 = 1,quad
                        beta3 = beta3 + weight(i1)*weight(j1)*(ux(i1,j1,:)**2 + uy(i1,j1,:)**2 + hx*hy*(uxx(i1,j1,:)**2 + uxy(i1,j1,:)**2 + uyy(i1,j1,:)**2))
                    end do
                end do
                
                beta3 = beta3*hx1*hy1
                
                ! 下边
                ux = 0
                uy = 0
                uxx = 0
                uxy = 0
                uyy = 0
                do m = 1,dim
                    do d = 1,dimpol
                        ux(:,:,m) = ux(:,:,m) + pD(d,m)*PolDx(:,:,d,m)
                        uy(:,:,m) = uy(:,:,m) + pD(d,m)*PolDy(:,:,d,m)
                        uxx(:,:,m) = uxx(:,:,m) + pD(d,m)*PolDxx(:,:,d,m)
                        uxy(:,:,m) = uxy(:,:,m) + pD(d,m)*PolDxy(:,:,d,m)
                        uyy(:,:,m) = uyy(:,:,m) + pD(d,m)*PolDyy(:,:,d,m)
                    end do
                end do
                
                do i1 = 1,quad
                    do j1 = 1,quad
                        beta4 = beta4 + weight(i1)*weight(j1)*(ux(i1,j1,:)**2 + uy(i1,j1,:)**2 + hx*hy*(uxx(i1,j1,:)**2 + uxy(i1,j1,:)**2 + uyy(i1,j1,:)**2))
                    end do
                end do
                
                beta4 = beta4*hx1*hy1
                        
                
                ! 计算权重
                w0 = gamma0/(beta0 + epsilon)**2
                w1 = gamma1/(beta1 + epsilon)**2
                w2 = gamma2/(beta2 + epsilon)**2
                w3 = gamma3/(beta3 + epsilon)**2
                w4 = gamma4/(beta4 + epsilon)**2
                
                S = w0 + w1 + w2 + w3 + w4
                
                w0 = w0/S
                w1 = w1/S
                w2 = w2/S
                w3 = w3/S
                w4 = w4/S
                
                ! 得到重构多项式
                do m = 1,dim
                    pnew(:,m) = w0(m)*pC(:,m) + w1(m)*pR(:,m) + w2(m)*pL(:,m) + w3(m)*pU(:,m) + w4(m)*pD(:,m)
                end do
                
                Polnew(i,j,:,:) = pnew
                
            end if
        end do
    end do
    
    do i = 1,Nx
        do j = 1,Ny
            if (Is_Trouble_Cell(i,j) == 1) then
                
                ! 代回
                Pol(i,j,:,:) = Polnew(i,j,:,:)
                
            end if
        end do
    end do
    
    end subroutine WENO_Limiter