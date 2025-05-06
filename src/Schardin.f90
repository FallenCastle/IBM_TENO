!!! Module setup 
module setup
    implicit none

    integer :: npx,npy,npz
    integer :: iter,maxiter,NFC,order
    integer :: nVar = 4, bound_type(4)
    double precision :: hx,hy,dt,ctime,tf,residglob,eta
    double precision,allocatable :: ub(:,:,:),flux1(:,:,:),&
            flux2(:,:,:),resid(:,:,:) ! res = -dF/dx-dG/dy
    double precision,allocatable :: dB(:,:),Flagp(:,:,:),Flagm(:,:,:)
    double precision,parameter :: gam = 1.4d0, mu = 1d-3
    double precision :: PI = 4.d0*DATAN(1.d0)
    
    double precision::timeStart, timeStart1,timeStop1, timeStart2,timeStop2,timeEnd
    double precision, allocatable:: dLossdBx(:,:), dLossdBy(:,:)
    double precision,dimension(4) :: uInlet,fpInlet,fmInlet
    double precision,dimension(4) :: ua1,ua2,fpa1,fma1,fpa2,fma2
    character*64 :: namerestart
    integer :: restart
    character(len=100) :: folder_path
    integer :: date_values(8), ierr

    !IBM
    character*64 :: grid_file
    double precision::xmin,xmax,ymin,ymax
    double precision,allocatable::xmsh(:),ymsh(:)
    double precision,allocatable::x(:,:,:),y(:,:,:),z(:,:,:),BI_Ainv(:,:,:),IP_Ainv(:,:,:)
    integer:: gridType,npSolid,npSobnd
    integer,allocatable::ij2solid(:,:),solid2ij(:,:),&
                         BI_bbox(:,:,:),IP_bbox(:,:,:),ij2iSolid(:,:),grid_mark(:,:)
    double precision,allocatable::BI2xy(:,:),IP2xy(:,:)
    double precision,allocatable::ub_BI(:,:),ub_IP(:,:)
    double precision:: delta_dis,eta_min

    double precision:: vertex(2,3) !coordinates of vertexes

end module

module setup_device
    implicit none

    integer, constant:: npx_d, npy_d, npxm1_d, npym1_d,order_d
    integer, constant:: npSolid_d
    integer, constant:: nVar_d
    integer, device, allocatable:: ij2solid_d(:,:), solid2ij_d(:,:),&
            BI_bbox_d(:,:,:), IP_bbox_d(:,:,:),ij2iSolid_d(:,:), grid_mark_d(:,:)
    
    double precision, constant:: gam_d,dt_d,eta_min_d,hx_d,hy_d,&
                                 delta_dis_d
    double precision, device:: residglob_d
    double precision, device, allocatable:: ub_d(:,:,:), resid_d(:,:,:),&
            ub0_d(:,:,:), uba_d(:,:,:),&
            uInlet_d(:), fpInlet_d(:), fmInlet_d(:), BI2xy_d(:,:), IP2xy_d(:,:),&
            ub_BI_d(:,:), ub_IP_d(:,:), BI_Ainv_d(:,:,:), IP_Ainv_d(:,:,:),&
            Flagp_d(:,:,:),Flagm_d(:,:,:)
    double precision, device, allocatable:: xmsh_d(:), ymsh_d(:)
    
end module setup_device

module kernel_function
    contains

    subroutine check_cuda_error(error_code,file,line)
        use iso_c_binding
        integer(c_int), value :: error_code
        character(len=*), intent(in) :: file
        integer(c_int), intent(in) :: line

        if (error_code /= 0) then
            write(*, *) 'CUDA Error: ', error_code, ' in file ', file, ' at line ', line
            stop
        endif
    
        
    end subroutine check_cuda_error
    
    subroutine SSPRK_time_stepping_gpu
        use setup
        use setup_device
        use cudafor
        use iso_c_binding
        use ieee_arithmetic
        implicit none
    
        integer:: tilex=32, tiley=8,tile = 256,&
                  bLkx,bLky,gridSolid
        integer:: istat,ipx,ipy,iSolid
        integer(c_int)::ierr_
    
    
        type(dim3)::grid0,gridx,gridy,block
        type(cudaEvent):: startEvent, stopEvent
    
        gridSolid = (npSolid+tile-1)/tile
    
        block = dim3(tilex,tiley,1)
        
        bLkx = (npx+tilex-1)/tilex
        bLky = (npy+tiley-1)/tiley
        grid0 = dim3(bLkx,bLky,1)
    
        bLkx = (npx+1+tilex-1)/tilex
        bLky = (npy  +tiley-1)/tiley
        gridx = dim3(bLkx,bLky,1)
    
        bLkx = (npx  +tilex-1)/tilex
        bLky = (npy+1+tiley-1)/tiley
        gridy = dim3(bLkx,bLky,1)
    
        istat = cudaEventCreate(startEvent)
        istat = cudaEventCreate(stopEvent)
        !!!=========================================
    
        ub0_d = ub_d
        
        
        !stage 1
        resid_d = 0.d0
        call compute_IBM_boundary_gpu
        !compute flux
        istat = cudaEventRecord(startEvent,0)
        call compute_fluxF<<<gridx,block>>>()
        call compute_fluxG<<<gridy,block>>>()
        istat = cudaEventRecord(stopEvent,0)
        istat = cudaEventSynchronize(stopEvent)
        !time advancing stage 1
        istat = cudaEventRecord(startEvent,0)
        call ssprk_stage_1<<<grid0,block>>>()
        istat = cudaEventRecord(stopEvent,0)
        istat = cudaEventSynchronize(stopEvent)
        
        !!--------------------------------------------
        !stage 2
        resid_d = 0.d0
        call compute_IBM_boundary_gpu
        !compute flux
        istat = cudaEventRecord(startEvent,0)
        call compute_fluxF<<<gridx,block>>>()
        call compute_fluxG<<<gridy,block>>>()
        istat = cudaEventRecord(stopEvent,0)
        istat = cudaEventSynchronize(stopEvent)
        
        !time advancing stage 2
        istat = cudaEventRecord(startEvent,0)
        call ssprk_stage_2<<<grid0,block>>>()
        istat = cudaEventRecord(stopEvent,0)
        istat = cudaEventSynchronize(stopEvent)
        !!--------------------------------------------
        uba_d = ub_d

        !stage 3
        resid_d = 0.d0
        call compute_IBM_boundary_gpu
        !compute flux
        istat = cudaEventRecord(startEvent,0)
        call compute_fluxF<<<gridx,block>>>()
        call compute_fluxG<<<gridy,block>>>()
        istat = cudaEventRecord(stopEvent,0)
        istat = cudaEventSynchronize(stopEvent)
        !time advancing stage 3
        istat = cudaEventRecord(startEvent,0)
        call ssprk_stage_3<<<grid0,block>>>()
        istat = cudaEventRecord(stopEvent,0)
        istat = cudaEventSynchronize(stopEvent)

        !!--------------------------------------------
        !stage 4
        resid_d = 0.d0
        call compute_IBM_boundary_gpu
        !compute flux
        istat = cudaEventRecord(startEvent,0)
        call compute_fluxF<<<gridx,block>>>()
        call compute_fluxG<<<gridy,block>>>()
        istat = cudaEventRecord(stopEvent,0)
        istat = cudaEventSynchronize(stopEvent)
        
        !time advancing stage 4
        istat = cudaEventRecord(startEvent,0)
        call ssprk_stage_4<<<grid0,block>>>()
        istat = cudaEventRecord(stopEvent,0)
        istat = cudaEventSynchronize(stopEvent)
        !!--------------------------------------------
        !stage 5
        resid_d = 0.d0
        call compute_IBM_boundary_gpu
        !compute flux
        istat = cudaEventRecord(startEvent,0)
        call compute_fluxF<<<gridx,block>>>()
        call compute_fluxG<<<gridy,block>>>()
        istat = cudaEventRecord(stopEvent,0)
        istat = cudaEventSynchronize(stopEvent)
        
        !time advancing stage 5
        istat = cudaEventRecord(startEvent,0)
        call ssprk_stage_5<<<grid0,block>>>()
        istat = cudaEventRecord(stopEvent,0)
        istat = cudaEventSynchronize(stopEvent)
        !!--------------------------------------------
        
        ctime = ctime + dt
        
    end subroutine SSPRK_time_stepping_gpu
    
    subroutine compute_IBM_boundary_gpu
        use setup
        use setup_device
        !use kernel_function
        use cudafor
        implicit none
        integer:: numBlock, numThread=256
        integer:: istat,iSolid,ipx,ipy
        integer(c_int)::ierr_
        type(cudaEvent):: startEvent, stopEvent
    
        numBlock = (npSolid+numThread-1)/numThread
    
        istat = cudaEventCreate(startEvent)
        istat = cudaEventCreate(stopEvent)
        
        !!!compute ub_BI_d
        istat = cudaEventRecord(startEvent,0)
        call compute_ub_BI_gpu<<<numBlock,numThread>>>()
        istat = cudaEventRecord(stopEvent,0)
        istat = cudaEventSynchronize(stopEvent)

       
        !!!compute ub_IP_d
        istat = cudaEventRecord(startEvent,0)
        call compute_ub_IP_gpu<<<numBlock,numThread>>>()
        istat = cudaEventRecord(stopEvent,0)
        istat = cudaEventSynchronize(stopEvent)
        
        !!!update solid grid points
        istat = cudaEventRecord(startEvent,0)
        call update_solid_points_gpu<<<numBlock,numThread>>>()
        istat = cudaEventRecord(stopEvent,0)
        istat = cudaEventSynchronize(stopEvent)
        
    end subroutine compute_IBM_boundary_gpu
    
    attributes(global) subroutine compute_ub_BI_gpu
        use setup_device
        implicit none
        
        integer iSolid, iGP, jGP
        integer:: i1,j1,i2,j2,i3,j3,i4,j4
        double precision:: xGP, yGP, xBI, yBI, xIP, yIP
        double precision:: x1,y1,x2,y2,x3,y3,x4,y4
        double precision:: phi(4,4)
        double precision:: rhoBI,uBI,vBI,prBI
        double precision:: unBI,utBI
        double precision:: dis,enx,eny,etx,ety,epsilon = 1d-20
    
        iSolid = threadidx%x+(blockidx%x-1)*blockdim%x
        if(iSolid<=npSolid_d)then
            iGP = solid2ij_d(1,iSolid)
            jGP = solid2ij_d(2,iSolid)
            
    
            xGP = xmsh_d(iGP)
            yGP = ymsh_d(jGP)
            xBI = BI2xy_d(1,iSolid)
            yBI = BI2xy_d(2,iSolid)
            xIP = IP2xy_d(1,iSolid)
            yIP = IP2xy_d(2,iSolid)
            
    
            !compute ub_BI_d
            i1 = BI_bbox_d(1,1,iSolid)
            j1 = BI_bbox_d(2,1,iSolid)
            i2 = BI_bbox_d(1,2,iSolid)
            j2 = BI_bbox_d(2,2,iSolid)
            i3 = BI_bbox_d(1,3,iSolid)
            j3 = BI_bbox_d(2,3,iSolid)
            i4 = BI_bbox_d(1,4,iSolid)
            j4 = BI_bbox_d(2,4,iSolid)
            
    
            phi(:,1) = ub_d(i1,j1,:)
            phi(:,2) = ub_d(i2,j2,:)
            phi(:,3) = ub_d(i3,j3,:)
            phi(:,4) = ub_d(i4,j4,:)
            
    
            call linear_interpolation_gpu(xBI,yBI,BI_Ainv_d(:,:,iSolid),phi(:,:),ub_BI_d(:,iSolid))
    
            dis = dsqrt((xIP-xGP)**2.d0 + (yIP-yGP)**2.d0)

            if(dis<=epsilon)then
                dis = dis*1.d5
            endif

            dis = 1.d0/dis
            
            !compute primitive variables
            rhoBI = ub_BI_d(1,iSolid)
            uBI = ub_BI_d(2,iSolid)/rhoBI
            vBI = ub_BI_d(3,iSolid)/rhoBI
            prBI = (ub_BI_d(4,iSolid) - 0.5d0*rhoBI*(uBI**2.d0+vBI**2.d0))*(gam_d - 1.d0)

            enx = (xIP - xGP)*dis
            eny = (yIP - yGP)*dis
            etx = -eny
            ety = enx

            unBI = uBI*enx + vBI*eny
            utBI = uBI*etx + vBI*ety

            unBI = 0.d0

            uBI = unBI*enx + utBI*etx
            vBI = unBI*eny + utBI*ety

            ub_BI_d(2,iSolid) = uBI*rhoBI
            ub_BI_d(3,iSolid) = vBI*rhoBI
            ub_BI_d(4,iSolid) = prBI/(gam_d-1.d0) + 0.5d0*rhoBI*(uBI*uBI+vBI*vBI)
            
        endif
    
        
    end subroutine compute_ub_BI_gpu 

    attributes(global) subroutine compute_ub_IP_gpu
        use setup_device
        implicit none

        integer iSolid, iGP, jGP
        integer:: i1,j1,i2,j2,i3,j3,i4,j4
        double precision:: xGP, yGP, xBI, yBI, xIP, yIP
        double precision:: x1,y1,x2,y2,x3,y3,x4,y4
        double precision:: phi(4,4)
        double precision:: rhoIP,uIP,vIP,prIP

        iSolid = threadidx%x+(blockidx%x-1)*blockdim%x
        if(iSolid<=npSolid_d)then
            xIP = IP2xy_d(1,iSolid)
            yIP = IP2xy_d(2,iSolid)
            !compute ub_IP_d
            i1 = IP_bbox_d(1,1,iSolid)
            j1 = IP_bbox_d(2,1,iSolid)
            i2 = IP_bbox_d(1,2,iSolid)
            j2 = IP_bbox_d(2,2,iSolid)
            i3 = IP_bbox_d(1,3,iSolid)
            j3 = IP_bbox_d(2,3,iSolid)
            i4 = IP_bbox_d(1,4,iSolid)
            j4 = IP_bbox_d(2,4,iSolid)

            if(ij2solid_d(i1,j1)==1)then
                phi(:,1) = ub_BI_d(:,ij2iSolid_d(i1,j1))
            else
                phi(:,1) = ub_d(i1,j1,:)
            endif
    
            if(ij2solid_d(i2,j2)==1)then
                phi(:,2) = ub_BI_d(:,ij2iSolid_d(i2,j2))
            else
                phi(:,2) = ub_d(i2,j2,:)
            endif

            if(ij2solid_d(i3,j3)==1)then
                phi(:,3) = ub_BI_d(:,ij2iSolid_d(i3,j3))
            else
                phi(:,3) = ub_d(i3,j3,:)
            endif
    
            if(ij2solid_d(i4,j4)==1)then
                phi(:,4) = ub_BI_d(:,ij2iSolid_d(i4,j4))
            else
                phi(:,4) = ub_d(i4,j4,:)
            endif
    
            call linear_interpolation_gpu(xIP,yIP,IP_Ainv_d(:,:,iSolid),phi(:,:),ub_IP_d(:,iSolid))

        endif
    
        
    end subroutine compute_ub_IP_gpu
    
    attributes(device) subroutine linear_interpolation_gpu(x,y,Ainv,phi,ub_out)
        use setup_device, only:nVar_d
        implicit none
    
        integer::i
        double precision, intent(in) :: x,y,Ainv(4,4),phi(4,4)
        double precision, intent(out) :: ub_out(4)
        double precision:: a, b, c, d
    
        do i=1,nVar_d
    
            a = Ainv(1,1)*phi(i,1) + Ainv(1,2)*phi(i,2) &
              + Ainv(1,3)*phi(i,3) + Ainv(1,4)*phi(i,4)
    
            b = Ainv(2,1)*phi(i,1) + Ainv(2,2)*phi(i,2) &
              + Ainv(2,3)*phi(i,3) + Ainv(2,4)*phi(i,4)
    
            c = Ainv(3,1)*phi(i,1) + Ainv(3,2)*phi(i,2) &
              + Ainv(3,3)*phi(i,3) + Ainv(3,4)*phi(i,4)
    
            d = Ainv(4,1)*phi(i,1) + Ainv(4,2)*phi(i,2) &
              + Ainv(4,3)*phi(i,3) + Ainv(4,4)*phi(i,4)
              
            ub_out(i) = a*x*y + b*x + c*y +d
    
        enddo
    
        
    end subroutine linear_interpolation_gpu
    
    attributes(global) subroutine update_solid_points_gpu()
        use setup_device
        implicit none
    
        integer:: iSolid,iGP,jGP
        double precision:: xGP,yGP,xBI,yBI,xIP,yIP
        double precision:: rhoBI,uBI,vBI,prBI,&
                           rhoIP,uIP,vIP,prIP,&
                           rhoGP,uGP,vGP,prGP
        double precision:: dis,enx,eny,etx,ety,epsilon = 1d-20
        double precision:: unIP,utIP,unGP,utGP,utBI
    
        
        iSolid = threadidx%x+(blockidx%x-1)*blockdim%x
        if(iSolid<=npSolid_d)then
    
            iGP = solid2ij_d(1,iSolid)
            jGP = solid2ij_d(2,iSolid)
            
            xGP = xmsh_d(iGP)
            yGP = ymsh_d(jGP)
            xBI = BI2xy_d(1,iSolid)
            yBI = BI2xy_d(2,iSolid)
            xIP = IP2xy_d(1,iSolid)
            yIP = IP2xy_d(2,iSolid)
    
            !compute original variables
            rhoBI = ub_BI_d(1,iSolid)
            uBI = ub_BI_d(2,iSolid)/rhoBI
            vBI = ub_BI_d(3,iSolid)/rhoBI
            prBI = (ub_BI_d(4,iSolid) - 0.5d0*rhoBI*(uBI**2.d0+vBI**2.d0))*(gam_d - 1.d0)
    
            rhoIP = ub_IP_d(1,iSolid)
            uIP = ub_IP_d(2,iSolid)/rhoIP
            vIP = ub_IP_d(3,iSolid)/rhoIP
            prIP = (ub_IP_d(4,iSolid) - 0.5d0*rhoIP*(uIP**2.d0+vIP**2.d0))*(gam_d - 1.d0)
    
            rhoGP = rhoIP
            prGP = prIP

            uGP = 2.d0*uBI-uIP
            vGP = 2.d0*vBI-vIP
    
            ub_d(iGP,jGP,1) = rhoGP
            ub_d(iGP,jGP,2) = rhoGP*uGP
            ub_d(iGP,jGP,3) = rhoGP*vGP
            ub_d(iGP,jGP,4) = prGP/(gam_d-1.d0) + 0.5d0*rhoGP*(uGP*uGP + vGP*vGP)
    
        endif
        
    end subroutine update_solid_points_gpu
    
    attributes(global) subroutine compute_fluxF()
        use setup_device
        implicit none
    
        integer::i,j,i1,im3,im2,im1,ip1,ip2,index,pt_index,mark
        double precision::ub_im3(4),ub_im2(4), ub_im1(4),&
                          ub_i(4), ub_ip1(4), ub_ip2(4)
        double precision::uavg(4), Rmat(4,4),&
                          Rinv(4,4), diag(4,4),&
                          fp_stencil(4,5), fm_stencil(4,5),&
                          fluxp(4), fluxm(4), num_flux(4),res,nf(2),&
                          ub_stencil(4,5),delta
    
        i=threadidx%x+(blockidx%x-1)*blockdim%x
        j=threadidx%y+(blockidx%y-1)*blockdim%y
    
        if(i<=npx_d+1 .and. j<=npy_d)then
            delta=xmsh_d(i)-xmsh_d(i-1)
            !if(i==npx_d+1)then
            !    mark=grid_mark_d(npx_d,j)
            !else
            !    mark=grid_mark_d(i,j)
            !endif
            mark=0
            
            nf(1) = 1.d0
            nf(2) = 0.d0
    
            if(i==1)then
                uavg(:) = uInlet_d(:)
            elseif(i==npx_d+1)then
                uavg(:) = ub_d(npx_d,j,:)
            else
                uavg(:) = 0.5d0*(ub_d(i,j,:) + ub_d(i-1,j,:))
            endif
    
            call compute_Rinv_gpu(uavg,nf,Rmat,Rinv,diag)
            
            !-----------------------------------------------
            
            !f+
            do index = 1,5
                pt_index = i-4+index
                if(pt_index.lt.1)then
                    ub_stencil(:,index) = uInlet_d(:)
                elseif(pt_index.gt.npx_d)then
                    ub_stencil(:,index) = ub_d(npx_d,j,:)
                else
                    ub_stencil(:,index) = ub_d(pt_index,j,:)
                endif
            enddo
            call flux_split_gpu(ub_stencil(:,1),nf,fp_stencil(:,1),1)
            call flux_split_gpu(ub_stencil(:,2),nf,fp_stencil(:,2),1)
            call flux_split_gpu(ub_stencil(:,3),nf,fp_stencil(:,3),1)
            call flux_split_gpu(ub_stencil(:,4),nf,fp_stencil(:,4),1)
            call flux_split_gpu(ub_stencil(:,5),nf,fp_stencil(:,5),1)

            call reconstruct_gpu_dis_ad(fp_stencil,Rmat,Rinv,fluxp,mark)
    
            !f-
            do index = 1,5
                pt_index = i+3-index
                if(pt_index.lt.1)then
                    ub_stencil(:,index) = uInlet_d(:)
                elseif(pt_index.gt.npx_d)then
                    ub_stencil(:,index) = ub_d(npx_d,j,:)
                else
                    ub_stencil(:,index) = ub_d(pt_index,j,:)
                endif
            enddo
            call flux_split_gpu(ub_stencil(:,1),nf,fm_stencil(:,1),-1)
            call flux_split_gpu(ub_stencil(:,2),nf,fm_stencil(:,2),-1)
            call flux_split_gpu(ub_stencil(:,3),nf,fm_stencil(:,3),-1)
            call flux_split_gpu(ub_stencil(:,4),nf,fm_stencil(:,4),-1)
            call flux_split_gpu(ub_stencil(:,5),nf,fm_stencil(:,5),-1)
            
            call reconstruct_gpu_dis_ad(fm_stencil,Rmat,Rinv,fluxm,mark)
            !----------------------------------------------------
            
            num_flux(:) = fluxp(:)+fluxm(:)
    
            res=atomicadd(resid_d(i-1,j,1), -num_flux(1)/delta)
            res=atomicadd(resid_d(i,  j,1),  num_flux(1)/delta)
    
            res=atomicadd(resid_d(i-1,j,2), -num_flux(2)/delta)
            res=atomicadd(resid_d(i,  j,2),  num_flux(2)/delta)
    
            res=atomicadd(resid_d(i-1,j,3), -num_flux(3)/delta)
            res=atomicadd(resid_d(i,  j,3),  num_flux(3)/delta)
    
            res=atomicadd(resid_d(i-1,j,4), -num_flux(4)/delta)
            res=atomicadd(resid_d(i,  j,4),  num_flux(4)/delta)
    
        endif
    
    end subroutine compute_fluxF
    
    attributes(global) subroutine compute_fluxG()
        use setup_device
        implicit none
    
        integer::i,j,j1,jm3,jm2,jm1,jp1,jp2,index,pt_index,mark
        double precision::ub_jm3(4),ub_jm2(4), ub_jm1(4),&
                          ub_j(4), ub_jp1(4), ub_jp2(4)
        double precision::uavg(4), Rmat(4,4),&
                          Rinv(4,4), diag(4,4),&
                          fp_stencil(4,5), fm_stencil(4,5),&
                          fluxp(4), fluxm(4), num_flux(4),res,nf(2),&
                          ub_stencil(4,5),delta
    
        i=threadidx%x+(blockidx%x-1)*blockdim%x
        j=threadidx%y+(blockidx%y-1)*blockdim%y
        if(i<=npx_d .and. j<=npy_d+1)then
            delta = ymsh_d(j)-ymsh_d(j-1)
            !if(j==npy_d+1)then
            !    mark = grid_mark_d(i,npy_d)
            !else
            !    mark = grid_mark_d(i,j)
            !endif
            mark=0

            nf(1) = 0.d0
            nf(2) = 1.d0
    
            if(j==1)then
                uavg(:) = ub_d(i,j,:)
                uavg(3) = 0.d0
            elseif(j==npy_d+1)then
                uavg(:) = ub_d(i,j-1,:)
                uavg(3) = 0.d0
            else
                uavg(:) = 0.5d0*(ub_d(i,j,:) + ub_d(i,j-1,:))
            endif
        
            call compute_Rinv_gpu(uavg,nf,Rmat,Rinv,diag)
    
            !-------------------------------------------
    
            !f+
            do index = 1,5
                pt_index = j-4+index
                if(pt_index.lt.1)then
                    ub_stencil(:,index) = ub_d(i,2-pt_index,:)
                    ub_stencil(3,index) = -ub_stencil(3,index)
                elseif(pt_index.gt.npy_d)then
                    ub_stencil(:,index) = ub_d(i,2*npy_d-pt_index,:)
                    ub_stencil(3,index) = -ub_stencil(3,index)
                else
                    ub_stencil(:,index) = ub_d(i,pt_index,:)
                endif
            enddo
            call flux_split_gpu(ub_stencil(:,1),nf,fp_stencil(:,1),1)
            call flux_split_gpu(ub_stencil(:,2),nf,fp_stencil(:,2),1)
            call flux_split_gpu(ub_stencil(:,3),nf,fp_stencil(:,3),1)
            call flux_split_gpu(ub_stencil(:,4),nf,fp_stencil(:,4),1)
            call flux_split_gpu(ub_stencil(:,5),nf,fp_stencil(:,5),1)

            call reconstruct_gpu_dis_ad(fp_stencil,Rmat,Rinv,fluxp,mark)
    
            !f-
            do index = 1,5
                pt_index = j+3-index
                if(pt_index.lt.1)then
                    ub_stencil(:,index) = ub_d(i,2-pt_index,:)
                    ub_stencil(3,index) = -ub_stencil(3,index)
                elseif(pt_index.gt.npy_d)then
                    ub_stencil(:,index) = ub_d(i,2*npy_d-pt_index,:)
                    ub_stencil(3,index) = -ub_stencil(3,index)
                else
                    ub_stencil(:,index) = ub_d(i,pt_index,:)
                endif
            enddo
            call flux_split_gpu(ub_stencil(:,1),nf,fm_stencil(:,1),-1)
            call flux_split_gpu(ub_stencil(:,2),nf,fm_stencil(:,2),-1)
            call flux_split_gpu(ub_stencil(:,3),nf,fm_stencil(:,3),-1)
            call flux_split_gpu(ub_stencil(:,4),nf,fm_stencil(:,4),-1)
            call flux_split_gpu(ub_stencil(:,5),nf,fm_stencil(:,5),-1)
            
            call reconstruct_gpu_dis_ad(fm_stencil,Rmat,Rinv,fluxm,mark)
    
            !----------------------------------------------------
            
            num_flux(:) = fluxp(:)+fluxm(:)
    
            res=atomicadd(resid_d(i,j-1,1), -num_flux(1)/delta)
            res=atomicadd(resid_d(i,  j,1),  num_flux(1)/delta)
    
            res=atomicadd(resid_d(i,j-1,2), -num_flux(2)/delta)
            res=atomicadd(resid_d(i,  j,2),  num_flux(2)/delta)
    
            res=atomicadd(resid_d(i,j-1,3), -num_flux(3)/delta)
            res=atomicadd(resid_d(i,  j,3),  num_flux(3)/delta)
    
            res=atomicadd(resid_d(i,j-1,4), -num_flux(4)/delta)
            res=atomicadd(resid_d(i,  j,4),  num_flux(4)/delta)
    
        endif
    
    end subroutine compute_fluxG
    
    attributes(device) subroutine flux_split_gpu(var,nf,flux,wind)
        use setup_device,only:gam_d,nVar_d
        implicit none
    
        double precision :: var(4),flux(4),nf(2)
        integer :: wind
        double precision :: c,rho,rhou,u,rhov,v,rhow,w,pr,unf,magnorm
    
        magnorm = dsqrt(nf(1)**2+nf(2)**2)
        nf(1) = nf(1)/magnorm
        nf(2) = nf(2)/magnorm
    
        rho = var(1)
        rhou= var(2)
        rhov= var(3)
        u   = rhou/rho
        v   = rhov/rho
        pr  = (var(4)-0.5*rho*(u**2+v**2))*(gam_d-1)
        c   = dsqrt(gam_d*pr/rho)
    
        unf = u*nf(1)+v*nf(2)
    
        flux(1) = rho*unf
        flux(2) = rho*u*unf + pr*nf(1)
        flux(3) = rho*v*unf + pr*nf(2)
        flux(4) = (var(4)+pr)*unf
    
        if (wind==1) then
            flux = 0.5*(flux + (dabs(unf)+c)*var)
        else if (wind==-1) then
            flux = 0.5*(flux - (dabs(unf)+c)*var)
        end if
    
    end subroutine flux_split_gpu
    
    attributes(device) subroutine compute_Rinv_gpu(ut,nf,Rmat,Rinv,diag)
        use setup_device,only:gam_d
        implicit none
        double precision,dimension(4),intent(in) :: ut
        double precision,dimension(2) :: nf
        double precision,dimension(4,4),intent(out) :: Rmat,diag
        double precision,dimension(4,4) :: Rtmp
        double precision,dimension(4,4),intent(out) :: Rinv
        double precision :: c,rho,rhou,rhov,u,v,pr,eH,ek,unf,magnorm
    
        magnorm = dsqrt(nf(1)**2+nf(2)**2)
        nf(1) = nf(1)/magnorm
        nf(2) = nf(2)/magnorm
    
        rho = ut(1)
        rhou= ut(2)
        rhov= ut(3)
        u   = rhou/rho
        v   = rhov/rho
        pr  = (ut(4)-0.5d0*rho*(u**2+v**2))*(gam_d-1)
        c   = dsqrt(gam_d*pr/rho) ! speed of sound
        eH  = (ut(4)+pr)/rho ! specific enthalpy
        unf = u*nf(1)+v*nf(2)
        ek = 0.5d0*(u**2+v**2)
    
        Rmat(1,1) = 1.d0
        Rmat(2,1) = u-c*nf(1)
        Rmat(3,1) = v-c*nf(2)
        Rmat(4,1) = eH-c*unf
    
        Rmat(1,2) = 1.d0
        Rmat(2,2) = u
        Rmat(3,2) = v
        Rmat(4,2) = ek
    
        Rmat(1,3) = 1.d0
        Rmat(2,3) = u+c*nf(1)
        Rmat(3,3) = v+c*nf(2)
        Rmat(4,3) = eH+c*unf
    
        Rmat(1,4) = 0.d0
        Rmat(2,4) = nf(2)
        Rmat(3,4) = -nf(1)
        Rmat(4,4) = u*nf(2)-v*nf(1)
    
        ! ===========================
        Rinv(1,1) = ((gam_d-1)*ek+c*unf)*0.5d0/(c**2.d0)
        Rinv(2,1) = (c**2-(gam_d-1)*ek)/(c**2.d0)
        Rinv(3,1) = ((gam_d-1)*ek-c*unf)*0.5d0/(c**2.d0)
        Rinv(4,1) = v*nf(1)-u*nf(2)
    
        Rinv(1,2) = ((1-gam_d)*u-c*nf(1))*0.5d0/(c**2.d0)
        Rinv(2,2) = (gam_d-1)*u/(c**2.d0)
        Rinv(3,2) = ((1-gam_d)*u+c*nf(1))*0.5d0/(c**2.d0)
        Rinv(4,2) = nf(2)
    
        Rinv(1,3) = ((1-gam_d)*v-c*nf(2))*0.5d0/(c**2.d0)
        Rinv(2,3) = (gam_d-1)*v/(c**2.d0)
        Rinv(3,3) = ((1-gam_d)*v+c*nf(2))*0.5d0/(c**2.d0)
        Rinv(4,3) = -nf(1)
    
        Rinv(1,4) = (gam_d-1)*0.5d0/(c**2.d0)
        Rinv(2,4) = (1-gam_d)/(c**2.d0)
        Rinv(3,4) = (gam_d-1)*0.5d0/(c**2.d0)
        Rinv(4,4) = 0.d0
    
        diag = 0.d0
        diag(1,1) = unf-c
        diag(2,2) = unf
        diag(3,3) = unf+c
        diag(4,4) = unf
    
    end subroutine compute_Rinv_gpu

    attributes(device) subroutine reconstruct_gpu_dis_ad(u_stencil,Rmat,Rinv,uf,mark)
    use setup_device,only:nVar_d
        implicit none
        double precision :: u_stencil(4,5),Rmat(4,4),Rinv(4,4),flux(4),&
                            v_stencil(4,5),vp0(4),vp1(4),vp2(4),vf(4)
        double precision:: beta0, beta1, beta2,&
                            omega0, omega1, omega2, omega_sum
        double precision,intent(out) :: uf(4)
        integer,intent(in):: mark
        double precision::tol = 1d-1
        integer :: i, d0, d1, d2
        double precision:: Cr = 0.2d0, xi = 0.001d0, epsilon, alpha1 = 10.5d0, alpha2 = 5.5d0
        double precision::g,beta,m
        double precision:: eta_im1, eta_i, eta_ip1, eta_id1d2
        double precision:: delF_imd3d2, delF_imd1d2, delF_ipd1d2, delF_ipd3d2
    
        epsilon = 0.9d0*Cr*xi**2.d0/(1.d0-0.9d0*Cr)
        !epsilon = 1d-40
    
        do i = 1,5
            v_stencil(:,i) = matmul(Rinv(:,:),u_stencil(:,i))
        end do
    
        vp0 = 1.d0/3*v_stencil(:,1)-7.d0/6*v_stencil(:,2)+11.d0/6*v_stencil(:,3)
        vp1 = -1.d0/6*v_stencil(:,2)+5.d0/6*v_stencil(:,3)+1.d0/3*v_stencil(:,4)
        vp2 = 1.d0/3*v_stencil(:,3)+5.d0/6*v_stencil(:,4)-1.d0/6*v_stencil(:,5)
    
        do i = 1,nVar_d
    
            delF_imd3d2 = v_stencil(i,2) - v_stencil(i,1)
            delF_imd1d2 = v_stencil(i,3) - v_stencil(i,2)
            delF_ipd1d2 = v_stencil(i,4) - v_stencil(i,3)
            delF_ipd3d2 = v_stencil(i,5) - v_stencil(i,4)
    
            eta_im1 = (dabs(2.d0*delF_imd1d2*delF_imd3d2) + epsilon)/&
                      (delF_imd1d2**2.d0+delF_imd3d2**2.d0 + epsilon)
    
            eta_i   = (dabs(2.d0*delF_ipd1d2*delF_imd1d2) + epsilon)/&
                      (delF_ipd1d2**2.d0+delF_imd1d2**2.d0 + epsilon)
    
            eta_ip1 = (dabs(2.d0*delF_ipd3d2*delF_ipd1d2) + epsilon)/&
                      (delF_ipd3d2**2.d0+delF_ipd1d2**2.d0 + epsilon)
    
            eta_id1d2 = min(eta_im1, eta_i, eta_ip1)


            m = 1.d0 - min(1.d0,eta_id1d2/Cr)
            g = (1.d0 - m)**4.d0*(1.d0 + 4.d0*m)
            beta = alpha1 - alpha2*(1.d0 - g)
            !beta = floor(beta)
            tol = 10**(-beta)

            
    
            beta0 = 13.d0/12*&
                    (v_stencil(i,1)-2.d0*v_stencil(i,2)+v_stencil(i,3))**2.d0+&
                    1.d0/4*(v_stencil(i,1)-4.d0*v_stencil(i,2)+3.d0*v_stencil(i,3))**2.d0
    
            beta1 = 13.d0/12*&
                    (v_stencil(i,2)-2.d0*v_stencil(i,3)+v_stencil(i,4))**2.d0+&
                    1.d0/4*(v_stencil(i,2)-v_stencil(i,4))**2.d0
    
            beta2 = 13.d0/12*&
                    (v_stencil(i,3)-2.d0*v_stencil(i,4)+v_stencil(i,5))**2.d0+&
                    1.d0/4*(3.d0*v_stencil(i,3)-4.d0*v_stencil(i,4)+v_stencil(i,5))**2.d0
    
            omega0 = (1.d0/(beta0+1d-40))**6.d0
            omega1 = (1.d0/(beta1+1d-40))**6.d0
            omega2 = (1.d0/(beta2+1d-40))**6.d0
    
            omega_sum = omega0 + omega1 + omega2
            omega0 = omega0/omega_sum
            omega1 = omega1/omega_sum
            omega2 = omega2/omega_sum
    
    
            if (omega0 < tol) then
                d0 = 0
            else
                d0 = 1
            end if
    
            if (omega1 < tol) then
                d1 = 0
            else
                d1 = 1
            end if
    
            if (omega2 < tol) then
                d2 = 0
            else
                d2 = 1
            end if
    
            if (d0==1 .and. d1==1 .and. d2==1)then
                vf(i) = 0.1d0*vp0(i) + 0.6d0*vp1(i) + 0.3d0*vp2(i)
            elseif(d0==1 .and. d1==1)then
                vf(i) = 0.25d0*vp0(i) + 0.75d0*vp1(i)
            elseif(d1==1 .and. d2==1)then
                vf(i) = 0.5d0*vp1(i) + 0.5d0*vp2(i)
            elseif(d0==1)then
                vf(i) = vp0(i)
            elseif(d1==1)then
                vf(i) = vp1(i)
            elseif(d2==1)then
                vf(i) = vp2(i)
            else
                write(*,*)'Error.'
                stop
            end if
    
        end do
    
        uf = matmul(Rmat,vf)
    end subroutine reconstruct_gpu_dis_ad
    
    attributes(global) subroutine ssprk_stage_1()
        use setup_device
        implicit none
        integer::i,j
        double precision, parameter :: b10=0.377268915331368d0
    
        i=threadidx%x+(blockidx%x-1)*blockdim%x
        j=threadidx%y+(blockidx%y-1)*blockdim%y
        if(i<=npx_d .and. j<=npy_d)then
    
            ub_d(i,j,:) = ub_d(i,j,:)+b10*dt_d*resid_d(i,j,:)
    
            !resid_d(i,j,:) = 0.d0
    
        endif
        
    end subroutine ssprk_stage_1
    
    attributes(global) subroutine ssprk_stage_2()
        use setup_device
        implicit none
        integer::i,j
        double precision, parameter :: b21=0.377268915331368d0
    
        i=threadidx%x+(blockidx%x-1)*blockdim%x
        j=threadidx%y+(blockidx%y-1)*blockdim%y
        if(i<=npx_d .and. j<=npy_d)then
    
            ub_d(i,j,:) = ub_d(i,j,:)+b21*dt_d*resid_d(i,j,:)
    
            !resid_d(i,j,:) = 0.d0
        endif
    
    end subroutine ssprk_stage_2
    
    attributes(global) subroutine ssprk_stage_3()
        use setup_device
        implicit none
        integer::i,j
        double precision, parameter :: a30=0.355909775063327d0
        double precision, parameter :: a32=0.644090224936674d0
        double precision, parameter :: b32=0.242995220537396d0
    
        i=threadidx%x+(blockidx%x-1)*blockdim%x
        j=threadidx%y+(blockidx%y-1)*blockdim%y
        if(i<=npx_d .and. j<=npy_d)then
    
            ub_d(i,j,:) = a30*ub0_d(i,j,:) + a32*ub_d(i,j,:) + b32*dt_d*resid_d(i,j,:)
    
            !resid_d(i,j,:) = 0.d0
    
            !uba_d(i,j,:) = ub_d(i,j,:)
    
        endif
        
    end subroutine ssprk_stage_3
    
    attributes(global) subroutine ssprk_stage_4()
        use setup_device
        implicit none
        integer::i,j
        double precision, parameter :: a40=0.367933791638137d0
        double precision, parameter :: a43=0.632066208361863d0
        double precision, parameter :: b43=0.238458932846290d0
    
        i=threadidx%x+(blockidx%x-1)*blockdim%x
        j=threadidx%y+(blockidx%y-1)*blockdim%y
        if(i<=npx_d .and. j<=npy_d)then
    
            ub_d(i,j,:) = a40*ub0_d(i,j,:) + a43*ub_d(i,j,:) + b43*dt_d*resid_d(i,j,:)
    
            !resid_d(i,j,:) = 0.d0
        endif
        
    end subroutine ssprk_stage_4
    
    attributes(global) subroutine ssprk_stage_5()
        use setup_device
        implicit none
        integer::i,j
        double precision, parameter :: a52=0.237593836598569d0
        double precision, parameter :: a54=0.762406163401431d0
        double precision, parameter :: b54=0.287632146308408d0
        
        i=threadidx%x+(blockidx%x-1)*blockdim%x
        j=threadidx%y+(blockidx%y-1)*blockdim%y
        if(i<=npx_d .and. j<=npy_d)then
    
            ub_d(i,j,:) = a52*uba_d(i,j,:) + a54*ub_d(i,j,:) + b54*dt_d*resid_d(i,j,:)
    
            !resid_d(i,j,:) = 0.d0
        endif
    
    end subroutine ssprk_stage_5
    
endmodule kernel_function

!!!Host Function
program Schardin
    use setup
    use setup_device
    use kernel_function
    use cudafor
    use, intrinsic::ieee_arithmetic
    implicit none

    integer :: ipx,ipy,argc,iostat,device_count,device_id,istat,case
    double precision :: rho, u, v, pr
    double precision,dimension(4,4) :: Rmat,Rinv,diag
    double precision,dimension(2) :: nf
    integer :: i, nCleaning, maxCleaning
    logical:: is_nan
    character*64 :: casename, filename, command, full_path,inputFile,device_id_c

    call cpu_time(timeStart)

    write(*,*)'Input GPU id: '
    read(*,*)device_id

    istat = cudaGetDeviceCount(device_count)
    if (istat /= cudaSuccess) then
        print *, 'Failed to get GPU device count:', cudaGetErrorString(istat)
        stop
    end if
    if(device_id>device_count-1)then
        device_id=device_count-1
    endif
    istat = cudaSetDevice(device_id)
    write(*,*)'Number of available GPU devices:', device_count
    write(*,*)'Successfully set GPU device to:', device_id

    ! 获取当前系统时间
    call date_and_time(values=date_values)

    !initialize parameters
    xmin = -5.d0
    xmax = 13.d0
    ymin = -5.d0
    ymax = 5.d0
    npx = 6000
    npy = 3000
    dt = 0.0001d0
    tf = 5.d0
    delta_dis=0.075
    restart = 0
    NFC=5000
    !
    vertex(1,1)=0.d0
    vertex(2,1)=0.d0
    vertex(1,2)=1.7d0
    vertex(2,2)=1.d0
    vertex(1,3)=1.7d0
    vertex(2,3)=-1.d0
    !!!===================!!!
 
    casename='Schardin'
    ! 格式化文件夹路径为 "月.日.时.分"
    write(folder_path, '(A,".",I2.2,".",I2.2,".",I2.2,".",I2.2)') &
        trim(casename), date_values(2), date_values(3), date_values(5), date_values(6)

    ! 创建文件夹命令
    write(command, '(A," ",A)') 'mkdir -p', trim(folder_path)
    call execute_command_line(command, wait=.true., exitstat=ierr)

    if (ierr /= 0) then
        print *, "Error creating folder: ", trim(folder_path)
        stop
    end if


    if(restart == 1)then
        write(*,*)'Input restart file name: '
        read(*,*)namerestart
    end if

    if(gridType==1)then
        write(*,*)'Input grid file name: '
        read(*,*)grid_file
    endif

    call initialize_grid(gridType)

    !geometry
    call compute_geometry

    maxiter = nint(tf/dt)

    allocate(ub(1:npx,1:npy,nVar)) !conserved variables
    allocate(resid(0:npx+1,0:npy+1,nVar))

    

    if (restart==0) then

        rho = 2.122d0
        u = 0.442d0
        v = 0.d0
        pr = 1.805d0
        uInlet(1) = rho
        uInlet(2) = rho*u
        uInlet(3) = rho*v
        uInlet(4) = pr/(gam-1.d0)+0.5*rho*(u**2+v**2)
        call Initialize

        ctime = 1d-10
        iter = 0
    else if (restart==1) then

        call read_data
        rho = 2.122d0
        u = 0.442d0
        v = 0.d0
        pr = 1.805d0
        uInlet(1) = rho
        uInlet(2) = rho*u
        uInlet(3) = rho*v
        uInlet(4) = pr/(gam-1.d0)+0.5*rho*(u**2+v**2)

    end if

    filename = "Info.dat"
    write(full_path,'(A,"/",A)') trim(folder_path), trim(filename)
    open(333,file=full_path)
    write(333,*)casename
    write(333,*)'npx = ',npx,',npy = ',npy,',dt = ', dt,&
    'grid type = ', gridType
    close(333)

    !!!---------------------------
    write(*,*)'allocating '
    call gpu_allocate
    write(*,*)'allocate end'
    residglob=0.d0
    call cpu_time(timeStart1)
    do while (ctime<=tf)
        call SSPRK_time_stepping_gpu
        call cpu_time(timeStart2)

        iter = iter+1

        if(iter==1 .or.mod(iter,50)==0 .or.iter==maxiter)then

            resid = resid_d
            residglob = norm2(resid(1:npx,1:npy,1:nVar))/dble(4*npx*npy)
            is_nan = ieee_is_nan(residglob)
            if(is_nan)then
                write(*,*)'Error: resid is nan'
                stop
            end if

        endif

        if(iter==1 .or.mod(iter,10)==0 .or.iter==maxiter)then
            print *,'iter=',iter,'time=',ctime,'resid=',residglob
        endif
        
        if (iter==maxiter &
            .or. dabs(ctime-1.50d0)<1d-8 &
            .or. dabs(ctime-2.10d0)<1d-8 &
            .or. dabs(ctime-2.50d0)<1d-8 &
            .or. dabs(ctime-3.09d0)<1d-8 &
            .or. dabs(ctime-4.26d0)<1d-8 &
            .or. mod(iter,NFC)==0) then
            
            is_nan = ieee_is_nan(residglob)
            if(is_nan)then
                write(*,*)'Error.'
                stop
            end if

            ub = ub_d


            call write_solu_tecplot_IBM_binary
            call write_data
        end if
    end do
    call cpu_time(timeEnd)
    print *, 'Total time = ',timeEnd-timeStart,'s'

end program Schardin

subroutine gpu_allocate
    use setup
    use setup_device
    implicit none
    !integer:: ierr

    gam_d = gam
    npx_d = npx
    npy_d = npy
    npSolid_d = npSolid
    nVar_d = nVar
    dt_d = dt
    eta_min_d = eta_min
    hx_d = hx
    hy_d = hy
    delta_dis_d = delta_dis
    allocate(ij2solid_d(npx,npy), stat=ierr); call gerror_report(ierr)
    ij2solid_d = ij2solid
    allocate(solid2ij_d(2,npSolid), stat=ierr); call gerror_report(ierr)
    solid2ij_d = solid2ij
    allocate(BI_bbox_d(2,4,npSolid), stat=ierr); call gerror_report(ierr)
    BI_bbox_d = BI_bbox
    allocate(IP_bbox_d(2,4,npSolid), stat=ierr); call gerror_report(ierr)
    IP_bbox_d = IP_bbox
    allocate(ij2iSolid_d(npx,npy),stat=ierr); call gerror_report(ierr)
    ij2iSolid_d = ij2iSolid

    allocate(ub_d(1:npx,1:npy,nVar),stat=ierr); call gerror_report(ierr)
    allocate(ub0_d(1:npx,1:npy,nVar),stat=ierr); call gerror_report(ierr)
    allocate(uba_d(1:npx,1:npy,nVar),stat=ierr); call gerror_report(ierr)
    ub_d = ub
    allocate(resid_d(0:npx+1,0:npy+1,nVar), stat=ierr); call gerror_report(ierr)
    resid_d = resid
    allocate(uInlet_d(nVar), stat=ierr); call gerror_report(ierr)
    uInlet_d = uInlet
    allocate(BI2xy_d(2,npSolid), stat=ierr); call gerror_report(ierr)
    BI2xy_d = BI2xy
    allocate(IP2xy_d(2,npSolid), stat=ierr); call gerror_report(ierr)
    IP2xy_d = IP2xy
    allocate(ub_BI_d(nVar,npSolid), stat=ierr); call gerror_report(ierr)
    ub_BI_d = ub_BI
    allocate(ub_IP_d(nVar,npSolid), stat=ierr); call gerror_report(ierr)
    ub_IP_d = ub_IP
    allocate(BI_Ainv_d(4,4,npSolid), stat=ierr); call gerror_report(ierr)
    BI_Ainv_d = BI_Ainv
    allocate(IP_Ainv_d(4,4,npSolid), stat=ierr); call gerror_report(ierr)
    IP_Ainv_d = IP_Ainv
    
    allocate(xmsh_d(0:npx+1), stat=ierr); call gerror_report(ierr)
    xmsh_d = xmsh
    allocate(ymsh_d(0:npy+1), stat=ierr); call gerror_report(ierr)
    ymsh_d = ymsh
    allocate(grid_mark_d(npx,npy),stat=ierr); call gerror_report(ierr)
    grid_mark_d = grid_mark

end subroutine gpu_allocate

subroutine initialize_grid(method)
    use setup
    implicit none
    integer, intent(in) :: method
    integer:: ipx, ipy, block_num_grid, itemp, jtemp, ktemp,&
            i,j,k,index
    double precision,allocatable::tmp(:)
    
    if(method==0)then
        allocate(xmsh(0:npx+1))
        allocate(ymsh(0:npy+1))
        
        hx = (xmax-xmin)/dble(npx-1)
        hy = (ymax-ymin)/dble(npy-1)

        do ipx = 1, npx
            xmsh(ipx) = xmin + hx*dble(ipx-1)
        enddo
        do ipy = 1,npy
            ymsh(ipy) = ymin + hy*dble(ipy-1)
        enddo

        xmsh(0) = 2.d0*xmsh(1) - xmsh(2)
        ymsh(0) = 2.d0*ymsh(1) - ymsh(2)
        xmsh(npx+1) = 2.d0*xmsh(npx) - xmsh(npx-1)
        ymsh(npy+1) = 2.d0*ymsh(npy) - ymsh(npy-1)

    elseif(method==1)then
        open(10,file=grid_file,form='binary')
        read(10)block_num_grid
        write(*,*)block_num_grid
        read(10)itemp, jtemp, ktemp
        write(*,*)itemp,jtemp, ktemp
        allocate(xmsh(1:itemp))
        allocate(ymsh(1:jtemp))
        allocate(x(1:itemp,1:jtemp,1:ktemp))
        allocate(y(1:itemp,1:jtemp,1:ktemp))
        allocate(z(1:itemp,1:jtemp,1:ktemp))
        read(10)(((x(i,j,k),i=1,itemp),j=1,jtemp),k=1,ktemp)
        read(10)(((y(i,j,k),i=1,itemp),j=1,jtemp),k=1,ktemp)
        read(10)(((z(i,j,k),i=1,itemp),j=1,jtemp),k=1,ktemp)
        close(10)
        npx = itemp
        npy = jtemp
        npz = ktemp
        xmsh(:) = x(:,1,1)
        ymsh(:) = y(1,:,1)
        deallocate(x,y,z)

        if(xmsh(1)>xmsh(npx))then
            write(*,*)'Reverse x-direction'
            allocate(tmp(itemp))
            tmp(:) = xmsh(:)

            do index = 1,npx
                xmsh(index) = tmp(npx-index+1)
            enddo
            
            deallocate(tmp)
        endif
        if(ymsh(1)>ymsh(npy))then
            write(*,*)'Reverse y-direction'
            allocate(tmp(jtemp))
            tmp(:) = ymsh(:)
            
            do index = 1,npy
                ymsh(index) = tmp(npy-index+1)
            enddo

            deallocate(tmp)
        endif

        xmin = xmsh(1)
        xmax = xmsh(npx)
        ymin = ymsh(1)
        ymax = ymsh(npy)

        xmsh(0) = 2.d0*xmsh(1) - xmsh(2)
        ymsh(0) = 2.d0*ymsh(1) - ymsh(2)
        xmsh(npx+1) = 2.d0*xmsh(npx) - xmsh(npx-1)
        ymsh(npy+1) = 2.d0*ymsh(npy) - ymsh(npy-1)

        write(*,*)"grid read in success."

    else
        write(*,*)'Wrong value of gridType'
        stop
    endif

end subroutine initialize_grid

subroutine compute_geometry
    use setup
    implicit none
    integer:: ipx, ipy, npFluid, iSolid
    integer::i1,j1,i2,j2,i3,j3,i4,j4
    double precision::x1,y1,x2,y2,x3,y3,x4,y4
    double precision:: tmp
    logical:: inside
    double precision::xGP, yGP, xIP, yIP, xBI, yBI,p(2)

    !!!compute npSolid
    npSolid = 0
    allocate(ij2solid(npx,npy))
    allocate(grid_mark(npx,npy))
    do ipy = 1,npy
        do ipx = 1,npx

            !tmp = dsqrt((xmsh(ipx)-x0)**2.d0+(ymsh(ipy)-y0)**2.d0)
            !if(tmp <= radius)then
            !    npSolid = npSolid + 1
            !    ij2solid(ipx,ipy) = 1
            !else
            !    ij2solid(ipx,ipy) = 0
            !endif
            p(1)=xmsh(ipx)
            p(2)=ymsh(ipy)
            call check_point_in_triangle(p,vertex,inside)
            if(inside)then
                npSolid = npSolid + 1
                ij2solid(ipx,ipy) = 1
            else
                ij2solid(ipx,ipy) = 0
            endif

            !if(tmp >= radius-delta_dis .and. tmp <= radius+delta_dis)then
            !    grid_mark(ipx,ipy) = 1
            !else
            !    grid_mark(ipx,ipy) = 0
            !endif


        enddo
    enddo
    npFluid = npx*npy-npSolid

    write(*,*)'Total node = ', npx*npy, ', fluid-node = ', npFluid,', solid-node = ', npSolid

    !!!allocate memory
    allocate(solid2ij(2,npSolid))
    allocate(BI_Ainv(4,4,npSolid))
    allocate(IP_Ainv(4,4,npSolid))
    allocate(BI_bbox(2,4,npSolid))
    allocate(IP_bbox(2,4,npSolid))
    allocate(BI2xy(2,npSolid))
    allocate(IP2xy(2,npSolid))
    allocate(ub_BI(nVar,npSolid))
    allocate(ub_IP(nVar,npSolid))
    allocate(ij2iSolid(npx,npy))
    ij2iSolid=-1
    
    iSolid = 0
    do ipy = 1,npy
        do ipx = 1,npx
            if(ij2solid(ipx,ipy)==1)then
                iSolid = iSolid + 1
                solid2ij(1,iSolid) = ipx
                solid2ij(2,iSolid) = ipy
                ij2iSolid(ipx,ipy) = iSolid
            endif
        enddo
    enddo

    !!!compute BI
    do iSolid = 1, npSolid
        ipx = solid2ij(1,iSolid)
        ipy = solid2ij(2,iSolid)
        xGP = xmsh(ipx)
        yGP = ymsh(ipy)

        call compute_IP_BI_Triangle(xGP,yGP,xIP,yIP,xBI,yBI)
        IP2xy(1,iSolid) = xIP
        IP2xy(2,iSolid) = yIP
        BI2xy(1,iSolid) = xBI
        BI2xy(2,iSolid) = yBI

        !!!find bounding box & compute Ainv
        !! BI
        call binary_search(xBI,yBI,i1,j1,i2,j2,i3,j3,i4,j4)

        BI_bbox(1,1,iSolid) = i1
        BI_bbox(1,2,iSolid) = i2
        BI_bbox(1,3,iSolid) = i3
        BI_bbox(1,4,iSolid) = i4
        BI_bbox(2,1,iSolid) = j1
        BI_bbox(2,2,iSolid) = j2
        BI_bbox(2,3,iSolid) = j3
        BI_bbox(2,4,iSolid) = j4
        x1 = xmsh(i1)
        x2 = xmsh(i2)
        x3 = xmsh(i3)
        x4 = xmsh(i4)
        y1 = ymsh(j1)
        y2 = ymsh(j2)
        y3 = ymsh(j3)
        y4 = ymsh(j4)

        call compute_Ainv(x1,y1,x2,y2,x3,y3,x4,y4,BI_Ainv(:,:,iSolid))
    enddo

    do iSolid=1,npSolid

        ipx = solid2ij(1,iSolid)
        ipy = solid2ij(2,iSolid)
        xGP = xmsh(ipx)
        yGP = ymsh(ipy)
        xIP = IP2xy(1,iSolid)
        yIP = IP2xy(2,iSolid)

        !!IP
        call binary_search(xIP,yIP,i1,j1,i2,j2,i3,j3,i4,j4)
        IP_bbox(1,1,iSolid) = i1
        IP_bbox(1,2,iSolid) = i2
        IP_bbox(1,3,iSolid) = i3
        IP_bbox(1,4,iSolid) = i4
        IP_bbox(2,1,iSolid) = j1
        IP_bbox(2,2,iSolid) = j2
        IP_bbox(2,3,iSolid) = j3
        IP_bbox(2,4,iSolid) = j4
        x1 = xmsh(i1)
        x2 = xmsh(i2)
        x3 = xmsh(i3)
        x4 = xmsh(i4)
        y1 = ymsh(j1)
        y2 = ymsh(j2)
        y3 = ymsh(j3)
        y4 = ymsh(j4)

        if(ij2solid(i1,j1)==1)then
            x1 = BI2xy(1,ij2iSolid(i1,j1))
            y1 = BI2xy(2,ij2iSolid(i1,j1))
        endif

        if(ij2solid(i2,j2)==1)then
            x2 = BI2xy(1,ij2iSolid(i2,j2))
            y2 = BI2xy(2,ij2iSolid(i2,j2))
        endif

        if(ij2solid(i3,j3)==1)then
            x3 = BI2xy(1,ij2iSolid(i3,j3))
            y3 = BI2xy(2,ij2iSolid(i3,j3))
        endif

        if(ij2solid(i4,j4)==1)then
            x4 = BI2xy(1,ij2iSolid(i4,j4))
            y4 = BI2xy(2,ij2iSolid(i4,j4))
        endif

        call compute_Ainv(x1,y1,x2,y2,x3,y3,x4,y4,IP_Ainv(:,:,iSolid))

    enddo

    !open(11,file='zsolid2ij_gpu.txt')
    !open(22,file='zBI_Ainv_gpu.txt')
    !open(33,file='zIP_Ainv_gpu.txt')
    !open(44,file='zBI_bbox_gpu.txt')
    !open(55,file='zIP_bbox_gpu.txt')
    !open(66,file='zBI2xy_gpu.txt')
    !open(77,file='zIP2xy_gpu.txt')
    !open(100,file='zij2iSolid_gpu.txt')
    !do iSolid=1,npSolid
    !    write(11,*),solid2ij(:,iSolid)
    !    write(22,*),BI_Ainv(:,:,iSolid)
    !    write(33,*),IP_Ainv(:,:,iSolid)
    !    write(44,*),BI_bbox(:,:,iSolid)
    !    write(55,*),IP_bbox(:,:,iSolid)
    !    write(66,*),BI2xy(:,iSolid)
    !    write(77,*),IP2xy(:,iSolid)
    !    write(100,*),ij2iSolid(solid2ij(1,iSolid),solid2ij(2,iSolid))
    !enddo
    !close(11)
    !close(22)
    !close(33)
    !close(44)
    !close(55)
    !close(66)
    !close(77)
    !close(100)
    !stop

    
end subroutine compute_geometry

subroutine check_left(p1, p2, p3, left)
    double precision, intent(in) :: p1(2), p2(2), p3(2)
    double precision, intent(out) :: left
    
    left = (p2(1) - p1(1)) * (p3(2) - p1(2)) - (p3(1) - p1(1)) * (p2(2) - p1(2))
end subroutine check_left

subroutine check_point_in_triangle(p, vertex, inside)
    double precision, intent(in) :: p(2)
    double precision, intent(in) :: vertex(2,3)  ! 2x3数组存储三个顶点坐标
    logical, intent(out) :: inside!0-不在三角形内，1-在三角形内
    
    double precision :: d1, d2, d3
    
    call check_left(vertex(:,1),vertex(:,2),p,d1)
    call check_left(vertex(:,2),vertex(:,3),p,d2)
    call check_left(vertex(:,3),vertex(:,1),p,d3)
    
    inside = (d1 >= 0.0d0 .and. d2 >= 0.0d0 .and. d3 >= 0.0d0) .or. &
            (d1 <= 0.0d0 .and. d2 <= 0.0d0 .and. d3 <= 0.0d0)
end subroutine check_point_in_triangle

subroutine compute_distance_symmetric_to_line(xGP,yGP,vertex1,vertex2,d,xIP,yIP,xBI,yBI)
    use, intrinsic:: ieee_arithmetic
    implicit none
    double precision, intent(in) :: xGP,yGP,vertex1(2),vertex2(2)
    double precision, intent(out) ::  d,xIP,yIP,xBI,yBI
    double precision:: a,b,c
    double precision::x_c,y_c,t
    logical:: is_nan

    a = vertex2(2) - vertex1(2)
    b = vertex1(1) - vertex2(1)
    c = vertex2(1)*vertex1(2) - vertex1(1)*vertex2(2)

    d = dabs(a*xGP+b*yGP+c)/dsqrt(a**2.d0+b**2.d0)

    t = (a*xGP+b*yGP+c)/(a**2.d0+b**2.d0)
    xBI = xGP - a*t
    yBI = yGP - b*t

    xIP = 2.d0*xBI - xGP
    yIP = 2.d0*yBI - yGP

    is_nan = ieee_is_nan(d)
    if(is_nan)then
        write(*,*)'Error: distance to line is NAN'
        stop
    endif
    
end subroutine compute_distance_symmetric_to_line

subroutine compute_IP_BI_Triangle(xGP,yGP,xIP,yIP,xBI,yBI)
    use setup, only: vertex
    implicit none
    double precision, intent(in) :: xGP, yGP
    double precision, intent(out) :: xIP, yIP, xBI, yBI
    integer:: selected(3),count
    double precision:: d1,d2,d3,dMin,epsilon=1d-20
    double precision:: xBI1,yBI1,xBI2,yBI2,xBI3,yBI3
    double precision:: xIP1,yIP1,xIP2,yIP2,xIP3,yIP3
    selected=0
    count=0

    call compute_distance_symmetric_to_line(xGP,yGP,vertex(:,1),vertex(:,2),d1,xIP1,yIP1,xBI1,yBI1)
    call compute_distance_symmetric_to_line(xGP,yGP,vertex(:,2),vertex(:,3),d2,xIP2,yIP2,xBI2,yBI2)
    call compute_distance_symmetric_to_line(xGP,yGP,vertex(:,3),vertex(:,1),d3,xIP3,yIP3,xBI3,yBI3)

    dMin=min(d1,d2,d3)

    if(dabs(dMin-d1)<epsilon)then
        xIP=xIP1
        yIP=yIP1
        xBI=xBI1
        yBI=yBI1
        count=count+1
        selected(1) = 1
    endif
    if(dabs(dMin-d2)<epsilon)then
        xIP=xIP2
        yIP=yIP2
        xBI=xBI2
        yBI=yBI2
        count=count+1
        selected(2) = 1
    endif
    if(dabs(dMin-d3)<epsilon)then
        xIP=xIP3
        yIP=yIP3
        xBI=xBI3
        yBI=yBI3
        count=count+1
        selected(3) = 1
    endif
    
    if(count<1)then
        write(*,*)'Error: "compite_IP_BI_Triangle", count<1'
        stop
    elseif(count==2)then
    elseif(count==3)then
    endif
    
end subroutine compute_IP_BI_Triangle

subroutine binary_search(xIP,yIP,i1,j1,i2,j2,i3,j3,i4,j4)
    use setup, only:xmsh,ymsh,npx,npy
    implicit none
    double precision, intent(in) :: xIP, yIP
    integer, intent(out) ::  i1,j1,i2,j2,i3,j3,i4,j4
    integer low, high, mid, iLow, iUp, jLow, jUp, iters, maxiters
    double precision::epsilon=1d-20

    if (xIP < xmsh(1) .or. xIP > xmsh(npx)) then
        write(*, *) "binary search error: x out of range!"
        stop
    endif
    if (yIP < ymsh(1) .or. yIP > ymsh(npy)) then
        write(*, *) "binary search error: y out of range!"
        stop
    endif

    !x-direction
    low = 1
    high = npx
    iters = 0
    maxiters = npx
    do while(low < high - 1)
        iters  = iters + 1
        if (iters > maxiters) then
            write(*, *) "Binary search iteration limit exceeded in x-direction!"
            stop
        endif
        mid = (low + high)/2
        if(xmsh(mid) <= xIP .and. xIP < xmsh(mid+1))then
            iLow = mid
            iUp = mid + 1
            exit
        elseif(xIP < xmsh(mid))then
            high = mid
        else
            low = mid
        endif
    enddo


    !y-direction
    low = 1
    high = npy
    iters = 0
    maxiters = npy
    do while(low < high - 1)
        iters = iters + 1
        if (iters > maxiters) then
            write(*, *) "Binary search iteration limit exceeded in y-direction!"
            stop
        endif
        mid = (low + high)/2
        if(ymsh(mid) <= yIP .and. yIP < ymsh(mid+1))then
            jLow = mid
            jUp = mid + 1
            exit
        elseif(yIP < ymsh(mid))then
            high = mid
        else
            low = mid
        endif
    enddo

    i1 = iLow
    i2 = iUp
    i3 = iUp
    i4 = iLow

    j1 = jLow
    j2 = jLow
    j3 = jUp
    j4 = jUp
    
end subroutine binary_search

subroutine compute_Ainv(x1,y1,x2,y2,x3,y3,x4,y4,Ainv)
    double precision, intent(in) :: x1,y1,x2,y2,x3,y3,x4,y4
    double precision, intent(out) ::Ainv(4,4)
    double precision:: d, epsilon = 1d-40

    !Compute the d
    d = x1*x2*y1*y3 - x1*x3*y1*y2 - x1*x2*y1*y4 &
      - x1*x2*y2*y3 + x1*x4*y1*y2 + x2*x3*y1*y2 &
      + x1*x2*y2*y4 + x1*x3*y1*y4 + x1*x3*y2*y3 &
      - x1*x4*y1*y3 - x2*x3*y1*y3 - x2*x4*y1*y2 &
      - x1*x3*y3*y4 - x1*x4*y2*y4 - x2*x3*y2*y4 &
      + x2*x4*y1*y4 + x2*x4*y2*y3 + x3*x4*y1*y3 &
      + x1*x4*y3*y4 + x2*x3*y3*y4 - x3*x4*y1*y4 &
      - x3*x4*y2*y3 - x2*x4*y3*y4 + x3*x4*y2*y4

    if(dabs(d)<epsilon)then
        write(*,*)'Warning: d is too small.'
    endif

    Ainv(1,1) = (x2*y3 - x3*y2 - x2*y4 + x4*y2 + x3*y4 - x4*y3)/(d)
    Ainv(1,2) = -(x1*y3 - x3*y1 - x1*y4 + x4*y1 + x3*y4 - x4*y3)/(d)
    Ainv(1,3) = (x1*y2 - x2*y1 - x1*y4 + x4*y1 + x2*y4 - x4*y2)/(d)
    Ainv(1,4) = -(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)/(d)

    Ainv(2,1) = -(x2*y2*y3 - x2*y2*y4 - x3*y2*y3 + x3*y3*y4 + x4*y2*y4 - x4*y3*y4)/(d)
    Ainv(2,2) = (x1*y1*y3 - x1*y1*y4 - x3*y1*y3 + x4*y1*y4 + x3*y3*y4 - x4*y3*y4)/(d)
    Ainv(2,3) = -(x1*y1*y2 - x2*y1*y2 - x1*y1*y4 + x2*y2*y4 + x4*y1*y4 - x4*y2*y4)/(d)
    Ainv(2,4) = (x1*y1*y2 - x1*y1*y3 - x2*y1*y2 + x2*y2*y3 + x3*y1*y3 - x3*y2*y3)/(d)

    Ainv(3,1) = (x2*x3*y2 - x2*x3*y3 - x2*x4*y2 + x2*x4*y4 + x3*x4*y3 - x3*x4*y4)/(d)
    Ainv(3,2) = -(x1*x3*y1 - x1*x4*y1 - x1*x3*y3 + x1*x4*y4 + x3*x4*y3 - x3*x4*y4)/(d)
    Ainv(3,3) = (x1*x2*y1 - x1*x2*y2 - x1*x4*y1 + x2*x4*y2 + x1*x4*y4 - x2*x4*y4)/(d)
    Ainv(3,4) = -(x1*x2*y1 - x1*x2*y2 - x1*x3*y1 + x1*x3*y3 + x2*x3*y2 - x2*x3*y3)/(d)

    Ainv(4,1) = -(x2*x3*y2*y4 - x2*x4*y2*y3 - x2*x3*y3*y4 + x3*x4*y2*y3 + x2*x4*y3*y4 - x3*x4*y2*y4)/(d)
    Ainv(4,2) = (x1*x3*y1*y4 - x1*x4*y1*y3 - x1*x3*y3*y4 + x3*x4*y1*y3 + x1*x4*y3*y4 - x3*x4*y1*y4)/(d)
    Ainv(4,3) = -(x1*x2*y1*y4 - x1*x4*y1*y2 - x1*x2*y2*y4 + x2*x4*y1*y2 + x1*x4*y2*y4 - x2*x4*y1*y4)/(d)
    Ainv(4,4) = (x1*x2*y1*y3 - x1*x3*y1*y2 - x1*x2*y2*y3 + x2*x3*y1*y2 + x1*x3*y2*y3 - x2*x3*y1*y3)/(d)
    
end subroutine compute_Ainv

subroutine Initialize
    use setup, only:gam, npx, npy, ub, xmsh, ymsh
    implicit none
    integer:: ipx, ipy
    double precision:: rho, u1, u2, pr, dis, x, y

    do ipy=1,npy
        do ipx=1,npx
            
            x=xmsh(ipx)
            y=ymsh(ipy)

            if(x<0.d0)then
                rho=2.122d0
                u1=0.442d0
                u2=0.d0
                pr=1.805d0
            else
                rho=1.4d0
                u1=0.d0
                u2=0.d0
                pr=1.d0
            endif

            ub(ipx,ipy,1) = rho
            ub(ipx,ipy,2) = rho*u1
            ub(ipx,ipy,3) = rho*u2
            ub(ipx,ipy,4) = pr/(gam-1.d0)+0.5*rho*(u1**2+u2**2)

        enddo
    enddo

    call write_solu_tecplot_IBM_binary
    
end subroutine Initialize

subroutine write_solu_tecplot_IBM_binary
    use setup
    implicit none

    character*6 :: decimal
    character*64 :: filename, full_path
    integer :: i,ipx,ipy,icx,icy,ipc
    double precision :: rho,u,v,pr,q(npx,npy,nVar),xx(npx,npy),yy(npx,npy)
    double precision rho_min,rho_max,u_min,u_max,&
                     v_min,v_max,pr_min,pr_max
    q = ub
    q(:,:,2) = q(:,:,2)/q(:,:,1)
    q(:,:,3) = q(:,:,3)/q(:,:,1)
    q(:,:,4) = (q(:,:,4) - 0.5*q(:,:,1)*(q(:,:,2)**2+q(:,:,3)**2))*(gam-1.0)
    rho_min = q(1,1,1)
    rho_max = rho_min
    u_min = q(1,1,2)
    u_max = u_min
    v_min = q(1,1,3)
    v_max = v_min
    pr_min = q(1,1,4)
    pr_max = pr_min

    do ipy=1,npy
        do ipx = 1,npx

            xx(ipx,ipy) = xmsh(ipx)
            yy(ipx,ipy) = ymsh(ipy)

            if(rho_min>q(ipx,ipy,1)) rho_min=q(ipx,ipy,1)

            if(rho_max<q(ipx,ipy,1)) rho_max=q(ipx,ipy,1)

            if(u_min>q(ipx,ipy,2)) u_min=q(ipx,ipy,2)

            if(u_max<q(ipx,ipy,2)) u_max=q(ipx,ipy,2)

            if(v_min>q(ipx,ipy,3)) v_min=q(ipx,ipy,3)

            if(v_max<q(ipx,ipy,3)) v_max=q(ipx,ipy,3)

            if(pr_min>q(ipx,ipy,4)) pr_min=q(ipx,ipy,4)

            if(pr_max<q(ipx,ipy,4)) pr_max=q(ipx,ipy,4)

        enddo
    enddo


    write(decimal,'(f6.3)') ctime-int(ctime)
    write(filename,'(a,i7.7,a4,a)')'TEC_IBM',int(ctime),decimal(3:6),'.plt'
    write(full_path,'(A,"/",A)') trim(folder_path), trim(filename)
    open(11,file=full_path,form='binary')

    write(11) "#!TDV108"            ! Magic number,version number 
	write(11) 1                     ! byte order
	write(11) ICHAR('a')            ! title
	write(11) 0                     
	write(11) 7                    ! Number of variables (NumVar) in the datafile
	write(11) ichar('x')            ! name of variables (NumVar) in the datafile
	write(11) 0
	write(11) ichar('y')
	write(11) 0
    write(11) ichar('f'),ichar('l'),ichar('a'),ichar('g')
	write(11) 0
	write(11) ichar('r'),ichar('h'),ichar('o')
	write(11) 0
	write(11) ichar('u')
	write(11) 0
	write(11) ichar('v')
	write(11) 0
	write(11) ichar('p')
	write(11) 0
	write(11) 299.0_4                ! Zone marker. Value = 299.0,single
	!write(zonename,'(a,i4.4)') 'zone',myid+1
    write(11) ichar('z'),ichar('o'),ichar('n'),ichar('e'),ichar('1')

    write(11) 0
	write(11) -1                    ! ParentZone: Zero based zone number within this datafile to which this zone is a child.
	write(11) -1                    ! StrandID   -1 static strandID
	write(11) ctime                ! Solution time,double
	write(11) -1                    ! Zone Color (set to -1 if you want tecplot to determine).
	write(11) 0                     ! ZoneType 0=ORDERED,1=FELINESEG,2=FETRIANGLE,
	                               ! 3=FEQUADRILATERAL,4=FETETRAHEDRON,5=FEBRICK
	write(11) 0                     ! DataPacking 0=Block, 1=Point
	write(11) 0                     ! Specify Var Location.  0 = Don't specify, 
	                               ! all data is located at the nodes.  1 = Specify
	write(11) 0                     ! Are raw local 1-to-1 face neighbors supplied? (0=FALSE 1=TRUE)
	write(11) 0                     ! Number of miscellaneous user defined face neighbor connections (value >= 0)
	write(11) npx,npy,1  ! if Ordered Zone:IMax,JMax,KMax
	write(11) 0                     ! 1=Auxiliary name/value pair to follow 0=No more Auxiliar name/value pairs
	                               ! 0=No more Auxiliar name/value pairs

	write(11) 357.0_4                 ! EOHMARKER,single

    !II. 数据区
	write(11) 299.0_4                 ! Zone marker  Value = 299.0,single
	do i=1,7
		write(11) 2                 ! variable data format, N=Total number of vars,1=Float, 2=Double, 3=LongInt, 4=ShortInt, 5=Byte, 6=Bit
	end do                         ! 1=Float, 2=Double, 3=LongInt, 4=ShortInt, 5=Byte, 6=Bit
	write(11) 0                     ! Has passive variables: 0 = no, 1 = yes.
	write(11) 0                     ! Has variable sharing 0 = no, 1 = yes.
	write(11) -1                    ! Zero based zone number to share connectivity list with (-1 = no sharing)

	write(11) xmin,xmax
	write(11) ymin,ymax
    write(11) 0.d0,1.d0
	write(11) rho_min,rho_max
	write(11) u_min,u_max
    write(11) v_min,v_max
    write(11) pr_min,pr_max

    write(11) xx(1:npx,1:npy)
    write(11) yy(1:npx,1:npy)
    write(11) dble(ij2solid(1:npx,1:npy))
    write(11) q(1:npx,1:npy,1)
    write(11) q(1:npx,1:npy,2)
    write(11) q(1:npx,1:npy,3)
    write(11) q(1:npx,1:npy,4)
		
    close(11)

    
end subroutine write_solu_tecplot_IBM_binary

subroutine write_data
    use setup
    implicit none

    integer :: io_status
    character*6 :: decimal
    character*64 :: filename, full_path

    write(decimal,'(f6.3)') ctime-int(ctime)
    write(filename,'(a,i7.7,a4,a)')'restart.',int(ctime),decimal(3:6),'.dat'
    write(full_path,'(A,"/",A)') trim(folder_path), trim(filename)

    open(unit=21,file=full_path,form="unformatted",&
            status="replace",action="write",iostat=io_status)
    write(21) ctime,iter,ub
    close(21)

end subroutine write_data

subroutine read_data
    use setup
    implicit none

    integer :: io_status
    character*64::full_path
    write(full_path,'(A,"/",A)') trim(folder_path), trim(namerestart)
    open(unit=22,file=full_path,form="unformatted",&
            status="old",action="read",iostat=io_status)
    read(22) ctime,iter,ub
    close(22)

end subroutine read_data

subroutine gerror_report(ierr)
    implicit none
        integer::ierr
    
        if(ierr/=0) then
            write(*,*)'Device cannot allocate memory!'
            stop
        end if
    
end subroutine gerror_report
