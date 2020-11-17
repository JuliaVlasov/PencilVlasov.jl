using FFTW

struct Vlasov4d

    f
    d_dx
    d_dy
    kx
    ky
    px
    py
    tmp_x
    tmp_y

    function Vlasov4d( maxwell )

        nx = maxwell.nx
        ny = maxwell.ny

        tmp_x = zeros(ComplexF64, nx)
        tmp_y = zeros(ComplexF64, nx)
        d_dx = similar(tmp_x)
        d_dy = similar(tmp_y)

        px = plan_fft( f, 1)
        py = plan_fft( f, 2)

        xmin = maxwell.xmin
        xmax = maxwell.xmax
        ymin = maxwell.ymin
        ymax = maxwell.ymax

        kx = 2π/(xmax-xmin)*[0:nx÷2-1;nx÷2-nx:-1]
        ky = 2π/(ymax-ymin)*[0:ny÷2-1;ny÷2-ny:-1]

        f = zeros(ComplexF64, (n_i, n_j, n_k, n_j))
        f̂ = zeros(ComplexF64, (n_i, n_j, n_k, n_j))

        new( d_dx, d_dy, kx, ky, px, py, tmp_x, tmp_y )

    end
   
end

# Euler explicite
#         this%tmp_x = this%tmp_x*(1._f64-cmplx(0.0_f64,this%kx,kind=f64)*vx)
# Euler implicite
#         this%tmp_x = this%tmp_x/(1._f64+cmplx(0.0_f64,this%kx,kind=f64)*vx)
# crank-nicolson
#         this%tmp_x = this%tmp_x/(1._f64+cmplx(0.0_f64,this%kx,kind=f64)*vx*0.5_f64)
#         this%tmp_x = this%tmp_x*(1._f64-cmplx(0.0_f64,this%kx,kind=f64)*vx*0.5_f64)
# Euler cn modified
#         this%tmp_x=this%tmp_x*(1._f64-0.5*vx*cmplx(0.0_f64,this%kx,kind=f64))
function advection_x1(this, dt)

  nc_x1    = this.nc_η1
  x3_min   = this.η3_min
  delta_x3 = this.delta_η3

  for l=1:nc_l
  for k=1:nc_k
     vx = (x3_min +(k-1)*delta_x3)*dt
     for j=1,nc_j
        fft_exec_r2c_1d(this%fwx, this%f(1:nc_x1,j,k,l),this%tmp_x)
        this%tmp_x = this%tmp_x*exp(-cmplx(0.0_f64,this%kx,kind=f64)*vx)
        call sll_s_fft_exec_c2r_1d(this%bwx, this%tmp_x, this%d_dx)
        this%f(1:nc_x1,j,k,l)= this%d_dx / nc_x1
     end
  end
  end

end 

#euler explicite
# tmp_y = this%tmp_y*(1._f64-cmplx(0.0_f64,this%ky,kind=f64)*vy)
#euler implicite
# tmp_y = this%tmp_y/(1._f64+cmplx(0.0_f64,this%ky,kind=f64)*vy)
#crank-nicolson
# tmp_y = this%tmp_y/(1._f64+cmplx(0.0_f64,this%ky,kind=f64)*vy*0.5_f64)
# tmp_y = this%tmp_y*(1._f64-cmplx(0.0_f64,this%ky,kind=f64)*vy*0.5_f64)
#Euler cn modified
# tmp_y = this%tmp_y*(1._f64-cmplx(0.0_f64,this%ky,kind=f64)*vy-0.5_f64*(this%ky*vy)**2)
function advection_x2(this,dt)

  for l=1,n_l
    vy = (x4_min +(gl-1)*delta_x4)*dt
    for k=1,n_k
    for i=1,n_i
       mul!(f̂, py, f)
       f̂ .= f̂ .* eky .* vy
       ldiv!(f, py, f̂)
       f[i,:,k,l] = d_dy 
    end
    end
  end

end

function advection_x3x4(this,dt)

  for i=1,n_i
  for j=1,n_j

     for k=1,n_k
     for l=1,n_l

        px = x3_min+(k-1)*delta_x3
        py = x4_min+(l-1)*delta_x4
        cthη = cos(bz[i,j]*dt)
        sthη = sin(bz[i,j]*dt)
        depvx  = 0.5dt * ex[i,j]
        depvy  = 0.5dt * ey[i,j]
        alpha_x[k,l] = - (px - (depvx+(px+depvx)*cthη-(py+depvy)*sthη))
        alpha_y[k,l] = - (py - (depvy+(px+depvx)*sthη+(py+depvy)*cthη))

     end
     end

     interpolate_disp(n_k,n_l, ft[i,j,:,:],alpha_x,alpha_y, ft[i,j,:,:])
  end
  end

end

@testset "Vlasov-Maxwell 2d2v" begin

    maxwell = PSTD( mesh )
    poisson = Poisson( mesh )
    for l=1,n_l, k=1,n_k, j=1,n_j, i=1,n_i
        x  = η1_min+(i-1)*delta_η1
        y  = η2_min+(j-1)*delta_η2
        vx = η3_min+(k-1)*delta_η3
        vy = η4_min+(l-1)*delta_η4
        v2 = vx*vx+vy*vy
        f[i,j,k,l] = ( 1 + ϵ * cos(kx * x) * cos( ky * y )) * exp(-v2)
    end


    compute_charge(vlasov4d)
    solve(poisson, ex, ey, rho)
    bz = zeros(nx,ny)
    
    transposevx!(vlasov4d)
    advection_x1!(vlasov4d,0.5dt)
    advection_x2!(vlasov4d,0.5dt)
    
    for istep = 1:nstep
    
         transposexv!(vlasov4d)
    
         call compute_current(vlasov4d)
    
         solve_ampere!(maxwell, dt)
    
         call advection_x3x4(vlasov4d, dt)
    
         call transposevx!(vlasov4d)
    
         call advection_x1(vlasov4d,dt)
         call advection_x2(vlasov4d,dt)
    
    end

end
