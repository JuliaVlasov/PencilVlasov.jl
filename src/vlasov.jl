struct Vlasov

    mesh :: Mesh
    fn
    ft
    kx
    ky

    function Vlasov4d( mesh )

        nx = mesh.nx1
        ny = mesh.nx2

        nvx = mesh.nx3
        nvy = mesh.nx4

        xmin = mesh.x1min 
        xmax = mesh.x1max
        ymin = mesh.x2min
        ymax = mesh.x2max

        kx = 2π/(xmax-xmin)*[0:nx÷2-1;nx÷2-nx:-1]
        ky = 2π/(ymax-ymin)*[0:ny÷2-1;ny÷2-ny:-1]


        fn = zeros(ComplexF64, (nx, ny, nvx, nvy))
        ft = zeros(ComplexF64, (nvx, nvy, nx, ny))
        new( mesh, fn, ft, kx, ky )

    end
   
end


function advection_x1(vlasov, dt)

  nx = vlasov.mesh.nx1
  vxmin  = vlasov.mesh.x3min
  dvx = vlasov.mesh.dx3

  for l=1:nvy, k=1:nvx
     vx = (vxmin +(k-1)*dvx)*dt
     for j=1:ny
        vlasov.f[:,j,k,l] .= ifft( exp(-1im .* vlasov.kx .* vx) .* fft(f[:,j,k,l]))
     end
  end

end 

#euler explicite
# tmp_y = vlasov%tmp_y*(1._f64-cmplx(0.0_f64,vlasov%ky,kind=f64)*vy)
#euler implicite
# tmp_y = vlasov%tmp_y/(1._f64+cmplx(0.0_f64,vlasov%ky,kind=f64)*vy)
#crank-nicolson
# tmp_y = vlasov%tmp_y/(1._f64+cmplx(0.0_f64,vlasov%ky,kind=f64)*vy*0.5_f64)
# tmp_y = vlasov%tmp_y*(1._f64-cmplx(0.0_f64,vlasov%ky,kind=f64)*vy*0.5_f64)
#Euler cn modified
# tmp_y = vlasov%tmp_y*(1._f64-cmplx(0.0_f64,vlasov%ky,kind=f64)*vy-0.5_f64*(vlasov%ky*vy)**2)
function advection_x2(vlasov,dt)

  for l=1:n_l
    vy = (x4_min +(l-1)*delta_x4)*dt
    for k=1:n_k, i=1:n_i
       mul!(f̂, py, f)
       f̂ .= f̂ .* eky .* vy
       ldiv!(f, py, f̂)
       f[i,:,k,l] = d_dy 
    end
  end

end

#function advection_x3x4(vlasov,dt)
#
#  for i=1:n_i, j=1:n_j
#     for k=1:n_k, l=1:n_l
#
#        px = x3_min+(k-1)*delta_x3
#        py = x4_min+(l-1)*delta_x4
#        cthη = cos(bz[i,j]*dt)
#        sthη = sin(bz[i,j]*dt)
#        depvx  = 0.5dt * ex[i,j]
#        depvy  = 0.5dt * ey[i,j]
#        alpha_x[k,l] = - (px - (depvx+(px+depvx)*cthη-(py+depvy)*sthη))
#        alpha_y[k,l] = - (py - (depvy+(px+depvx)*sthη+(py+depvy)*cthη))
#     end
#
#     interpolate_disp(n_k,n_l, ft[i,j,:,:],alpha_x,alpha_y, ft[i,j,:,:])
#  end
#
#end
