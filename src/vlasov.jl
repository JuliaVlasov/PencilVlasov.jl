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

        vxmin = mesh.x3min 
        vxmax = mesh.x3max
        vymin = mesh.x4min
        vymax = mesh.x4max

        kx = 2π/(xmax-xmin)*[0:nx÷2-1;nx÷2-nx:-1]
        ky = 2π/(ymax-ymin)*[0:ny÷2-1;ny÷2-ny:-1]
        kvx = 2π/(vxmax-vxmin)*[0:nvx÷2-1;nvx÷2-nvx:-1]
        kvy = 2π/(vymax-vymin)*[0:nvy÷2-1;nvy÷2-nvy:-1]


        new( mesh, fn, ft, kx, ky )

    end
   
end


function advection_x1(fn, vlasov, dt)

  nx  = vlasov.mesh.nx1
  ny  = vlasov.mesh.nx2
  nvx = vlasov.mesh.nx3
  nvy = vlasov.mesh.nx4

  xmin   = vlasov.mesh.x1min
  ymin   = vlasov.mesh.x2min
  vxmin  = vlasov.mesh.x3min
  vymin  = vlasov.mesh.x4min

  dx  = vlasov.mesh.dx1
  dy  = vlasov.mesh.dx2
  dvx = vlasov.mesh.dx3
  dvy = vlasov.mesh.dx4

  for l=1:nvy, k=1:nvx
     vx = (vxmin +(k-1)*dvx)*dt
     for j=1:ny
        vlasov.f[:,j,k,l] .= ifft(exp(-1im .* vlasov.kx .* vx) .* fft(f[:,j,k,l]))
     end
  end

end 

function advection_x2( fn, vlasov,dt)

  nx  = vlasov.mesh.nx1
  ny  = vlasov.mesh.nx2
  nvx = vlasov.mesh.nx3
  nvy = vlasov.mesh.nx4

  xmin  = vlasov.mesh.x1min
  ymin  = vlasov.mesh.x2min
  vxmin = vlasov.mesh.x3min
  vymin = vlasov.mesh.x4min

  dx  = vlasov.mesh.dx1
  dy  = vlasov.mesh.dx2
  dvx = vlasov.mesh.dx3
  dvy = vlasov.mesh.dx4

  for l=1:nvy
    vy = (vymin +(l-1)*dvy)*dt
    for k=1:nvx, i=1:nx
       vlasov.f[i,:,k,l] .= ifft(exp(-1im .* vlasov.kx .* vy) .* fft(vlasov.f[i,:,k,l]))
    end
  end

end

function advection_x3(ft, vlasov, ex, dt)

  nx  = vlasov.mesh.nx1
  ny  = vlasov.mesh.nx2
  nvx = vlasov.mesh.nx3
  nvy = vlasov.mesh.nx4

  xmin   = vlasov.mesh.x1min
  ymin   = vlasov.mesh.x2min
  vxmin  = vlasov.mesh.x3min
  vymin  = vlasov.mesh.x4min

  dx  = vlasov.mesh.dx1
  dy  = vlasov.mesh.dx2
  dvx = vlasov.mesh.dx3
  dvy = vlasov.mesh.dx4

  for j=1:ny, i=1:nx
  for l=1:nvy
       vlasov.ft[i,j,:,l] .= ifft(exp(-1im .* vlasov.kvx .* ex) .* fft(vlasov.f[i,j,:,l]))
    end
  end

end

function advection_x4(ft, vlasov, ey, dt)

  nx  = vlasov.mesh.nx1
  ny  = vlasov.mesh.nx2
  nvx = vlasov.mesh.nx3
  nvy = vlasov.mesh.nx4

  xmin   = vlasov.mesh.x1min
  ymin   = vlasov.mesh.x2min
  vxmin  = vlasov.mesh.x3min
  vymin  = vlasov.mesh.x4min

  dx  = vlasov.mesh.dx1
  dy  = vlasov.mesh.dx2
  dvx = vlasov.mesh.dx3
  dvy = vlasov.mesh.dx4

  for k=1:nvy
    for j=1:ny, i=1:nx
       vlasov.f[i,j,k,:] .= ifft(exp(-1im .* vlasov.kvy .* ey) .* fft(vlasov.f[i,j,k,:]))
    end
  end

end
