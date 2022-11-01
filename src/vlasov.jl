export Vlasov

struct Vlasov

    nx :: Int
    ny :: Int
    nvx :: Int
    nvy :: Int
    fn :: Array{Float64, 4}
    ft :: Array{Float64, 4}
    kx :: Vector{Float64}
    ky :: Vector{Float64}
    dvx :: Float64
    dvy :: Float64

    function Vlasov( mesh )

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

        kx = 2π / (xmax-xmin) * fftfreq(nx, nx)
        ky = 2π / (ymax-ymin) * fftfreq(ny, ny)
        kvx = 2π / (vxmax-vxmin) * fftfreq(nvx, nvx)
        kvy = 2π / (vymax-vymin) * fftfreq(nvy, nvy)

        dvx = (vxmax - vxmin) / nvx
        dvy = (vymax - vymin) / nvy

        fn = zeros(nx, ny, nvx, nvy)
        ft = zeros(nvx, nvy, nx, ny)

        new( nx, ny, nvx, nvy, fn, ft, kx, ky, dvx, dvy )

    end
   
end

export transposevx!

function transposevx!(vlasov :: Vlasov)

     p = (3, 4, 1, 2) # permutation
     permutedims!(vlasov.fn, vlasov.ft, p)

end

export transposexv!

function transposexv!(vlasov :: Vlasov)

     p = (3, 4, 1, 2) # permutation
     permutedims!(vlasov.ft, vlasov.fn, p)

end

function advection_x1(vlasov, dt)

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
        vlasov.fn[:,j,k,l] .= ifft(exp(-1im .* vlasov.kx .* vx) .* fft(f[:,j,k,l]))
     end
  end

end 

function advection_x2( vlasov, dt)

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

function advection_x3(vlasov, ex, dt)

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

function advection_x4(vlasov, ey, dt)

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
