export compute_charge

function compute_charge(vlasov)

    dvx, dvy = vlasov.dvx, vlasov.dvy
    nvx, nvy = vlasov.nvx, vlasov.nvy
    nx, ny = vlasov.nx, vlasov.ny

    dvxvy = dvx * dvy

    rho = zeros(nx, ny)

    for j=1:nx, i=1:ny
        rho[i,j] = sum(view(vlasov.ft,1:nvx,1:nvy,i,j)) * dvxvy
    end

    return rho

end 
