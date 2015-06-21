function rhoouter=MantleDensity(M,rcorei,rhocorei,router)

rhoouter=-(3*M-4*pi*(rcorei.^3).*rhocorei)./(4*pi*(rcorei.^3)-4*pi*(router^3));