allocate(rs(0:nr))
allocate(thetas(0:ntheta))

allocate(dr3inv(0:nr-1))
allocate(drinv(0:nr))
allocate(dsinv(0:ntheta-1))

allocate(ur(0:nr,0:ntheta-1))
allocate(uw(0:nr-1,0:ntheta))

allocate(xs(0:nr-1,0:ntheta-1))
allocate(ys(0:nr-1,0:ntheta-1))
allocate(c(-1:nr,0:ntheta-1))
allocate(advis(-1:nr,0:ntheta-1))
