function R = rotmat(phi,theta,psi)

cphi = cos(phi);
sphi = sin(phi);
ctheta  = cos(theta);
stheta = sin(theta);
cpsi  = cos(psi);
spsi = sin(psi);

R = [ctheta*cpsi, sphi*stheta*cpsi - cphi*spsi, cphi*stheta*cpsi + sphi*spsi;
     ctheta*spsi, sphi*stheta*spsi + cphi*cpsi, cphi*stheta*spsi - sphi*cpsi;
      -stheta,    sphi*ctheta,              cphi*ctheta];

end