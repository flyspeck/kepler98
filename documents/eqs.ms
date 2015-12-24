{VERSION 2 3 "SUN SPARC SOLARIS" "2.3" }
{USTYLETAB {CSTYLE "Maple Input" -1 0 "Courier" 0 1 255 0 0 1 0 1 0 0 
1 0 0 0 0 }{PSTYLE "Normal" -1 0 1 {CSTYLE "" -1 -1 "Times" 1 14 0 0 
0 0 2 2 2 0 0 0 0 0 0 }0 0 0 -1 -1 -1 0 0 0 0 0 0 -1 0 }{PSTYLE "R3 Fo
nt 0" -1 256 1 {CSTYLE "" -1 -1 "Helvetica" 1 14 0 0 0 0 2 1 2 0 0 0 
0 0 0 }0 0 0 -1 -1 -1 0 0 0 0 0 0 -1 0 }{PSTYLE "R3 Font 2" -1 257 1 
{CSTYLE "" -1 -1 "Courier" 1 14 0 0 0 0 2 2 2 0 0 0 0 0 0 }0 0 0 -1 
-1 -1 0 0 0 0 0 0 -1 0 }}
{SECT 0 {EXCHG {PARA 0 "> " 0 "" {MPLTEXT 1 0 1486 "delta:=(a,b,c,d,e,
f)->a^2*d^2*(-a^2+b^2+c^2-d^2+e^2+f^2)+\n            b^2*e^2*(a^2-b^2+
c^2+d^2-e^2+f^2)+\n            c^2*f^2*(a^2+b^2-c^2+d^2+e^2-f^2)\n    \+
       -b^2*c^2*d^2-a^2*c^2*e^2-a^2*b^2*f^2-d^2*e^2*f^2:\n\nu:=(i,j,k)
->-i^4-j^4-k^4+2*i^2*k^2+2*i^2*j^2+2*j^2*k^2:\n\nrho:=(a,b,c,d,e,f)->-
a^4*d^4-b^4*e^4-c^4*f^4+2*a^2*b^2*d^2*e^2+2*a^2*c^2*d^2*f^2+2*b^2*c^2*
e^2*f^2:\n\nfacecircumradius:=(i,j,k)->i*j*k/sqrt(u(i,j,k)):\n\ncircum
sphereradius:=(a,b,c,d,e,f)->1/2*sqrt(rho(a,b,c,d,e,f)/delta(a,b,c,d,e
,f)):\n\nF:=(x,y,z)->sqrt(facecircumradius(x,y,z)^2-(z/2)^2):\n\nH:=(a
,b,c,d,e,f,i,j,k)->sqrt(circumsphereradius(a,b,c,d,e,f)^2-facecircumra
dius(i,j,k)^2):\n\nVorovol:=(a,b,c,d,e,f)->evalf(1/3*1/2*(\nF(f,b,a)*H
(a,b,c,d,e,f,f,b,a)*a/2+\nF(f,a,b)*H(a,b,c,d,e,f,f,a,b)*b/2+\nF(e,c,a)
*H(a,b,c,d,e,f,e,c,a)*a/2+\nF(e,a,c)*H(a,b,c,d,e,f,e,a,c)*c/2+\nF(d,b,
c)*H(a,b,c,d,e,f,d,b,c)*c/2+\nF(d,c,b)*H(a,b,c,d,e,f,d,c,b)*b/2)):\n\n
phi:=(x,y,z)->arccos((x^2+y^2-z^2)/(2*x*y)):\ntheta:=(i,j,k)->arccos((
cos(i)-cos(j)*cos(k))/(sin(j)*sin(k))):\nalpha:=(a,b,c,d,e,f)->theta(p
hi(b,c,d),phi(a,c,e),phi(a,b,f)):\nbeta:=(a,b,c,d,e,f)->theta(phi(a,c,
e),phi(b,c,d),phi(a,b,f)):\nGamma:=(a,b,c,d,e,f)->theta(phi(a,b,f),phi
(b,c,d),phi(a,c,e)):\n\nDihedral:=(a,b,c,d,e,f)->evalf(alpha(a,b,c,d,e
,f)):\n\nSolid:=(a,b,c,d,e,f)->evalf(alpha(a,b,c,d,e,f)+beta(a,b,c,d,e
,f)+Gamma(a,b,c,d,e,f)-Pi):\nDigits:=16:\n\naslope222TTT:=Vorovol(2,2,
2,2,2,2)/Solid(2,2,2,2,2,2):\n\nnewf:=(a,b,c,d,e,f)->Vorovol(a,b,c,d,e
,f)-aslope222TTT*Solid(a,b,c,d,e,f):" }}}}{MARK "0 0 0" 1477 }
{VIEWOPTS 1 1 0 1 1 1803 }
