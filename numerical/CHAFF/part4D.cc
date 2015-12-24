
//////////
// u135M = - 4 doct u135 D[quoi135+quoin315,x5]
//
static double u135M(double y1,double y3,double y5)
    {
    static const double doc = 0.72090294951746509284124;
    double x1 = y1*y1, x3 = y3*y3, x5=y5*y5;
    double a = y1/2.0, ap= y3/2.0, b = radf(y1,y3,y5), c = 1.255;
    double c2b2 = c*c-b*b;
    if (c2b2<0) return 0;
    double root = -c2b2*sqrt(c2b2);
    double Dquo135 = a*root/(3.0*b*sqrt(b*b-a*a));
    double Dquo315 = ap*root/(3.0*b*sqrt(b*b-ap*ap));
    double u135 = U(x1,x3,x5);
    return (Dquo135+Dquo315)*x1*x3*(x5+x3-x1)*(x5+x1-x3)*(-4.0*doc)/
                (2.0*b*u135);
    }

// The function f(y1,y2) from Section A14 of Sphere Packings IV.
static double crossTerm(double y1,double y2,double lambda)
    {
    double lon=3.2;
    double zp = global::zetapt;
    double t0=1.255;
    double alpha1=dihedraly(y1,t0,y2,lambda,lon,lambda);
    double alpha2=dihedraly(y2,t0,y1,lambda,lon,lambda);
    double sol=solid(y2,t0,y1,lambda,lon,lambda);
    double phi0=phi(t0,t0);
    double phi1=phi(y1/2.0,t0);
    double phi2=phi(y2/2.0,t0);
    double pi=global::pi;
    return 2.0*(sol*(zp-phi0)+
                alpha1*(1-y1/(2*t0))*(-phi1+phi0)+
                alpha2*(1-y2/(2*t0))*(-phi2+phi0))+
        (pi-2.0*alpha2)*cortauRoundLong(y2,lambda)/(2.0*pi)+
        (pi-2.0*alpha1)*cortauRoundLong(y1,lambda)/(2.0*pi);
    }

static double pitau(int k0,int k1,int k2) // Section A17 of Sphere PackingsIV.
	{
	if (k2==0) return 0;
	if ((k0==1)&&(k1==1)&&(k2==1)) return 0.0254;
	return 0.0463 + (double(k0)+2.*double(k2)-3.)*0.008/3. +
			double (k2)*0.0066;
	}
static double pisigma(int k0,int k1,int k2) // Section A18 of SPIV.
	{
	if (k2=0) return 0;
	if ((k0==2)&&(k2==1)) return 1;
	return double(k0+2*k2)*0.008/3. + 0.009*double(k2);
	}

// ALL THE INEQUALITIES FOR PART 4:

		// Section I.A1:
        case 757995764 : *ret = dihedraly(x[2],x[1],x[0],x[5],x[4],x[3])-
            beta(acos(x[0]/2.77),x[0],x[2],x[4]); break;
        case 735258244: *ret =dihedraly(x[2],x[1],x[0],x[5],x[4],x[3]) -
                beta(acos(x[0]/2.51),x[0],x[2],x[4]); break;
		case 343330051:
			*ret=dihedraly(x[1],x[2],x[0],x[4],x[5],x[3])-
				beta(arc(x[0],1.255,1.6),x[0],x[1],x[5]); break;
		case 49446087:
			*ret=dihedraly(x[1],x[2],x[0],x[4],x[5],x[3])-
				beta(arc(x[0],1.255,1.6),x[0],x[1],x[5]); break;
		case 799187442L
			*ret=dihedraly(x[1],x[2],x[0],x[4],x[5],x[3])-
				dihR(y2/2.,radf(y1,y2,y6),y1/(2.*cos(arc(y1,1.255,1.6))));
        case 275706375 : *ret = 0.00005 -
            vorVc(x[0],x[1],x[2],x[3],x[4],x[5], 1.385); break;
        case 324536936 : *ret = 0.00005 -
            vorVc(x[0],x[1],x[2],x[3],x[4],x[5], 1.385); break;
        case 983547118 : *ret = tauVc(x[0],x[1],x[2],x[3],x[4],x[5]) -
			0.0682; break;
        case 206278009 : *ret =  tauVc(x[0],x[1],x[2],x[3],x[4],x[5]) -
            0.0682; break;

		// Section I.A2:
		case 413688580+0 : *ret =
            -4.3223 +4.10113*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-
            gammaNu(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
		case 413688580+1 : *ret =
            -4.3223 +4.10113*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-
            vorNu(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 805296510+0 : *ret =
            -0.9871 +0.80449 *dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-
            gammaNu(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 805296510+1 : *ret =
            -0.9871 +0.80449 *dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-
            vorNu(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
 
        case 136610219+0 : *ret =
            -0.8756 + 0.70186 *dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-
            gammaNu(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 136610219+0 : *ret =
            -0.8756 + 0.70186 *dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-
            vorNu(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
 
        case 379204810+0 : *ret =
            -0.3404 +0.24573*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-
            gammaNu(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 379204810+1 : *ret =
            -0.3404 +0.24573*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-
            vorNu(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
 
        case 878731435+0 : *ret =
            -0.0024+0.00154*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-
            gammaNu(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 878731435+1 : *ret =
            -0.0024+0.00154*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-
            vorNu(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
 
        case 891740103+0 : *ret =
            0.1196- 0.07611*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-
            gammaNu(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 891740103+1 : *ret =
            0.1196- 0.07611*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-
            vorNu(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;



		// Section I.A3:
        case 334002329+0 : *ret =
            -4.42873+4.16523*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+
            tauGnu(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 334002329+1 : *ret =
            -4.42873+4.16523*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+
            tauVnu(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
 
        case 883139937+0 : *ret =
            -1.01104+0.78701*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+
            tauGnu(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 883139937+1 : *ret =
            -1.01104+0.78701*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+
            tauVnu(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
 
        case 50798176+0 : *ret =
            -0.99937+0.77627*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+
            tauGnu(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 50798176+1 : *ret =
            -0.99937+0.77627*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+
            tauVnu(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
 
        case 244435805+0 : *ret =
            -0.34877+0.21916*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+
            tauGnu(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 244435805+1 : *ret =
            -0.34877+0.21916*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+
            tauVnu(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
 
        case 930176500+0 : *ret =
            -0.11434+ 0.05107*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+
            tauGnu(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 930176500+1 : *ret =
            -0.11434+ 0.05107*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+
            tauVnu(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
 
        case 815681339+0 : *ret =
            0.07749 - 0.07106*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+
            tauGnu(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 815681339+1 : *ret =
            0.07749 - 0.07106*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+
            tauVnu(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;

		// Section I.A4:

		case 649592321 : case 649592322: *ret =
            -3.421 + 2.28501*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-
            vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 649592323 : case 649592324: *ret =
            -3.421 + 2.28501*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-
            vorVc(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 600996944 : case 600996945 : *ret =
            -2.616+1.67382*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-
            vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 600996946 : case 600996947 : *ret =
            -2.616+1.67382*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-
            vorVc(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 70667639 : case 70667640: *ret =
            -1.4486+0.8285*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-
            vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 70667641 : case 70667642: *ret =
            -1.4486+0.8285*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-
            vorVc(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 99182343 : case 99182344 : *ret =
            -0.79+0.390925*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-
            vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 99182345 : case 99182346 : *ret =
            -0.79+0.390925*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-
            vorVc(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 578762805 : case 578762806 : *ret =
            -0.3088+0.12012*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-
            vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 578762807 : case 578762808 : *ret =
            -0.3088+0.12012*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-
            vorVc(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 557125557 : case 557125558 : *ret =
            -0.1558 + 0.0501*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-
            vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 557125559 : case 557125560 : *ret =
            -0.1558 + 0.0501*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-
            vorVc(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;

		// Section I.A5:
        case 719735900: case 719735901: *ret =
            -3.3407+2.1747*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+
            tau_analytic(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 719735902: case 719735903: *ret =
            -3.3407+2.1747*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+
            tauVc(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 359616783: case 359616784: *ret =
            -2.945+1.87427*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+
            tau_analytic(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 359616785: case 359616786: *ret =
            -2.945+1.87427*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+
            tauVc(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 440833181 : case 440833182 : *ret =
            -1.5035+0.83046*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+
            tau_analytic(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 440833183 : case 440833184 : *ret =
            -1.5035+0.83046*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+
            tauVc(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 578578364 : case 578578365 : *ret =
            -1.0009 + 0.48263*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+
            tau_analytic(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 578578366 : case 578578367 : *ret =
            -1.0009 + 0.48263*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+
            tauVc(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 327398152 : case 327398153 : *ret =
            -0.7787+0.34833*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+
            tau_analytic(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 327398154 : case 327398155 : *ret =
            -0.7787+0.34833*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+
            tauVc(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 314861952 : case 314861953: *ret =
            -0.4475+0.1694*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+
            tau_analytic(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 314861954 : case 314861955: *ret =
            -0.4475+0.1694*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+
            tauVc(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 234753056 : case 234753057 : *ret =
            -0.2568+0.0822*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+
            tau_analytic(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;
        case 234753058 : case 234753059 : *ret =
            -0.2568+0.0822*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+
            tauVc(x[0],x[1],x[2],x[3],x[4],x[5]);
            break;


		// Section I.A6: inequalities
		case 555481748:
			*ret = -3.58+2.28501*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 615152889:
			*ret = -2.715+1.67382*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 647971645:
			*ret = -1.517+0.8285*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 516606403:
			*ret = -0.858+0.390925*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 690552204:
			*ret = -0.358+0.12012*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 852763473:
			*ret = -0.186+0.0501*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;

		// Section I.A7: inequalities
		case 679673664:
			*ret = -3.48+2.1747*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				+tauVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 926514235:
			*ret = -3.06+1.87427*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				+tauVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 459744700:
			*ret = -1.58+0.83046*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				+tauVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 79400832:
			*ret = -1.06+0.48263*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				+tauVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 277388353:
			*ret = -0.83+0.34833*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				+tauVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 839852751:
			*ret = -0.50+0.1694*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				+tauVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 787458652:
			*ret = -0.29+0.0822*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				+tauVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;

		// Section I.A8: inequalities
		case 499014780:
			*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.23; break;
		case 901845849:
			*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.4167; break;
		case 410091263:
			*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.65; break;
		case 125103581:
			*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-0.956; break;
		case 504968542:
			*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-0.28; break;
		case 770716154:
			*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.714; break;
		case 666090270:
			*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.714; break;
		case 971555266:
			*ret = -dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+2.184; break;

		// Section I.A9: inequalities
		case 956875054:
			*ret = -0.003521 - kappa(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 664200787: 
			*ret = -0.017 - kappa(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 390273147: 
			*ret = -0.017 - kappa(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 654422246: 
			*ret = -0.02274 - kappa(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 366536370: 
			*ret = -0.029 - kappa(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 62532125: 
			*ret = -0.03883 - kappa(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 370631902: 
			*ret = -0.0325 - kappa(x[0],x[1],x[2],x[3],x[4],x[5]); break;

		// Section I.A10: inequalities
		case 214637273: *ret=
			octavorVc(x[0],x[1],x[2],x[3],x[4],x[5])-
			gamma(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 751772680: *ret=
			octavorVc(x[0],x[1],x[2],x[3],x[4],x[5])+0.01561-
			gamma(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 366146051: *ret=
			octavorVc(x[0],x[1],x[2],x[3],x[4],x[5])+0.00935 -
			gamma(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 675766140: *ret=
			octavorVc(x[0],x[1],x[2],x[3],x[4],x[5])+0.00928 -
			gamma(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 520734758: *ret=
			octavorVc(x[0],x[1],x[2],x[3],x[4],x[5]) -
			gamma(x[0],x[1],x[2],x[3],x[4],x[5]); break;

		// Section I.A11: inequalities
		case 378432183: *ret=
			octavorVc(x[0],x[1],x[2],x[3],x[4],x[5]) -
			octavor(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 572206659: *ret=
			octavorVc(x[0],x[1],x[2],x[3],x[4],x[5]) -
			octavor(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 310679005: *ret=
			vorVc(x[0],x[1],x[2],x[3],x[4],x[5])+0.003521 -
			vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 284970880: *ret=
			vorVc(x[0],x[1],x[2],x[3],x[4],x[5])-0.003521 -
			vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 972111620: *ret=
			vorVc(x[0],x[1],x[2],x[3],x[4],x[5])-0.009 -
			vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 875762896: *ret=
			octavorVc(x[0],x[1],x[2],x[3],x[4],x[5]) -
			octavor(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 385332676: *ret=
			octavorVc(x[0],x[1],x[2],x[3],x[4],x[5])- 0.004131 -
			octavor(x[0],x[1],x[2],x[3],x[4],x[5]); break;

		// Section I.A12: inequalities
		case 630927837: *ret=
			tau_analytic(x[0],x[1],x[2],x[3],x[4],x[5])-0.13
				-0.2*(dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-global::pi/2.);
			break;
		case 825907374: *ret=
			tauVc(x[0],x[1],x[2],x[3],x[4],x[5],global::sqrt2)-0.13
				-0.2*(dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-global::pi/2.);
			break;
		case 812894433+1: *ret=
			-0.349+0.024573*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-
			vorNu(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 812894433+2: *ret=
			-0.349+0.024573*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-
            gammaNu(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 404793781+1: *ret=
			-0.0571-vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 404793781+2: *ret=
			-0.0571-vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 404703781+3: *ret=
			-0.0571-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 404703781+4: *ret=
			-0.0571-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;

		// Section I.A13: inequalities
		case 705592875+1: *ret=
			tauVNu(x[0],x[1],x[2],x[3],x[4],x[5]) - 0.033; break;
		case 705592875+2: *ret=
			tauGNu(x[0],x[1],x[2],x[3],x[4],x[5]) - 0.033; break;
		case 747727191: *ret=
			tauVc(x[0],x[1],x[2],x[3],x[4],x[5]) - 0.06585+0.0066;  break;
		case 474496219: *ret=
			0.009-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 649551700: *ret=
			0.0461-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 74657942: *ret=
			-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 897129160: *ret=
			-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 760840103: *ret=
			tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-0.014; break;
		case 675901554: *ret=
			tauVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 712696695: *ret=
			tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-0.06585; break;

		// Section I.A14: inequalities
		case 424011442: *ret=
			-Vee0(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 140881233: *ret=
			-Vee1(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 601456709+0: *ret=
			-Vee0(x[0],x[1],x[2],x[3],x[4],x[5])-0.82*sqrt(421.); break;
		case 601456709+1: *ret=
			-Vee1(x[0],x[1],x[2],x[3],x[4],x[5])-0.82*sqrt(421.); break;
		case 292977281+0: *ret=
			-Vee0(x[0],x[1],x[2],x[3],x[4],x[5])-0.82*sqrt(421.); break;
		case 292977281+1: *ret=
			-Vee1(x[0],x[1],x[2],x[3],x[4],x[5])-0.82*sqrt(421.); break;
		case 927286061+0: *ret=
			-Vee0(x[0],x[1],x[2],x[3],x[4],x[5])-0.5*sqrt(421.); break;
		case 927286061+1: *ret=
			-Vee1(x[0],x[1],x[2],x[3],x[4],x[5])-0.5*sqrt(421.); break;
		case 340409511+0: *ret=
			-Vee0(x[0],x[1],x[2],x[3],x[4],x[5])-0.5*sqrt(421.); break;
		case 340409511+1: *ret=
			-Vee1(x[0],x[1],x[2],x[3],x[4],x[5])-0.5*sqrt(421.); break;
		case 727498658: *ret=
			421.-Deltay(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 484314425: *ret=
			0.82-u135M(x[0],x[2],x[4]); break;
		case 440223030: *ret=
			0.5-u135M(x[0],x[2],x[4]); break;
		case 115756648: *ret=
			crossTerm(x[0],x[1],1.945)-0.887;  break;

		// Section I.A15: inequalities
		case 329882546+0:
		case 427688691+0:
		case 564506426+0:
		case 562103670+0: 
		case 288224597+0: *ret=
			D2Vee0(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 329882546+1:
		case 427688691+1:
		case 564506426+1:
		case 562103670+1:
		case 288224597+1: *ret=
			D2Vee1(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 979916330: *ret=
			D2Vee0(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 749968927: *ret=
			D2Vee1(x[0],x[1],x[2],x[3],x[4],x[5]); break;

		// Section I.A16: inequalities
		case 695180203+1: *ret=
			tau(x[0],x[1],x[2],x[3],x[4],x[5]) - 0.06585; break;
		case 695180203+2: *ret=
			tau_analytic(x[0],x[1],x[2],x[3],x[4],x[5]) - 0.06585; break;
		case 695180203+3: *ret=
			tau_analytic(x[0],x[1],x[2],x[3],x[4],x[5]) - 0.06585; break;
		case 695180203+4: *ret=
			tauVc(x[0],x[1],x[2],x[3],x[4],x[5])+0.0063 - 0.06585; break;
		case 695180203+5: *ret=
			tauVc(x[0],x[1],x[2],x[3],x[4],x[5])+0.0114 - 0.06585; break;
		case 695180203+6: *ret=
			tauGnu(x[0],x[1],x[2],x[3],x[4],x[5]) - 0.06585/2.; break;
		case 695180203+7: *ret=
			tauVnu(x[0],x[1],x[2],x[3],x[4],x[5]) - 0.06585/2.; break;
		case 695180203+8: *ret=
			tauVc(x[0],x[1],x[2],x[3],x[4],x[5],1.385) - 0.06585; break;
		case 695180203+9: *ret=
			tauVc(x[0],x[1],x[2],x[3],x[4],x[5]) - 0.06585; break;

		case 690626704+1: *ret=
			-gamma(x[0],x[1],x[2],x[3],x[4],x[5]) ; break;
		case 690626704+2: *ret=
			-vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]) ; break;
		case 690626704+3: *ret=
			-vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]) ; break;
		case 690626704+4: *ret=
			-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])+0.0063 ; break;
		case 690626704+5: *ret=
			-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])+0.0114 ; break;
		case 690626704+6: *ret=
			-gamma(x[0],x[1],x[2],x[3],x[4],x[5]) ; break;
		case 690626704+7: *ret=
			-vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]) ; break;
		case 690626704+9: *ret=
			-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]) ; break;

		case 807023313: *ret=
			-0.05714-vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]) ; break;
		case 590577214: *ret=
			tau_analytic(x[0],x[1],x[2],x[3],x[4],x[5])-0.13943 ; break;
		case 949210508+1: *ret=
			-0.05714-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]) ; break;
		case 949210508+2: *ret=
			-0.05714-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]) ; break;
		case 671961774+1: *ret=
			tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-0.13943 ; break;
		case 671961774+2: *ret=
			tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-0.13943 ; break;
		case 173606284: *ret=
			tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-0.0592; break;
		case 568664009: *ret=
			0.009-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;

		// Section I.A17: inequalities
		case 645264496+1: *ret=
			tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-pitau(1,1,1)
			-0.13943; break;
		case 645264496+2: *ret=
			tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-pitau(1,0,2)
			-0.13943; break;
		case 645264496+3: *ret=
			tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-pitau(0,3,0)
			-0.21301; break;
		case 645264496+4: *ret=
			tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-pitau(0,2,1)
			-0.21301; break;
		case 645264496+5: *ret=
			tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-pitau(0,1,2)
			-0.21301; break;
		case 645264496+5: *ret=
			tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-pitau(0,0,3)
			-0.21301; break;
		case 910154674: *ret=
			tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-0.034052
			-0.13943; break;
		case 877743345: *ret=
			tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-0.034052-0.0066
			-0.13943; break;


		// Section I.A18: inequalities
		case 612259047+1: *ret=
			-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])-pisigma(1,1,1)-0.05714;
			break;
		case 612259047+2: *ret=
			-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])-pisigma(1,0,2)-0.05714;
			break;
		case 612259047+3: *ret=
			-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])-pisigma(0,3,0)-0.11423;
			break;
		case 612259047+4: *ret=
			-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])-pisigma(0,2,1)-0.11423;
			break;
		case 612259047+5: *ret=
			-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])-pisigma(0,1,2)-0.11423;
			break;
		case 612259047+6: *ret=
			-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])-pisigma(0,0,3)-0.11423;
			break;

		// Section I.A19: inequalities
		case 357477295+1: *ret=
			tauVc(x[0],x[1],x[2],x[3],x[4],x[5])+
			tauVc(x[6],x[1],x[2],x[3],x[7],x[8])-0.235;
			break;
		case 357477295+2: *ret=
			-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])+
			-vorVc(x[6],x[1],x[2],x[3],x[7],x[8])-0.075;
			break;
		case 357477295+3: *ret=
			tauVc(x[0],x[1],x[2],x[3],x[4],x[5])+
			tauVc(x[6],x[1],x[2],x[3],x[7],x[8])-0.3109;
			break;
		case 357477295+4: *ret=
			-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])+
			-vorVc(x[6],x[1],x[2],x[3],x[7],x[8])-0.137;
			break;

		// Section I.A20: inequalities

		case 193776341+1: 
		case 193776341+2: 
		case 193776341+3: 
		case 193776341+4: 
		case 193776341+5: 
			{
			double k0,k1,k2;
			double Z4k;
			switch (INEQ_NUMBER-193776341)
				{
				case 1 : k0=3; k1=1; k2=0; Z4k= -0.05709; break; 
				case 2 : k0=3; k1=0; k2=1; Z4k= -0.05709; break;
				case 3 : k0=2; k1=2; k2=0; Z4k= -0.11418; break;
				case 4 : k0=2; k1=1; k2=1; Z4k= -0.11418; break;
				case 5 : k0=2; k1=0; k2=2; Z4k= -0.11418; break;
				default: cout << "ERROR!! in 193776341" << endl;
				}
				*ret=
			-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
			+Z4k-0.009*k2-(k0+2.*k2)*0.008/3.; break;
			}

		case 898647773+1:
		case 898647773+2:
		case 898647773+3:
		case 898647773+4:
		case 898647773+5:
			{
			double k0,k1,k2;
			double D4k;
			switch (INEQ_NUMBER-898647773)
				{
				case 1 : k0=3; k1=1; k2=0; D4k= 0.20528; break;
				case 2 : k0=3; k1=0; k2=1; D4k= 0.20528; break;
				case 3 : k0=2; k1=2; k2=0; D4k= 0.27886; break;
				case 4 : k0=2; k1=1; k2=1; D4k= 0.27886; break;
				case 5 : k0=2; k1=0; k2=2; D4k= 0.27886; break;
				default: cout << "ERROR!! in 898647773" << endl;
				}
				*ret=
			tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
			-D4k-0.04683-(k0+2.*k2-3.)*0.008/3. -0.0066*k2; break;
			}

		// Section I.A21: inequalities
		case 275286804:*ret=
			-0.008
			-vorVc(2,2,2,x[0],2,2)
			-vorVc(2,2,2,x[1],2,2)
			-vorVc(2,2,2,x[0],x[1],2); break;
		case 627654828:*ret=
			tauVc(2,2,2,x[0],2,2)+tauVc(2,2,2,x[1],2,2)+
			tauVc(2,2,2,x[0],x[1],2)-0.008-0.27113; break;
		case 995177961:*ret=
			-2.*0.008 -0.11408-3.*0.0461-vorVc(2,2,2,x[3],x[4],x[5]); break;
		case 735892048:*ret=
			tauVc(2,2,2,x[3],x[4],x[5])-0.41056-0.06688; break;

			
		// Section I.A22: inequalities
		case 53502142:*ret=
			-3.58-2.28501*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
			-vorVc(x[6],x[1],x[2],x[3],x[7],x[8]); break;
		case 134398524:*ret=
			-2.715+1.67382*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
			-vorVc(x[6],x[1],x[2],x[3],x[7],x[8]); break;
		case 371491817:*ret=
			-1.517+0.8285*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
			-vorVc(x[6],x[1],x[2],x[3],x[7],x[8]); break;
		case 832922998:*ret=
			-0.858+0.390925*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
			-vorVc(x[6],x[1],x[2],x[3],x[7],x[8]); break;
		case 724796759:*ret=
			-0.358+0.009+0.12012*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
			-vorVc(x[6],x[1],x[2],x[3],x[7],x[8]); break;
		case 431940343:*ret=
			-0.186+0.009+0.0501*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
			-vorVc(x[6],x[1],x[2],x[3],x[7],x[8]); break;
		case 980721294:*ret=
			-3.58/2.+2.28501*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 989564937:*ret=
			-2.715/2.+1.67382*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 263355808:*ret=
			-1.517/2.+0.8285*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 445132132:*ret=
			-0.858/2.+0.390925*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 806767374:*ret=
			(-0.358+0.009)/2.+0.12012*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 511038592:*ret=
			(-0.186+0.009)/2.+0.0501*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;

		// Section I.A23: inequalities
		case 4591018:*ret=
			-3.48+2.1747*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+tauVc(x[6],x[1],x[2],x[3],x[7],x[8])
			+tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-0.06585; break;
		case 193728878:*ret=
			-3.06+1.87427*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+tauVc(x[6],x[1],x[2],x[3],x[7],x[8])
			+tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-0.06585; break;
		case 2724096:*ret=
			-1.58+0.83046*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+tauVc(x[6],x[1],x[2],x[3],x[7],x[8])
			+tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-0.06585; break;
		case 213514158:*ret=
			-1.06+0.48263*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+tauVc(x[6],x[1],x[2],x[3],x[7],x[8])
			+tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-0.06585; break;
		case 750768322:*ret=
			-0.83+0.34833*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+tauVc(x[6],x[1],x[2],x[3],x[7],x[8])
			+tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-0.06585; break;
		case 371464244:*ret=
			-0.50+0.1694*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+tauVc(x[6],x[1],x[2],x[3],x[7],x[8])
			+tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-0.06585; break;
		case 657011065:*ret=
			-0.29+0.0014+0.0822*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+tauVc(x[6],x[1],x[2],x[3],x[7],x[8])
			+tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-0.06585; break;
		case 953023504:*ret=
			-3.48/2.+2.1747*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-0.06585/2.; break;
		case 887276655:*ret=
			-3.06/2.+1.87427*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-0.06585/2.; break;
		case 246315515:*ret=
			-1.58/2.+0.83046*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-0.06585/2.; break;
		case 784421604:*ret=
			-1.06/2.+0.48263*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-0.06585/2.; break;
		case 258632246:*ret=
			-0.83/2.+0.34833*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-0.06585/2.; break;
		case 404164527:*ret=
			-0.50/2.+0.1694*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-0.06585/2.; break;
		case 163088471:*ret=
			(-0.29+0.0014)/2.+0.0822*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-0.06585/2.; break;

		case 0000:

// CONSTRAINTS
		// Section C.A1:
		case 275706375 : case 983547118 :
            *ret = global::sqrt2 - radf(x[3],x[4],x[5]); break;
		case 324536936 : case 206278009 : 
            *ret = global::sqrt2 - radf(x[1],x[2],x[3]); break;


		// Section C.A2,C.A3:
		case 413688580: case 805296510: case 136610219:
		case 379204810: case 878731435: case 891740103:
		case 334002329: case 883139937: case 507989176:
		case 244435805: case 930176500: case 815681339:
	            switch (whichFn)
                {
                case 1 : *ret= radf(x[0],x[1],x[5])- global::sqrt2; break;
                case 2 : *ret= radf(x[0],x[2],x[4])- global::sqrt2; break;
                }
            break;
		case 413688580+1: case 805296510+1: case 136610219+1:
		case 379204810+1: case 878731435+1: case 891740103+1:
		case 334002329+1: case 883139937+1: case 507989176+1:
		case 244435805+1: case 930176500+1: case 815681339+1:
            *ret = global::sqrt2 - radf(x[0],x[1],x[5]); break;


		// Section C.A4:
		// Section C.A5:
		        case 649592321 : case 600996944 :
        case 70667639 : case 99182343 : case 578762805 :
        case 557125557 : case 719735900 : case 359616783:
        case 440833181 : case 578578364 : case 327398152:
        case 314861952:  case 234753056 :
            switch(whichFn)
                {
                case 1: *ret= dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-2.46;
                    break;
                case 2: *ret= radf(x[1],x[2],x[3])-global::sqrt2; break;
                case 3: *ret= radf(x[3],x[4],x[5])-global::sqrt2; break;
                }
            break;
        case 649592322 : case 600996945 :
        case 70667640 : case 99182344 : case 578762806 :
        case 557125558 : case 719735901 : case 359616784:
        case 440833182 : case 578578365 : case 327398153:
        case 314861953:  case 234753057 :
        case 555481748 : case 615152889 : case 647971645 :
        case 516606403 : case 690552204 : case 852763473 :
        case 679673664 : case 926514235 : case 459744700 :
        case 79400832 : case 277388353 :
        case 787458652:  case 839852751 :
                *ret= dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-2.46;
            break;
 
        case 649592323 : case 600996946 :
        case 70667641 : case 99182345 : case 578762807 :
        case 557125559 : case 719735902 : case 359616785:
        case 440833183 : case 578578366 : case 327398154:
        case 314861954:  case 234753058 :
            switch(whichFn)
                {
                case 1: *ret= dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-2.46;
                    break;
                case 2: *ret= -radf(x[1],x[2],x[3])+global::sqrt2; break;
                }
            break;
 
        case 649592324 : case 600996947 :
        case 70667642 : case 99182346 : case 578762808 :
        case 557125560 : case 719735903 : case 359616786:
        case 440833184 : case 578578367 : case 327398155:
        case 314861955:  case 234753059 :
            switch(whichFn)
                {
                case 1: *ret= dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-2.46;
                    break;
                case 2: *ret= -radf(x[3],x[4],x[5])+global::sqrt2; break;
                }
            break;
 


		// Section C.A6:
		case 555481748: case 615152889: case 647971645:
		case 516606403: case 690552204: case 852763473:
			*ret=dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-2.46; 
			break;

		// Section C.A7:
		case 679673664: case 926514235: case 459744700:
		case 79400832 : case 277388353: case 839852751:
		case 787458652:
			*ret=dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-2.46; break;

		// Section C.A8:
		// Section C.A9:
		case 664200787 :
			*ret=global::sqrt2-radf(x[1],x[2],x[3]); break;	
		case 390273147:
			*ret=global::sqrt2-radf(x[3],x[4],x[5]); break;
		case 654422246: case 366536370: case 62532125:
		case 370631902:
			*ret=-Deltay(x[0],x[1],x[2],x[3],x[4],x[5]); break;

		// Section C.A10:
		// Section C.A11:
		case 875762896:
			*ret=global::sqrt2-radf(x[0],x[1],x[5]); break;
		case 385332676:
			switch(whichFn) {
			case 1: *ret=radf(x[0],x[1],x[5])-global::sqrt2; break;
			case 2: *ret=-radf(x[0],x[2],x[4])+global::sqrt2; break;
			}
			break;

		// Section C.A12:
		case 630927837:
			*ret=rady(x[0],x[1],x[2],x[3],x[4],x[5])-global::sqrt2;
			break;
		case 825907374:
			*ret=-rady(x[0],x[1],x[2],x[3],x[4],x[5])+global::sqrt2;
			break;
		case 812894433+1:
			*ret=-radf(x[0],x[1],x[5])+global::sqrt2;
			break;
		case 812894433+2:
			switch(whichFn)
				{
				case 1 : *ret=radf(x[0],x[1],x[5])-global::sqrt2; break;
				case 2 : *ret=radf(x[0],x[2],x[4])-global::sqrt2; break;
				} break;
		case 404793781+1:
			*ret=dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-2.2; break;
		case 404793781+2:
			switch(whichFn)
				{
				case 1: *ret=dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-2.2; break;
				case 2 : *ret=radf(x[3],x[4],x[5])-global::sqrt2; break;
				case 3 : *ret=radf(x[3],x[2],x[1])-global::sqrt2; break;
				} break;
		case 404793781+3:
			switch(whichFn)
				{
				case 1: *ret=dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-2.2; break;
				case 2 : *ret=radf(x[3],x[4],x[5])-global::sqrt2; break;
				} break;
		case 404793781+4:
			switch(whichFn)
				{
				case 1: *ret=dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-2.2; break;
				case 2 : *ret=radf(x[3],x[2],x[1])-global::sqrt2; break;
				} break;
			break;

		// Section C.A13:
		case 705592875+2:
			*ret=-radf(x[0],x[2],x[4])+global::sqrt2; break;
			break;

		// Section C.A14:
		case 424011442:
			switch(whichFn) {
			case 1: *ret=x[4]-x[5]; break;
			case 2: *ret=-deltay(x[0],x[1],x[2],x[3],x[4],x[5]); break;
			case 3: *ret=x[3]-x[2]-x[1]; break;
			} break;
		case 140881233:
			switch(whichFn) {
			case 1: *ret=x[4]-x[5]; break;
			case 2: *ret=-deltay(x[0],x[1],x[2],x[3],x[4],x[5]); break;
			case 3: *ret=x[3]-x[2]-x[1]; break;
			} break;
		case 601456709+0:
		case 601456709+1:
			*ret=-deltay(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 292977281+0:
		case 292977281+1:
			switch(whichFn) {
			case 1: *ret=-deltay(x[0],x[1],x[2],x[3],x[4],x[5]); break;
			case 2: *ret=x[3]-x[2]-x[1]; break;
			} break;
		case 927286061+0:
		case 927286061+1:
			*ret=-deltay(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 340409511+0:
		case 340409511+1:
			switch(whichFn) {
			case 1: *ret=-deltay(x[0],x[1],x[2],x[3],x[4],x[5]); break;
			case 2: *ret=x[3]-x[2]-x[1]; break;
			} break;
		case 727498658:
			switch(whichFn) {
			case 1: *ret=radf(x[0],x[2],x[4])-1.255; break;
			case 2: *ret=x[3]-x[2]-x[1]; break;
			} break;

		// Section C.A15:

		case 329882546+0:
		case 329882546+1:
		case 427688691+0:
		case 427688691+1:
		case 564506426+0:
		case 564506426+1:
		case 562103670+0:
		case 562103670+1:
		case 288224597+0:
		case 288224597+1:
			switch(whichFn) {
			case 1: *ret=-deltay(x[0],x[1],x[2],x[3],x[4],x[5]); break;
			case 2: *ret=x[3]-x[2]-x[1]; break;
			case 3: *ret=x[3]-x[4]-x[5]; break;
			} break;
		case 979916330:
			switch(whichFn) {
			case 1: *ret=-deltay(x[0],x[1],x[2],x[3],x[4],x[5]); break;
			case 2: *ret=x[3]-x[2]-x[1]; break;
			case 3: *ret=D1Vee0(x[0],x[1],x[2],x[3],x[4],x[5]); break;
			} break;
		case 749968927:
			switch(whichFn) {
			case 1: *ret=-deltay(x[0],x[1],x[2],x[3],x[4],x[5]); break;
			case 2: *ret=x[3]-x[2]-x[1]; break;
			case 3: *ret=D1Vee1(x[0],x[1],x[2],x[3],x[4],x[5]); break;
			} break;
			
		// Section C.A16:
		case 695180203+1:
		case 690626704+1:
			break; // gamma, try it w/o constraints.
		case 695180203+2:
		case 690626704+2:
			*ret=global::sqrt2-radf(x[3],x[4],x[5]); break;
			break; // vor etatop>sqrt2
		case 695180203+3:
		case 690626704+3:
			*ret=global::sqrt2-radf(x[0],x[1],x[2]); break;
			break; // vor etaside>sqrt2
		case 695180203+9:
		case 690626704+9:
			*ret=global::sqrt2-radf(x[3],x[4],x[5]); break;
			break; // 

		// Section C.A17:
		// Section C.A18:

		// Section C.A19:
		case 357477295+1:
		case 357477295+2:
		case 357477295+3:
		case 357477295+4:
			*ret=-crossdiag+global::sqrt8; break;
			break;

		// Section C.A20:
		case 193776341+1:
		case 898647773+1: 
		case 193776341+2:
		case 898647773+2: 
		case 193776341+3:
		case 898647773+3: 
		case 193776341+4:
		case 898647773+4: 
		case 193776341+5:
		case 898647773+5: 
			*ret=-crossdiag(x)+3.2; break;

		// Section C.A21:

		// Section C.A22:
		case 53502142: case 134398524: case 371491817:
		case 832922998: case 724796759: case 431940343:
			*ret=dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-2.46;
			break;

		// Section C.A23:
		case 4591018: case 193728878: case 2724096:
		case 213514168: case 750768322: case 371464244:
		case 657011065:
			*ret=dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-2.46;
			break;


// BOUNDS:

		// Section B.A1:
		case 757995764 : xmax[1]=xmax[2]=2.23; xmin[3]=2.77;
            xmax[3]=global::sqrt8; break;
		case 735258244 : 
			xmin[3]=xmax[3]=3.2;
			xmin[4]=xmax[4]=2.51;
			xmin[5]=xmax[5]=2.;
            break;
		case 34330051:
			xmax[3]=3.2;
			xmin[4]=xmax[4]=2.51;
			xmin[5]=xmax[5]=2.51;
			break;
		case 49446087:
			xmax[3]=3.2;
			xmin[4]=xmax[4]=3.2;
			xmax[5]=2.;
			break;
		case 799187442:
			xmax[0]=2.2;
			xmin[2]=xmax[2]=2.51;
			xmin[3]=xmax[3]=3.2;
			xmin[4]=xmax[4]=3.2;
			break;
        case 275706375 : case 324536936 :
        case 983547118 : case 206278009 :
			xmin[3]=2.77; xmax[3]=global::sqrt8;
            nconstr = 1;  break;

		// Section B.A2:
		case 413688580: case 805296510: case 136610219:
		case 379204810: case 878731435: case 891740103:
			xmin[0]=2.51; xmax[0]=global::sqrt8;
            nconstr=2; 
            break;
		case 413688580+1: case 805296510+1: case 136610219+1:
		case 379204810+1: case 878731435+1: case 891740103+1:
			xmin[0]=2.51; xmax[0]=global::sqrt8;
            nconstr=1; 
            break;

		// Section B.A3:
		case 334002329: case 883139937: case 507989176:
		case 244435805: case 930176500: case 815681339:
            xmin[0]=2.51; xmax[0]=global::sqrt8;
            nconstr=2; 
            break;
		case 334002329+1: case 883139937+1: case 507989176+1:
		case 244435805+1: case 930176500+1: case 815681339+1:
            xmin[0]=2.51; xmax[0]=global::sqrt8;
            nconstr=1; 
            break;


		// Section B.A4:
		// Section B.A5:
        case 649592321 : case 600996944 :
        case 70667639 : case 99182343 : case 578762805 :
        case 557125557 : case 719735900 : case 359616783:
        case 440833181 : case 578578364 : case 327398152:
        case 314861952:  case 234753056 :
            xmin[0]=2.51; xmax[0]=global::sqrt8;
            xmin[3]=2.51; xmax[3]=global::sqrt8;
            nconstr=3; 
            break;
        case 649592322 : case 600996945 :
        case 70667640 : case 99182344 : case 578762806 :
        case 557125558 : case 719735901 : case 359616784:
        case 440833182 : case 578578365 : case 327398153:
        case 314861953:  case 234753057 :
            xmin[0]=2.51; xmax[0]=global::sqrt8;
            xmin[3]=2.51; xmax[3]=2.77;
            nconstr=1; 
            break;
        case 649592323 : case 600996946 :
        case 70667641 : case 99182345 : case 578762807 :
        case 557125559 : case 719735902 : case 359616785:
        case 440833183 : case 578578366 : case 327398154:
        case 314861954:  case 234753058 :
            xmin[0]=2.51; xmax[0]=global::sqrt8;
            xmin[3]=2.77; xmax[3]=global::sqrt8;
            nconstr=2; 
            break;
        case 649592324 : case 600996947 :
        case 70667642 : case 99182346 : case 578762808 :
        case 557125560 : case 719735903 : case 359616786:
        case 440833184 : case 578578367 : case 327398155:
        case 314861955:  case 234753059 :
            xmin[0]=2.51; xmax[0]=global::sqrt8;
            xmin[3]=2.77; xmax[3]=global::sqrt8;
            nconstr=2; 
            break;
 

		// Section B.A6:
		case 555481748: case 615152889: case 647971645:
		case 516606403: case 690552204: case 852763473:
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			nconstr=1; //dih<2.46.
			break;

		// Section B.A7:
		case 679673664: case 926514235: case 459744700:
		case 79400832 : case 277388353: case 839852751:
		case 787458652:
			xmin[0]= 2.51; xmax[0]=global::sqrt8;
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			nconstr=1;
			break;

		// Section B.A8:
		case 499014780:
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			xmin[3]=xmax[3]=2.51;
			break;
		case 901845849:
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			xmin[3]=xmax[3]=global::sqrt8;
			break;
		case 410091263:
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			xmin[3]=xmax[3]=3.2;
			break;
		case 125103581:
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			xmin[3]=xmax[3]=2.;
			break;
		case 504968542:
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			xmin[3]=xmax[3]=2.;
			xmin[4]=2.; xmax[4]=global::sqrt8;
			break;
		case 770716154:
			xmin[0]=2.7; xmax[0]=global::sqrt8;
			xmin[3]=xmax[3]=3.2;
			break;
		case 666090270:
			xmin[0]=2.51; xmax[0]=2.7;
			xmin[3]=xmax[3]=3.2;
			xmax[1]=2.25;
			break;
		case 971555266:
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			break;

		// Section B.A9:
		case 956875054:
			xmin[0]=2.696; xmax[0]=global::sqrt8;
			xmin[1]=2.45; 
			xmin[5]=2.45;
			xmin[3]=xmax[3]=2.77;
			break;
		case 664200787: case 390273147:
			xmin[0]=2.51; xmax[0]=2.696;
			xmin[3]=2.77; xmax[3]=global::sqrt8;
			nconstr=1;
			break;
		case 654422246:
			xmin[0]=2.57; xmax[0]=global::sqrt8;
			xmin[3]=xmax[3]=3.2;
			nconstr=1;
			break;
		case 366536370:
			xmin[0]=2.51; xmax[0]=2.57;
			xmin[3]=xmax[3]=3.2;
			nconstr=1;
			break;
		case 62532125:
			xmin[0]=2.51; xmax[0]=2.57;
			xmin[3]=xmax[3]=3.2;
			xmax[1]=xmax[2]=xmax[4]=xmax[5]=2.25;
			nconstr=1;
			break;
		case 370631902:
			xmin[0]=2.51; xmax[0]=2.57;
			xmin[3]=xmax[3]=3.2;
			xmax[1]=xmax[2]=xmax[4]=2.25;
			nconstr=1;
			break;

		// Section B.A10:
		case 214637273:
			xmin[0]=2.696; xmax[0]=global::sqrt8;
			break;
		case 751772680:
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			break;
		case 366146051:
			xmin[0]=2.57; xmax[0]=global::sqrt8;
			break;
		case 675766140:
			xmin[0]=2.51; xmax[0]=2.57;
			xmin[1]=2.25;
			break;
		case 520734758:
			xmin[0]=2.51; xmax[0]=2.57;
			xmin[1]=2.25;
			xmin[5]=2.25;
			break;


		// Section B.A11:
		case 378432183:
			xmin[0]=2.696; xmax[0]=global::sqrt8;
			xmax[1]=xmax[2]=2.45;
			break;
		case 572206659:
			xmin[0]=2.696; xmax[0]=global::sqrt8;
			xmin[1]=xmin[4]=2.45;
			break;
		case 310679005:
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			break;
		case 284970880:
			xmin[0]=2.696; xmax[0]=global::sqrt8;
			xmin[1]=xmax[5]=2.45;
			xmin[3]=2.51; xmax[3]=2.77;
			break;
		case 972111620:
			xmin[0]=2.51; xmax[0]=2.696;
			xmin[3]=2.51; xmax[3]=global::sqrt8;
			break;
		case 875762896:
			xmin[0]=2.51; xmax[0]=2.57;
			nconstr=1;
			break;
		case 385332676:
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			xmax[2]=2.2;
			nconstr=2;
			break;

		// Section B.A12:
		case 630927837:
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			xmin[1]=2.51; xmax[1]=global::sqrt8;
			nconstr=1;
			break;
		case 825907374:
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			xmin[1]=2.51; xmax[1]=global::sqrt8;
			nconstr=1;
			break;
		case 812894433+1:
			xmin[0]=2.75; xmax[0]=global::sqrt8;
			nconstr=1;
			break;
		case 812894433+2:
			xmin[0]=2.75; xmax[0]=global::sqrt8;
			nconstr=2;//gamma type
			break;
		case 404793781+1:
			xmin[0]=2.51; xmax[0]=2.75;
			xmin[3]=2.51; xmax[3]=2.77;
			nconstr=1;//
			break;
		case 404793781+2:
			xmin[0]=2.51; xmax[0]=2.75;
			xmin[3]=2.77; xmax[3]=global::sqrt2;
			nconstr=3;// both eta<sqrt2
			break;
		case 404793781+3:
			xmin[0]=2.51; xmax[0]=2.75;
			xmin[3]=2.77; xmax[3]=global::sqrt2;
			nconstr=2;// etatop>sqrt2
			break;
		case 404793781+4:
			xmin[0]=2.51; xmax[0]=2.75;
			xmin[3]=2.77; xmax[3]=global::sqrt2;
			nconstr=2;// etaside>sqrt2
			break;

		// Section B.A13:
		case 705592875+1:
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			break;
		case 705592875+2:
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			nconstr=1;// eta(0,2,4)>sqrt2.
			break;
		case 747727191:
			xmin[3]=xmax[3]=global::sqrt8;
			break;
		case 474496219:
			xmin[3]=xmax[3]=global::sqrt8;
			break;
		case 649551700:
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			xmax[0]=xmax[1]=xmax[4]=xmax[5]=2.;
			break;
		case 74657942:
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			xmin[0]=xmax[0]=2.51;
			xmax[1]=xmax[4]=xmax[5]=2.;
			break;
		case 897129160:
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			xmin[2]=2.51;
			xmax[4]=xmax[5]=2.;
			break;
		case 7606840103:
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			xmax[0]=xmax[1]=xmax[4]=xmax[5]=2.;
			break;
		case 675901554:
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			xmax[0]=xmax[1]=xmax[2]=xmax[4]=xmax[5]=2.;
			break;
		case 712696695:
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			xmin[2]=2.51;
			xmax[4]=xmax[5]=2.;
			break;

		// Section B.A14:
		case 424011442:
			xmin[3]=2.; xmax[3]=2.51+2.51;
			xmin[4]=2.; xmax[4]=3.2;
			xmin[5]=2.; xmax[5]=3.2;
			nconstr=3; // y6>y5, del>0, y4<y2+y3.
			break;
		case 140881233:
			xmin[3]=2.; xmax[3]=2.51+2.51;
			xmin[4]=2.; xmax[4]=3.2;
			xmin[5]=2.; xmax[5]=3.2;
			nconstr=3; // y6>y5, del>0, y4<y2+y3.
			break;
		case 601456709+0:
		case 601456709+1:
			xmin[4]=2.; xmax[4]=2.189;
			xmin[3]=global::sqrt8;  xmax[3]=3.2;
			nconstr=1; // del>0
			break;
		case 292977281+0:
		case 292977281+1:
			xmin[4]=2.; xmax[4]=2.189;
			xmin[3]=3.2; xmax[3]=2.51+2.51;
			xmin[5]=2.; xmax[5]=3.2;
			nconstr=2; // del>0, y4<y2+y3;
			break;
		case 927286061+0:
		case 927286061+1:
			xmin[4]=2.189; xmax[4]=2.51;
			xmin[3]=global::sqrt2; xmax[3]=3.2;
			nconstr=1; // del>0, 
			break;
		case 340409511+0:
		case 340409511+1:
			xmin[4]=2.189; xmax[4]=3.2;
			xmin[3]=3.2;  xmax[3]=2.51+2.51;
			xmax[4]=xmax[5]=3.2;
			nconstr=2; // del>0, y4<y2+y3;
			break;
		case 727498658:
			xmin[3]=global::sqrt8; xmax[3]=2.51+2.51;
			xmax[4]=xmax[5]=3.2;
			nconstr=2; // eta(1,3,5)<t0, y4<y2+y3;
			break;
		case 484314425:
			xmax[1]=xmax[3]=xmax[5]=2.;
			break;
		case 440223030:
			xmax[1]=xmax[3]=xmax[5]=2.;
			xmin[4]=2.189;
			break;
		case 115756648:
			xmax[3]=xmax[4]=xmax[5]=2.;
			break;

		// Section B.A15:
		case 329882546+0:
		case 329882546+1:
			xmin[4]=xmax[4]=2.;
			xmin[5]=xmax[5]=2.;
			xmax[3]=2.51+2.51;
			nconstr=3; // del>0, y4<y2+y3,y5+y6.
			break;
		case 427688691+0:
		case 427688691+1:
			xmin[4]=xmax[4]=2.;
			xmin[5]=xmax[5]=2.51;
			xmax[3]=2.51+2.51;
			nconstr=3; // del>0, y4<y2+y3,y5+y6.
			break;
		case 564506426+0:
		case 564506426+1:
			xmin[4]=xmax[4]=2.51;
			xmin[5]=xmax[5]=2.51;
			xmax[3]=2.51+2.51;
			nconstr=3; // del>0, y4<y2+y3,y5+y6.
			break;
		case 562103670+0:
		case 562103670+1:
			xmin[4]=xmax[4]=2.;
			xmin[5]=xmax[5]=global::sqrt8;
			xmax[3]=2.51+2.51;
			nconstr=3; // del>0, y4<y2+y3,y5+y6.
			break;
		case 288224597+0:
		case 288224597+1:
			xmin[4]=xmax[4]=2.51;
			xmin[5]=xmax[5]=global::sqrt8;
			xmax[3]=2.51+2.51;
			nconstr=3; // del>0, y4<y2+y3,y5+y6.
			break;
		case 979916330:
			xmin[4]=xmax[4]=global::sqrt8;
			xmin[5]=xmax[5]=global::sqrt8;
			xmax[3]=2.51+2.51;
			nconstr=3; // del>0, y4<y2+y3, D1f<0..
			break;
		case 749968927:
			xmin[4]=xmax[4]=global::sqrt8;
			xmin[5]=xmax[5]=global::sqrt8;
			xmax[3]=2.51+2.51;
			nconstr=3; // del>0, y4<y2+y3, D1f1<0..
			break;

		// Section B.A16:
		case 695180203+1:
		case 690626704+1:
			xmin[3]=2.51; xmax[3]=global::sqrt8;
			break; // gamma, try it w/o constraints.
		case 695180203+2:
		case 690626704+2:
			xmin[3]=2.51; xmax[3]=global::sqrt8;
			nconstr=1;
			break; // vor etatop>sqrt2
		case 695180203+3:
		case 690626704+3:
			xmin[3]=2.51; xmax[3]=global::sqrt8;
			nconstr=1;
			break; // vor etaside>sqrt2
		case 695180203+4:
		case 690626704+4:
			xmin[3]=2.6; xmax[3]=global::sqrt8;
			xmin[0]=2.2;
			break; // 
		case 695180203+5:
		case 690626704+5:
			xmin[3]=2.7; xmax[3]=global::sqrt8;
			xmax[0]=2.2;
			break; // 
		case 695180203+6:
		case 690626704+6:
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			break; // gamma, try w/0 constraints.
		case 695180203+7:
		case 690626704+7:
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			break; // vor, try w/o constraints.
		case 695180203+9:
		case 690626704+9:
			xmin[3]=2.6; xmax[3]=global::sqrt8;
			xmin[0]=2.2;
			nconstr=1; // vor, top>sqrt2.
			break; // 

		// Section B.A17:
		// Section B.A18:
		case 645264496+1: 
		case 612259047+1:
			xmin[3]=2.; xmax[3]=2.51;
			xmin[4]=2.51;  xmax[4]=global::sqrt8;
			xmin[5]=global::sqrt8; xmax[5]=3.2;
			break;
		case 645264496+2: 
		case 612259047+2:
			xmin[3]=2.; xmax[3]=2.51;
			xmin[4]=global::sqrt8; xmax[4]=3.2;
			xmin[5]=global::sqrt8; xmax[5]=3.2;
			break;
		case 645264496+3: 
		case 612259047+3:
			xmin[3]=2.51; xmax[3]=global::sqrt8;
			xmin[4]=2.51; xmax[4]=global::sqrt8;
			xmin[5]=2.51; xmax[5]=global::sqrt8;
			break;
		case 645264496+4: 
		case 612259047+4:
			xmin[3]=2.51; xmax[3]=global::sqrt8;
			xmin[4]=2.51; xmax[4]=global::sqrt8;
			xmin[5]=global::sqrt8; xmax[5]=3.2;
			break;
		case 645264496+5: 
		case 612259047+5:
			xmin[3]=2.51; xmax[3]=global::sqrt8;
			xmin[4]=global::sqrt8; xmax[4]=3.2;
			xmin[5]=global::sqrt8; xmax[5]=3.2;
			break;
		case 645264496+6: 
		case 612259047+6:
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			xmin[4]=global::sqrt8; xmax[4]=3.2;
			xmin[5]=global::sqrt8; xmax[5]=3.2;
			break;
		case 910154674:
			xmin[4]=2.6; xmax[4]=global::sqrt8;
			xmin[5]=global::sqrt8; xmax[5]=3.2;
			break;
		case 877743345:
			xmin[3]=xmax[3]=2.;
			xmin[4]=xmax[4]=2.51;
			xmin[5]=xmax[5]=3.2;
			break;
		

		// Section B.A19:
		case 357477295+1:
		case 357477295+2:
			//BUGBUG set number of variables to 9.
			xmax[1]=xmax[5]=xmax[6]=xmax[7]=xmax[8]=2.;
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			xmin[4]=2.51; xmax[4]=global::sqrt8;
			nconstr=1;
			break;
		case 357477295+3:
		case 357477295+4:
			//BUGBUG set number of variables to 9.
			xmax[1]=xmax[5]=xmax[6]=xmax[7]=xmax[8]=2.;
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			xmin[4]=global::sqrt8; xmax[4]=3.2;
			nconstr=1;
			break;

		// Section B.A20:
			/* copied from above.
			case 1 : k0=3; k1=1; k2=0; Z4k= -0.05709; break; 
			case 2 : k0=3; k1=0; k2=1; Z4k= -0.05709; break;
			case 3 : k0=2; k1=2; k2=0; Z4k= -0.11418; break;
			case 4 : k0=2; k1=1; k2=1; Z4k= -0.11418; break;
			case 5 : k0=2; k1=0; k2=2; Z4k= -0.11418; break;
			*/
		case 193776341+1:
		case 898647773+1: 
		case 193776341+2:
		case 898647773+2: 
		case 193776341+3:
		case 898647773+3: 
		case 193776341+4:
		case 898647773+4: 
		case 193776341+5:
		case 898647773+5: 
			//BUGBUG set number of variables to 9.
			{
			double b3,b4;
			switch(INEQ_NUMBER) {
				case 193776341+1:
				case 898647773+1: b3=2.; b4=2.51; break;
				case 193776341+2: 
				case 898647773+2: b3=2.; b4=global::sqrt8; break;
				case 193776341+3:
				case 898647773+3: b3=2.51; b4=2.51; break;
				case 193776341+4:
				case 898647773+4: b3=2.51; b4=global::sqrt8; break;
				case 193776341+5:
				case 898647773+5: b3=global::sqrt8; b4=global::sqrt8; break;
				}
			xmin[0]=2.; xmax[0]=2.51;  // a2
			xmin[1]=xmax[1]=2.;
			xmin[2]=xmax[2]=2.;
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			xmin[4]=xmax[4]=2.;
			xmin[5]=xmax[5]=2.;
			xmin[6]=2.; xmax[6]=2.51; // a4
			xmin[7]=xmax[7]=b3;   // b3
			xmin[8]=xmax[8]=b4; // b4
			nconstr=1; // cross diag>3.2.
			}
			break;

		// Section B.A21:
		case 275286804: 
		case 627654828: 
			xmin[0]=global::sqrt8; xmax[0]=3.2;
			xmin[1]=global::sqrt8; xmax[1]=3.2;
			xmax[2]=xmax[3]=xmax[4]=xmax[5]=2.;
			break;
		case 995177961:
		case 735892048:
			xmax[0]=xmax[1]=xmax[2]=2.;
			xmin[3]=xmin[4]=xmin[5]=global::sqrt8;
			xmax[3]=xmax[4]=xmax[5]=3.2;
			break;

		// Section B.A22:
		case 53502142: case 134398524: case 371491817:
		case 832922998: case 724796759: case 431940343:
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			nconstr=1;
			break;
		case 980721294: case 989564937: case 263355808:
		case 445132132: case 806767374: case 511038592:
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			xmin[4]=xmax[4]=2.51;
			break;

		// Section B.A23:
		case 4591018: case 193728878: case 2724096:
		case 213514168: case 750768322: case 371464244:
		case 657011065:
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			nconstr=1;
			break;
		case 953023504: case 887276655: case 246315515:
		case 784421604: case 258632246: case 404164527:
		case 163088471:
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			xmin[4]=xmax[4]=2.51;
			break;

