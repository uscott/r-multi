
#ifdef FUCKOFF

double MargrabeFcn::Pearson(double S1, double S2, double K, 
                            double vol1, double vol2, double rho, 
                            double r, double q1, double q2, double tau, 
                            LPCTSTR callput, long GrdPts) 
{
  /*clock_t start, finish;
    double duration;
    start = clock();*/
	
  double h, valG, eps = 0.01;
  long N = GrdPts;
  CString ctemp = callput;

  h = 5.0*vol1*sqrt(tau)/N;
  OneStepGridSpread2(N, S1, S2, vol1, vol2, rho, r, q1, q2, K, tau, valG, h);

  if (ctemp == "c")
    return valG;
  else
    return valG-(S2*exp(-q2*tau)-S1*exp(-q1*tau)-K*exp(-r*tau));

  /*finish = clock();
    duration = (double)(finish-start)/CLOCKS_PER_SEC;
    return duration;*/
}

//////////////////////////////////////////////////////////////////////
void OneStepGridSpread2(long N, double S10, double S20, 
                        double sig1, double sig2, 
                        double rho, double r, double q1, double q2, double X, double T, 
                        double &value, double h)
{
  double omega;
  double *mat1, *mat2, *z;
  double m, diffExp, *G1, *G2;
  double a_0, b_0, a_N, b_N;
  long i, N2 = 2*N;

  mat1 = new double [N2];
  mat2 = new double [N2];
  z = new double [N2+1];

  omega = sig1*sqrt(T)/h;

  double *n1 = new double [N2+1];
  double *n2 = new double [N2+1];
  double siggy = sig1*sqrt(T), wind;

  for (i = 0; i <= N; i++) //rows and columns of underlying Toeplitz matrix
    {
      wind = (i-N)/omega;
      n1[i] = NDist(wind-siggy);
      n2[i] = NDist(wind);
    }

  for (i = N+1; i <= N2; i++) 
    {
      n1[i] = NDist((i-N)/omega-siggy);
      n2[i] = 1-n2[N2-i];
    }

  for (i = 0; i < N2; i++) 
    {
      mat1[i] = n1[i+1]-n1[i];
      mat2[i] = n2[i+1]-n2[i];
    }

  delete [] n1;
  delete [] n2;
	
  G1 = new double [N2];
  G2 = new double [N2];
	
  m = log(S10)+(r-q1-0.5*sig1*sig1)*T;

  double *art = new double [N2+1];

  for (i = 0; i <= N2; i++)
    {
      art[i] = exp(m+h*(i-N)); //grid points
      z[i] = Inner(art[i], S10, S20, sig1, sig2, rho, r, q1, q2, T, X);
    }

  for (i = 0; i < N2; i++) //linear interpolation between grid points
    {
      diffExp = art[i+1]-art[i];
      G1[i] = (z[i+1]-z[i])/diffExp;
      G2[i] = (art[i+1]*z[i]-art[i]*z[i+1])/diffExp;
    }

  /*
  //boundary conditions (end terms)
  if (G1[0] > 0.0) //flat for call
  a_0 = 0.0;
  else //smooth join for put
  a_0 = G1[0];

  b_0 = z[0]-a_0*art[0];

  if (G1[N2-1] > 0.0) //smooth join for call
  a_N = G1[N2-1];
  else //flat for put
  a_N = 0.0;

  b_N = z[N2]-a_N*art[N2];
  */

  double A = (r*(1.0-rho*sig2/sig1)-(q2-q1*rho*sig2/sig1)+0.5*rho*sig2*(sig1-rho*sig2))*T;
  a_N = exp(A)*S20/pow(S10, rho*sig2/sig1);
  b_N = z[N2]-a_N*pow(art[N2], rho*sig2/sig1);
  a_0 = exp(A)*S20/pow(S10, rho*sig2/sig1);
  b_0 = z[0]-a_0*pow(art[0], rho*sig2/sig1);

  value = 0.0;

  for (i = 0; i < N2; i++)
    value += mat1[i]*G1[i];

  value *= exp(m+0.5*sig1*sig1*T);

  for (i = 0; i < N2; i++)
    value += mat2[i]*G2[i];

  //boundary terms
  //value += a_0*NDist(-N/omega-siggy)*exp(m+0.5*sig1*sig1*T)+b_0*NDist(-N/omega);
  //value += a_N*(1.0-NDist(N/omega-siggy))*exp(m+0.5*sig1*sig1*T)+b_N*(1.0-NDist(N/omega));
  value += a_0*NDist(-N/omega-rho*sig2*sqrt(T))*exp(rho*sig2*m/sig1+0.5*rho*rho*sig2*sig2*T)+b_0*NDist(-N/omega);
  value += a_N*(1.0-NDist(N/omega-rho*sig2*sqrt(T)))*exp(rho*sig2*m/sig1+0.5*rho*rho*sig2*sig2*T)+b_N*(1.0-NDist(N/omega));

  value *= exp(-r*T);

  delete [] mat1;
  delete [] mat2;
  delete [] z;
  delete [] G1;
  delete [] G2;
  delete [] art;
}

//Inner integration for spread options (from Pearson 1995)
double Inner(double S, double S10, double S20, double sig1, double sig2, double rho, 
             double r, double q1, double q2, double tau, double X)
{
  double A, x1, x2, M2, m1, m2, sig;

  A = (r*(1.0-rho*sig2/sig1)-(q2-q1*rho*sig2/sig1)+0.5*rho*sig2*(sig1-rho*sig2))*tau;
  sig = sig2*sqrt((1.0-rho*rho)*tau);
  m1 = log(S10)+(r-q1-0.5*sig1*sig1)*tau;
  m2 = log(S20)+(r-q2-0.5*sig2*sig2)*tau;
  M2 = m2+(rho*sig2/sig1)*(log(S)-m1);

  if (S+X > 0)
    {
      if (rho*rho == 1)
        {
          if (M2 > log(S+X))
            return exp(A)*S20*pow(S/S10, rho*sig2/sig1)-(S+X);
          else
            return 0;
        }
      else
        {
          x1 = (M2+sig*sig-log(S+X))/sig;
          x2 = x1-sig;

          return exp(A)*S20*pow(S/S10, rho*sig2/sig1)*NDist(x1)-(S+X)*NDist(x2);
        }
    }
  else
    return exp(A)*S20*pow(S/S10, rho*sig2/sig1)-(S+X);
}



#endif
