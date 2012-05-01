#ifdef FUCKOFF
double MargrabeFcn::PearsonGreeks(double S1, double S2, double K, double vol1, 
								  double vol2, double rho, double r, double q1, 
								  double q2, double tau, LPCTSTR callput, 
								  LPCTSTR output, long GrdPts) 
{
	double h, valG, eps = 0.0001;
	long N = GrdPts;
	CString ctemp = callput, otemp = output;

	h = 5.0*vol1*sqrt(tau)/N;
	
	if (otemp == "p")
	{
		OneStepGridSpread2(N, S1, S2, vol1, vol2, rho, r, q1, q2, K, tau, valG, h);

		if (ctemp == "c")
			return valG;
		else
			return valG-(S2*exp(-q2*tau)-S1*exp(-q1*tau)-K*exp(-r*tau));
	}
	else if (otemp == "d1")
	{
		double valGup, valGdn;

		OneStepGridSpread2(N, S1+eps, S2, vol1, vol2, rho, r, q1, q2, K, tau, valGup, h);
		OneStepGridSpread2(N, S1-eps, S2, vol1, vol2, rho, r, q1, q2, K, tau, valGdn, h);

		if (ctemp == "c")
			return (valGup-valGdn)/(2*eps);
		else
			return (valGup-valGdn)/(2*eps)+exp(-q1*tau);
	}
	else if (otemp == "d2")
	{
		double valGup, valGdn;

		OneStepGridSpread2(N, S1, S2+eps, vol1, vol2, rho, r, q1, q2, K, tau, valGup, h);
		OneStepGridSpread2(N, S1, S2-eps, vol1, vol2, rho, r, q1, q2, K, tau, valGdn, h);

		if (ctemp == "c")
			return (valGup-valGdn)/(2*eps);
		else
			return (valGup-valGdn)/(2*eps)-exp(-q2*tau);
	}
	else if (otemp == "g1")
	{
		double valGup, valGdn;

		OneStepGridSpread2(N, S1+eps, S2, vol1, vol2, rho, r, q1, q2, K, tau, valGup, h);
		OneStepGridSpread2(N, S1-eps, S2, vol1, vol2, rho, r, q1, q2, K, tau, valGdn, h);
		OneStepGridSpread2(N, S1, S2, vol1, vol2, rho, r, q1, q2, K, tau, valG, h);

		return (valGup-2*valG+valGdn)/(eps*eps);
	}
	else if (otemp == "g2")
	{
		double valGup, valGdn;

		OneStepGridSpread2(N, S1, S2+eps, vol1, vol2, rho, r, q1, q2, K, tau, valGup, h);
		OneStepGridSpread2(N, S1, S2-eps, vol1, vol2, rho, r, q1, q2, K, tau, valGdn, h);
		OneStepGridSpread2(N, S1, S2, vol1, vol2, rho, r, q1, q2, K, tau, valG, h);

		return (valGup-2*valG+valGdn)/(eps*eps);
	}
	else if (otemp == "g3")
	{
		double val1, val2, val3, val4;

		OneStepGridSpread2(N, S1+eps, S2+eps, vol1, vol2, rho, r, q1, q2, K, tau, val1, h);
		OneStepGridSpread2(N, S1+eps, S2-eps, vol1, vol2, rho, r, q1, q2, K, tau, val2, h);
		OneStepGridSpread2(N, S1-eps, S2+eps, vol1, vol2, rho, r, q1, q2, K, tau, val3, h);
		OneStepGridSpread2(N, S1-eps, S2-eps, vol1, vol2, rho, r, q1, q2, K, tau, val4, h);

		return (val1-val2-val3+val4)/(4*eps*eps);
	}
	else if (otemp == "th")
	{
		double valGup = 0, valGdn = 0;
		eps = 1.0/365.0;

		//OneStepGridSpread2(N, S1, S2, vol1, vol2, rho, r, q1, q2, K, tau+eps, valGup, h);
		OneStepGridSpread2(N, S1, S2, vol1, vol2, rho, r, q1, q2, K, tau-eps, valGdn, h);
		OneStepGridSpread2(N, S1, S2, vol1, vol2, rho, r, q1, q2, K, tau, valG, h);

		//return (valGup-2*valG+valGdn)/(eps*eps);
		return valGdn-valG;
	}
	else if (otemp == "v1")
	{
		double valGup, valGdn;

		OneStepGridSpread2(N, S1, S2, vol1+eps, vol2, rho, r, q1, q2, K, tau, valGup, h);
		OneStepGridSpread2(N, S1, S2, vol1-eps, vol2, rho, r, q1, q2, K, tau, valGdn, h);

		return (valGup-valGdn)/(2*eps);
	}
	else if (otemp == "v2")
	{
		double valGup, valGdn;

		OneStepGridSpread2(N, S1, S2, vol1, vol2+eps, rho, r, q1, q2, K, tau, valGup, h);
		OneStepGridSpread2(N, S1, S2, vol1, vol2-eps, rho, r, q1, q2, K, tau, valGdn, h);

		return (valGup-valGdn)/(2*eps);
	}
	else if (otemp == "e") //eta
	{
		double valGup, valGdn;

		eps = MIN(eps, (1-rho)/100);

		OneStepGridSpread2(N, S1, S2, vol1, vol2, rho+eps, r, q1, q2, K, tau, valGup, h);
		OneStepGridSpread2(N, S1, S2, vol1, vol2, rho-eps, r, q1, q2, K, tau, valGdn, h);

		return (valGup-valGdn)/(2*eps);
	}
	else
		return -69;
}

#endif
