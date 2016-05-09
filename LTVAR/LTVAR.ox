/*----------------------------------------------------------
                        LTVAR.ox
----------------------------------------------------------*/
											 
/*
**  MCMC estimation for
**  Latent Threshold Vector AutoRegression (LTVAR) model 
*/

#include<oxstd.h>
#include<oxprob.h>
#include<oxfloat.h>
#include<oxdraw.h>
#import <maximize>

class LTVAR
{
	decl m_irs, m_nl, m_fli, m_flLTb, m_flLTa, m_asvar;
	decl m_vimp, m_vcum, m_flib, m_vidb, m_vpd, m_io;
	decl m_my, m_mimpm, m_nf, m_mf;
	
	LTVAR();				//constructer

	SetRanseed(const irs);			//set ranseed

	SetData(const my, const nl);	//set data and lag
	SetIntCept(const fli);			//set intercept
	SetLTb(const flLT);				//set threshold for b
	SetLTa(const flLT);				//set threshold	for a
	SetImpulse(const vimp, const vcum, const flib);
									//set impulse option
	SetFcst(const nf);				//set forcast periods
	GetFcst();						//return forcasts
	SetVarName(const asvar);		//set variable name
	SetPeriod(const iy, const ip, const np);		  
									//set data period
	SetOutput(const io);			//set output
	SetOutputB(const vid);			//set output of b
	
	MCMC(const nsim);				//MCMC estimation
	
	fSSmoother(const my, const amZ, const amG2, const vWb,
			   const mT, const mH2, const va0, const mH02);
			   						 //simulation smoother
	fSampLTBv(const my, const amX, const mPhi, const vmu,
			  const mSig, const amG2inv, const vd,
			  const mbo);  //sampling beta
	fSampSPM(const mp, const mSig, const mPhi, const vmu,
			 const dnu, const dV0, const da0, const db0,
			 const dm0, const ds0, const vd, const dk);
			 						//sampling (Sig,Phi,mu)
	fSampSPG(const my, const mh, const mPhi, const vgam,
			 const dnu, const dV0, const da0, const db0,
			 const dg0, const dG0);//sampling (Sig,Phi,gam)

	fSampLT(const my, const amX, const mp, const amG2inv,
			const vd, const vdc); //sampling threshold

	fSVSamp(const vy, const vh, const dphi, const dsig2,
			const dh00, const dsig02, const nK);
						   //sampling stochastic volatility
	
	ImpRes(const ns, const nk, const nl, const nlen,
		   const mb, const ma, const mh, const vout,
		   const vcum);
								  //comput impulse response

	fTsvar(const vx, const iBm);  //function for 
	fGeweke(const vx, const iBm); //convergence diagnostics
};

LTVAR::LTVAR()
{												   
	/*--- default options ---*/

	m_irs = 1;		//ranseed
	m_fli = 0;		//intersept
	m_flLTb = 1;	//Latent threshold for b (def: on)
	m_flLTa = 1;	//Latent threshold for a (def: on)
	m_vimp = <4, 8, 12>;	//time point of impulse
	m_vcum = -1;			//columns of cumulative
	m_flib = 0;				//impulse band
	m_nf = 0;				//periods to forecast
	m_asvar = <>;			//variable name
	m_vpd = ones(3, 1);		//sample period
	m_io = 1;				//set output (def: on)
	m_vidb = <>;			//output series of b
}

/*
** SetRanseed set ranseed for random generater
*/
LTVAR::SetRanseed(const irs)
{
	m_irs = irs;
}

/*
** SetData sets data
**
** [arguments]
**	 (my, nlag)
**		my:		n*k matrix
**		nlag:	integer
**
*/
LTVAR::SetData(const my, const nlag)
{
	m_my = my;
	m_nl = nlag;
}

/*
** SetIntCept sets intercept
**
** [arguments]
**	 fli  0:without intercept(default), 1:with intercept
**
*/
LTVAR::SetIntCept(const fli)
{
	m_fli = fli;
}

/*
** SetLT sets latent threshold
**
** [arguments]
**	 flLT  0:no-threshold, 1: with threshold (default)
**
*/
LTVAR::SetLTb(const flLT)
{
	m_flLTb = flLT;
}
LTVAR::SetLTa(const flLT)
{
	m_flLTa = flLT;
}

/*
** SetImpulse sets impluse response (random-walk only)
**
** [arguments]
**	 (nimp, vcum)
**		vimp:   time point to draw impulse response
**		vcum:   column No. of cumulative response
**		flib:	compute bands of impulse
*/
LTVAR::SetImpulse(const vimp, const vcum, const flib)
{
	m_vimp = vimp;
	m_vcum = vcum;
	m_flib = flib;
}

/*
** SetFcst sets periods to forecast
*/
LTVAR::SetFcst(const nf)
{
	m_nf = nf;
}

/*
** GetFcst returns forecasts
*/
LTVAR::GetFcst()
{
	return m_mf;
}

/*
** SetVarName sets variable name
*/
LTVAR::SetVarName(const asvar)
{
	m_asvar = asvar;
}

/*
** SetPeriod sets sample period
**
** [arguments]
**	  iy:  year of the first observation
**	  up:  month/quarter(etc.) of the first observation
**	  np:  frequency of the period	
*/
LTVAR::SetPeriod(const iy, const ip, const np)
{
	m_vpd = iy | ip | np;
}

/*
** SetOutput set output
**
** [arguments]
**	 io  0:off, 1:on (default), 2:including dat file
**
*/
LTVAR::SetOutput(const io)
{
	m_io = io;
}

/*
** SetOutputB set output of b
*/
LTVAR::SetOutputB(const vid)
{
	m_vidb = vid;
}

/*
**  fXt sets matrix "X_t" from y(t-1)...y(t-p)
**
**  [output]
**	  mXt	nk*(nk^2*ns) matrix
*/
fXt(const my, const fli)
{											 
	return unit(columns(my))
		   ** ((fli .? 1 .: <>) ~ vecr(reversec(my))');
}

/*
**  fAt sets lower-triangular matrix with free
**  elements inserted from 1*na vector "a"
**
**  [output]
**	  mAt	nk*nk lower-triangular matrix with diagonal
**			elements equal to one.
*/

fAt(const va, const nk)
{
	decl i, mAt = unit(nk);
	
	for(i=1 ; i<nk ; i++)
		mAt[i][:i-1] = va[i*(i-1)/2:(i+1)*i/2-1];

	return mAt;
}

/*
**  fXh sets matrix "X_hat" from a 1*nk vector "y_hat"
**
**  [output]
**	  mXh	nk*na matrix of "X_hat"
*/
fXh(const vyh, const nk, const na)
{
	decl i, mXh = zeros(nk, na);
	
	for(i=1 ; i<nk ; i++)
		mXh[i][i*(i-1)/2:(i+1)*i/2-1] = -vyh[:i-1];

	return mXh;
}

/*
** fK0 computes upper bound of threshold
*/
fK0(const dmu, const dsig2, const dphi, const dk0)
{
	return fabs(dmu)+ dk0 * sqrt(dsig2 / (1 - dphi^2));
}

/*
**  fcor computes correlation vector from covariance matrix
**
**  [output]
**	  vcor	1*na correlation matrix
*/
fcor(const mS, const vid)
{
	decl vd = diagonal(mS);
	
	return vec(mS ./ sqrt(vd ' vd))[vid]';
}

/*
**  franWis returns random value from wishart
**  distribution W(nu, S)
**
**	[NOTE]
**	  "ranwishart" in Ox is miscoded
**
**  [output]
**	  mW	np*np matrix from wishart distribution
*/

franWis(const dnu, const mS)
{
	decl np, mT, i, mW;

	np = rows(mS);
	mT = zeros(np, np);
	
	if(dnu >= np){
		for(i=0 ; i<np ; i++)
			mT[i][i] = sqrt(ranchi(1, 1, dnu - i));

		mT = choleski(mS) * setlower(mT, rann(np, np));
	}
	else
		mT = choleski(mS) * rann(np, dnu);
	
	mW = mT * mT';

	return mW;
}

/*
** MCMC computes MCMC estimation
*/
LTVAR::MCMC(const nsim)
{
	decl ns, nk, nl, nb, na, mb, ma, mh, i, j, k;
	decl mSigb, mSiga, mSigh, myh, mya, amX, amXh, mOmsq;
	decl amOm, amOminv, amG2, amG2inv, vmub, vmua, vgam;
	decl mPhib, mPhia, mPhih, vidi, dnub, dnua, dnuh, vW;
	decl dvb0, dva0, dvh0, dVb0, dVa0, dVh0, dm0, ds0;
	decl da0, db0, dg0, dG0, mSb0, mSa0, mSh0, dk0;
	decl npmt, msamp, msampb, msampa, msamph, msampi;
	decl msampsb, msampsa, msampf, msampdb, msampda;
	decl nf, vym, vbf, vaf, vhf, myf, vbfd, vafd, vyf;
	decl mA, mAinv, dSh, vh, vac, ca, vi, ni, mm, ml, mu;
	decl vdn, vdo, mdn, mdo, dln, dlo, vyh, dbn, mSn;
	decl vdb, vda, mbd, mad, vdbc, vdac, fldo, vidB;
	decl nburn, nK, iBm, iacf, aspar1, aspar2, mpost, ascf;
	decl amimp, vline, id, mimp, mimpi, mimps, sname;
	decl vsamp, nimp, nimax, time = timer();

	ranseed(m_irs);

	
	/*--- set variables ---*/

	ns = sizer(m_my);	//# of sample
	nk = sizec(m_my);	//# of series
	nl = m_nl;			//# of lag
	
	nb = nk * (nk*nl+m_fli);	//# of coefficients in beta
	na = nk * (nk-1) / 2;		//# of parameters in a

	nf = m_nf;			//# of forecasts

	vym = m_fli .? zeros(1, nk) .: meanc(m_my);
	m_my -= vym;
	
	myh = mya = zeros(ns, nk);
	amX = amXh = amOm = amOminv
			   = amG2 = amG2inv = new array[ns];
	for(i=nl ; i<ns ; i++){
		if(nl) amX[i] = fXt(m_my[i-nl:i-1][], m_fli);
		else if(m_fli) amX[i] = unit(nk); 
		if(m_flLTb)
			amOminv[i] = unit(nk);
		else
			amOm[i] = unit(nk);
	}
	
	mb = mbd = zeros(ns, nb);
	ma = mad = zeros(ns, na);
	mh = zeros(ns, nk);

	mSigb = unit(nb) * 0.01;
	mSiga = unit(na) * 0.01;
	mSigh =	unit(nk) * 0.01;

	//for AR case
	vmub = zeros(nb, 1);
	vmua = zeros(na, 1);
	vgam = ones(1, nk);
	mPhib = unit(nb) * 0.95;
	mPhia = unit(na) * 0.95;
	mPhih = unit(nk) * 0.95;
				
	vdb = zeros(nb, 1);		//threshold for b
	vda = zeros(na, 1);		//threshold for a

	vidB = range(0, nb-1);
	if(m_fli){
		vidi = range(0, nk-1) * (nk*nl+1);
		vidB = dropc(vidB, vidi);
	}
	
								 
	/*--- prior ---*/

	dvb0 = 40;		//sigb_i^2 ~ IG(vb0/2, Vb0/2)
	dVb0 = 0.02;
	dva0 = 4;		//siga_i^2 ~ IG(va0/2, Vb0/2)
	dVa0 = 0.02;
	dvh0 = 4;		//sigh_i^2 ~ IG(vh0/2, Vh0/2)
	dVh0 = 0.02;

	dm0 = 0;		//mu ~ N(m0, s0^2)
	ds0 = 1;
	da0 = 20;		//(phi+1)/2 ~ Beta(a0, b0)
	db0 = 1.5;
	dg0 = 6;		//gamma_sv ~ InvGamma(g0s/2, G0s/2)
	dG0 = 0.06;

	dnub = dvb0 + ns - nl;
 	dnua = dva0 + ns - nl;
	dnuh = dvh0 + ns - nl;

	dk0 = 3;		//latent threshold prior
	

	/*--- set sampling option ---*/

	nburn = 0.1 * nsim;			//burn-in period
	npmt = 3 * (nb>0) + 6	   //# of parameter
		 + (m_flLTb + m_flLTa)*2;

	msamp = zeros(nsim, npmt);	//sample box
	if(m_io){
	msampb = zeros(nsim, ns*sizec(m_vidb));
	msampa = zeros(nsim, ns*na);
	msamph = zeros(nsim, ns*nk);
	if(m_fli) msampi = zeros(nsim, ns*nk);
	if(m_flLTb) msampsb = zeros(ns, sizec(m_vidb));
	if(m_flLTa) msampsa = zeros(ns, na);
	}
	if(nl){
	  nimp = sizec(m_vimp);
	  nimax = max(m_vimp) + 1;
	  if(m_flib){
	    m_mimpm = zeros(nsim, ns*nimp*nk^2);
		mimpi = zeros(ns*nimp, nk^2);
	  }
	  else m_mimpm = zeros(ns*nimp, nk^2);
	}

	msampf = zeros(nsim, nf*nk);

	vac = zeros(4, 1);
	nK = int(ns/30);		//# of blocks for sampling h


	
	/*--- MCMC sampling ---*/
	
	println("\nIteration:");
							  
	/*----------- S A M P L I N G   S T A R T -----------*/

	for(k=-nburn ; k<nsim ; k++){
				  
	if(nb){

	/*--- sampling beta ---*/
	

	if(m_flLTb){	//threshold case

	  [mb[nl:][], ca]
	    = fSampLTBv(m_my[nl:][], amX[nl:], mPhib, vmub,
	  				mSigb, amOminv[nl:], vdb, mb[nl:][]);

	  if(k>=0) vac[0] += ca;
	}

	else{	//no-threshold case

	  vW = (unit(nb) - mPhib) * vmub;
	  mSb0 = mSigb * invert(unit(nb) - mPhib.^2);

	  mb[nl:][]
	   = fSSmoother(m_my[nl:][], amX[nl:], amOm[nl:], vW,
					mPhib, mSigb, vmub, mSb0);
	}

	
	/*--- sampling (Sigb, phib, mub) ---*/

	[mSigb, mPhib, vmub]
	 = fSampSPM(mb[nl:][], mSigb, mPhib, vmub, dnub, dVb0,
	 			da0, db0, dm0, ds0, vdb, dk0);



	/*--- sampling threshold of beta ---*/

	if(m_flLTb){
	
	  vdbc = fabs(vmub)'
		   + dk0 * sqrt( diagonal( mSigb
					   * invert( unit(nb) - mPhib.^2 ) ) );

	  [vdb, ca] = fSampLT(m_my[nl:][], amX[nl:], mb[nl:][],
						  amOminv[nl:], vdb, vdbc);
						  
	  if(k>=0) vac[1] += ca;
	}

	mbd = fabs(mb) .< vdb' .? 0 .: mb;

	}	//if(nb) end

	
	/*--- sampling a ---*/

	for(i=nl ; i<ns ; i++){
		if(nb) myh[i][] = m_my[i][] - mbd[i][] * amX[i]';
		else   myh[i][] = m_my[i][];
		amXh[i] = fXh(myh[i][], nk, na);
		if(m_flLTa)
			amG2inv[i] = diag(exp(-mh[i][])./vgam);
		else
			amG2[i] = diag(exp(mh[i][]) .* vgam);
	}

	if(m_flLTa){	//threshold case

	  [ma[nl:][], ca]
		= fSampLTBv(myh[nl:][], amXh[nl:], mPhia, vmua,
	 				mSiga, amG2inv[nl:], vda, ma[nl:][]);

	  if(k>=0) vac[2] += ca;
	}
	
	else{	//no-threshold case
	
	  vW = (unit(na) - mPhia) * vmua;
	  mSa0 = mSiga * invert(unit(na) - mPhia.^2);

	  ma[nl:][]
	   = fSSmoother(myh[nl:][], amXh[nl:], amG2[nl:], vW,
	 				mPhia, mSiga, vmua, mSa0);
	}
	
	/*--- sampling (Siga, phia, mua) ---*/

	[mSiga, mPhia, vmua]
	 = fSampSPM(ma[nl:][], mSiga, mPhia, vmua, dnua, dVa0,
	 			da0, db0, dm0, ds0, vda, dk0);

		
	/*--- sampling threshold of a ---*/

	if(m_flLTa){

	  vdac = fabs(vmua)'
		   + dk0 * sqrt( diagonal( mSiga
					   * invert( unit(na) - mPhia.^2 ) ) );

	  [vda, ca] = fSampLT(myh[nl:][], amXh[nl:], ma[nl:][],
						  amG2inv[nl:], vda, vdac);

	  if(k>=0) vac[3] += ca;
	}
										
	mad = fabs(ma) .< vda' .? 0 .: ma;


	/*--- sampling h ---*/

	for(i=nl ; i<ns ; i++)
		mya[i][] = myh[i][] * fAt(mad[i][], nk)';

	mSh0 =	mSigh * invert(unit(nk) - mPhih.^2);
	
	for(i=0 ; i<nk ; i++)
		mh[nl:][i]
		 = fSVSamp(mya[nl:][i]/sqrt(vgam[i]), mh[nl:][i],
		 		   mPhih[i][i],	mSigh[i][i], 0,
				   mSh0[i][i], nK);

				   
	/*--- sampling (Sigh, phih, gam) ---*/

	[mSigh, mPhih, vgam]
	 = fSampSPG(mya[nl:][], mh[nl:][], mPhih, vgam,	dnuh,
	 			dVh0, da0, db0, dg0, dG0);


	for(i=nl ; i<ns ; i++){
		mA = fAt(mad[i][], nk);
		mAinv = invert(mA);
		if(m_flLTb)
		  amOminv[i] = mA ' diag(exp(-mh[i][]) ./ vgam)
				 	 * mA;
		else
		  amOm[i] = mAinv * diag(exp(mh[i][]) .* vgam)
				  * mAinv';
	}

						
	/*------------- sampling parameters end -------------*/

	
	if(k >= 0){

	/*--- storing sample ---*/

	vsamp = <>;
	if(nb){
		vsamp ~= vmub[0] ~ mPhib[0][0] ~ sqrt(mSigb[0][0]);
	}
	vsamp ~= vmua[0] ~ mPhia[0][0] ~ sqrt(mSiga[0][0])
		   ~ vgam[0] ~ mPhih[0][0] ~ sqrt(mSigh[0][0]);
	if(m_flLTb) vsamp ~= vdb[<0,1>]';
	if(m_flLTa) vsamp ~= vda[<0,1>]';

	msamp[k][] = vsamp;

	if(m_io){
	  msampa[k][] = vec(ma)';
	  msamph[k][] = vec(mh + log(vgam))';
	  if(sizec(m_vidb))
		msampb[k][] = vec(mb[][m_vidb])';
	  if(m_fli) msampi[k][] = vec(mb[][vidi])';	
	  if(m_flLTb) msampsb += fabs(mb[][m_vidb])
						  .< vdb[m_vidb]';
	  if(m_flLTa) msampsa += fabs(ma) .< vda';
	}
	
	if(nl){
	
	/*--- impulse response ---*/
	
	amimp = ImpRes(ns, nk, nl, nimax, mbd[][vidB], mad,
				   mh, m_vimp, m_vcum);
	 
	for(i=0 ; i<nk ; i++){
	  for(j=0 ; j<nk ; j++){
	    if(m_flib)
		  mimpi[][i*nk+j] = vec(amimp[i*nk+j]);
		else
		  m_mimpm[][i*nk+j] += vec(amimp[i*nk+j]);
	  }
	}
	if(m_flib) m_mimpm[k][] = vec(mimpi)';


	if(nf){
							   
	/*--- forecasting ---*/

	vbf = mb[ns-1][]';
	vaf = ma[ns-1][]';
	vhf = mh[ns-1][]';
	myf = m_my[ns-nl:ns-1][];
	
	for(i=0 ; i<nf ; i++){
		vbf = vmub + mPhib*(vbf - vmub)
			+ choleski(mSigb) * rann(nb, 1);
		vbfd = fabs(vbf) .< vdb .? 0 .: vbf;

		vaf = vmua + mPhia*(vaf - vmua)
			+ choleski(mSiga) * rann(na, 1);
		vhf = mPhih * vhf
			+ choleski(mSigh) * rann(nk, 1);
		vafd = fabs(vaf) .< vda .? 0 .: vaf;
		mOmsq = invert(fAt(vafd, nk))
			  * diag(exp(vhf'/2) .* sqrt(vgam));
	  
		vyf = fXt(myf, m_fli) * vbfd + mOmsq * rann(nk, 1);
		msampf[k][range(0,nk-1)+i*nk] = vyf';
		myf = (myf | vyf')[1:][];
	}
	
	}	//if(m_nf) end
	
	}	//if(nl) end
	
	}	//if(k>=0) end
	
	
	if(!imod(k,10000)) println(k);	//print counter
										  
	}		 

    /*----------- I T E R A T I O N   E N D -------------*/
	

	/*--- output result ---*/
	
	iBm = iacf = min(500, nsim/2);	//bandwidth

	aspar1 = aspar2 = {};
	if(nb){
		aspar1 ~= {"mub", "phib", "sigb"};
		aspar2 ~= {"$\\mu_b$", "$\\phi_b$", "$\\sigma_b$"};
	}
	aspar1 ~= {"mua", "phia", "siga",
			   "gam", "phih", "sigh"};
	aspar2 ~= {"$\\mu_a$", "$\\phi_a$", "$\\sigma_a$",
			   "$\\gamma$", "$\\phi_h$","$\\sigma_h$"};
	if(m_flLTb){
		aspar1 ~= {"db1", "db2"};
		aspar2 ~= {"$d_{b1}$", "$d_{b2}$"};
	}
	if(m_flLTa){
		aspar1 ~= {"da1", "da2"};
		aspar2 ~= {"$d_{a1}$", "$d_{a2}$"};
	}
							   
	print("\n================================"
		  "=================================="
		  "\n                        ",
		  "ESTIMATION RESULT"
		  "\n--------------------------------"
		  "----------------------------------"
		  "\nParameter    Mean     Stdev     "
		  " 95%L      95%U   Geweke     Inef."  
		  "\n--------------------------------"
		  "----------------------------------");

	mpost = zeros(npmt, 6);
	for(i=0 ; i<npmt ; i++)
		mpost[i][]
			= meanc(msamp[][i])
			~ sqrt(varc(msamp[][i]))
			~ quantilec(msamp[][i], <.025, .975>)'
			~ fGeweke(msamp[][i],iBm)
			~ fTsvar(msamp[][i], iBm) / varc(msamp[][i]);

	ascf = {"%11.4f", "%10.4f", "%10.4f",
			"%10.4f", "%9.3f", "%9.2f"};
	print("%r", aspar1, "%cf", ascf, mpost);

	print("================================="
		  "=================================\n",
		  "Intercept: ", m_fli .? "ON " .: "OFF",
		  "        Lag: ", nl, "\nThreshold(b): ");
	if(m_flLTb) print("ON      AccRate(b): ",
					 "%5.1f", vac[0]/(nsim*ns)*100);
	else print("OFF");
	if(m_flLTb) print("\n                      ",
			 "AccRate(d): ","%5.1f", vac[1]/(nsim*nb)*100);
	print("\nThreshold(a): ");
	if(m_flLTa) print("ON      AccRate(a): ",
					 "%5.1f", vac[2]/(nsim*ns)*100);
	else print("OFF");
	if(m_flLTa) print("\n                      ",
			 "AccRate(d): ","%5.1f", vac[3]/(nsim*na)*100);

	/*--- draw MCMC result ---*/

	if(nb){
	
	vi = range(0, 5+(nb>0)*3);
	
	ni = sizec(vi);
	DrawAdjust(ADJ_AREAMATRIX, 3, ni);
	for(i=0 ; i<ni ; i++){
 		DrawCorrelogram(i, msamp[][vi[i]]', {""}, iacf);
    	DrawTitle(i, aspar2[vi[i]]);
    	DrawMatrix(i+ni, thinr(msamp[][vi[i]], 200)',
				   {""}, 1, nsim/200);
    	DrawTitle(i+ni, aspar2[vi[i]]);
		DrawDensity(i+ni*2, msamp[][vi[i]]', {""}, 1);
    	DrawTitle(i+ni*2, aspar2[vi[i]]);
	} 
	SaveDrawWindow("samp.gwg");
  	CloseDrawWindow();

	}  //if(nb) end

	//draw threshold
	if(m_flLTb || m_flLTa){
	vi = range(npmt-(1+m_flLTb*m_flLTa)*2, npmt-1);
	ni = sizec(vi);
	DrawAdjust(ADJ_AREAMATRIX, 3, ni);
	for(i=0 ; i<ni ; i++){
 		DrawCorrelogram(i, msamp[][vi[i]]', {""}, iacf);
    	DrawTitle(i, aspar2[vi[i]]);
    	DrawMatrix(i+ni, thinr(msamp[][vi[i]], 200)',
				   {""}, 1, nsim/200);
    	DrawTitle(i+ni, aspar2[vi[i]]);
		DrawDensity(i+ni*2, msamp[][vi[i]]', {""}, 0, 1);
    	DrawTitle(i+ni*2, aspar2[vi[i]]);
	} 
	SaveDrawWindow("samp_d.gwg");
  	CloseDrawWindow();
	}
	
	if(m_io){
	
	if(sizec(m_vidb)){
	
	//draw beta
	ca = sizec(m_vidb);
	mm = shape(meanc(msampb), ns, ca); 
	ml = shape(quantilec(msampb, .05), ns, ca);
	mu = shape(quantilec(msampb, .95), ns, ca);
	if(m_io == 2) savemat("samp_b.xls", mm ~ ml ~ mu);
	for(i=0 ; i<ca ; i++){
	  DrawTMatrix(i, (mm[nl:][i]~ml[nl:][i]~mu[nl:][i])',
				  i .? "" .: {"pos.mean", "5%", "95%"},
				  m_vpd[0], m_vpd[1]+nl, m_vpd[2], 0, 4);
	  DrawTitle(i, sprint("$\\beta_{", m_vidb[i]+1,"t}$"));
	}					   

	if(m_flLTb){
	msampsb /= nsim;
	if(m_io == 2)
	  savemat("samp_sb.xls", msampsb);
	for(i=0 ; i<ca ; i++){
	  DrawHistogram(i+ca, msampsb[nl:][i]', 1, 1, 0, 14);
	  DrawHistogram(i+ca, 1~zeros(1,ns-nl-1), 1, 1, 0, 0);
	  DrawTitle(i+ca, 
	  		i .? sprint("($\\beta_{", m_vidb[i]+1, "t}$)")
			  .: sprint("Pos.Prob. of $s_t=0$ ($\\beta_{",
					    m_vidb[i]+1, "t}$)"));
	}					   
	}
	SaveDrawWindow("samp_b.gwg");
  	CloseDrawWindow();
	
	}	//if(sizec(m_vidb)) end
	
	//draw a
	mm = shape(meanc(msampa), ns, na);
	ml = shape(quantilec(msampa, .05), ns, na);
	mu = shape(quantilec(msampa, .95), ns, na);
	if(m_io == 2)
	  savemat("samp_a.xls", mm ~ ml ~ mu);
	for(i=0 ; i<na ; i++){
	  DrawTMatrix(i, (mm[nl:][i]~ml[nl:][i]~mu[nl:][i])',
				  i .? "" .: {"pos.mean", "5%", "95%"},
				  m_vpd[0], m_vpd[1]+nl, m_vpd[2], 0, 4);
	  DrawTitle(i, sprint("$a_{", i+1, "t}$"));
	}					   

	if(m_flLTa){
	msampsa /= nsim;
	if(m_io == 2) savemat("samp_sa.xls", msampsa);
	for(i=0 ; i<na ; i++){
	  DrawHistogram(i+na, msampsa[nl:][i]', 1, 1, 0, 14);
	  DrawHistogram(i+na, 1~zeros(1,ns-nl-1), 1, 1, 0, 0);
	  DrawTitle(i+na, i .? sprint("($a_{", i+1, "t}$)")
						.: sprint("Pos.Prob. of $s_t=0$ ",
					    "($a_{", i+1, "t}$)"));
	} }					   
	SaveDrawWindow("samp_a.gwg");
  	CloseDrawWindow();

	//draw h
	mm = shape(meanc(msamph), ns, nk);
	ml = shape(quantilec(msamph, .05), ns, nk);
	mu = shape(quantilec(msamph, .95), ns, nk);
	for(i=0 ; i<nk ; i++){
	  DrawTMatrix(i*2, (mm[nl:][i]~ml[nl:][i]~mu[nl:][i])',
				  i .? "" .: {"pos.mean", "5%", "95%"},
				  m_vpd[0], m_vpd[1]+nl, m_vpd[2], 0, 4);
	  DrawTMatrix(i*2+1,
	  			exp((mm[nl:][i]~ml[nl:][i]~mu[nl:][i])/2)',
				i .? "" .: {"pos.mean", "5%", "95%"},
				m_vpd[0], m_vpd[1]+nl, m_vpd[2], 0, 4);
	  DrawTitle(i*2, sprint("$h_{", i+1, "t}$"));
	  DrawTitle(i*2+1, sprint("$\\exp(h_{", i+1,"t}/2)$"));
	}					   
	SaveDrawWindow("samp_h.gwg");
  	CloseDrawWindow();

	//draw intercept
	if(m_fli){
	mm = shape(meanc(msampi), ns, nk);
	ml = shape(quantilec(msampi, .05), ns, nk);
	mu = shape(quantilec(msampi, .95), ns, nk);
	for(i=0 ; i<nk ; i++){
	  DrawTMatrix(i, (mm[nl:][i]~ml[nl:][i]~mu[nl:][i])',
				  i .? "" .: {"pos.mean", "5%", "95%"},
				  m_vpd[0], m_vpd[1]+nl, m_vpd[2], 0, 4);
	  DrawTitle(i, sprint("$c_{", i+1, "t}$"));
	}					   
	SaveDrawWindow("samp_i.gwg");
  	CloseDrawWindow();
	if(m_io == 2) savemat("samp_i.xls", mm ~ ml ~ mu);
	}

	}	//if(m_io) end
	
	if(nl){

	//draw time-varying impulse response

	if(m_flib){
		vline = <5, 16, 84, 95> / 100;
		for(i=0 ; i<4 ; i++)
		  savemat(sprint("impulse", i+1, ".xls"),
		  		  shape(quantilec(m_mimpm, vline[i]),
		  				ns*nimp, nk^2));
		m_mimpm = shape(quantilec(m_mimpm, 0.5),
						ns*nimp, nk^2);
		savemat("impulse0.xls", m_mimpm);
	}
	else
	  m_mimpm /= nsim;

	vline = <2, 8, 6, 12>;

	if(!isarray(m_asvar)) m_asvar = range(1, nk);

	mimps = <>;
	DrawAdjust(ADJ_AREAMATRIX, nk, nk);
	for(i=0 ; i<nk ; i++){
	for(j=0 ; j<nk ; j++){
	  id = i*nk + j;
	  mimp = shape(m_mimpm[][id], ns, nimp);
				
	  for(k=0 ; k<nimp ; k++){
		DrawTMatrix(id, mimp[nl:][k]',
					i+j .? ""
						.: {sprint(m_vimp[k], "-period")},
					m_vpd[0], m_vpd[1]+nl, m_vpd[2], 0,
					vline[k]);
	  }
  	  sname = sprint("$\\varepsilon_{", m_asvar[i], "}{"
					 "\\small\\uparrow}\\ \\rightarrow "
					 "\\ ", m_asvar[j], "$");
	  DrawTitle(i*nk+j, sname);

	  mimps ~= mimp;	
	}
	}
	SaveDrawWindow("impulse.gwg");
  	CloseDrawWindow();
	if(!m_flib) savemat("impulse.xls", mimps);
		
	}	//if(nl) end

	if(nf){

	//draw forecast
	m_mf = reshape(meanc(msampf), nf, nk);
	mm = m_my[ns-4:ns-1][] | m_mf;
	ml = m_my[ns-4:ns-1][]
	   | reshape(quantilec(msampf, .05), nf, nk);
	mu = m_my[ns-4:ns-1][]
	   | reshape(quantilec(msampf, .95), nf, nk);
	for(i=0 ; i<nk ; i++){
	  DrawTMatrix(i, (mm[][i]~ml[][i]~mu[][i])' + vym[i],
				  i .? "" .: {"pos.mean", "5%", "95%"},
				  1, -3, 1, 0, 4);
	  DrawTitle(i, sprint("$y_{", i+1, ",n+t}^f$"));
	}					   
	SaveDrawWindow("samp_f.gwg");
  	CloseDrawWindow();

	m_mf += vym;

	}	//if(nf) end
	
	println("\n\nRanseed: ", m_irs, "\nIteration: ", nsim,
			"\nTime: ", timespan(time));
}
  

/*
**  fSSmoother implements simulation smoother
**  by de Jong & Shephard (1995)
**
**  [NOTE]
**    y_t = Z_t*a_t + G_t*u_t,
**    a_{t+1} = Wb + T_t*a_t + H_t*u_t,  u_t ~ N(0,I)
**    a_1 ~ N(va0, H02)
**
**  [output]
**    malpha:	ns*np sampled state variable
*/

LTVAR::fSSmoother(const my, const amZ, const amG2,
				  const vWb, const mT, const mH2,
				  const va0, const mH02)
{
	decl ns, nk, np, va, vr, mP, mU;
	decl me, amDinv, mK, amL, meta, i;
	decl mC, mCinv, veps, mV, veta0, malpha;
	
	ns = rows(my);		//# of observation
	nk = columns(my);	//# of series
	np = columns(amZ[0]);	//# of state
		 
	va = va0;
	vr = zeros(np, 1);
	mP = mH02;
	mU = zeros(np, np);

	me = zeros(nk, ns);
	amDinv = amL = new array[ns];
	meta = zeros(np, ns);
	
	for(i=0 ; i<ns ; i++){		//Kalman filter
		me[][i] = my[i][]' - amZ[i] * va;
		amDinv[i] = invert(amZ[i]*mP*amZ[i]' + amG2[i]);
		mK = mT * mP * amZ[i] ' amDinv[i];
		amL[i] = mT - mK * amZ[i];
				  
		va = vWb + mT * va + mK * me[][i];
		mP = mT * mP * amL[i]' + mH2;
	}
	
	for(i=ns-1 ; i>=0 ; i--){	//Simulation smoother
		mC = mH2 - mH2 * mU * mH2;
		veps = choleski(mC) * rann(np, 1);
		meta[][i] = mH2 * vr + veps;
		mV = mH2 * mU * amL[i];

		mCinv = invert(mC);
		vr = amZ[i] ' amDinv[i] * me[][i] + amL[i] ' vr
		   - mV ' mCinv * veps;
		mU = amZ[i] ' amDinv[i] * amZ[i]
		   + amL[i] ' mU * amL[i]
		   + mV ' mCinv * mV;
	}

	mC = mH02 - mH02 * mU * mH02;
	veta0 = mH02 * vr + choleski(mC) * rann(np, 1);

	malpha = zeros(np, ns);
	malpha[][0] = va0 + veta0;
	for(i=0 ; i<ns-1 ; i++)
		malpha[][i+1] = vWb + mT * malpha[][i] + meta[][i];

	return malpha';
}
	
/*
**  fSampLTBv samples beta with latent threshold
**  in multi-variate model by single-move sampler
**
**  [NOTE]
**    y_t = X_t*b_t + e_t,    e_t ~ N(0, G2)
**    b_t = B_t*s_t,          s_t = I(B_t>=d)
**
**    B_{t+1} = mu + Phi*(B_t - mu) + eta_t,
**    eta_t ~ N(0, Sig),  eta_0 ~ N(0, Sig/(I-Phi^2))
**
*/
LTVAR::fSampLTBv(const my, const amX, const mPhi,
				 const vmu, const mSig, const amG2inv,
				 const vd, const mbo)
{
	decl ns, nk, np, mba, mSigi, mSigp, vpm1, vpm2;
	decl mSig0i, mSighi, mxh, mSigh, vbh, vbn, vbo, vba_1;
	decl dhn, dho, dln, dlo, dfrac, ca, i;

	ns = rows(my);			//# of sample
	nk = columns(my);		//# of series
	np = columns(amX[0]);	//# of states
	ca = 0;
							   
	mba = mbo;
	
	mSigi = invert(mSig);		//Sig^(-1)
	mSigp = mSigi * (unit(np) + mPhi.^2);
	vpm1 = (unit(np) - mPhi) * vmu;
	vpm2 = (unit(np) - 2 * mPhi + mPhi.^2) * vmu;
		
	mSig0i = (unit(np) - mPhi.^2) * mSigi;  //Sig0^(-1)
		
	for(i=0 ; i<ns ; i++){
					   
	if(!i){

	mSigh = invert( amX[i] ' amG2inv[i] * amX[i]
				    + mSig0i + mSigi * mPhi.^2 );
	vbh = mSigh
		* ( amX[i] ' amG2inv[i] * my[i][]' + mSig0i * vmu
			+ mSigi * mPhi * (mba[i+1][]' - vpm1) );
	}

	else if(i < ns-1){
	 
	mSigh = invert(amX[i] ' amG2inv[i] * amX[i] + mSigp);
	vbh = mSigh
		* ( amX[i] ' amG2inv[i] * my[i][]' 
		 	+ mSigi
			  * ( mPhi * (vba_1 + mba[i+1][]') + vpm2 ) );
	}

	else{

	mSigh = invert(amX[i] ' amG2inv[i] * amX[i]  + mSigi);
	vbh = mSigh
		* ( amX[i] ' amG2inv[i] * my[i][]' 
		 	+ mSigi * ( mPhi * vba_1 + vpm1 ) );
	}

	vbn = vbh + choleski(mSigh) * rann(np, 1);  //candidate
	vbo = mba[i][]';
							
	mSighi = invert(mSigh);
	dhn = dho = -0.5 * log(determinant(mSigh));
	dhn += -0.5 * (vbn - vbh) ' mSighi * (vbn - vbh);
	dho += -0.5 * (vbo - vbh) ' mSighi * (vbo - vbh);

	if(sumc(fabs(vbn) .< vd)){
								
	  mxh = amX[i] .* (fabs(vbn) .>= vd)';
 
	  if(!i){

	  mSigh = invert( mxh ' amG2inv[i] * mxh
				      + mSig0i + mSigi * mPhi.^2 );
	  vbh = mSigh
		  * ( mxh ' amG2inv[i] * my[i][]' + mSig0i * vmu
			  + mSigi * mPhi * (mba[i+1][]' - vpm1) );
	  }
	  
	  else if(i < ns-1){
	  
	  mSigh = invert( mxh ' amG2inv[i] * mxh + mSigp );
	  vbh = mSigh
		  * ( mxh ' amG2inv[i] * my[i][]'
		 	  + mSigi
				* ( mPhi * (vba_1 + mba[i+1][]') + vpm2) );
	  }

	  else{

	  mSigh = invert( mxh ' amG2inv[i] * mxh + mSigi );
	  vbh = mSigh
		  * ( mxh ' amG2inv[i] * my[i][]'
		  	  + mSigi * ( mPhi * vba_1 + vpm1 ) );
	  }

	  dln =	-0.5 * log(determinant(mSigh))
	  		-0.5 * (vbn-vbh) ' invert(mSigh) * (vbn-vbh);
	}
	else dln = dhn;
	
	if(sumc(fabs(vbo) .< vd)){
	
	  mxh = amX[i] .* (fabs(vbo) .>= vd)';
						  
	  if(!i){
	  
	  mSigh = invert( mxh ' amG2inv[i] * mxh
				      + mSig0i + mSigi * mPhi.^2 );
	  vbh = mSigh
		  * ( mxh ' amG2inv[i] * my[i][]' + mSig0i * vmu
			  + mSigi * mPhi * (mba[i+1][]' - vpm1));

	  }
	  
	  else if(i < ns-1){
	  
	  mSigh = invert( mxh ' amG2inv[i] * mxh + mSigp );
	  vbh = mSigh
		  * ( mxh ' amG2inv[i] * my[i][]'
		 	  + mSigi
				* ( mPhi * (vba_1 + mba[i+1][]') + vpm2) );
	  }
	  
 	  else{

	  mSigh = invert( mxh ' amG2inv[i] * mxh  + mSigi );
	  vbh = mSigh
		  * ( mxh ' amG2inv[i] * my[i][]'
		  	  + mSigi * ( mPhi * vba_1 + vpm1 ) );
	  }

	  dlo =	-0.5 * log(determinant(mSigh))
	  		-0.5 * (vbo-vbh) ' invert(mSigh) * (vbo-vbh);
	}
	else dlo = dho;
		
	dfrac = exp(dln - dhn - dlo + dho);

	if(ranu(1, 1) < dfrac){
		mba[i][] = vbn';
		vba_1 = vbn;
		ca++;
	}
	else
		vba_1 = mba[i][]';

	}

	return {mba, ca};
}

/*
**  fSampSPM samples (Sig, Phi, mu)
**  conditional on time-varying parameter
**
**  [NOTE]
**   B_{t+1} = mu + Phi*(B_t - mu) + eta_t,
**   eta_t ~ N(0, Sig),  eta_0 ~ N(0, Sig/(I-Phi^2))
**
*/
LTVAR::fSampSPM(const mp, const mSig, const mPhi, const vmu,
				const dnu, const dV0, const da0, const db0,
				const dm0, const ds0, const vd, const dk0)
{
	decl ns, np, mSign, mPhin, vmun, i, fldo, dup;
	decl vV, vsum, vphii, vsigi, dphin, dphio, dfrac, vmui;
	decl dsign, dmun, mSig0;

	ns = rows(mp);
	np = columns(mp);

	mSign = mSig;
	mPhin = mPhi;
	vmun  = vmu;
					
	/*--- Sampling Sig ---*/

	vV = dV0
	   + (mp[0][] - vmu').^2 * (unit(np) - mPhi.^2)
	   + sumsqrc( mp[1:][] - mp[:ns-2][] * mPhi
	   			  - vmu' * (unit(np) - mPhi) );

	for(i=0 ; i<np ; i++){
	  fldo = 0;
	  do{
	    dsign = 1 / rangamma(1, 1, dnu/2, vV[i]/2);
		dup = fK0(vmu[i], dsign, mPhi[i][i], dk0);
	  }while((fldo < 100) && (vd[i] > dup));

	  if(fldo < 100){
		dfrac = fK0(vmu[i], mSig[i][i], mPhi[i][i], dk0)
			  / dup;
		if(ranu(1, 1) < dfrac)
		  mSign[i][i] = dsign;
	  }
	  
	}


	/*--- sampling Phi ---*/

	vsum = sumsqrc(mp[1:ns-2][] - vmu');
	vphii = sumc((mp[1:][] - vmu') .* (mp[:ns-2][] - vmu'))
		  ./ vsum;
	vsigi = diagonal(mSign) ./ vsum;

	for(i=0 ; i<np ; i++){
		fldo = 0;
		do{
		  dphin = vphii[i] + rann(1, 1) * sqrt(vsigi[i]);
		  dup = fK0(vmu[i], mSign[i][i], dphin, dk0);
		  fldo++;
		}while((fabs(dphin) >= 1 && fldo < 100)
			    || (vd[i] > dup));

		if(fabs(dphin) < 1){
		  dphio = mPhi[i][i];				
		  dfrac = densbeta((dphin+1)/2, da0, db0)
				/ densbeta((dphio+1)/2, da0, db0)
				* sqrt(1 - dphin^2) / sqrt(1 - dphio^2)
				* fK0(vmu[i], mSign[i][i], dphio, dk0)
				/ dup;
		  if(ranu(1, 1) < dfrac)
			mPhin[i][i] = dphin;
		}
	}

	mSig0 = sqrt(mSign * invert(unit(np)-mPhin.^2));

		
	/*--- sampling mu ---*/

	vsigi = 1 ./ (1 / ds0^2
			+ (1 - diagonal(mPhin).^2
		  	   + (ns-1) * (1 - diagonal(mPhin)).^2)
		    			  		./ diagonal(mSign) );
	vmui = vsigi .* 
	  	   (dm0 / ds0^2
	   		+ (mp[0][] .* (1 - diagonal(mPhin).^2)
			 + (1 - diagonal(mPhin))
	   	   	   .* sumc(mp[1:][]
					   - diagonal(mPhin) .* mp[:ns-2][]))
			 ./ diagonal(mSign) );

	for(i=0 ; i<np ; i++){
	  fldo = 0;
	  do{								
	    dmun = (vmui[i] + rann(1, 1) .* sqrt(vsigi[i]))';
	    dup = fabs(dmun) + dk0 * mSig0[i][i];
	    fldo++;
	  }while((fldo < 100) && sumc(vd[i] .> dup));

	  if(fldo < 100){
	    dfrac = (fabs(vmu[i]) + dk0 * mSig0[i][i]) / dup;
		if(ranu(1, 1) < dfrac)
	      vmun[i] = dmun;
	  }
	}

	return {mSign, mPhin, vmun};
}

/*
**  fSampSPG samples (Sig, Phi, gamma)
**  conditional on time-varying parameter
**
**  [NOTE]
**    y_t = sqrt(gam)*exp(h_t/2)eps_t,  eps_t ~ N(0,1)
**    h_{t+1} = Phi*h_t + eta_t,
**    eta_t ~ N(0, Sig),  eta_0 ~ N(0, Sig/(I-Phi^2))
*/
LTVAR::fSampSPG(const my, const mh, const mPhi, const vgam,
				const dnu, const dV0, const da0, const db0,
				const dg0, const dG0)
{
	decl ns, nk, mSign, mPhin, vgamn, i, fldo;
	decl vV, vG, vsum, vphii, vsigi, dphin, dphio, dfrac;

	ns = rows(my);
	nk = columns(my);

	mSign = unit(nk);
	mPhin = mPhi;
	vgamn = vgam;
	
	/*--- Sampling Sig ---*/

	vV = dV0
	   + mh[0][].^2 * (unit(nk) - mPhi.^2)
	   + sumsqrc( mh[1:][] - mh[:ns-2][] * mPhi );

	for(i=0 ; i<nk ; i++)
		mSign[i][i] = 1 / rangamma(1, 1, dnu/2, vV[i]/2);

					   
	/*--- sampling gamma ---*/

	vG = dG0 + sumc(my.^2 ./ exp(mh));

	for(i=0 ; i<nk ; i++)
		vgamn[i] = 1 / rangamma(1, 1, (dg0+ns)/2, vG[i]/2);


	/*--- sampling Phi ---*/

	vsum = sumsqrc(mh[1:ns-2][]);
	vphii = sumc(mh[1:][] .* mh[:ns-2][]) ./ vsum;
	vsigi = diagonal(mSign) ./ vsum;

	for(i=0 ; i<nk ; i++){
		fldo = 0;
		do{
		  dphin = vphii[i] + rann(1, 1) * sqrt(vsigi[i]);
		  fldo++;
		}while(fabs(dphin) >= 1 && fldo < 100);
		if(fabs(dphin) < 1){
		  dphio = mPhi[i][i];				
		  dfrac = densbeta((dphin+1)/2, da0, db0)
				/ densbeta((dphio+1)/2, da0, db0)
				* sqrt(1 - dphin^2) / sqrt(1 - dphio^2);
		  if(ranu(1, 1) < dfrac)
			mPhin[i][i] = dphin;
		}
	}

	return {mSign, mPhin, vgamn};
}

/*
**  fSampLT samples latent threshold
**  conditional on time-varying parameter
**
**  [NOTE]
**    y_t = X_t*b_t + e_t,    e_t ~ N(0, G2)
*/
LTVAR::fSampLT(const my, const amX, const mp,
			   const amG2inv, const vd, const vdc)
{
	decl ns, nk, vdn, vdo, mdn, mdo, dln, dlo;
	decl ve, ca, i, t;

	ns = rows(mp);
	nk = columns(mp);

	vdo = vd;
	ca = 0;
	
	for(i=0 ; i<nk ; i++){
		vdn = vdo;
		vdn[i] = ranu(1, 1) * vdc[i];

		mdn = fabs(mp) .< vdn' .? 0 .: mp;
		mdo = fabs(mp) .< vdo' .? 0 .: mp;
		dln = dlo = 0;
	
		for(t=0 ; t<ns ; t++){
		  ve = my[t][] - mdn[t][] * amX[t]';
		  dln += ve * amG2inv[t] * ve';
		  ve = my[t][] - mdo[t][] * amX[t]';
		  dlo += ve * amG2inv[t] * ve';
		}

		if(ranu(1, 1) < exp(-0.5 * (dln - dlo))){
			vdo = vdn;
			ca++;
		}
	}

	return {vdo, ca};
}

/*
**  fSVSamp implements multi-move sampler for SV model
**  by Shephard & Pitt (1997) and Omori & Watanabe (2004)
**
**  [NOTE]
**    y_t = exp(h_t/2)eps_t,  eps_t ~ N(0,1)
**    h_{t+1} = phi*h_t + eta_t,  eta_t ~ N(0,sig^2)
**    eta_0 ~ N(0, sig0^2)
**
**  [output]
**    vhs:	ns*1 sampled stochastic volatility
*/
LTVAR::fSVSamp(const vy, const vh, const dphi, const dsig2,
			   const dh00, const dsig02, const nK)
{
	decl ns, nite, vhs, vk, ir, id, ird, vyi, i, j, t;
	decl vho, vhn, vhh, vgder1, vgder2, vsiga2, vha, dhrd1;
	decl da, dh0, dr, dU, dH2, dH20, dP;
	decl ve, vDinv, vK, vL, vu, dC, deps, dV, dCinv, du0;
	decl dpron, dposn, dproo, dposo, dfrac, fldo;

	ns = rows(vy);
	nite = 5;		//# of iteration for h_hat
						
	vhs = vh;		//current point
	
	do{
	  vk = 0
		 ~ floor(ns * (range(1, nK)+ranu(1, nK))/(nK+2))
		 ~ ns;		//stochastic knots
	}while(sumc(diff0(vk', 1)[1:] .< 2));

	for(i=0 ; i<nK+1 ; i++){
								
	ir  = vk[i];
	id  = vk[i+1] - vk[i]; 
	ird = vk[i+1] - 1;

	vyi = vy[ir:ird];
	vho = vh[ir:ird];	//current (old) point
	vhn = zeros(id, 1);	//new point		

	if(i<nK) dhrd1 = vh[ird+1];	//h(r+d+1)

	/*--- finding mode & draw candidate ---*/
	
	for(j=0 ; j<=nite ; j++){
	   
	vhh = j ? vhn : vho;		//h_hat

	vgder2 = -0.5 * vyi.^2 ./ exp(vhh); //g''(h)	
	vgder1 = -0.5 - vgder2;				//g'(h)
	vsiga2 = -1 ./ vgder2;				//sig2_ast
	vha = vhh + vsiga2 .* vgder1;		//h_ast

	if(i<nK){
		vsiga2[id-1] = 1 / (-vgder2[id-1] + dphi^2/dsig2);
		vha[id-1]
		 = vsiga2[id-1]
		   * (vgder1[id-1] - vgder2[id-1] * vhh[id-1]
			   			   + dphi*dhrd1/dsig2);
	}		

	/*- simulation smoother	-*/
	
	da = dh0 = i ? dphi*vhs[ir-1] : dh00;
	dH2 = dsig2;						
	dH20 = dP = i ? dsig2 : dsig02;	   
	ve = vDinv = vK = vL = vu = zeros(id, 1);
					  
	for(t=0 ; t<id ; t++){		//Kalman filter
		ve[t] = vha[t] - da;
		vDinv[t] = 1 / (dP + vsiga2[t]);
		vK[t] = dphi * dP * vDinv[t];
		vL[t] = dphi - vK[t];

		da = dphi * da + vK[t] * ve[t];
		dP = dphi * dP * vL[t] + dH2;	
	}
		   
	if(j < nite){	//finding mode
	
	dr = dU = 0;

	for(t=id-1 ; t>=0 ; t--){	//simulation smoother
		dC = dH2 * (1 - dU * dH2);
		deps = 0;
		vu[t] = dH2 * dr + deps;
		dV = dH2 * dU * vL[t];

		dCinv = 1 / dC;
		dr = vDinv[t] * ve[t] + vL[t] * dr
		   - dV * dCinv * deps;
		dU = vDinv[t] + vL[t]^2 * dU + dV^2 * dCinv;
	}
	
	du0 = dH20 * dr;

	vhn[0] = dh0 + du0;
	for(t=0 ; t<id-1 ; t++)
		vhn[t+1] = dphi * vhn[t] + vu[t];

	}	//end finding mode

	else{			//draw candidate

	fldo = 0;

	do{

	dr = dU = 0;

	for(t=id-1 ; t>=0 ; t--){	//simulation smoother
		dC = dH2 * (1 - dU * dH2);
		deps = sqrt(dC) * rann(1, 1);
		vu[t] = dH2 * dr + deps;
		dV = dH2 * dU * vL[t];

		dCinv = 1 / dC;
		dr = vDinv[t] * ve[t] + vL[t] * dr
		   - dV * dCinv * deps;
		dU = vDinv[t] + vL[t]^2 * dU + dV^2 * dCinv;
	}
	
	dC = dH20 * (1 - dU * dH20);
	du0 = dH20 * dr + sqrt(dC) * rann(1, 1);

	vhn[0] = dh0 + du0;
	for(t=0 ; t<id-1 ; t++)
		vhn[t+1] = dphi * vhn[t] + vu[t];
					
	/*- AR step -*/
	
	dpron
	 = sumc(
	 	-0.5 * (vhh + vyi.^2 ./ exp(vhh))	//g(h_hat)
		+ vgder1 .* (vhn - vhh)				//g'(h_hat)
		+ 0.5 * vgder2 .* (vhn - vhh).^2);	//g''(h_hat)

	dposn = sumc(-0.5 * (vhn + vyi.^2 ./ exp(vhn)));
		   
	fldo++;
	}while((ranu(1, 1) >= exp(dposn - dpron))
			&& (fldo < 100));
		   
	}	//end draw candidate

	}	//end finding mode & draw candidate	


	if(fldo < 100){

	/*- MH step	-*/
	
	dproo
	 = sumc(
	 	-0.5 * (vhh + vyi.^2 ./exp(vhh))	//g(h_hat)
		+ vgder1 .* (vho - vhh)				//g'(h_hat)
		+ 0.5 * vgder2 .* (vho - vhh).^2);	//g''(h_hat)

	dposo = sumc(-0.5 * (vho + vyi.^2 ./ exp(vho)));
					 
	dfrac = exp(  dposn + min(dposo|dproo)
				- dposo - min(dposn|dpron) );
					
	if(ranu(1, 1) < dfrac)
		vhs[ir:ird] = vhn;

	}
	
	}	//end sampling all blocks

	return vhs;
}

/*
**  ImpRes computes impulse response function
**  (using average shock size)
**
**  [output]
**    amimp:  (nk^2)*ns*(nlen+1) array matrix for impulse
*/
LTVAR::ImpRes(const ns, const nk, const nl, const nlen,
			  const mb, const ma, const mh, const vout,
			  const vcum)
{
	decl amimp, amOmsq, my, mbs, vh, ve, vimp, nout;
	decl i, j, t;

	amimp = new array[nk^2];
	amOmsq = new array[ns];
	my = zeros(nl+nlen, nk);
	mbs = mb | (ones(nlen, 1) ** mb[ns-1][]);
	nout = sizec(vout);

	vh = meanc(mh[nl:][]);		//average shock size
	for(t=nl ; t<ns ; t++)
	  amOmsq[t] = invert(fAt(ma[t][], nk))
	  			* diag(exp(vh/2));
	  
	for(i=0 ; i<nk ; i++){

	  for(j=0 ; j<nk ; j++)
	  	amimp[i*nk+j] = zeros(ns, nout);
	
	  ve = zeros(nk, 1);
	  ve[i] = 1;		//unit shock
	
	  for(t=nl ; t<ns ; t++){
		my[nl][] = (amOmsq[t] * ve)';

		for(j=nl+1 ; j<nl+nlen ; j++)
		  my[j][] = mbs[t+j-nl][]
		  		  * fXt(my[j-nl:j-1][], 0)';
		  		  
		for(j=0 ; j<nk ; j++){
		  vimp = sumr(vcum.==j) ? cumulate(my[nl:][j])
		  						: my[nl:][j];
		  amimp[i*nk+j][t][] = vimp[vout]';
		}
		
	  }
	}		
		
	return amimp;
}

/*
**  fTsvar & fGeweke return the diagnostics for convergence
*/
LTVAR::fTsvar(const vx, const iBm)
{
	decl ns, vsp;

	ns = rows(vx);
	vsp = periodogram(vx, iBm, 2, 1) / ns;
	
	return M_2PI * vsp[0];
}

LTVAR::fGeweke(const vx, const iBm)
{
	decl ns, n1, n2, vx1, vx2, dx1bar, dx2bar;
	decl dvar1, dvar2, dz;

	ns = rows(vx);
	n1 = floor(0.1*ns);
	n2 = floor(0.5*ns);
	vx1 = vx[:n1-1];
	vx2 = vx[ns-n2:];
	dx1bar = meanc(vx1);
	dx2bar = meanc(vx2);
	dvar1 = fTsvar(vx1, iBm);
	dvar2 = fTsvar(vx2, iBm);
	dz = (dx1bar - dx2bar) / sqrt(dvar1 / n1 + dvar2 / n2);
	
	return 2 * tailn(fabs(dz));
}