/*
**  MCMC estimation for
**  Latent Threshold Vector AutoRegression (LTVAR) model 
*/

#include<oxstd.h>
#include<oxprob.h>
#include<oxfloat.h>
#include<oxdraw.h>
#include<LTVAR.ox>	//LTVAR class

main()
{										
	decl nlag, my, asvar;


	nlag = 2;	//# of lag


	/*--- data load ---*/

	my = loadmat("usdata.xls");	//data

	/*--- some options ---*/
	
	asvar = {"p", "x", "i"};	//variable name
	
	
	/*--- MCMC estimation ---*/

	decl Ltvar = new LTVAR();	//LTVAR class

	Ltvar.SetData(my, nlag);	//set data and lag
	Ltvar.SetVarName(asvar);	//set variable name(*)

	Ltvar.SetIntCept(1);		//set intercept
	Ltvar.SetLTb(1);			//set threshold
	Ltvar.SetLTa(1);			//set threshold

	Ltvar.SetPeriod(1977, 3, 4);//set initial period (*)
								//(year, period, frequency)
	Ltvar.SetImpulse(<4, 8, 12>, -1, 1);
								//set impulse option

	Ltvar.MCMC(50000);			//MCMC estimation
	
	delete Ltvar;
}