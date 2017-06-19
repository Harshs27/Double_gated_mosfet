#include<stdio.h>
//#include<conio.h>
//#include<unistd.h>
#include<stdlib.h>
#include<math.h>
//#include<iostream>
#include<time.h>


void main_part(double*,double*,double* ,double* ,double* ,double* ,double ,double ,double , double ,double ,double ,double ,double ,double ,double );
void my_psigzp(double* , double, double,double);
void root_find_trig (double* ,double* , double* ,double* ,double* ,double* ,double* ,double* ,double ,double ,double , double ,double ,double ,double ,double ,double);
void trig_func(double *,double ,double ,double , double ,double ,double ,double);
void trig_limit (double* ,double* ,double* ,double ,double ,double ,double ,double ,double , double );
void hyp_limit(double*,double* ,double * ,double , double ,double ,double ,double ,double ,double ,double ,double );
void root_find_hyp(double *,double *,double *, double *, double *, double *, double *, double *,double , double , double ,double , double , double , double , double , double , double , double , double );
void hyp_func (double *, double,double ,double ,double ,double ,double ,double ,double ,double ,double );
void hyp_func_der (double *,double ,double ,double ,double ,double ,double ,double ,double ,double ,double ,double ,double );

double abs1(double nos)
{
	if((nos >0)||(nos==0))
	return nos;
	else
	return -nos;
}

double min(double x,double y)
{
	if(x<y)
		return x;
	else if(y<x)
		return y;
	else
		return x;
}

double sign(double x)
{
	if (x>0)
	return 1;
	else if(x<0)
	return -1;
	else
	return 0;

}

double coth(double x)
{
	return (1/tanh(x));
}

double cot(double y)
{
	return (1/tan(y));
}

double asinh (double x)
{
	double r;
	r=log(x+sqrt(x*x+1));
	return r;
}


double B=1.1496658940635037000000000e+006, realmin=2.2250738585072014000e-308,  eps=2.2204460492503131000e-016;
double pi=3.14159265358979310e+000, realmax=1.797693134862315700000000000000e+308, q=1.602176462000000100e-019, ni=1.45e16, epsi=1.044794419999999900000000000000e-010, epox=3.453134099999999500000000000000e-011, beta=3.868172675200000300000000000000e+001, mu = 2.999999999999999900000000000000e-002, W=9.999999999999999500000000000000e-007, L =9.999999999999999500000000000000e-007, WI_LT_H=4.000000000000000200000000000000e-001, ERC=1e-2, ERW=1e-8, ER=9.999999999999999500000000000000e-007, Neps_H=100, Neps_T=100, Neps_P=10, Neps_TL=10, T_bound=1e110, H_bound=1e110, NR_thr=9.999999999999999500000000000000e-008, THETA_thr=270, pi_bound=1e-2, ERT, ERH,               /*note  ERH and ERT have no values assiged yet*/
accu_level=1,

max_itr_NR=100,
max_itr_rid=100,
max_itr_bis=500,
max_itr_t=1,

tsi=20e-9,
tox1=2.000000000000000100000000000000e-009,
tox2=1.000000000000000100000000000000e-009,




del2 =0.02,
step=0,

psi1_D=0,
psi2=0,
psi2_old=0,

psi2_D=0,
I_d=0,

cox1, cox2, csi, r1, r2,UT,

del1=0.1,
vgs1=4,
vgs2=12,
V=0,

min_limit=0,
VD;

int main()
{
	double  VD=2,z=0,  psi1_old=0, psi1, psi1_D_old=0, psi2_D_old=0, I_d_old=0, A=0, G=0, theta=0, colour=0, fid1=11, A_D=0, G_D=0, theta_D=0, colour1=0, fid2=12, Qifd=0, Qibd=0, Qifs=0, Qibs=0, I_f=0, I_r=0;
	FILE *fp1;
	fp1=fopen("final_c.txt","w");
	double cox1=epox/tox1, cox2=epox/tox2, csi=epsi/tsi, r1=csi/cox1, r2=csi/cox2,UT=1.0/beta,c_t=0;
	//printf("%18.20e %18.20e %18.20e %18.20e %18.20e %18.20e\n\n",cox1,cox2,csi,r1,r2,UT);
//double vgs1=20,vgs2=14;
	//double i;
	//for(i=1;i<10;i++)
	//{
		for(VD=0; VD <=20 ; VD=VD+1)
		{
			for(vgs1=0;vgs1<=20;vgs1=vgs1+1)
			{
				for(vgs2=0;vgs2<=20;vgs2=vgs2+1)
				{
					c_t++;
					
			//vgs1 = 1 + (int)( 3.0 * rand() / ( RAND_MAX + 1.0 ) );
			//vgs1 = 1 + (int)( 3.0 * rand() / ( RAND_MAX + 1.0 ) );
//vgs2=vgs1;
         if ((vgs1 <= 1.5)||(vgs2 <= 1.5))
         {//printf("1st enter\n");
             ERH=1e-13;
             ERT=1e-17;
         }
         else if  ((vgs1 <= 3)||(vgs2 <= 3))
         {//printf("2st enter\n");
             ERH=1e-25;
             ERT=1e-25;
         }
         else
         {//printf("3st enter\n");
             ERH=0;
             ERT=0;
         }
		//	printf("%e %e\n",ERH,ERT);
			step=step+1;
			z=0;            /* some stuff abt clear*/
			z++;
			B=(2*q*ni/(beta*epsi));
			psi1_old=psi1;
   		psi1_D_old=psi1_D;
    		psi2_old=psi2;
    		psi2_D_old=psi2_D;
    		I_d_old=I_d;
    		
    		main_part(&psi1,&psi2, &A, &G, &theta, &colour, vgs1, vgs2, V, tsi, tox1, tox2,  csi, cox1, cox2, fid1);
    		//printf("back\n");
    		main_part(&psi1_D, &psi2_D,  &A_D, &G_D, &theta_D, &colour1,vgs1, vgs2, VD, tsi, tox1, tox2, csi, cox1, cox2, fid2);
    		
    		Qifd= cox1*(vgs1 -psi1_D);
    		Qibd= cox2*(vgs2 -psi2_D);
    		
    		Qifs= cox1*( vgs1 - psi1);
    		Qibs= cox2*( vgs2 - psi2);
    		
    		I_f=(((Qifs*Qifs)/(2*cox1))+((Qibs*Qibs)/(2*cox2))+((2/beta)*(Qifs+Qibs))+(epsi*tsi/2)*G);
    		I_r=(((Qifd*Qifd)/(2*cox1))+((Qibd*Qibd)/(2*cox2))+((2/beta)*(Qifd+Qibd))+(epsi*tsi/2)*G_D);
  
    		I_d=mu*(W/L)*(I_f-I_r);
    	fprintf(fp1,"%18.20e %18.20e %18.20e %18.20e %18.20e\n",vgs1,vgs2,VD,psi1,psi2);
    		
    		printf("a %18.20e %18.20e %18.20e %18.20e %18.20e %18.20e %18.20e \n\n",psi1,psi2, A, G, theta, colour,vgs2);
    		//printf("b %18.20e %18.20e %18.20e %18.20e %18.20e %18.20e \n\n",psi1_D, psi2_D,  A_D, G_D, theta_D, colour1);
    		
    	}}}
 	fclose(fp1);
printf("\n%e\n",c_t);
 return 0;
 }
 
 
 
 void main_part(double*psi1,double* psi2,double*A,double*G,double* theta,double*colour,double vgs1,double vgs2,double V, double tsi,double tox1,double tox2,double csi,double cox1,double cox2,double fid)
{
	double a1=0,a2=0,a3=0,a4=0,a5=0,a6=0,a7=0,a8=0,a9=0, f=0,aq=0,d=0,w=0,u=0,y=0,v=0,r=0,t=0,h1=0,h2=0,h3=0,h4=0,h5=0,h6=0,vgs1crit=0, error_seen, err_trig=0, err_trig1=0, indx=0, x=0 ,x1=0 ,x0=0, exp1=0, exp2=0, method_overhead=0, method=0,  i=0, psigzp=0, psi1_r=0, psi2_r=0, psi_wi=0, psigzp1=0, psigzp2=0, psi_r=0, psi_limit=0, psi1_wi, psi2_wi, fac=0, time_overhead_old=0, time_overhead=0, iter_overhead=0, iteration=1, FLAG_err=0, vgscrit=0, vgs2crit=0, S1=0, S2=0, FLAG_pi=0 , kappa=0, fg0=0, fg1=0 , fg2=0;
	char *mode1,*mode;
	double w1=0,w2=0,w3=0,w4=0,w5=0,w6=0,w7=0,b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,b7=0,b8=0,b9=0,u1=0,u2=0,u3=0,u4=0,u5=0,u6=0,u7=0,u8=0,u9=0,u10=0,u11=0,t1=0,t2=0,t3=0,t4=0;
	double dr1=0,dr2=0,dr3=0,dr4=0,dr5=0,rt1=0,rt2=0,op=0,c1=0,c2=0;
	*psi1=0;*psi2=0;*A=0;*G=0;* theta=0;*colour=0;
	mode1=(char*)malloc(100*sizeof(char));
	mode=(char*)malloc(100*sizeof(char));
	clock_t tstart_overhead=0,tstop_overhead;
	tstart_overhead=clock();
	mode="S";
	method=0;

	for (i=1;i<=max_itr_t;i++)
	{//printf("e1\n");
		*A=0;
		*A=(B*exp(-beta*V));
		if (vgs2==vgs1)
		{
			psi1_wi=vgs1;
			psi2_wi=vgs2;
		}
		else
		{//printf("e2\n");
			//t1=vgs1*(cox1*csi+cox1*cox2)/(cox1*csi+cox1*cox2+cox2*csi);
			//t2=vgs2*(cox2*csi)/(cox1*csi+cox1*cox2+cox2*csi);
			//t3=vgs2*(cox2*csi+cox1*cox2)/(cox2*csi+cox1*cox2+cox1*csi);
			//t4=vgs1*(cox1*csi)/(cox2*csi+cox1*cox2+cox1*csi);
			psi1_wi=vgs1*(cox1*csi+cox1*cox2)/(cox1*csi+cox1*cox2+cox2*csi)+vgs2*(cox2*csi)/(cox1*csi+cox1*cox2+cox2*csi);
			psi2_wi=vgs2*(cox2*csi+cox1*cox2)/(cox2*csi+cox1*cox2+cox1*csi)+vgs1*(cox1*csi)/(cox2*csi+cox1*cox2+cox1*csi);
		}

		if(vgs2 <= vgs1)
		{
			my_psigzp(&psigzp1,vgs1, V,cox1);

			if(psigzp1 > V)
				vgs2crit=((-2/beta)*(-beta*V*1/2+log((beta*sqrt(B)*tsi*0.5)+exp(-beta*(psigzp1-V)*0.5)))-((epsi)/(cox2*((beta*tsi*0.5)+(1/sqrt(B))*exp(-beta*(psigzp1-V)*0.5)))));
			else 
				vgs2crit=((-2/beta)*(-beta*V*1/2+log(exp(beta*(psigzp1-V)*0.5)*(beta*sqrt(B)*tsi*0.5)+1)-(beta*(psigzp1-V)*0.5))-((epsi)*exp(beta*(psigzp1-V)*0.5)/(cox2*((beta*tsi*0.5)*exp(beta*(psigzp1-V)*0.5)+(1/sqrt(B))))));

			psigzp=psigzp1;
			vgscrit=vgs2crit;
			

			if (abs1(vgs1-psigzp1)< ERW)
			{//printf("entry1\n");
				*psi1=psi1_wi;
				*psi2=psi2_wi;
				
				if(vgs2crit <= vgs2)
				{	
					*colour=1;
					*G=-(((cox1/epsi)*(vgs1-*psi1))*((cox1/epsi)*(vgs1-*psi1))-B*exp(beta*(*psi1-V)));
					if (*G <= 0)
						*G=realmin;

					fac=sqrt(B*exp(beta*(*psi1-V)));

					if(*G==realmin)
						*theta=0;
					else
						*theta=(beta*tsi*0.5*sqrt(*G))+asin(sqrt(*G)/fac);
				}

				else
				{	
					*colour=2;
					*G=(((cox1/epsi)*(vgs1-*psi1))*((cox1/epsi)*(vgs1-*psi1))-B*exp(beta*(*psi1-V)));

					if (*G<=0)
						*G=realmin;

					fac=sqrt(B*exp(beta*(*psi1-V)));

					if (fac <=1e15*realmin)
						*theta=(beta*tsi*0.5*sqrt(*G))+asinh(realmax);
					else 
						*theta=(beta*tsi*0.5*sqrt(*G))+asinh(sqrt(*G)/fac);
				}
				//h1=((cox1/epsi)*(vgs1-*psi1))*((cox1/epsi)*(vgs1-*psi1));
				//h2=B*exp(beta*(*psi1-V));
				//*G=h1-h2;
				*G=(((cox1/epsi)*(vgs1-*psi1))*((cox1/epsi)*(vgs1-*psi1))-B*exp(beta*(*psi1-V)));
			
			}

			else if (abs1(vgs2-vgs2crit) <ER)
			{
				x= psigzp1;
				
				exp1=sqrt(B*exp(beta*(x-V)));
				exp2=beta*exp1*tsi;
				S1=((12*cox2*epsi*epsi*exp1*exp1*(exp2 + 2)*(exp2 + 2)*(exp2 + 2))/(exp1*exp1*(16*cox2*(epsi*epsi*(9*exp2 + 6) + beta*cox1*cox1*tsi*tsi*(vgs1 - x)*(5*exp2 + 9)) + 48*beta*cox1*cox1*epsi*tsi*(vgs1 - x)*(2*exp2 + 3)) + 96*cox1*cox1*exp1*(vgs1 - x)*(epsi + cox2*tsi) + 2*beta*exp1*exp1*exp1*exp2*(epsi*epsi*epsi*(24*exp2 + 24) + 4*cox2*epsi*epsi*tsi*(5*exp2 + 12) + beta*cox1*cox1*cox2*tsi*tsi*tsi*(vgs1 - x)*(exp2 + 10) + 2*beta*cox1*cox1*epsi*tsi*tsi*(vgs1 - x)*(exp2 + 8)) + beta*epsi*epsi*exp1*exp1*exp1*exp2*exp2*exp2*(16*epsi + 10*cox2*tsi + 2*exp2*epsi + cox2*exp2*tsi)));
				psi1_r= psigzp1+S1*(vgs2-vgs2crit);
				*psi1=psi1_r;
				
				if(*psi1==psi1_wi)
					*psi2=psi2_wi;
				else if (*psi1==psigzp1)
				{
					psigzp2=V-(2/beta)*log(0.5*sqrt(B)*tsi*beta+exp(0.5*beta*(V-psigzp1)));
					*psi2=psigzp2;
				}
				else 
				{
					if(vgs2crit <= vgs2)
					{
						*colour=1;
						*G=-(((cox1/epsi)*(vgs1-*psi1))*((cox1/epsi)*(vgs1-*psi1))-B*exp(beta*(*psi1-V)));

						if(*G <= 0)
							*G=realmin;

						fac=sqrt(B*exp(beta*(*psi1-V)));

						if(*G==realmin)
							*theta=0;
						else
							*theta=(beta*tsi*0.5*sqrt(*G))+asin(sqrt(*G)/fac);
					}
					else
					{
						*colour=2;
						*G=(((cox1/epsi)*(vgs1-*psi1))*((cox1/epsi)*(vgs1-*psi1))-B*exp(beta*(*psi1-V)));

						if(*G <= 0)
							*G=realmin;

						fac=sqrt(B*exp(beta*(*psi1-V)));

						if (fac <=1e15*realmin)
							*theta=(beta*tsi*0.5*sqrt(*G))+asinh(realmax);
						else
							*theta=(beta*tsi*0.5*sqrt(*G))+asinh(sqrt(*G)/fac);


						*psi2=V+(-2/beta)*log(sqrt(B/(*G))*sinh(*theta));
					}
				}
				
				*G=(((cox1/epsi)*(vgs1-*psi1))*((cox1/epsi)*(vgs1-*psi1))-B*exp(beta*(*psi1-V)));

				if(vgs2crit <= vgs2)
					*colour=1;
				else
					*colour=2;

			}
			else
			{ 
				
				if(vgs2crit<= vgs2)
				{//printf("entry3b\n");
					*colour=1;
					x= psigzp1;
				//	printf("x==%18.18e\n\n",x);
					exp1=sqrt(B*exp(beta*(x-V)));
					exp2=beta*exp1*tsi;
				
		
			
	f=9*exp2;
	aq=epsi*epsi;
	a9=(vgs1 - x);
	b1=(5*exp2 + 9);
	b2=(vgs1 - x);
	b3=(2*exp2 + 3);
	b4=(vgs1 - x);
	b5=(epsi + cox2*tsi);
	b6=(vgs1 - x);
	b7=(exp2 + 8);
	b8=(exp2 + 10);
	b9=(5*exp2 + 12);
	w1=(24*exp2 + 24);
	u=beta*cox1*cox1*tsi*tsi*a9*b1;
	y=48*beta*cox1*cox1*epsi*tsi*b2*b3;
	v=96*cox1*cox1*exp1*b4*b5;
	a1=2*beta*cox1*cox1*epsi*tsi*tsi*b6*b7;
	a2=beta*cox1*cox1*cox2*tsi*tsi*tsi*b6*b8 ;
	a3=4*cox2*epsi*epsi*tsi*b9 ;
	a4=epsi*epsi*epsi*w1;
	a5=16*epsi;
	a6=10*cox2*tsi;
	a7=2*exp2*epsi;
	a8=cox2*exp2*tsi;
	
	r= beta*epsi*epsi*exp1*exp1*exp1*exp2*exp2*exp2*( a5+ a6 + a7 + a8);
	t= 2*beta*exp1*exp1*exp1*exp2*(a4 + a3+ a2+ a1);
	w2=v + t +r;
	w3=(f + 6);
	w4=(aq*w3 + u);
	w5=(16*cox2*w4 + y);
	w=(exp1*exp1*w5 + w2);
	d=12*cox2*epsi*epsi*exp1*exp1*(exp2 + 2)*(exp2 + 2)*(exp2 + 2);
	S1=(d/w);

					S1=((12*cox2*epsi*epsi*exp1*exp1*(exp2 + 2)*(exp2 + 2)*(exp2 + 2))/(exp1*exp1*(16*cox2*(epsi*epsi*(9*exp2 + 6) + beta*cox1*cox1*tsi*tsi*(vgs1 - x)*(5*exp2 + 9)) + 48*beta*cox1*cox1*epsi*tsi*(vgs1 - x)*(2*exp2 + 3)) + 96*cox1*cox1*exp1*(vgs1 - x)*(epsi + cox2*tsi) + 2*beta*exp1*exp1*exp1*exp2*(epsi*epsi*epsi*(24*exp2 + 24) + 4*cox2*epsi*epsi*tsi*(5*exp2 + 12) + beta*cox1*cox1*cox2*tsi*tsi*tsi*(vgs1 - x)*(exp2 + 10) + 2*beta*cox1*cox1*epsi*tsi*tsi*(vgs1 - x)*(exp2 + 8)) + beta*epsi*epsi*exp1*exp1*exp1*exp2*exp2*exp2*(16*epsi + 10*cox2*tsi + 2*exp2*epsi + cox2*exp2*tsi)));
		//	printf("S1==%18.20e\n",S1);
					psi1_r= psigzp1+S1*(vgs2-vgs2crit);

					root_find_trig(psi1, &psi_limit, &method, &method_overhead, &iteration, &iter_overhead,  &time_overhead, &FLAG_err,vgs1, vgs2, V, tsi, cox1, cox2, psigzp1, psi1_r, fid);
//printf("%18.16Le intrufgy\n\n",*psi1);
					//h5=((cox1/epsi)*(vgs1-*psi1))*((cox1/epsi)*(vgs1-*psi1));
					//h6=B*exp(beta*(*psi1-V));
					//*G=0;
					//*G=-(h5-h6);
					*G=-(((cox1/epsi)*(vgs1-*psi1))*((cox1/epsi)*(vgs1-*psi1))-B*exp(beta*(*psi1-V)));
//printf("G==%18.20e %18.20e\n\n",*G,*psi1-V);
					if(*G <= 0)
						*G=realmin;

					fac=sqrt(B*exp(beta*(*psi1-V)));
//printf("faccchk==%18.16Le *psi1==%18.16Le G==%18.16Le\n\n",fac,*psi1,*G);
					if (*G==realmin)
						*theta=0;
					else
						*theta=(beta*tsi*0.5*sqrt(*G))+asin(sqrt(*G)/fac);
//						printf("theta1==%18.16Le G1==%18.16Le sqrt==%18.16Le fac==%18.16Le\n\n",(beta*tsi*0.5*sqrt(*G)),asin(sqrt(*G)/fac),sqrt(*G),fac);   /*point of load*/
//printf("theta==%18.16Le G==%18.16Le\n\n",*theta,*G);
					FLAG_pi=1;

					if (abs1(*theta-pi) < pi_bound)
						FLAG_pi=1;

					if(*psi1==psi_wi)
						*psi2=psi2_wi;
					else if (*psi1==psigzp1)
					{
						 psigzp2=V-(2/beta)*log(0.5*sqrt(B)*tsi*beta+exp(0.5*beta*(V-psigzp1)));
                   *psi2=psigzp2;
					}
					else 
					{//printf("entry3b-a\n");
						if(FLAG_pi==0)
							*psi2=V+(-2/beta)*log(sqrt(B/(*G))*sin(*theta));
						else
						{
							kappa=(epsi/cox2)*sqrt(B);
							if (vgs2 > 0)
								x0=min(vgs2, ((2 / beta) * (log ((vgs2 / kappa)) + beta * V * 0.5)));
							else
								x0 = vgs2;

							for (indx=1;indx<= max_itr_NR; indx++)
							{
								exp1=B*exp(beta*(x0-V));
								fg0=*G+(cox2/epsi)*(cox2/epsi)*((vgs2-x0)*(vgs2-x0))-exp1;
								fg1=-2*(cox2/epsi)*(cox2/epsi)*(vgs2-x0)-exp1*beta;
								fg2=2*(cox2/epsi)*(cox2/epsi)-exp1*beta*beta;
                            
                           		x1=x0 - 2*fg0*fg1/(2*fg1*fg1-fg0*fg2);
                            
								if (( abs1(x0 - x1)<= 1e-10))
									break;
								else
									x0=x1;
							}

							*psi2=x1;
						}
					}
					//h3=((cox1/epsi)*(vgs1-*psi1))*((cox1/epsi)*(vgs1-*psi1));
					//h4=B*exp(beta*(*psi1-V));
					//*G=(h3-h4);
					*G=(((cox1/epsi)*(vgs1-*psi1))*((cox1/epsi)*(vgs1-*psi1))-B*exp(beta*(*psi1-V)));

					trig_func(&err_trig,vgs1, vgs2, tsi, cox1, cox2, V, *psi1);

					error_seen=abs1(err_trig);
				}
				else 
				{ 
					*colour=2;
					root_find_hyp(psi1, &psi_limit, &method, &method_overhead, &iteration, &iter_overhead,  &time_overhead, &FLAG_err, vgs1, vgs2, V,  tsi, csi, cox1, cox2, psi1_wi, psigzp1, vgs2crit, *A, fid);
					//*G=0;
					*G=(((cox1/epsi)*(vgs1-*psi1))*((cox1/epsi)*(vgs1-*psi1))-B*exp(beta*(*psi1-V)));

					if(*G <= 0)
						*G=realmin;

					fac=sqrt(B*exp(beta*(*psi1-V)));

					if (fac <=1e15*realmin)
						*theta=(beta*tsi*.5*sqrt(*G))+asinh(realmax);
					else
						*theta=(beta*tsi*.5*sqrt(*G))+asinh(sqrt(*G)/sqrt(B*exp(beta*(*psi1-V))));

					if (*psi1==psi1_wi)
						*psi2=psi2_wi;
					else if (*psi1==psigzp1)
					{
						 psigzp2=V-(2/beta)*log(0.5*sqrt(B)*tsi*beta+exp(0.5*beta*(V-psigzp1)));
                         *psi2=psigzp2;
					}
					else 
					{
						if(*theta > 700)
							*psi2=V+(-2/beta)*(log(sqrt(B/(4*(*G))))+(*theta));
						else 
							*psi2=V+(-2/beta)*log(sqrt(B/(*G))*sinh(*theta));
					}

				}
			}

		error_seen=1;
		}
		else
		{//printf("e3\n");
			psigzp2=0;
			my_psigzp(&psigzp2, vgs2, V, cox2);

			if (psigzp > V)
				vgs1crit=((-2/beta)*(-beta*V*1/2+log((beta*sqrt(B)*tsi*0.5)+exp(-beta*(psigzp2-V)*0.5)))-((epsi)/(cox1*((beta*tsi*0.5)+(1/sqrt(B))*exp(-beta*(psigzp2-V)*0.5)))));
			else
            vgs1crit=((-2/beta)*(-beta*V*1/2+log(exp(beta*(psigzp2-V)*0.5)*(beta*sqrt(B)*tsi*0.5)+1)-(beta*(psigzp2-V)*0.5))-((epsi)*exp(beta*(psigzp2-V)*0.5)/(cox1*((beta*tsi*0.5)*exp(beta*(psigzp2-V)*0.5)+(1/sqrt(B))))));
			psigzp=0;vgscrit=0;
			psigzp=psigzp2;
			vgscrit=vgs1crit;
//printf("vgscrit==%18.18e  psigzp==%18.18e\n",vgscrit,psigzp);
			if(abs1(vgs2-psigzp2)<ERW)
			{
				*psi2=psi2_wi;
				*psi1=psi1_wi;

				if (vgs1crit <= vgs1)
				{
					*colour=1;
					u1=((cox2/epsi)*(vgs2-*psi2))*((cox2/epsi)*(vgs2-*psi2));
					u2=B*exp(beta*(*psi2-V));
					*G=0;
					*G=-(u1-u2);
					//*G=-(((cox2/epsi)*(vgs2-*psi2))*((cox2/epsi)*(vgs2-*psi2))-B*exp(beta*(*psi2-V)));
					if (*G<=0)
						*G=realmin;

					fac=sqrt(B*exp(beta*(*psi2-V)));

					if(*G==realmin)
						*theta=0;
					else
						*theta=(beta*tsi*0.5*sqrt(*G))+asin(sqrt(*G)/fac);
				}

				else 
				{
					*colour=2;
					*G=(((cox2/epsi)*(vgs2-*psi2))*((cox2/epsi)*(vgs2-*psi2))-B*exp(beta*(*psi2-V)));

					if(*G <= 0)
						*G=realmin;

					fac=sqrt(B*exp(beta*(*psi2-V)));

					if (fac <=1e15*realmin)
						*theta=(beta*tsi*0.5*sqrt(*G))+asinh(realmax);
					else
						*theta=(beta*tsi*0.5*sqrt(*G))+asinh(sqrt(*G)/fac);
				}
			}
			else if ((abs1(vgs1-vgs1crit) <ER))
			{
				x= psigzp2;
				exp1=sqrt(B*exp(beta*(x-V)));
				exp2=beta*exp1*tsi;
				S2=((12*cox1*epsi*epsi*exp1*exp1*(exp2 + 2)*(exp2 + 2)*(exp2 + 2))/(exp1*exp1*(16*cox1*(epsi*epsi*(9*exp2 + 6) + beta*cox2*cox2*tsi*tsi*(vgs2 - x)*(5*exp2 + 9)) + 48*beta*cox2*cox2*epsi*tsi*(vgs2 - x)*(2*exp2 + 3)) + 96*cox2*cox2*exp1*(vgs2 - x)*(epsi + cox1*tsi) + 2*beta*exp1*exp1*exp1*exp2*(epsi*epsi*epsi*(24*exp2 + 24) + 4*cox1*epsi*epsi*tsi*(5*exp2 + 12) + beta*cox2*cox2*cox1*tsi*tsi*tsi*(vgs2 - x)*(exp2 + 10) + 2*beta*cox2*cox2*epsi*tsi*tsi*(vgs2 - x)*(exp2 + 8)) + beta*epsi*epsi*exp1*exp1*exp1*exp2*exp2*exp2*(16*epsi + 10*cox1*tsi + 2*exp2*epsi + cox1*exp2*tsi)));
				psi2_r= psigzp2+S2*(vgs1-vgs1crit);
				*psi2=psi2_r;

				if (*psi2==psi2_wi)
					*psi1=psi1_wi;
				else if (*psi2==psigzp2)
				{
					 psigzp1=V-(2/beta)*log(0.5*sqrt(B)*tsi*beta+exp(0.5*beta*(V-psigzp2)));
                     *psi1=psigzp1;
				}
				else 
				{
					if(vgs1crit <= vgs1)
					{
						*G=-(((cox2/epsi)*(vgs2-*psi2))*((cox2/epsi)*(vgs2-*psi2))-B*exp(beta*(*psi2-V)));

						if(*G <= 0)
							*G=realmin;

						fac=sqrt(B*exp(beta*(*psi2-V)));

						if (*G==realmin)
							*theta=0;
						else
							*theta=(beta*tsi*0.5*sqrt(*G))+asin(sqrt(*G)/fac);

						*psi1=V+(-2/beta)*log(sqrt(B/(*G))*sin(*theta));
					}
					else
					{
						*G=(((cox2/epsi)*(vgs2-*psi2))*((cox2/epsi)*(vgs2-*psi2))-B*exp(beta*(*psi2-V)));

						if(*G <= 0)
							*G=realmin;

						fac=sqrt(B*exp(beta*(*psi2-V)));

						if (fac <=1e15*realmin)
							*theta=(beta*tsi*0.5*sqrt(*G))+asinh(realmax);
						else
							*theta=(beta*tsi*0.5*sqrt(*G))+asinh(sqrt(*G)/fac);

						*psi1=V+(-2/beta)*log(sqrt(B/(*G))*sinh(*theta));
					}
				}
				*G=0;
				*G=(((cox2/epsi)*(vgs2-*psi2))*((cox2/epsi)*(vgs2-*psi2))-B*exp(beta*(*psi2-V)));

				if(vgs1crit <= vgs1)
					*colour=1;
				else
					*colour=2;
			}
			
			else 
			{//printf("e4\n");
				if( vgs1crit <= vgs1 )
				{
					*colour=1;
					x= psigzp2;
					exp1=sqrt(B*exp(beta*(x-V)));
					exp2=beta*exp1*tsi;
					S2=((12*cox1*epsi*epsi*exp1*exp1*(exp2 + 2)*(exp2 + 2)*(exp2 + 2))/(exp1*exp1*(16*cox1*(epsi*epsi*(9*exp2 + 6) + beta*cox2*cox2*tsi*tsi*(vgs2 - x)*(5*exp2 + 9)) + 48*beta*cox2*cox2*epsi*tsi*(vgs2 - x)*(2*exp2 + 3)) + 96*cox2*cox2*exp1*(vgs2 - x)*(epsi + cox1*tsi) + 2*beta*exp1*exp1*exp1*exp2*(epsi*epsi*epsi*(24*exp2 + 24) + 4*cox1*epsi*epsi*tsi*(5*exp2 + 12) + beta*cox2*cox2*cox1*tsi*tsi*tsi*(vgs2 - x)*(exp2 + 10) + 2*beta*cox2*cox2*epsi*tsi*tsi*(vgs2 - x)*(exp2 + 8)) + beta*epsi*epsi*exp1*exp1*exp1*exp2*exp2*exp2*(16*epsi + 10*cox1*tsi + 
2*exp2*epsi + cox1*exp2*tsi)));
					
					psi2_r= psigzp2+S2*(vgs1-vgs1crit);
               // printf("S2==%18.18e %18.20e \n\n",S2, psi2_r);
              	root_find_trig(psi2, &psi_limit, &method, &method_overhead, &iteration, &iter_overhead,  &time_overhead, &FLAG_err,vgs2, vgs1, V, tsi, cox2, cox1, psigzp2, psi2_r, fid);
              //*psi2=7.7367898968720006e-001;
               dr3=(vgs2-*psi2);dr4=(*psi2-V);op=((cox2/epsi)*dr3);
              	dr1=(op*op);
              	dr2=B*exp(beta*dr4);
               //printf("%18.20e %18.20e psi2==%18.20e\n",(cox2/epsi),dr3,*psi2);
              *G=-(dr1-dr2);
					 *G=-(((cox2/epsi)*(vgs2-*psi2))*((cox2/epsi)*(vgs2-*psi2))-B*exp(beta*(*psi2-V)));
               // printf("G_test==%18.20e %18.20e %18.20e\n",*G,((cox2/epsi)*(vgs2-*psi2))*((cox2/epsi)*(vgs2-*psi2)),B*exp(beta*(*psi2-V)));
					if (*G<=0)
						*G=realmin;

					fac=sqrt(B*exp(beta*(*psi2-V)));
//printf("fac==%18.20e\n",fac);
					if (*G==realmin)
						*theta=0;
					else
						{
							c1=asin(sqrt(*G)/fac);
							c2=(beta*tsi*0.5*sqrt(*G));
							*theta=c2+c1;
						}		
					FLAG_pi=0;
				//printf("theta==%18.20e\n",*theta);
					if(abs1(*theta-pi) <= pi_bound)
						FLAG_pi=1;

					if (*psi2==psi2_wi)
						*psi1=psi1_wi;
						
					
					else if (*psi2==psigzp2)
					{
						psigzp1=V-(2/beta)*log(0.5*sqrt(B)*tsi*beta+exp(0.5*beta*(V-psigzp2)));
                        *psi1=psigzp1;
					}
					else
					{
						if(FLAG_pi==0)
							*psi1=V+(-2/beta)*log(sqrt(B/(*G))*sin(*theta));
						else
						{
							kappa=(epsi/cox1)*sqrt(B);
							if (vgs1 > 0)
								x0=min(vgs1, ((2 / beta) * (log ((vgs1 / kappa)) + beta * V * 0.5)));
							else
								x0 = vgs1;

							for (indx=1; indx <= max_itr_NR; indx++)
							{
								exp1=B*exp(beta*(x0-V));
								fg0=*G+(cox1/epsi)*(cox1/epsi)*((vgs1-x0)*(vgs1-x0))-exp1;
								fg1=-2*(cox1/epsi)*(cox1/epsi)*(vgs1-x0)-exp1*beta;
								fg2=2*(cox1/epsi)*(cox1/epsi)-exp1*beta*beta;
                            
                                x1=x0 - 2*fg0*fg1/(2*fg1*fg1-fg0*fg2);
                            
								if (( abs1(x0 - x1)<= 1e-10))
									break;
								else
									x0=x1;
							}
							*psi1=x1;
						}
					}
				//	rt1=((cox2/epsi)*(vgs2-*psi2))*((cox2/epsi)*(vgs2-*psi2));
			//		rt2=B*exp(beta*(*psi2-V));
				//	*G=0;
				//	*G=(rt1-rt2);
					*G=(((cox2/epsi)*(vgs2-*psi2))*((cox2/epsi)*(vgs2-*psi2))-B*exp(beta*(*psi2-V)));

					trig_func(&err_trig1,vgs2, vgs1, tsi, cox2, cox1, V, *psi2);
					error_seen=abs1(err_trig1);
//printf("psi_temp==%18.20e\n",*psi1);
				}

				else
				{//printf("e5\n");
					*colour=2;
					root_find_hyp(psi2, &psi_limit, &method, &method_overhead, &iteration, &iter_overhead, &time_overhead, &FLAG_err, vgs2, vgs1, V,  tsi, csi, cox2, cox1, psi2_wi, psigzp2, vgs1crit, *A, fid);
					//u9=(vgs2-*psi2);u10=(vgs2-*psi2);u11=(*psi2-V);
					//u7=((cox2/epsi)*u9)*((cox2/epsi)*u10);
					//u8=B*exp(beta*u11);
					//*G=0;
					//*G=(u7-u8);
					*G=(((cox2/epsi)*(vgs2-*psi2))*((cox2/epsi)*(vgs2-*psi2))-B*exp(beta*(*psi2-V)));
                //printf("u %18.20e %18.20e\n\n ",u11,B*exp(beta*u11));
					if (*G<=0)
						*G=realmin;
//printf("g==%18.18e psi2==%18.18e\n",*G,*psi2);
					fac=sqrt(B*exp(beta*(*psi2-V)));

					if (fac <=1e15*realmin)
						*theta=(beta*tsi*0.5*sqrt(*G))+asinh(realmax);
					else
						*theta=(beta*tsi*0.5*sqrt(*G))+asinh(sqrt(*G)/fac);
//printf("theta==%18.18e\n",*theta);
					if (*psi2==psi2_wi)
						*psi1=psi1_wi;
					else if (*psi2==psigzp2)
					{
						psigzp1=V-(2/beta)*log(0.5*sqrt(B)*tsi*beta+exp(0.5*beta*(V-psigzp2)));
                  *psi1=psigzp1;
					}
					else
					{
						if(*theta > 700)
							*psi1=V+(-2/beta)*(log(sqrt(B/(4*(*G))))+(*theta));
						else 
							*psi1=V+(-2/beta)*log(sqrt(B/(*G))*sinh(*theta));
					}
				}
			}
		error_seen=1;
		}
	//time_overhead=time_overhead_old+time_overhead;  /*r dey imp*/
    //time_overhead_old=time_overhead;
	}
	tstop_overhead=clock();
	time_overhead = ((double)( tstop_overhead - tstart_overhead )) / CLOCKS_PER_SEC;
	/* some lines related to time deleted*/

	if (*colour==1)
		mode1="trig";
	else
		mode1="hyp";
	
printf("%s\n",mode1);
	return;
}


void root_find_hyp(double *psi,double *psi_limit,double *method, double *method_overhead, double *iteration, double *iter_overhead, double *time_overhead, double *FLAG_err,double vgs1, double vgs2, double V,double tsi, double csi, double cox1, double cox2, double psi_wi, double psigzp, double vgscrit, double A, double fid)
{
	double FLAG_WI=0, flag_approx=0, FLAG, FLAG_theta, psi_approx=0, x0, x_max, xnew, x_min, xm, x1, j, iter_rid, G_0, G_min, G_max, G_new, G_m, exp1, exp1_0, exp4_0, exp1_min, exp1_new, exp4_min, exp4_new, theta_new, exp1_m, exp1_max, exp4_max, exp4_m, theta_m, theta_min, theta_max, theta_0, fg0, fg1, f_min, f_max, fm, fnew, A_0, A_min, A_max, A_m ,A_new , s;
	
	*psi_limit=0;
	*method_overhead=0;
	*iter_overhead=0;
	*FLAG_err=0;
	clock_t tstart_overhead=0,tstop_overhead;
	tstart_overhead=clock();
	FLAG_theta=0;
	
	if((abs1(vgs1-V)<=WI_LT_H)||((vgscrit-vgs2)<ERC))
	{
		FLAG_WI=1;
		psi_approx=min(psigzp-(1e-7),psi_wi);
	}
	else 
	{
		hyp_limit(psi_limit, method_overhead,iter_overhead,vgs1, vgs2, V,  csi, cox1, cox2, psi_wi, psigzp, flag_approx);
		psi_approx=*psi_limit;
	}
	
	tstop_overhead=clock();
	*time_overhead = ((double)( tstop_overhead - tstart_overhead )) / CLOCKS_PER_SEC;
	
	x0=psi_approx;
	FLAG=0;
	j=1;
	
	x_max=min(psi_wi+10*eps, psigzp-10*eps);
	x_min=x0;
	x1=x_max;
	
	iter_rid=0;
	
	for(j=1;j<=max_itr_NR;j++)
	{
		if(FLAG==0)
		{
			*method=1;
			G_0=((((cox1/epsi)*(vgs1-x0))*((cox1/epsi)*(vgs1-x0)))-(B*exp(beta*(x0-V))));
			
			if (G_0<=0)
			{
				G_0=realmin;
			}	
			A_0=(B*exp(beta*(x0-V)));
			if (A_0<=0)
				A_0=realmin;

			exp1_0=sqrt(G_0);
			exp4_0=(exp1_0/sqrt(A_0));
			theta_0=(beta*tsi*.5*exp1_0)+asinh(exp4_0);

			if ((theta_0 > THETA_thr))
			{
				FLAG_theta=1;
				FLAG=1;
			}

			if(FLAG_theta==0)
			{
				hyp_func(&fg0,vgs1, vgs2, tsi, cox1, cox2, V, x0, FLAG_theta, exp1_0, theta_0);
				hyp_func_der(&fg1,vgs1, vgs2, tsi, cox1, cox2, V, x0, G_0, exp1_0, A_0, exp4_0, theta_0 );

				if(abs1(fg0)<=ERH)
				{
					x1=x0;
					break;
				}
				else 
				{
					exp1=fg0/fg1;

					if ((j==1) && (abs1(exp1) < NR_thr))
						FLAG=1;

					x1=x0-exp1;
				}
			}

			if (((x1 > x_max)||(j > 9))||(FLAG_theta==1))
			{
				tstart_overhead=clock();
				if(FLAG_WI==1)
				{
					hyp_limit(psi_limit, method_overhead,iter_overhead ,vgs1, vgs2, V, csi, cox1, cox2, psi_wi, psigzp, flag_approx);
				}
				tstop_overhead=clock();
				*time_overhead =*time_overhead+((double)( tstop_overhead - tstart_overhead )) / CLOCKS_PER_SEC;

				FLAG=1;

				if(j>2)
				{
					x_min=*psi_limit;

					x_max=min(psi_wi+10*eps, psigzp-10*eps);
					
					G_min=((((cox1/epsi)*(vgs1-x_min))*((cox1/epsi)*(vgs1-x_min)))-(B*exp(beta*(x_min-V))));
					 if (G_min<=0)
						G_min=realmin;

					 A_min=(B*exp(beta*(x_min-V)));

					 if (A_min<=0)
						A_min=realmin;

					 exp1_min=sqrt(G_min);
					 exp4_min=(exp1_min/sqrt(A_min));
					 theta_min=(beta*tsi*.5*exp1_min)+asinh(exp4_min);
                
					 G_max=((((cox1/epsi)*(vgs1-x_max))*((cox1/epsi)*(vgs1-x_max)))-(B*exp(beta*(x_max-V))));

					 if (G_max<=0)
						G_max=realmin;

					 A_max=(B*exp(beta*(x_max-V)));
					 if (A_max<=0)
						 A_max=realmin;

					 exp1_max=sqrt(G_max);
					 exp4_max=(exp1_max/sqrt(A_max));
					 theta_max=(beta*tsi*.5*exp1_max)+asinh(exp4_max);

					 if ((theta_min>THETA_thr)||(theta_max>THETA_thr))
						FLAG_theta=1;

					 hyp_func(&f_min,vgs1, vgs2, tsi, cox1, cox2, V, x_min, FLAG_theta, exp1_min, theta_min);
					 hyp_func(&f_max,vgs1, vgs2, tsi, cox1, cox2, V, x_max, FLAG_theta, exp1_max, theta_max);

					 if(x1 < x_max)
					 {
						 fg0=fg1;
						 x0=x1;
					 }

					  if (FLAG_theta==1)
						 hyp_func(&fg0,vgs1, vgs2, tsi, cox1, cox2, V, x0, FLAG_theta, exp1_0, theta_0);

					  if (fg0*f_min > 0)
					  {
						  x_min=x0;
						  f_min=fg0;
					  }
					  else
					  {
						   x_max=x0;
						   f_max=fg0;
					  }
				}
			}

			if(FLAG==0)
			{
				if((( abs1(x0 - x1)<=(Neps_H+1)*eps)))
					break;
				else
					x0=x1;
			}
		}

		else
		{
			*method=2;
			tstart_overhead=clock();

			if((j<=2) &&(iter_rid==0))
			{
				hyp_limit(psi_limit, method_overhead,iter_overhead,vgs1, vgs2, V, csi, cox1, cox2, psi_wi,  psigzp, flag_approx);
				x_min=*psi_limit;
				x_max=min(psi_wi+10*eps, psigzp-10*eps);

				G_min=((((cox1/epsi)*(vgs1-x_min))*((cox1/epsi)*(vgs1-x_min)))-(B*exp(beta*(x_min-V))));
				if (G_min<=0)
					G_min=realmin;

				A_min=(B*exp(beta*(x_min-V)));
				if (A_min<=0)
					A_min=realmin;

				exp1_min=sqrt(G_min);
				exp4_min=(exp1_min/sqrt(A_min));
				theta_min=(beta*tsi*.5*exp1_min)+asinh(exp4_min);
            
				G_max=((((cox1/epsi)*(vgs1-x_max))*((cox1/epsi)*(vgs1-x_max)))-(B*exp(beta*(x_max-V))));
				if (G_max<=0)
					G_max=realmin;

				A_max=(B*exp(beta*(x_max-V)));
				if (A_max<=0)
					A_max=realmin;

				exp1_max=sqrt(G_max);
				exp4_max=(exp1_max/sqrt(A_max));
				theta_max=(beta*tsi*.5*exp1_max)+asinh(exp4_max);
            
            
				if ((theta_min > THETA_thr)||(theta_max>THETA_thr))
					FLAG_theta=1;

				hyp_func(&f_min,vgs1, vgs2, tsi, cox1, cox2, V, x_min, FLAG_theta, exp1_min, theta_min);
				hyp_func(&f_max,vgs1, vgs2, tsi, cox1, cox2, V, x_max, FLAG_theta, exp1_max, theta_max);
			}

			tstop_overhead=clock();
			*time_overhead =((double)( tstop_overhead - tstart_overhead )) / CLOCKS_PER_SEC;

			iter_rid=iter_rid+1;
			xm=0.5*(x_min+x_max);

			if ((abs1(x_min-x_max)<=(Neps_H+1)*eps)||(abs1(x_min-xm)<=(Neps_H)*eps)||(abs1(x_max-xm)<=(Neps_H)*eps))
			{
				x1=xm;
				break;
			}
			if (accu_level==1)
			{
				if (abs1(f_min) <= ERH)
				{
					x1=x_min;
		            break;
				}
				else if (abs1(f_max) <= ERH)
				{
					x1=x_max;
					break;
				}
			}

			if ((sign(f_max)*sign(f_min)<0) && ((abs1(f_min)>H_bound*realmin)&&abs1(f_max)>H_bound*realmin))
			{
				G_m=((((cox1/epsi)*(vgs1-xm))*((cox1/epsi)*(vgs1-xm)))-(B*exp(beta*(xm-V))));
				if (G_m<=0)
					G_m=realmin;

				A_m=(B*exp(beta*(xm-V)));
				if (A_m<=0)
					A_m=realmin;

				exp1_m=sqrt(G_m);
				exp4_m=(exp1_m/sqrt(A_m));
				theta_m=(beta*tsi*.5*exp1_m)+asinh(exp4_m);

				hyp_func(&fm,vgs1, vgs2, tsi, cox1, cox2, V, xm, FLAG_theta, exp1_m, theta_m);
				if (abs1(fm) <= ERH)
				{
					x1=xm;
					break;
				}

				s=abs1(fm*fm-f_min*f_max);

				if ((s==0)||(1-abs1(fm/sqrt(s))) < eps )
				{
					*method=2;
					if (sign(f_min)*sign(fm) > 0)
					{
						x_min=xm;
						f_min=fm;
					}

					else if (sign(f_min)*sign(fm) <0)
					{
						x_max=xm;
						f_max=fm;
					}
					else 
					{
						x1=xm;
						break;
					}
				}
				else
				{
					if (f_max<=f_min)
					{
						xnew=xm+(xm-x_min)*fm/sqrt(s);
					}
					else
					{
						xnew=xm-(xm-x_min)*fm/sqrt(s);
					}
					
					G_new=((((cox1/epsi)*(vgs1-xnew))*((cox1/epsi)*(vgs1-xnew)))-(B*exp(beta*(xnew-V))));
					
					if (G_new<=0)
						G_new=realmin;

					A_new=(B*exp(beta*(xnew-V)));
					if (A_new<=0)
						A_new=realmin;

					exp1_new=sqrt(G_new);
					exp4_new=(exp1_new/sqrt(A_new));
					theta_new=(beta*tsi*.5*exp1_new)+asinh(exp4_new);
                
                	hyp_func(&fnew,vgs1, vgs2, tsi, cox1, cox2, V, xnew, FLAG_theta, exp1_new, theta_new);

					if (abs1(fnew) <= ERH)
					{
						x1=xnew;
                        break;
					}

					if ((sign(fm)*sign(fnew) < 0) && (xnew > xm))
					{
						x_max=xnew;
						f_max=fnew;
						x_min=xm;
						f_min=fm;
					}

					else if ((sign(fm)*sign(fnew) < 0) && (xnew < xm))
					{
						x_min=xnew;
						f_min=fnew;
						x_max=xm;
						f_max=fm;
					}

					else if (xnew==xm)
					{
						if((sign(f_max)*sign(fnew) >0))
						{
							x_max=xm;
							f_max=fm;
						}
						else if ((sign(fnew)==0))
						{
							x1=xnew;
							break;
						}

						else 
						{
							x_min=xm;
							f_min=fm;
						}
					}
					else if (sign(fnew==0))
					{
						x1=xnew;
						break;
					}
					else if ((sign(f_min)*sign(fnew) < 0))
					{
						x_max=xnew;
						f_max=fnew;
					}
					else if((sign(f_max)*sign(fnew) <= 0))
					{
						x_min = xnew;
						f_min=fnew;
					}
					else 
					{
						*FLAG_err=1;
						break;
					}
				}
			}
			else 
			{
				if((abs1(f_min)<=H_bound*realmin))
				{
					x1=x_min;
					break;
				}
				else if ((abs1(f_max)<=H_bound*realmin))
				{
					x1=x_max;
					break;
				}
				else if (((iter_rid==1) &&(FLAG_theta==0)&&(flag_approx==1)))
				{
					flag_approx=0;
					hyp_limit(psi_limit, method_overhead,iter_overhead ,vgs1, vgs2, V, csi, cox1, cox2, psi_wi, psigzp, flag_approx);

					x_min=*psi_limit;

					x_max=min(psi_wi+10*eps, psigzp-10*eps);
					G_min=((((cox1/epsi)*(vgs1-x_min))*((cox1/epsi)*(vgs1-x_min)))-(B*exp(beta*(x_min-V))));
					if (G_min<=0)
						G_min=realmin;

					 A_min=(B*exp(beta*(x_min-V)));
                
					if (A_min<=0)
						A_min=realmin;

					exp1_min=sqrt(G_min);
					exp4_min=(exp1_min/sqrt(A_min));
					theta_min=(beta*tsi*.5*exp1_min)+asinh(exp4_min);
                
					G_max=((((cox1/epsi)*(vgs1-x_max))*((cox1/epsi)*(vgs1-x_max)))-(B*exp(beta*(x_max-V))));
					if (G_max<=0)
						G_max=realmin;

					A_max=(B*exp(beta*(x_max-V)));
					if (A_max<=0)
						A_max=realmin;

					exp1_max=sqrt(G_max);
					exp4_max=(exp1_max/sqrt(A_max));
					theta_max=(beta*tsi*.5*exp1_max)+asinh(exp4_max);

					if ((theta_min>700)||(theta_max>700))
						FLAG_theta=1;
						
					hyp_func(&f_min,vgs1, vgs2, tsi, cox1, cox2, V, x_min, FLAG_theta, exp1_min, theta_min);
					hyp_func(&f_max,vgs1, vgs2, tsi, cox1, cox2, V, x_max, FLAG_theta, exp1_max, theta_max);
				}
				else if (iter_rid<=2)
				{
					flag_approx=0;
					hyp_limit(psi_limit, method_overhead,iter_overhead,vgs1, vgs2, V, csi, cox1, cox2, psi_wi, psigzp, flag_approx);
					x_min=*psi_limit;
                
					x_max=min(psi_wi+10*eps, psigzp-10*eps);
					G_min=((((cox1/epsi)*(vgs1-x_min))*((cox1/epsi)*(vgs1-x_min)))-(B*exp(beta*(x_min-V))));
					if (G_min<=0)
						G_min=realmin;

					A_min=(B*exp(beta*(x_min-V)));
                
					if (A_min<=0)
						A_min=realmin;

					exp1_min=sqrt(G_min);
					exp4_min=(exp1_min/sqrt(A_min));
					theta_min=(beta*tsi*.5*exp1_min)+asinh(exp4_min);
                
					G_max=((((cox1/epsi)*(vgs1-x_max))*((cox1/epsi)*(vgs1-x_max)))-(B*exp(beta*(x_max-V))));
					if (G_max<=0)
						G_max=realmin;

					A_max=(B*exp(beta*(x_max-V)));
					if (A_max<=0)
						A_max=realmin;

					exp1_max=sqrt(G_max);
					exp4_max=(exp1_max/sqrt(A_max));
					theta_max=(beta*tsi*.5*exp1_max)+asinh(exp4_max);
                
					FLAG_theta=1;
                
                	hyp_func(&f_min,vgs1, vgs2, tsi, cox1, cox2, V, x_min, FLAG_theta, exp1_min, theta_min);
					hyp_func(&f_max,vgs1, vgs2, tsi, cox1, cox2, V, x_max, FLAG_theta, exp1_max, theta_max);
				}

				else
				{
					*FLAG_err=1;
					break;
				}
			}
		}
	}

	*iteration=j;
	*psi=x1;

	if (j==max_itr_NR)
		*FLAG_err=1;

	if(*FLAG_err==1)
	{
		*method=9;
		if(j<5)
		{
			*psi=x_min;
		}
		else
		{
			if ((abs1(f_min)< sqrt(H_bound)*realmin))
				*psi=x_min;
            
			else if ((abs1(f_max)<sqrt(H_bound)*realmax))
				*psi=x_max;
            
			else
				*psi=((x_max+x_min)/2);
		}
	}

	if(*psi > psi_wi)
	{
		*psi=psi_wi;
	}
	return;
}

void hyp_limit(double*psi_limit,double* method_overhead,double * iter_overhead,double vgs1, double vgs2,double V,double csi,double cox1,double cox2,double psi_wi,double psigzp,double FLAG_approx)
{
	double x0=0 , x1=0, fx0=0 ,fx1=0, fx2=0, fg0=0, fg1=0, fg2=0, fg3=0, h=0,i=0;
	if(FLAG_approx==1)
	{
		*method_overhead=1;
		*iter_overhead=1;
		
		x0=min(psigzp,psi_wi);
		
		fx0=B*exp(beta*(x0-V));
    	fx1=(cox1/epsi)*(cox1/epsi)*(vgs1-x0);
   	fx2=((cox2/epsi)*(cox2/epsi)*(x0-vgs2)*(1/((1+(cox2/csi))*(1+(cox2/csi)))));   
    	fg0=fx0-fx1*(vgs1-x0)+fx2*(x0-vgs2);
    	fg1=beta*fx0+2*fx1+2*fx2;
    	fg2=beta*beta*fx0- (2*cox1*cox1)/(epsi*epsi)+(2*cox2*cox2)/(epsi*epsi*(cox2/csi + 1)*(cox2/csi + 1));
    	fg3=beta*beta*beta*fx0;
    
    	h=-(fg0/fg1)*(1+((fg0*fg2)/(2*fg1*fg1))+((fg0*fg0)*(3*fg2*fg2-fg1*fg3)/(6*fg1*fg1*fg1*fg1)));
    
    	*psi_limit=x0+h;
    	
    }
    
    else
    {
	    *method_overhead=2;
	    x0=min(psigzp,psi_wi);
	    
	    for(i=0;i<=max_itr_NR;i++)
	    {
		    fx0=B*exp(beta*(x0-V));
        	 fx1=(cox1/epsi)*(cox1/epsi)*(vgs1-x0);
			 fx2=((cox2/epsi)*(cox2/epsi)*(x0-vgs2)*(1/((1+(cox2/csi))*(1+(cox2/csi)))));
        
        	 fg0=fx0-fx1*(vgs1-x0)+fx2*(x0-vgs2);
        	 fg1=beta*fx0+2*fx1+2*fx2;
		    
			 x1=x0- fg0/fg1;
		
		 	 if(( abs1(x0 - x1) <= 10*eps))
            break;
        	 else
            x0=x1;
       }
       
       *psi_limit=x0+eps;
       *iter_overhead=i;
   return;
	}
}
void hyp_func_der (double *value,double vgs1,double vgs2,double tsi,double cox1,double cox2,double V,double psi,double G,double exp1,double exp0,double alpha,double theta)
{
	double A=0, exp11=0, exp2=0, exp3=0, exp7=0,exp5=0,exp6=0;
	
	A=B*exp(beta*(-V));
	exp11=exp(-0.5*beta*psi);
//printf("%18.16le\n",exp11);
	if (exp0 <= (1e25*realmin))
     exp0=1e25*realmin;
     
   exp2=1/sqrt(1 + G/exp0);
	exp3=beta*exp0 - (cox1*cox1)/(epsi*epsi)*(2*psi - 2*vgs1);
	exp7= sqrt(A*1.0/G);
//	printf("chk %18.16Le %18.20e %18.20en\n",sqrt(A/G),A,G);
	if(exp7 <= pow(realmin,1.0/3.0))
		exp7=pow(realmin,1.0/3.0);
		
	exp5=exp2*(.5*beta*exp11/exp7 + (A*exp11*(exp3))/(2*(exp7)*(exp7)*(exp7)*(G)*G));
	exp6=(exp5 + .25*beta*tsi/exp1*(exp3));
//	printf("%18.16le %18.16le %18.16le %18.16le %18.16le %18.16le\n\n",exp11,exp2,exp3,exp5,exp6,exp7);
	*value=(exp6*exp7*cosh(theta) + (beta*((epsi*exp3*coth(theta))/(2*cox2*exp1) - (epsi*exp1*exp6*(coth(theta)*coth(theta) - 1))/cox2))/(2*exp((beta*(vgs2 + (epsi*exp1*coth(theta))/cox2))*0.5)) - (A*exp3*sinh(theta))/(2*G*G*exp7));
	return;
}

void hyp_func (double *value, double vgs1,double vgs2,double tsi,double cox1,double cox2,double V,double psi,double FLAG_theta,double exp1,double theta)
{
	double a=0,b=0;
	if(FLAG_theta==1)
	{
		if(theta>THETA_thr)
		{
			*value=exp1-(cox2/epsi)*((-2/beta)*log(0.5*(sqrt(B)/exp1))+(-2/beta)*theta-vgs2+V);
		}
		else
		{
			*value=exp1*(1/tanh(theta))+(cox2/epsi)*(((2/beta)*(log((sqrt(B)/exp1)*sinh(theta))))+vgs2-V);
		}
	}
	
	else
	{
		a=(exp((beta*.5)*(-(epsi/cox2)*exp1*(1/tanh(theta))-vgs2)));
		 b=((sqrt(B*exp(beta*(-V)))*sinh(theta))/exp1);
		 *value=a-b;
		// *value=exp((beta*.5)*(-(epsi/cox2)*exp1*(1/tanh(theta))-vgs2))-(sqrt(B*exp(beta*(-V)))*sinh(theta))/exp1;
	}
	return;
}




                

void root_find_trig (double*psi1, double*psi_pi_orig, double*method, double* method_overhead, double* iteration,double* iter_overhead,double* time_overhead,double* FLAG_err,double vgs1,double vgs2,double V, double tsi,double cox1,double cox2,double psigzp,double psi_r,double fid)
{
	double x_min=0 ,x_max=0, xnew=0, x1=0, xm=0, f_min=0, f_max=0, fm=0, fnew=0, j=0, s=0,a_v=0;
	*FLAG_err=0;
	clock_t tstart_overhead=0,tstop_overhead;
	//tstart_overhead=clock();
	
	trig_limit(psi_pi_orig, method_overhead, iter_overhead,vgs1, vgs2, V, tsi, cox1,  psigzp, fid);
	//printf("psi_pi==%18.20e\n",*psi_pi_orig);
	//tstop_overhead=clock();
	//*time_overhead = ((double)( tstop_overhead - tstart_overhead )) / CLOCKS_PER_SEC;

	x_min=psigzp+50*eps;
	
	x_max=*psi_pi_orig-50*eps;
//printf("x_testing %18.20e %18.20e %18.20e\n",psigzp,*psi_pi_orig,eps);
	//printf("x series start %18.20e %18.20e\n\n",x_max,x_min);
	trig_func(&f_min, vgs1, vgs2, tsi, cox1, cox2, V, x_min);
	
	trig_func(&f_max, vgs1, vgs2, tsi, cox1, cox2, V, x_max);
	//printf("f series start %18.20e %18.20e\n\n",f_max,f_min);
	x1=psi_r;

	for(j=1;j<=max_itr_bis;j++)
	{
		*method=2;
		xm=0.5*(x_min+x_max);
		
		if ((abs1(x_min-x_max)<=(Neps_T+1)*eps)||(abs1(x_min-xm)<=(Neps_T)*eps)||(abs1(x_max-xm)<=(Neps_T)*eps))
		{//printf("q1\n");
			x1=xm;
			break;
		}
		
		if(accu_level==1)
		{	// printf("q2\n");
			if(abs1(f_min)<= ERT)
			{
				x1=x_min;
				break;
			}
			else if ((abs1(f_max) <= ERT))
			{
				x1=x_max;
				break;
			}
		}
		//printf("%18.20e %18.20e \n",f_max,f_min);
		if ((sign(f_max)*sign(f_min)<0) && ((abs1(f_min)>T_bound*realmin)&&abs1(f_max)>T_bound*realmin))
		{	//printf("q3\n");
			trig_func(&fm,vgs1, vgs2, tsi, cox1, cox2, V, xm);
			//printf("fm==%18.20e\n",fm);
			if((abs1(fm) <= ERT))
			{
				x1=xm;
				break;
			}
a_v=fm*fm-f_min*f_max;
			s=abs1(a_v);
			//printf("%18.20e\n",s);
			if((s==0)||(1-abs1(fm/sqrt(s))) < eps )
			{
				*method=2;
				if(sign(f_min)*sign(fm) > 0)
				{
					x_min=xm;
					f_min=fm;
				}
				else if (sign(f_min)*sign(fm) <0)
				{
					x_max=xm;
					f_max=fm;
					
				}
				else
				{
					x1=xm;
					break;
				}
			}

			else
			{
				if(f_max<=f_min)
					{
						xnew=xm+(xm-x_min)*fm/sqrt(s);
					}
				else
					xnew=xm-(xm-x_min)*fm/sqrt(s);
				
		/*	if((xnew==x_min)||(xnew==x_max)||(xnew==xm))
				{
					*method=2;
					if(sign(f_min)*sign(fm) > 0)
					{
					x_min=xm;
					f_min=fm;
					}
					else if (sign(f_min)*sign(fm) <0)
					{
					x_max=xm;
					f_max=fm;
					
					}
					else
					{
					x1=xm;
					break;
					}
				}
				else
				{*/					
					
				trig_func(&fnew,vgs1, vgs2, tsi, cox1, cox2, V, xnew);
			     //	printf("x series %18.20Le %18.20Le %18.20Le %18.20Le\n\n",x_max,x_min,xnew,xm);
				//printf("f series %18.20Le %18.20Le %18.20Le %18.20Le\n\n",f_max,f_min,fnew,fm);

				if(abs1(fnew) <= ERT)
				{
					x1=xnew;
					//printf( "1st f_max track_2  %18.20e   \n\n",f_max);             /* the prob lies here...BEWARE...!!!    */
					break;
				}

				if((sign(fm)*sign(fnew) < 0) && (xnew > xm))
				{
					x_max=xnew;
					f_max=fnew;
					x_min=xm;
					f_min=fm;
				}

				else if ((sign(fm)*sign(fnew) < 0) && (xnew < xm))
				{
					x_min=xnew;
					f_min=fnew;
					x_max=xm;
					f_max=fm;
				}

				else if(xnew==xm)
				{
					
					if((sign(f_max)*sign(fnew) >0))
					{
						x_max=xm;
						f_max=fm;
					}
					else if (sign(fnew)==0)
					{
						x1=xnew;
						break;
					}
					else
					{
						x_min=xm;
						f_min=fm;
					}
				}
				else if (sign(fnew)==0)
				{
					x1=xnew;
					
					break;
				}
				
				else if (sign(f_min)*sign(fnew) < 0)
				{
					x_max=xnew;
					f_max=fnew;
					
				}

				else if ((sign(f_max)*sign(fnew) <= 0))
				{
					x_min=xnew;
					f_min=fnew;
					//printf( "2nd f_max track_2  %18.20e  %18.20e \n\n",f_max,fnew);
				}

				else 
				{
					*FLAG_err=1;       /*File format of fid not considered , chk it*/
					//fprintf(fid, 'Error in ridder loop of psi vgs1=%18.20f vgs2=%18.20f V=%18.20f iter=%d \n', vgs1, vgs2, V, j);
					//printf("in 1st case %Le",*FLAG_err);
					break;
					
				}
			}
		}
	
		else
		{
			if(abs1(f_min)<=T_bound*realmin)
			{
				x1=x_min;
				break;
			}
			else if ((abs1(f_max)<=T_bound*realmin))
			{
				x1=x_max;
				break;
			}
			else         /*File pointer fid ,chk out*/
			{
				//fprintf( fid, 'Error entering ridder due to Sign error.. root not bracketed ...j=%d vgs1= %18.20f  vgs2=%18.20f V=%18.20f  \n', j, vgs1, vgs2, V);
				*FLAG_err=1;
				//printf("in 2nd case %Le\n\n",*FLAG_err);
				break;
			}
		}
	}
	//printf("fseries %18.20Le %18.20Le %18.20Le %18.20Le\n\n",f_max,f_min,fnew,fm);
	*iteration=j;
	*psi1=x1;
//	printf("Flag err =%18.20Le\n\n\n",*FLAG_err);
	if((j==max_itr_bis))
	{
		*FLAG_err=1;
	}
	//printf("Flag err =%18.20e\n\n\n",*FLAG_err);

	if(*FLAG_err==1)
	{
		
		*method=9;
		if(j<2)
			*psi1=psi_r;
		else
		{
			if((abs1(f_min)< sqrt(T_bound)*realmin))
				*psi1=x_min;
			else if ((abs1(f_max)<sqrt(T_bound)*realmax))
				*psi1=x_max;
			else
				*psi1=((x_max+x_min)/2.0);
		}
	}
	
return;
}

void trig_limit (double* psi_pi,double* method,double* iteration_overhead,double vgs1,double vgs2,double V,double tsi,double cox1,double psigzp, double fid)  /* see fid type as file open der*/
{
	double x0, x1, x_min, x_max, xm=0, xnew, xans, FLAG, FLAG_err, psi_star=0, G, exp1, exp2, exp3, f_min=0, f_max=0, fm=0, fnew=0, x_guess, fg0=0, fg1, fg2,term1,s,j;
	
	x0=psigzp+1*eps;
	//printf("%18.20e\n",x0);
	x_min=x0;
	xans=x_min;
	FLAG=0;
	FLAG_err=0;
	
	G=-((cox1/epsi)*(vgs1-x_min))*((cox1/epsi)*(vgs1-x_min))+B*exp(beta*(x_min-V));
	if(G<=0)
		G=realmin;
		
	exp2=(B*exp(beta*(x_min-V)));
	f_min= pi-beta*tsi*0.5*sqrt(G)-asin(sqrt(G/exp2));

	if((vgs1-V)<18.2)
	{
		x_max=vgs1;
		G=(B*exp(beta*(x_max-V)));
		
		if (G<=0)
			G=realmin;
			
		f_max= pi/2.0-beta*tsi*0.5*sqrt(G);
		
		if((f_min*f_max) > 0)
		{
			*psi_pi=x_max+eps;
			*method=0;
			*iteration_overhead=0;
			return;
		}
	}
	
	x_guess=(1/beta)*log(((2*pi)/(beta*tsi))*((2*pi)/(beta*tsi))+(cox1/epsi)*(cox1/epsi)*(vgs1-psigzp)*(vgs1-psigzp))-(1/beta)*log(B)+V;
	x0=x_guess;
	//printf("raw %18.20e %18.20e %18.20e %18.20e\n",x0,
	for(j=1;j<=100;j++)
	{
		exp1=((4*pi*pi)/(beta*beta*tsi*tsi) + (cox1*cox1*(vgs1 - x0)*(vgs1 - x0))/(epsi*epsi));
		
		fg0=(((log(exp1)-log(B)+(beta*V))/beta) -x0);
    	fg1= (- (cox1*cox1*(2*vgs1 - 2*x0))/(beta*epsi*epsi*exp1) - 1);
    	fg2=((2*cox1*cox1)/(beta*epsi*epsi*(exp1)))-beta*(1+fg1)*(1+fg1);
    	
    	x1 = x0 - 2*fg0*fg1/(2*fg1*fg1-fg0*fg2);
    	//printf("raw [%e] ==%18.20e\n",j,x1);
    	if(abs1(x1-x0)<10*eps)
    	{//printf("break\n");
	    	break;
	   }
	   
	   else
	   	x0=x1;
   }
	
	x_max=x1;
	
	psi_star=x_max;
	//printf("star %18.20e %18.20e %18.20e\n",psi_star,x1,x_max);
	G=(B*exp(beta*(x_max-V)));
	
	if(G<=0)
	{G=realmin;}
	f_max= pi/2.0-beta*tsi*0.5*sqrt(G);
	FLAG=0;
	
	for (j=1;j<=max_itr_NR;j++)
	{
		if(FLAG==0)
		{//printf("go1\n");
			*method=1;
	
			x0=x1;
			
			G=-((cox1/epsi)*(vgs1-x0))*((cox1/epsi)*(vgs1-x0))+B*exp(beta*(x0-V));
			exp2=(B*exp(beta*(x0-V)));
        
			fg0= pi-beta*tsi*0.5*sqrt(G)-asin(sqrt(G/exp2));
         
			exp1=(cox1*cox1*(2*x0 - 2*vgs1))/(epsi*epsi);
			exp3=(-((cox1/epsi)*(vgs1-x0))*((cox1/epsi)*(vgs1-x0))+exp2);
			
			if (abs1(exp2-exp3) < 2*eps)
            fg1=- (beta*tsi*(exp2*beta))/(4*sqrt(exp3));
			else
         	fg1 = (- ((exp2*beta - exp1)/(exp2) - (beta*(exp3))/(exp2))/(2*sqrt((1 - (exp3)/(exp2)))*sqrt(((exp3)/(exp2)))) - (beta*tsi*(exp2*beta - exp1))/(4*sqrt(exp3)));
         	
			term1=fg0/fg1;
         
			if((j==1)&&(abs1(term1)<NR_thr))
         		{//printf("go1-a\n");
				FLAG=1;}
         					
			if (FLAG==0)
			{
				x1=x0-fg0/fg1;
				
				if((x1>x_max)||(x1 < x_min)|| (j >4))
				{
					FLAG=1;
					if((j>1)||(j<1))
					{
						if((x1<x_max)&&(x1> x_min))
						{
							fg0=fg1;
							x0=x1;
						}
						if(fg0*f_min>0)
						{
							x_min=x0;
							f_min=fg0;
						}
						else
						{
							x_max=x0;
							f_max=fg0;
						}
					}
				}
				
				if(( abs1(x0 - x1)<= Neps_TL*eps))
					break;
				else
					x0=x1;
			}
		}
			
			else
			{
				*method=2;
				xm=0.5*(x_min+x_max);
				
				if ((abs1(x_min-x_max)<=(Neps_TL+1)*eps)||(abs1(x_min-xm)<=(Neps_TL)*eps)||(abs1(x_max-xm)<=(Neps_TL)*eps))
				{//printf("go2\n");
					x1=xm;
					break;
				}
				
				
				 if ((sign(f_max)*sign(f_min)<0) && (abs1(f_min)>T_bound*realmin)&&(abs1(f_max)>T_bound*realmin))
				 {//printf("go3\n");
					 
					 G=-((cox1/epsi)*(vgs1-xm))*((cox1/epsi)*(vgs1-xm))+B*exp(beta*(xm-V));
					 exp2=(B*exp(beta*(xm-V)));
					 fm= pi-beta*tsi*0.5*sqrt(G)-asin(sqrt(G/exp2));
					 
					 s=abs1(fm*fm-f_min*f_max);
					 
					 if((s==0)||((1-abs1(fm/sqrt(s)))<eps))
					 {
						 *method=2;
						 if(sign(f_min)*sign(fm) > 0)
						 {
							 x_min=xm;
							 f_min=fm;
						 }
						 else if (sign(f_min)*sign(fm) > 0)
						 {
							 x_max=xm;
						     f_max=fm;
						 }
						 else
						 {
							 x1=xm;
							 break;
						 }
					 }
					
					 else
					 {
						 s=sqrt(s);
						 if(f_min>=f_max)
							 xnew=xm+(xm-x_min)*fm/s;
						 else
							 xnew=xm-(xm-x_min)*fm/s;
							
	/*if((xnew==x_min)||(xnew==x_max)||(xnew==xm))
				{
					*method=2;
					if(sign(f_min)*sign(fm) > 0)
					{
					x_min=xm;
					f_min=fm;
					}
					else if (sign(f_min)*sign(fm) <0)
					{
					x_max=xm;
					f_max=fm;
					
					}
					else
					{
					x1=xm;
					break;
					}
				}
				else
				{*/
					



						 G=-((cox1/epsi)*(vgs1-xnew))*((cox1/epsi)*(vgs1-xnew))+B*exp(beta*(xnew-V));
				       	 exp2=(B*exp(beta*(xnew-V)));
						 fnew= pi-beta*tsi*0.5*sqrt(G)-asin(sqrt(G/exp2));

						 if ((sign(fm)*sign(fnew) < 0) && (xnew > xm))
						 {
							x_max=xnew;
							f_max=fnew;
							x_min=xm;
							f_min=fm;
						 }

						 else if ((sign(fm)*sign(fnew) < 0) && (xnew < xm))
						 {
							 x_min=xnew;
							 f_min=fnew;
							 x_max=xm;
							 f_max=fm;
						 }

						 else if (xnew==xm)
						 {
							 if(sign(f_max)*sign(fnew) >0)
							 {
								 x_max=xm;
								 f_max=fm;
							 }
							 else if (sign(fnew)==0)
							 {
								 x1=xnew;
								 break;
							 }

							 else
							 {
								 x_min=xm;
								 f_min=fm;
							 }
						 }

						 else if (sign(fnew)==0)
						 {
							 x1=xnew;
							 break;
						 }

						 else if((sign(f_min)*sign(fnew) < 0))
						 {
							 x_max=xnew;
							 f_max=fnew;
						 }

						 else if ((sign(f_max)*sign(fnew) <= 0))
						 {
							 x_min=xnew;
							 f_min=fnew;
						 }

						 else 
						 {
							 FLAG_err=1;
							//fid=fopen("fid","w");                                                                          /*f statement to b added*/
							/* fprintf(fid,"Error in ridder loop of psi_pi_orig vgs1=%18.20f vgs2=%18.20f V=%18.20f iter=%d \n", &vgs1, &vgs2, &V, &j);
							 fclose(fid);*/
							 break;
						 }
					}
				}
					
					else
						{
							 if ((abs1(f_min)<=T_bound*realmin))
							 {
								  x1=x_min;
								  break;
							 }

							 else if ((abs1(f_max)<=T_bound*realmin))
							 {
								  x1=x_max;
							     break;
							 }

							 else 
							 {
								  /*fid=fopen("fid","w");
						        fprintf( fid, "Error entering ridder due to Sign error in psi_pi_orig.. root not bracketed ...j=%d vgs1= %18.20f  vgs2=%18.20f V=%18.20f  \n", &j, &vgs1, &vgs2, &V);
						        fclose(fid);*/
								  FLAG_err=1;                            /*fprintf statement to b added*/
								  break;
							 }
					   }
				}
			
		}					
	
	*iteration_overhead=j;
	
	if(j==max_itr_NR)
	{
		//printf("hmmm.1..\n");
		*psi_pi=x_max;
	}
	else
	{		
		//printf("x1==%18.20e psi_star=%18.20e\n",x1,psi_star);
		*psi_pi=x1;
	}
	
	if (FLAG_err==1)
		{//printf("wrong place\n");
		*psi_pi=psi_star;}

	if (vgs1<=*psi_pi)
		*psi_pi=vgs1;
	 
	return;
}

void trig_func(double *value,double vgs1,double vgs2,double tsi, double cox1,double cox2,double V,double psi)
{
	double A=0,G=0,t1=0,t3=0,theta=0,interim=0,exp_part=0,t4=0,sign1=0,d=0,i1=0,r=0,y=0,er=0,p1=0,p2=0,p3=0,p4=0;
	
	A=(B*exp(-beta*V));
	
	if(A<=0)
	{
		A=realmin;
	}
	//i1=(cox1/epsi) * (vgs1-psi);
	//r=(i1*i1);
	//er= exp ( beta * ( psi-V));
	//y=B *er;
	//G=-(r-y);
	G=-((cox1/epsi)*(cox1/epsi)*(vgs1-psi)*(vgs1-psi)-B * exp ( beta * ( psi-V ) ) ) ;
	//printf("trig_func1 G==%18.20e psi==%18.20e %18.20e %18.20e\n\n",G,psi,r,y );
	//G=2.9552640000000000e+006;
	if (G<=0)
    	G=realmin;
    	
   t1=sqrt(G);  //p3=B*exp(beta*(psi-V)); p4=sqrt(p3);
	t3=(t1/sqrt(B*exp(beta*(psi-V))));
	
	if((t3-1) >0 )
    	t3=1;
    	
   //p1=beta*t1*tsi*.5;
  // p2=asin(t3);
   theta=beta*t1*tsi*.5 +asin(t3) ;  
   //printf("%18.20le %18.20le\n",beta*t1*tsi*.5,asin(t3));
   //theta= 3.139735e+000;   /* value given by asin not appropriate*/
   d=theta-pi;
  //printf("%Le %le\n",abs1(d),t3);
//printf("%18.20e\n",theta);
   if (abs1(theta-pi) <= 100*realmin)
   {//printf("yo1\n");
	   t4=exp((beta*.5)*((epsi/cox2)*t1*cot(theta)));
	   interim = (sin(theta)/t1*t4-1/sqrt(A)*exp((-beta*.5)*vgs2));
	   sign1=sign(interim);
	   *value=sign1*realmax;
	   
	}
	else 
	{//printf("abs===%e\n",abs1(-100));
		if(abs1(theta) <= 100*realmin)
		{//printf("yo2-a\n");
			t4=exp((beta*.5)*((epsi/cox2)*t1*cot(theta)));
    		interim = (sin(theta)/t1*t4-1/sqrt(A)*exp((-beta*.5)*vgs2));
    		sign1=sign(interim);
    		*value= sign1*realmax;
    		
    	}
    	else
    	{//printf("yo3\n");
	   	exp_part=(beta*.5)*((epsi/cox2)*t1*cot(theta));
				
			if (exp_part >300)
				exp_part=300;
				
			t4=exp(exp_part);
    
    		*value = (sin(theta)/t1*t4-1/sqrt(A)*exp((-beta*.5)*vgs2));
    		//printf("%18.20e %18.20e %18.20e %18.20e %18.20e\n",*value,t4,cot(theta),t1,theta);
    	}
    }
    
	return;
}



void my_psigzp(double*psigzp , double vgs1, double V,double cox1)       /* psigzp chk type */
{
	double kappa=0 ,x0=0,x1=0,test=0,j=0,g=0,dg=0,dg2=0;
	kappa = (epsi/cox1)*sqrt(B);
	test=((2.0 / beta) * (log ((vgs1 / kappa)) + beta * V * 0.5));
/**/	//printf("%Le\t%Le\n",test,kappa);
	if (vgs1 > 0)
	{
		if(test<vgs1)
		x0=test;
		else
		x0=vgs1;
	}
	else x0=vgs1;
/*printf("%Le\t%Le\n",x0,kappa);	*/
	for (j=0;j<=max_itr_NR;j++)        /*no use of loop as j not inside the loop*/
	{
		g = vgs1 - x0 -kappa*exp(0.5*beta*(x0-V));
		
    	dg  = -1+0.5*beta*(g - vgs1 +x0);
    	dg2 = (beta*beta)*0.25*(g - vgs1 +x0);
    	x1 = x0 - 2*g*dg/(2*dg*dg-g*dg2);
      
    //  /**/printf("%Le\t%Le\t%Le\t%Le\t%Le\n",g,dg,dg2,x1,x0);
      
      if ((x0>x1)&&((x0-x1)<Neps_P*eps))
      break;
      
      if ((x1>x0)&&((x1-x0)<Neps_P*eps))
      break;
      
      x0=x1;
      
     
   }
*psigzp = x0;
return;
}
