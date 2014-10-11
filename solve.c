#include<stdio.h>
#include<math.h>
#define PI 3.141592653589793

double EoS(double K, double P, double gamma);
double TOV(double rho, double P, double m, double r);

int main(int argc, char* argv[])
{
  /* 
   * TOV equations with G=c=1
   * dm/dr = 4*pi(r)*r^2
   * dP/dr = -(rho+P)(m+4*pi*r^3*P)/(r^2-2mr)
   *
   */

  /* Define variables */
  FILE *ifp;
  double P,rho,K=7.32,gamma=5.0/3.0, rho_c, P_c;
  double dPdr, r, m, delta_r;
  int i,j=0, plot=0;
  char path[]="./";
  
  ifp=fopen("TOV_output","w");
  plot=(int)atoi(argv[1]);

  for(j=0;j<1000;j++){

    /* central values */
    rho_c=7.42e-5;                   // This is 10^14 g/cm^3 in units of km^-2
    rho_c+=(double)(j*0.0000742);
    //   printf("%d %lf\n",j,rho_c);
    P_c=K*pow(rho_c,1.67);
    rho=rho_c;
    P=P_c;
    r=0.0;
    m=0.0;
    dPdr=0.0;
    //  printf("%.10lf %.10lf %.10lf %.10lf %.10lf\n",rho,P,r,m,dPdr);
    delta_r=0.001;
    r=delta_r;
    
    /* first step 
     * define small radius st P and rho still have central values */
    m=4*PI*pow(delta_r,3.0)*rho_c;
    dPdr=TOV(rho,P,m,r);
    //printf("%.10lf %.10lf %.10lf %.10lf %.10lf\n",rho,P,r,m,dPdr);
    
    for(i=1;i<100000;i++){
      r+=delta_r;
      P+=(dPdr*delta_r);
      rho=EoS(K,P,gamma)+(P/(gamma-1.0));
      m+=(4*PI*pow(r,2.0)*rho*delta_r);
      dPdr=TOV(rho,P,m,r);
      if (P<=1.0e-9){
      	fprintf(ifp,"%.10lf %.10lf %.10lf\n",rho_c,m,r);
      	break;
      }
      //printf("%.10lf %.10lf %.10lf %.10lf %.10lf\n",rho,P,r,m,dPdr);
    }
  }

  rewind(ifp);
  if (plot == 1){
    rewind(ifp);
    FILE *gnuplot = popen("/usr/bin/gnuplot -persist","w");
    //    fprintf(gnuplot,"set title \"Solutions of TOV Equations for n=3/2 polytrope.\"\n");
    fprintf(gnuplot,"set multiplot\n");
    fprintf(gnuplot,"set origin 0.0,0.0\n");
    fprintf(gnuplot,"set size 0.98,0.48\n");
    fprintf(gnuplot,"set logscale x\n");
    fprintf(gnuplot,"set ylabel \"M (km)\"\n");
    fprintf(gnuplot,"set xlabel \"{/Symbol r} (km)\"\n");
    fprintf(gnuplot,"plot \"%sTOV_output\" using 1:2 notitle lt 3\n",path);
    fprintf(gnuplot,"set origin 0.00,0.50\n");
    fprintf(gnuplot,"set size 0.98,0.48\n");
    fprintf(gnuplot,"set ylabel \"M (km)\"\n");
    fprintf(gnuplot,"set xlabel \"R (km)\"\n");
    fprintf(gnuplot,"plot \"%sTOV_output\" using 3:2 notitle lt 3\n",path);
    fprintf(gnuplot,"exit \n");
    pclose(gnuplot);
  }

  fclose(ifp);
  return(0);
}

/* Equation of state
 * P=K*rho^gamma
 * rho=(P/K)^(1/gamma)*/
double EoS(double K, double P, double gamma)
{
  return(pow((P/K),(1.0/gamma)));
}

double TOV(double rho, double P, double m, double r)
{
  double dPdr;
  //  dPdr = (rho+P)*(m+(4*PI*pow(r,3.0)*P))/(pow(r,2.0)-(2*m*r));
  //  printf("%f val1\n",dPdr);
  dPdr = -(rho+P)*(m+(4*PI*pow(r,3.0)*P))/(pow(r,2.0)-(2*m*r));
  //  printf("%f val1\n",dPdr);
  return(dPdr);
}
