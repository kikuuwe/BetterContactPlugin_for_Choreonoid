/* 
 This file is a part of BetterContactPlugin.
 
 Author: Shin'ichiro Nakaoka
 Author: Ryo Kikuuwe
 
 Copyright (c) 2007-2015 Shin'ichiro Nakaoka
 Copyright (c) 2014-2015 Ryo Kikuuwe
 Copyright (c) 2007-2015 National Institute of Advanced Industrial
                         Science and Technology (AIST)
 Copyright (c) 2014-2015 Kyushu University

 BetterContactPlugin is a plugin for better simulation of frictional contacts.
 
 BetterContactPlugin is a free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 BetterContactPlugin is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with BetterContactPlugin; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 
 Contact: Ryo Kikuuwe, kikuuwe@ieee.org
*/

#include <cnoid/EigenUtil>
#include <fstream>
#include <iomanip>


#include "BCCoreQMR.h"


using namespace cnoid;


BCCoreQMR::BCCoreQMR(int maxNumGaussSeidelIteration, double gaussSeidelErrorCriterion):
	    randomgen(boost::mt19937(), boost::uniform_real<>(-1.0, 1.0))
{
	thebuf=0;
	SZ = 0;
	setGaussSeidelMaxNumIterations(maxNumGaussSeidelIteration);
	setGaussSeidelErrorCriterion  (gaussSeidelErrorCriterion);
    randomgen.engine().seed();
}
void BCCoreQMR::setGaussSeidelMaxNumIterations(int n)
{
	MAXITE = n;
}
void BCCoreQMR::setGaussSeidelErrorCriterion(double e)
{
}

BCCoreQMR::~BCCoreQMR()
{
  DeleteBuffer();
}

void BCCoreQMR::NewBuffer(int NC3)
{
  if(NC3<=0){SZ=0;DeleteBuffer();return;}
  SZ= NC3;
  thebuf = new double[SZ*12];
  int i = 0;
  x .elm = &thebuf[i];i+= SZ;
  b .elm = &thebuf[i];i+= SZ;
  r0.elm = &thebuf[i];i+= SZ;
  p .elm = &thebuf[i];i+= SZ;
  q .elm = &thebuf[i];i+= SZ;
  d .elm = &thebuf[i];i+= SZ;
  Ap.elm = &thebuf[i];i+= SZ;
  Aq.elm = &thebuf[i];i+= SZ;
  v .elm = &thebuf[i];i+= SZ;
  w .elm = &thebuf[i];i+= SZ;
  vt.elm = &thebuf[i];i+= SZ;
  wt.elm = &thebuf[i];
}

void BCCoreQMR::DeleteBuffer()
{
  if(thebuf!=0)  delete [] thebuf;
  thebuf=0;
}


bool BCCoreQMR::callSolver(const MatrixX& A, const VectorX& ab, VectorX& ax, const VectorX& contactIndexToMu, ofstream& os)    
{
	const double EPSTHRESH = 1.0e-20;
	ini_copy(&x, ax);
	for(int i=0;i<SZ;i++){b(i)=-ab(i);}
	/*********************/
	int cause_of_termination=0; 
	iniV_zero(&p);
	iniV_zero(&q);
	iniV_zero(&d);
	iniV_mulMOVO(&Ap ,A, x );
	iniV_minVOVO(&r0 ,b, Ap);
	double rho;
	iniS_squVTVO(&rho,r0    ); rho = sqrt(rho);
	if(fabs(rho)>EPSTHRESH)
	{
		double c2 = 1, eps = 1, xi = 1,  theta2= 0, eta =-1;
		iniV_mulVOS(&v, r0, 1./rho );
#if 1
		ini_copy   (&w, v)          ;
#else
		for(int i=0;i<SZ;i++){w(i)= v(i);}
		double w_abs; iniS_squVTVO(&w_abs,w    );
		w_abs = 1./sqrt(w_abs);
		for(int i=0;i<SZ;i++){w(i)*=w_abs;}
#endif
		for(int iteration=0;iteration<MAXITE ;iteration++)
		{
			double delta; iniS_mulVTVO (&delta, w, v ); 
			if(fabs(eps)<EPSTHRESH){cause_of_termination=1;break; }
			mulV_S_plusVO(&p  , -xi * delta / eps,  v);
			mulV_S_plusVO(&q  , -rho* delta / eps,  w); 
			iniV_mulMOVO (&Ap ,  A ,  p); 
			iniV_mulMTVO (&Aq ,  A ,  q); 
			iniS_mulVTVO (&eps,  q , Ap);
			if(fabs(delta)<EPSTHRESH){cause_of_termination=1;break;}
			double beta = eps / delta ;
			double rho_old = rho;
			iniV_AmBC(&vt , Ap, v, beta ) ;
			iniV_AmBC(&wt , Aq, w, beta ) ; 
			iniS_squVTVO(&rho, vt ); rho = sqrt(rho ) ;
			iniS_squVTVO(&xi , wt ); xi  = sqrt(xi  ) ;
			double theta2_old = theta2;
			double c2_old      = c2;
			theta2 = rho*rho/((c2_old)*(beta*beta));
			c2     = 1./(1.+theta2);
			eta   *= - (rho_old*c2)/(beta*c2_old); 
			mulV_S_plusVOS(&d  , theta2_old*c2, p , eta ) ;
			addV_VO       (&x  , d ) ; 
			if(fabs(rho)<EPSTHRESH){cause_of_termination=1;break; }
			if(fabs(xi )<EPSTHRESH){cause_of_termination=1;break; }
			iniV_mulVOS(&v, vt, 1./rho);
			iniV_mulVOS(&w, wt, 1./xi );
		}
	}
	int NC = SZ/3;
	if(NC*3!=SZ) exit(19);
	for(int j=0;j<1;j++)
	{
		double   an = x(j);//contactIndexToMu[j] * 
		double&  a0= x(NC+2*j+0);
		double&  a1= x(NC+2*j+1);
		double at2 = a0 * a0 + a1 * a1;
		//if(at2>an*an && an>EPSTHRESH)
		//{		
		//	at2 = fabs(an)/sqrt(at2);
		//	x(NC+2*j+0) *= 0;//at2;
		//	x(NC+2*j+1) *= 0;//at2;
		//}
		//	if(x(j)<0)x(j)=0;
		//x(NC+2*j+0) = 0;
		//x(NC+2*j+1) = 0;
		
	}
	ini_copy(&ax,x);
//	if(cause_of_termination==1) return false; 
	return true;
}

