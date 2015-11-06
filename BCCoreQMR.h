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


#ifndef CNOID_BCPLUGIN_BCCOREQMR_H
#define CNOID_BCPLUGIN_BCCOREQMR_H

#include <boost/random.hpp>

using namespace std;


namespace cnoid
{

class BCCoreQMR
{
  public:
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixX;
    typedef VectorXd VectorX;
    void NewBuffer   (int aNC) ;
    void DeleteBuffer();

  
    BCCoreQMR(int maxNumGaussSeidelIteration, double gaussSeidelErrorCriterion);
    ~BCCoreQMR();
    bool   callSolver(const MatrixX& Mlcp, const VectorX& b, VectorX& solution, const VectorX& contactIndexToMu,ofstream& os);    
	void setGaussSeidelErrorCriterion(double e);
	void setGaussSeidelMaxNumIterations(int n);
	double * thebuf;
    int SZ;
    int    MAXITE;
  private:
    boost::variate_generator<boost::mt19937, boost::uniform_real<> > randomgen;
	class KKVector
	{
	  public:
		double * elm;
		double& operator()(int i)const {return *(&(elm[i]));}
	};
    KKVector x ;/*01*/
    KKVector b ;/*02*/
    KKVector r0;/*03*/
    KKVector p ;/*04*/
    KKVector q ;/*05*/
    KKVector d ;/*06*/
    KKVector Ap;/*07*/
    KKVector Aq;/*08*/
    KKVector v ;/*09*/
    KKVector w ;/*10*/
    KKVector vt;/*11*/
    KKVector wt;/*12*/
    void iniV_mulMOVO(KKVector* x, const MatrixX& A, const KKVector& b)
    {
		for(int i=0;i<SZ;i++){(*x)(i)=0;for(int j=0;j<SZ;j++)(*x)(i)+=A(i,j)*b(j);}
	} 
    void iniV_mulMTVO(KKVector* x, const MatrixX& A, const KKVector& b)
    {
		for(int i=0;i<SZ;i++){(*x)(i)=0;for(int j=0;j<SZ;j++)(*x)(i)+=A(j,i)*b(j);}
	} 
	void iniV_zero   (KKVector* x){for(int i=0;i<SZ;i++){(*x)(i)=0;}}
	void iniS_squVTVO(double  * x, const KKVector& a                   ){(*x)=0;for(int i=0;i<SZ;i++){(*x)+=a(i)*a(i);}}
	void iniS_mulVTVO(double  * x, const KKVector& a, const KKVector& b){(*x)=0;for(int i=0;i<SZ;i++){(*x)+=a(i)*b(i);}}
    void iniV_minVOVO(KKVector* x, const KKVector& a, const KKVector& b){for(int i=0;i<SZ;i++){(*x)(i)=a(i)-b(i);}} 
	void iniV_mulVOS (KKVector* x, const KKVector& a, const double& b  ){for(int i=0;i<SZ;i++){(*x)(i)=a(i)*b;}}
	void ini_copy  (KKVector* x, const KKVector& a                   ){for(int i=0;i<SZ;i++){(*x)(i)=a(i);}} 
	void ini_copy  (KKVector* x, const VectorX & a                   ){for(int i=0;i<SZ;i++){(*x)(i)=a(i);}} 
	void ini_copy  (VectorX* x , const KKVector & a                  ){for(int i=0;i<SZ;i++){(*x)(i)=a(i);}} 
	void iniV_AmBC (KKVector* x, const KKVector& a, const KKVector& b, const double& c){for(int i=0;i<SZ;i++){(*x)(i)=a(i)-b(i)*c;}} 
	void mulV_S_plusVOS(KKVector* x, const double& a, const KKVector& b, const double& c){for(int i=0;i<SZ;i++){(*x)(i)=(*x)(i)*a+b(i)*c;}} 
	void mulV_S_plusVO (KKVector* x, const double& a, const KKVector& b                 ){for(int i=0;i<SZ;i++){(*x)(i)=(*x)(i)*a+b(i)  ;}} 
	void addV_VO(KKVector* x, const KKVector& a                   ){for(int i=0;i<SZ;i++){(*x)(i)+=a(i);}} 
};

};

#endif
