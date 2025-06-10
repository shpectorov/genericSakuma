//
// 
// Construction of the universal 2-generated axial algebra of Monster type (al,bt), with al<>4*bt
//
// (c) 2020 C. Franchi, M. Mainardis, S. Shpectorov
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


ring F=(0,al,bt),( lm, lmf, lm2, lm2f), dp;
freemodule(9);
vector am2=gen(1);
vector am1=gen(2);
vector a0=gen(3);
vector a1=gen(4);
vector a2=gen(5);
vector s1=gen(6);
vector s2=gen(7);
vector s2f=gen(8);
vector pp=gen(9);
list a=am2,am1,a0,a1,a2,s1,s2,s2f;
list r=pp,pp,pp,pp,pp,pp,pp,pp;
list Mult=r,r,r,r,r,r,r,r;

Mult[1][1]=am2; Mult[2][2]=am1; Mult[3][3]=a0; Mult[4][4]=a1; Mult[5][5]=a2; 

Mult[1][2]=s1+bt*(a[1]+a[2]); Mult[2][1]=Mult[1][2];
Mult[2][3]=s1+bt*(a[2]+a[3]); Mult[3][2]=Mult[2][3];
Mult[3][4]=s1+bt*(a[3]+a[4]); Mult[4][3]=Mult[3][4];
Mult[4][5]=s1+bt*(a[4]+a[5]); Mult[5][4]=Mult[4][5];
Mult[1][3]=s2+bt*(a[1]+a[3]); Mult[3][1]=Mult[1][3];
Mult[3][5]=s2+bt*(a[3]+a[5]); Mult[5][3]=Mult[3][5];
Mult[2][4]=s2f+bt*(a[2]+a[4]); Mult[4][2]=Mult[2][4];


Mult[2][6]=(al-bt)*a[6]+(lmf*(1-al)+bt*(al-bt-1))*a[2]+
  (1/2)*(al-bt)*bt*(a[1]+a[3]);
Mult[3][6]=(al-bt)*a[6]+(lm*(1-al)+bt*(al-bt-1))*a[3]+
  (1/2)*(al-bt)*bt*(a[2]+a[4]);
Mult[4][6]=(al-bt)*a[6]+(lmf*(1-al)+bt*(al-bt-1))*a[4]+
  (1/2)*(al-bt)*bt*(a[3]+a[5]);
Mult[3][7]=(al-bt)*a[7]+(lm2*(1-al)+bt*(al-bt-1))*a[3]+
  (1/2)*(al-bt)*bt*(a[1]+a[5]);
  
Mult[6][2]=Mult[2][6]; 
Mult[6][3]=Mult[3][6]; 
Mult[6][4]=Mult[4][6]; 
Mult[7][3]=Mult[3][7]; 
 
// Times returns the product of the vectors u and v

proc Times (vector u, vector v)
{
vector ans=[0];
list rem=list();
 for(int i=1; i<=8; i=i+1)
  { if (u[i]<>0)
     {
        for(int j=1; j<=8; j=j+1)
          {
              if (v[j]<>0)
               { if (Mult[i][j]<>pp)
                   {
                   ans=ans+u[i]*v[j]*Mult[i][j];
                   }
                   else
                  { list l=list(i,j), u[i]*v[j];
                      rem=rem+l;
                   }             
          }  
     }
  }
  }
    if (size(rem)==0)
       {   return(ans);
           }
   else
       {list k=ans, rem;
      return(k);}
}

//
// Symmetries: 
// tau0 is the Miyamoto automorphism associated to the axis a0   
// 

proc tau0 (vector v)
  {
    return(v[5]*a[1]+v[4]*a[2]+v[3]*a[3]+v[2]*a[4]+v[1]*a[5]+v[6]*a[6]+v[7]*a[7]+v[8]*a[8]);        
  } 

// ff is the automorphism of the ring F swapping lm with lmf and lm2 with lm2f 
 
map ff=F, lmf,lm, lm2f, lm2;

// fold is the semi-automorphism of the algebra swapping a0 and a1

proc fold (vector v)
{
poly v1=v[1];
poly v2=v[2];
poly v3=v[3];
poly v4=v[4];
poly v5=v[5];
poly v6=v[6];
poly v7=v[7];
poly v8=v[8];
poly ffv2=ff(v2); 
poly ffv3=ff(v3);
poly ffv4= ff(v4); 
poly ffv5=ff(v5); 
poly ffv6=ff(v6); 
poly ffv7=ff(v7); 
poly ffv8=ff(v8);

if (v1==0)
       {vector s=ffv5*a[2]+ffv4*a[3]+ffv3*a[4]+ffv2*a[5]+ffv6*a[6]+ffv8*a[7]+ffv7*a[8];
        return(s);}
        else
          {return("fail");}
          } 

// Function lambda
          
proc lambda (vector v)
{
 if (v[8]==0) 
  { poly p= lm2*v[1]+lm*v[2]+v[3]+lm*v[4]+lm2*v[5]+(lm-bt-bt*lm)*v[6]+(lm2-bt-bt*lm2)*v[7];
  return(p); }
else
          {return("fail");}
}        

//
// 
// 
 
vector v1= (s1+(bt-lm)*a0+(bt/2)*(a1+am1))/al; // al-eigenvector for a0
vector u1=a1-lm*a0-v1-(a1-am1)/2;              // 0-eigenvector for a0

vector uu=Times(u1, u1)[1];  
vector uv=Times(u1, v1)[1]; 
vector vv=Times(v1, v1)[1]; 

// (uu+uv) is a 0-eigenvector for a0 and contains s2f: thus lambda(uu+uv) must be 0 and we deduce lambda(s2f)

vector w=(uu+uv)-(uu+uv)[8]*s2f;
 
poly ls2f=-lambda(w)/((uu+uv)[8]); // ls2f is the value of lambda(s2f)

// update lambda

proc lambda (vector v)
{
   poly p= lm2*v[1]+lm*v[2]+v[3]+lm*v[4]+lm2*v[5]+(lm-bt-bt*lm)*v[6]+(lm2-bt-bt*lm2)*v[7]+ls2f*v[8];
  return(p); 
  } 
  
// u1*u1 is a 0-eigenvector for a0, hence it has 0-projection on a0, and it contains 1/(al^2)*s1*s1  
// Hence we can obtain the value of lambda(s1*s1)

poly ls1s1=-(al^2)*lambda(uu);; // lambda(s1*s1)

//
// Express a0*s2f
//
vector w=uu-vv -lambda(uu-vv)*a0; // this is a 0-eigenvector for a0

// Times(a0, w)  this vector must be zero
// it involves a0*s2f with coefficient (al-2*bt)/(2*al), hence we can express a0*s2f

Mult[3][8]=-2*al*Times(a0, w)[1]/(al-2*bt);
Mult[8][3]=Mult[3][8];

//
// Express s1*s1
//

// a0*(u1*u1+u1*v1)=al*u1*v1, and  al*u1*v1 contains -s1*s1/al

Mult[6][6]=-al*(Times(a0,uu+uv)-al*uv);

//
// Express a3 and am3
//

// s1*s1 contains bt*(al-bt)^2*(al-4*bt)/(4*(al-2*bt)) * am2
//
// so fold(s1*s1)-s1*s1=0 and it contains bt*(al-bt)^2*(al-4*bt)/(4*(al-2*bt))* a3 

vector A3=-4*(al-2*bt)*(fold(Times(s1,s1)-bt*((al-bt)^2)*(al-4*bt)*am2/(4*(al-2*bt)))-Times(s1,s1))/(bt*((al-bt)^2)*(al-4*bt));
     
// Update function f

proc f (vector v)
{
poly v1=v[1];
poly v2=v[2];
poly v3=v[3];
poly v4=v[4];
poly v5=v[5];
poly v6=v[6];
poly v7=v[7];
poly v8=v[8];
poly ffv1=ff(v1);
poly ffv2=ff(v2); 
poly ffv3=ff(v3);
poly ffv4= ff(v4); 
poly ffv5=ff(v5); 
poly ffv6=ff(v6); 
poly ffv7=ff(v7); 
poly ffv8=ff(v8);
       vector s=ffv1*A3+ffv5*a[2]+ffv4*a[3]+ffv3*a[4]+ffv2*a[5]+ffv6*a[6]+ffv8*a[7]+ffv7*a[8];
        return(s);
          } 
     
// Express am3 and a4

vector Am3=tau0(A3);
vector A4=f(Am3);

// Define tau1 (Miyamoto involution associated to a1)

proc tau1(vector v)
{
vector s=v[1]*A4+v[2]*A3+v[3]*a2+v[4]*a1+v[5]*a0+v[6]*s1+v[7]*s2+v[8]*s2f;
 return(s); 
}

// Add the products am2*s1 and a2*s1

Mult[1][6]=(al-bt)*s1+(lm*(1-al)+bt*(al-bt-1))*a[1]+(1/2)*(al-bt)*bt*(Am3+a[2]);
Mult[6][1]=Mult[1][6];
Mult[5][6]=(al-bt)*s1+(lm*(1-al)+bt*(al-bt-1))*a[5]+(1/2)*(al-bt)*bt*(a[4]+A3);
Mult[6][5]=Mult[5][6];


//
// Express S03, S13, S23
//

vector S03=Times(a0,A3)-bt*(a0+A3);
vector S13=f(S03);
vector S23=tau1(S03);

// D=<tau0,f> should induce D_6 on {S03,S13,S23}
// hence the following must be relations

vector diff0=S23-tau0(S13);
vector diff1=f(diff0); // =-diff0
vector diff2=tau1(diff0);

vector doublediff0=tau0(diff0)-diff0;;
vector doublediff1=tau0(diff1)-diff1;; // equal to -doublediff0
vector doublediff2=tau0(diff2)-diff2;;

// here are some polynomial relations

poly poly1=lambda(diff0); // p1 in the paper
poly poly2=lambda(diff2); // p2 in the paper

// function lambdb: lambda returns the value of lambda_{a_1}(v)

proc lambdb(vector v)
{
poly l=lambda(f(v));
 return(ff(l)); 
}

poly  poly3=lambdb(diff0); // p3 in the paper


// complete the multiplication table

Mult[4][7]=f(Mult[3][8]);;
Mult[7][4]=Mult[7][4];;

Mult[2][7]=tau0(Mult[4][7]);;
Mult[7][2]=Mult[2][7];;

Mult[5][8]=f(Mult[2][7]);;
Mult[8][5]=Mult[5][8];;

Mult[1][8]=tau0(Mult[5][8]);;
Mult[8][1]=Mult[1][8];;

Mult[4][8]=f(Mult[3][7]);;
Mult[8][4]=Mult[4][8];;

Mult[2][8]=tau0(Mult[4][8]);;
Mult[8][2]=Mult[2][8];;

Mult[5][7]=tau1(Mult[3][7]);;
Mult[7][5]=Mult[5][7];;

Mult[1][7]=tau0(Mult[5][7]);;
Mult[7][1]=Mult[1][7];;

Mult[4][7]=tau0(Mult[2][7]);;
Mult[7][4]=Mult[4][7];;

Mult[1][4]=f(Times(a0, A3));;
Mult[4][1]=Mult[1][4];;

Mult[2][5]=tau0(Mult[1][4]);;
Mult[5][2]=Mult[2][5];

Mult[1][5]=f(tau0(f(Times(a0, A4))));;
Mult[5][1]=Mult[1][5];;

vector u2=-s2+(lm2-bt-al*lm2)*a0+((al-bt)/2)*(a2+am2);;  
vector v2=s2+(bt-lm2)*a0+bt*(a2+am2)/2;;  

vector u1u2=Times(u1, u2)[1]; // contains [6][7] with coefficient 1/al
vector v1u2=Times(v1, u2)[1]; // contains [6][7] with coefficient -1/al

// we obtain the products s1*s2,  s1*s2f, and s2*s2 with the method described in Proposition 5.10

vector p12=u1u2+v1u2;;
Mult[6][7]=al*Times(v1, u2)[1]-Times(a0, p12);;
Mult[7][6]=Mult[6][7];; 
Mult[6][8]=f(Mult[6][7]);;
Mult[8][6]=Mult[6][8];; 

vector u2u2=Times(u2, u2)[1];; // contains [7][7] with coefficient 1
vector v2u2=Times(v2, u2)[1];; // contains [7][7] with coefficient -1
vector p22=u2u2+v2u2;;
Mult[7][7]=Times(v2, u2)[1]-Times(a0, p22)/al;;
Mult[8][8]=f(Mult[7][7]);;

// finally we compute s2*s2f. A3*A3-A3 must be the zero vector and it contains s2*s2f

Mult[7][8]=(A3-Times(A3,A3)[1])/(2*Times(A3,A3)[2][2]);;
Mult[8][7]=Mult[7][8];;

//////// two more polynomials

poly poly4=lambda(Times(A4,A4)-A4); // p4 in the paper
poly poly5=ff(poly4); 

proc P(int n)
{
list Po=poly1, poly2, poly3, poly4;
return(Po[n]);
}

// map val returns the polynomial p(lm, lmf, lm2, lm2f) evaluated in (lm, lm, lm2, lm2)
map val=F, lm, lm, lm2, lm2;
poly q1=val(poly1);
poly q2=val(poly2);
poly q3=val(poly3);
poly q4=val(poly5);

proc Q(int n)
{
list Qo=q1, q2, q3, q4;
return(Qo[n]);
}

