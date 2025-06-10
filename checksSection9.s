//
//  Computations for Section 9 of the paper “2-generated axial algebras of Monster type”
// by Clara Franchi, Mario Mainardis, and Sergey Spectorov
//
/////////////////////////////////////////////////////
//
// Construct the universal algebra by running genericsakuma.s to produce the polynomials poly1, poly2, poly3, poly4, q1, q2.


// Check for the proof of Lemma 9.1 

factorize(resultant(q1,q2,lm2)); // we get the five values for lm in the field: 1, 0, bt/2, al/2, (3*al^2+3*al*bt-al-2*bt)/(8*al-4)

// for each of the previous values for lm we compute the respective value of lm2

map val1=F, 1,1, lm2,lm2;
factorize(val1(q1)); 
factorize(val1(q2)); // we get lm2=1 

map val2=F, 0,0, lm2,lm2; 
factorize(val2(q1));  
factorize(val2(q2));  // we get lm2=1

map val3=F, bt/2,bt/2, lm2,lm2; 
factorize(val3(q1));  
factorize(val3(q2));  // we get lm2=bt/2

map val4=F, al/2,al/2, lm2,lm2; 
factorize(val4(q1));  
factorize(val4(q2));  // we get lm2=1

map val5=F, (3*al^2+3*al*bt-al-2*bt)/(8*al-4), (3*al^2+3*al*bt-al-2*bt)/(8*al-4), lm2,lm2; 
factorize(val5(q1));  
factorize(val5(q2));  // we get lm2=(3*al^2+3*al*bt-al-2*bt)/(8*al-4)


// Check for the proof of Lemma 9.2

list L= 1, 0, bt/2, al/2, (3*al^2+3*al*bt-al-2*bt)/(8*al-4);

map vval1=F,1,1,1,1;
vval1(poly1);
vval1( poly2);
vval1(poly3);
vval1(poly4);

// we define some auxiliary functions
// mapL(i,j)  is the endomorphism of the ring F mapping  lm2 to L[i] and lm2f=L[j]

proc mapL(int i, int j) 
 { 
    map m = F, lm, lmf, L[i],L[j];
    return(m);
    }

// given two lists of polynomials (of degree 1) l1 and l2 , samep returns the  polynomials that are, up to a scalar, in the intersection of l1 with l2
    
proc  samep(list l1, list l2) 
{ 
list good=list();
 for (int i=1; i<=size(l1); i=i+1)
     {
	for (int j=1; j<=size(l2); j=j+1)	  	
	   {
               if (resultant(l1[i], l2[j], lm)==0*al) 
		{   
		list LL = l1[i];
	           good=good+LL;
		}
	    }
      }
      return(good);
      }

// given the polynomials pol1 and pol2, degone returns the list of the factors of degree 1 of the resultant of mapL(i,j)(pol1) and mapL(i,j)(pol2) with respect to lmf
      
proc degone(int i, int j, poly pol1, poly pol2) 
{
   list Fact=list();
   list fact = factorize(resultant(mapL(i,j)(pol1), mapL(i,j)(pol2), lmf))[1];
      for (int k = 2; k<= size(fact[1]); k=k+1)
	  {
	       if ( deg(fact[1][k])==1) 
			{
			     list l=fact[1][k];
			     Fact=Fact+l;
			     }
			     }
   return(Fact);
   }

// commonf returns the common factors of degree 1 of the resultants, with respect to lmf, of  mapL(i,j)(pol1) and mapL(i,j)(pol2), and of mapL(i,j)(pol1) and mapL(i,j)(pol3)

 
proc commonf(int i, int j, poly pol1, poly pol2, poly pol3) {   
   return(samep(degone(i,j, pol1, pol2), degone(i,j, pol1, pol3)));
   }

// for 1<=i,j<=5 we compute the common factors of degree 1 of the resultants, with respect to lmf, of the polynomials poly1, poly4, and of poly1 and  poly3, evaluated for lm2=L[i], lm2f=L[j].

    
  for (int i = 1; i <= 5; i=i+1)
       { 
        for (int j = 1; j<= 5; j=j+1)
         {
            list l=L[i],L[j],commonf(i,j,poly1, poly2, poly3);
            print( l );
            } 
        }

// we get 9 pairs (i,j) such that l[3] is not the empty list:
   
// for each such pair (i,j), we evaluate poly1 and/or  poly 3 in ( x0, lmf , l[I], l[j] ), where x0 is one of the solutions of the polynomial in l[3]   
// and solve the equation poly1=0 or poly3=0 in the indeterminate  lmf

// case (1,lmf,1,1)
map m1=F, 1, lmf, 1,1; 
factorize(m1(poly1)); // get lmf=1

// case (0,lmf,1,1)
map m2=F, 0, lmf, 1,1; 
factorize(m2(poly3)); // get lmf=1

// case (al/2,lmf,1,1)
map m3=F, al/2, lmf, 1,1; 
factorize(m3(poly1));  // get lmf=al/2

// case (0,lmf,1,0)
map m4=F, 0, lmf, 1,0; 
factorize(m4(poly3));  // get no lmf in the field

// case (0,lmf,1,bt/2)
map m5=F, 0, lmf, 1,bt/2; 
factorize(m5(poly3));  // get no lmf in the field

// case (0,lmf,1,al/2)
map m6=F, 0, lmf, 1,al/2; 
factorize(m6(poly3));  // get no lmf in the field

// case (0,lmf,1,(3*al^2+3*al*bt-al-2*bt)/(8*al-4))
map m7=F, 0, lmf, 1,(3*al^2+3*al*bt-al-2*bt)/(8*al-4); 
factorize(m7(poly3));  // get no lmf in the field

// case (bt/2,lmf,bt/2,bt/2)
map m8=F, bt/2, lmf, bt/2,bt/2; 
factorize(m8(poly1));  // get lmf=bt/2

// case ((3*al^2+3*al*bt-al-2*bt)/(8*al-4),lmf,(3*al^2+3*al*bt-al-2*bt)/(8*al-4),(3*al^2+3*al*bt-al-2*bt)/(8*al-4))
map m9=F, (3*al^2+3*al*bt-al-2*bt)/(8*al-4), lmf, (3*al^2+3*al*bt-al-2*bt)/(8*al-4),(3*al^2+3*al*bt-al-2*bt)/(8*al-4); 
factorize(m9(poly1)); // get lmf=(3*al^2+3*al*bt-al-2*bt)/(8*al-4)

/////////////////////////////////////////////////////////
// proof of Theorem 1.3

// value returns the image of the vector v when we set lm=p1, lmf=p2,  lm2=p3, lm2f=p4

proc value(poly p1, poly p2, poly p3, poly p4, vector v)
 { 
 poly v1=v[1];
 poly v2=v[2];
 poly v3=v[3];
 poly v4=v[4];
 poly v5=v[5];
 poly v6=v[6];
 poly v7=v[7];
 poly v8=v[8];
     map m=F, p1, p2, p3, p4;
     return( m(v1)*am2+m(v2)*am1+m(v3)*a0+m(v4)*a1+m(v5)*a2+m(v6)*s1+m(v7)*s2+m(v8)*s2f);
  }

// 
// Determination of the algebras U_P for the different choices of P
// case P=(bt/2, bt/2, bt/2, bt/2)

value(bt/2, bt/2, bt/2, bt/2, doublediff0); // the coefficient of am2 is <>0

vector Am2=-(value(bt/2, bt/2, bt/2, bt/2, doublediff0)-value(bt/2, bt/2, bt/2, bt/2, doublediff0)[1]*am2)/value(bt/2, bt/2, bt/2, bt/2, doublediff0)[1]; // we get Am2=a2+a1-am1

proc update(vector v)
  { 
  return( v-v[1]*am2+v[1]*Am2);
  }

vector S2f=-(value(bt/2, bt/2, bt/2, bt/2, update(diff0))-value(bt/2, bt/2, bt/2, bt/2, update(diff0))[8]*s2f)/value(bt/2, bt/2, bt/2, bt/2, update(diff0))[8]; 

proc update(vector v)
  { 
  return( v-v[1]*am2+v[1]*Am2-v[8]*s2f+v[8]*S2f);
  }
  
vector A2=-(value(bt/2, bt/2, bt/2, bt/2, update(diff2))-value(bt/2, bt/2, bt/2, bt/2, update(diff2))[5]*a2)/value(bt/2, bt/2, bt/2, bt/2, update(diff2))[5]; // A2=am1 whence Am2=a1 and s2=s1, whence s2f=s2

proc update(vector v)
  { 
  return( v-v[1]*am2+v[1]*a1-v[8]*s2f+v[8]*s1-v[5]*a2+v[5]*am1-v[7]*s2+v[7]*s1);
  }

// the following vector C must be 0  
vector C=update(value(bt/2, bt/2, bt/2, bt/2,update(Times(s2,s1)-Times(s1,s1))));  
// from here we get s1=-bt/2*(am1+a0+a1)


// Case P=((3*al^2+3*al*bt-al-2*bt)/(8*al-4), (3*al^2+3*al*bt-al-2*bt)/(8*al-4), (3*al^2+3*al*bt-al-2*bt)/(8*al-4), (3*al^2+3*al*bt-al-2*bt)/(8*al-4))

poly nu=(3*al^2+3*al*bt-al-2*bt)/(8*al-4);

value(nu, nu, nu, nu, doublediff0); // the coefficient of am2 is <>0

vector Am2=-(value(nu, nu, nu, nu, doublediff0)-value(nu, nu, nu, nu, doublediff0)[1]*am2)/value(nu, nu, nu, nu, doublediff0)[1]; // we get Am2=a2+a1-am1

proc update(vector v)
  { 
  return( v-v[1]*am2+v[1]*Am2);
  }

vector S2f=-(value(nu, nu, nu, nu, update(diff0))-value(nu, nu, nu, nu, update(diff0))[8]*s2f)/value(nu, nu, nu, nu, update(diff0))[8]; 

proc update(vector v)
  { 
  return( v-v[1]*am2+v[1]*Am2-v[8]*s2f+v[8]*S2f);
  }
  
vector A2=-(value(nu, nu, nu, nu, update(diff2))-value(nu, nu, nu, nu, update(diff2))[5]*a2)/value(nu, nu, nu, nu, update(diff2))[5]; // A2=am1

// case P=(al/2, al/2, 1,1)

value(al/2, al/2, 1, 1, doublediff0); // the coefficient of am2 is <>0

vector Am2=-(value(al/2, al/2, 1, 1, doublediff0)-value(al/2, al/2, 1, 1, doublediff0)[1]*am2)/value(al/2, al/2, 1, 1, doublediff0)[1]; // 

proc update(vector v)
  { 
  return( v-v[1]*am2+v[1]*Am2);
  }

vector S2f=-(value(al/2, al/2, 1,1, update(diff0))-value(al/2, al/2, 1, 1, update(diff0))[8]*s2f)/value(al/2, al/2, 1, 1, update(diff0))[8]; 

proc update(vector v)
  { 
  return( v-v[1]*am2+v[1]*Am2-v[8]*s2f+v[8]*S2f);
  }
 
value(al/2, al/2, 1,1, update(doublediff2))/value(al/2, al/2, 1,1, update(doublediff2))[2];  // am1-a1
value(al/2, al/2, 1,1, update(diff2))/value(al/2, al/2, 1,1, update(diff2))[2];  // am1+a0-a1-a2

// case P=(0, 0, 1,1)

value(0, 0, 1, 1, doublediff0); // the coefficient of am2 is <>0

vector Am2=-(value(0, 0, 1, 1, doublediff0)-value(0, 0, 1, 1, doublediff0)[1]*am2)/value(0, 0, 1, 1, doublediff0)[1]; // 

proc update(vector v)
  { 
  return( v-v[1]*am2+v[1]*Am2);
  }

vector S2f=-(value(0,0, 1,1, update(diff0))-value(0,0,1,1, update(diff0))[8]*s2f)/value(0, 0, 1, 1, update(diff0))[8]; 

proc update(vector v)
  { 
  return( v-v[1]*am2+v[1]*Am2-v[8]*s2f+v[8]*S2f);
  }
 
value(0,0, 1,1, update(doublediff2))/value(0,0, 1,1, update(doublediff2))[2];  // am1-a1
value(0,0, 1,1, update(diff2))/value(0,0, 1,1, update(diff2))[2];  // am1+a0-a1-a2

proc update(vector v)
  { 
  return( v-v[1]*am2+v[1]*a0-v[2]*am1+v[2]*a1-v[5]*a2+v[5]*a0-v[7]*s2+v[7]*(1-2*bt)*a0-v[8]*s2f+v[8]*(1-2*bt)*a1);
  }

vector C1=value(0,0,1,1, Times(a0, (1-2*bt)*a1))-value(0,0,1,1,update(Times(a0, s2f)));
vector S1=-(C1-C1[6]*s1)/C1[6];

// case P=(1,1, 1,1)

value(1,1, 1, 1, doublediff0); // the coefficient of am2 is <>0

vector Am2=-(value(1,1, 1, 1, doublediff0)-value(1,1, 1, 1, doublediff0)[1]*am2)/value(1,1, 1, 1, doublediff0)[1]; // 

proc update(vector v)
  { 
  return( v-v[1]*am2+v[1]*Am2);
  }

vector S2f=-(value(1,1, 1,1, update(diff0))-value(1,1,1,1, update(diff0))[8]*s2f)/value(1,1, 1, 1, update(diff0))[8]; 

proc update(vector v)
  { 
  return( v-v[1]*am2+v[1]*Am2-v[8]*s2f+v[8]*S2f);
  }

// vector d2 
value(1,1, 1,1, update(doublediff2))/value(1,1, 1,1, update(doublediff2))[2];  // am1-a1, hence s2f=(1-2*bt)*a1 and s2=(1-2*bt)*a0

proc update(vector v)
  { 
  return( v-v[1]*am2+v[1]*Am2-v[8]*s2f-v[8]*(1-2*bt)*a1-v[2]*am1+v[2]*a1-v[7]*s2+v[7]*(1-2*bt)*a0);
  }
  
update(Am2); // a2  

proc update(vector v)
  { 
  return( v-v[1]*am2+v[1]*a2-v[8]*s2f+v[8]*(1-2*bt)*a1-v[2]*am1+v[2]*a1-v[7]*s2+v[7]*(1-2*bt)*a0);
  }
  
  vector C1=value(1,1,1,1, Times(a0, (1-2*bt)*a1))-value(1,1,1,1,update(Times(a0, s2f)));
  vector S1=-(C1-C1[6]*s1)/C1[6];

// compute diff2-diff2^f

(value(1,1, 1,1, update(diff2))-update(tau1(value(1,1, 1,1, update(diff2)))))/(value(1,1, 1,1, update(diff2))-update(tau1(value(1,1, 1,1, update(diff2)))))[5]; //a2-a0

proc update(vector v)
  { 
  return( v-v[1]*am2+v[1]*a0-v[8]*s2f+v[8]*(1-2*bt)*a1-v[2]*am1+v[2]*a1-v[7]*s2+v[7]*(1-2*bt)*a0-v[5]*a2+v[5]*a0);
  }

// vector d2  
value(1,1, 1,1, update(diff2))/value(1,1, 1,1, update(diff2))[3]; // a0-a1
