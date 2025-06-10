//
// Universal 2-generated axial algebra of Monster type (4*bt, bt), bt<>1/2.
//  Computations of the formulas in Section 7 of the paper "2-generated axial algebras of Monster type" by 
// Clara Franchi, Mario Mainardis, Sergey Shpectorov
// 
//
// (c) 2024 C. Franchi, M. Mainardis, S. Shpectorov 
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Construction of the algebra
//
ring F=(0, bt), (lm, lmf, lm2, lm2f, lm3, lm3f, lm4, lm4f), dp;
freemodule(18);
vector am4=gen(1);  
vector am3=gen(2);  
vector am2=gen(3);  
vector am1=gen(4);  
vector a0=gen(5);   
vector a1=gen(6);   
vector a2=gen(7);   
vector a3=gen(8);   
vector a4=gen(9);  
vector s1=gen(10); 
vector s2=gen(11);  
vector s2f=gen(12);  
vector s3=gen(13);  
vector s3f=gen(14);  
vector s3t=gen(15);  
vector s4=gen(16);   
vector s42=gen(17);   
vector pp=gen(18);  // extra generator needed to write the function Time properly
poly al=4*bt;

list a=am4, am3, am2,am1,a0,a1,a2,a3,a4, s1,s2,s2f, s3, s3f, s3t, s4, s42, pp;



// Construct the multiplication table
 
list r=pp,pp,pp,pp,pp,  pp,pp,pp,pp,pp,  pp,pp,pp,pp,pp,  pp,pp,pp;
list Mult=r,r,r,r,r, r,r,r,r,r, r,r,r,r,r,  r,r,r ;

// Known products

Mult[1][1]=am4; Mult[2][2]=am3; Mult[3][3]=am2; Mult[4][4]=am1; Mult[5][5]=a0; Mult[6][6]=a1; Mult[7][7]=a2; Mult[8][8]=a3; Mult[9][9]=a4;

Mult[1][2]=s1+bt*(am4+am3);     Mult[2][1]=Mult[1][2];
Mult[2][3]=s1+bt*(am3+am2);     Mult[3][2]=Mult[2][3];
Mult[3][4]=s1+bt*(am2+am1);     Mult[4][3]=Mult[3][4];
Mult[4][5]=s1+bt*(am1+a0);        Mult[5][4]=Mult[4][5];
Mult[5][6]=s1+bt*(a0+a1);           Mult[6][5]=Mult[5][6];
Mult[6][7]=s1+bt*(a1+a2);           Mult[7][6]=Mult[6][7];
Mult[7][8]=s1+bt*(a2+a3);           Mult[8][7]=Mult[7][8];
Mult[8][9]=s1+bt*(a3+a4);           Mult[9][8]=Mult[8][9];

Mult[1][3]=s2+bt*(am4+am2);     Mult[3][1]=Mult[1][3];
Mult[2][4]=s2f+bt*(am3+am1);    Mult[4][2]=Mult[2][4];
Mult[3][5]=s2+bt*(am2+a0);        Mult[5][3]=Mult[3][5];
Mult[4][6]=s2f+bt*(am1+a1);       Mult[6][4]=Mult[4][6];
Mult[5][7]=s2+bt*(a0+a2);           Mult[7][5]=Mult[5][7];
Mult[6][8]=s2f+bt*(a1+a3);          Mult[8][6]=Mult[6][8];
Mult[7][9]=s2+bt*(a2+a4);           Mult[9][7]=Mult[7][9];

Mult[1][4]=s3t+bt*(am4+am1);      Mult[4][1]=Mult[1][4];
Mult[2][5]=s3+bt*(am3+a0);          Mult[5][2]=Mult[2][5];
Mult[3][6]=s3f+bt*(am2+a1);         Mult[6][3]=Mult[3][6];
Mult[4][7]=s3t+bt*(am1+a2);         Mult[7][4]=Mult[4][7];
Mult[5][8]=s3+bt*(a0+a3);             Mult[8][5]=Mult[5][8];
Mult[6][9]=s3f+bt*(a1+a4);            Mult[9][6]=Mult[6][9];

Mult[1][5]=s4+bt*(am4+a0);        Mult[5][1]=Mult[1][5];
Mult[3][7]=s42+bt*(am2+a2);       Mult[7][3]=Mult[3][7];
Mult[5][9]=s4+bt*(a0+a4);         Mult[9][5]=Mult[5][9];

Mult[5][10]=(al-bt)*s1+(lm*(1-al)+bt*(al-bt-1))*a0+
  (1/2)*(al-bt)*bt*(am1+a1);                                            Mult[10][5]=Mult[5][10];
Mult[6][10]=(al-bt)*s1+(lmf*(1-al)+bt*(al-bt-1))*a1+
  (1/2)*(al-bt)*bt*(a0+a2);                                               Mult[10][6]=Mult[6][10]; 
Mult[4][10]=(al-bt)*s1+(lmf*(1-al)+bt*(al-bt-1))*am1+
  (1/2)*(al-bt)*bt*(a0+am2);                                            Mult[10][4]=Mult[4][10];    
Mult[7][10]=(al-bt)*s1+(lm*(1-al)+bt*(al-bt-1))*a2+
  (1/2)*(al-bt)*bt*(a1+a3);                                               Mult[10][7]=Mult[7][10];  
Mult[3][10]=(al-bt)*s1+(lm*(1-al)+bt*(al-bt-1))*am2+
  (1/2)*(al-bt)*bt*(am1+am3);                                         Mult[10][3]=Mult[3][10];    
  
Mult[5][11]=(al-bt)*s2+(lm2*(1-al)+bt*(al-bt-1))*a0+
  (1/2)*(al-bt)*bt*(am2+a2);                                            Mult[11][5]=Mult[5][11];
Mult[6][12]=(al-bt)*s2f+(lm2f*(1-al)+bt*(al-bt-1))*a1+
  (1/2)*(al-bt)*bt*(a3+am1);                                            Mult[12][6]=Mult[6][12];  
Mult[4][12]=(al-bt)*s2f+(lm2f*(1-al)+bt*(al-bt-1))*am1+
  (1/2)*(al-bt)*bt*(am3+a1);                                            Mult[12][4]=Mult[4][12];  
Mult[7][11]=(al-bt)*s2+(lm2*(1-al)+bt*(al-bt-1))*a2+
  (1/2)*(al-bt)*bt*(a4+a0);                                               Mult[11][7]=Mult[7][11];  
Mult[3][11]=(al-bt)*s2+(lm2*(1-al)+bt*(al-bt-1))*am2+
  (1/2)*(al-bt)*bt*(am4+a0);                                               Mult[11][3]=Mult[3][11]; 


Mult[5][13]=(al-bt)*s3+(lm3*(1-al)+bt*(al-bt-1))*a0+
  (1/2)*(al-bt)*bt*(am3+a3);                                           Mult[13][5]=Mult[5][13];

Mult[5][16]=(al-bt)*s4+(lm4*(1-al)+bt*(al-bt-1))*a0+
  (1/2)*(al-bt)*bt*(am4+a4);                                           Mult[16][5]=Mult[5][16];   


// Times returns the product of the vectors u and v if computable. Otherwise it returns a list where the first element is the part of the product of u and v which can be calculated and the second element is the list of the pairs of base elements whose product is still unknown and which is needed to obtain the product of u and v, each one with the right coefficient 
 
proc Times (vector u, vector v)
{
vector ans=[0];
list rem=list();
 for(int i=1; i<=18; i=i+1)
  { if (u[i]<>0)
     {
        for(int j=1; j<=18; j=j+1)
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

// Symmetries: tau0 is the Miyamoto involution associated to the axis a0
// tau1 is the Miyamoto involution associated to the axis a1

proc tau0 (vector v)
  {
 return(v[9]*am4+v[8]*am3+v[7]*am2+v[6]*am1+v[5]*a0+v[4]*a1+v[3]*a2+v[2]*a3+v[1]*a4+v[10]*s1+v[11]*s2+v[12]*s2f+v[13]*s3+v[14]*s3t+v[15]*s3f+v[16]*s4+v[17]*s42);        
}

proc tau1 (vector v)
  {
 if (v[1]==0 and v[2]==0)
{
return(v[3]*a4+v[4]*a3+v[5]*a2+v[6]*a1+v[7]*a0+v[8]*am1+v[9]*am2+v[10]*s1+v[11]*s2+v[12]*s2f+v[13]*s3t+v[14]*s3f+v[15]*s3+v[16]*s42+v[17]*s4);        
  } 
  else {return ("fail");}
  } 


// ff is an endomorphism of F swapping lm with lmf, lm2 with lm2f, lm3 with lm3f, lm4 with lm4f 

map ff=F, lmf, lm, lm2f, lm2, lm3f, lm3, lm4f, lm4; 

// f is the semiautomorphism of the algebra swapping a0 and a1

 proc f (vector v)
  {
poly v1=v[1];
poly v4=v[4];
poly v3=v[3];
poly v2=v[2];
poly v5=v[5];
poly v6=v[6];
poly v7=v[7];
poly v8=v[8];
poly v9=v[9];
poly v10=v[10];
poly v11=v[11];
poly v12=v[12];
poly v13=v[13];
poly v14=v[14];
poly v15=v[15];
poly v16=v[16];
poly v17=v[17];

 if (v[1]==0 and v[16]==0 and v[17]==0)
{
return(ff(v9)*am3+ff(v8)*am2+ff(v7)*am1+ff(v6)*a0+ff(v5)*a1+ff(v4)*a2+ff(v3)*a3+ff(v2)*a4+ff(v10)*s1+ff(v11)*s2f+ff(v12)*s2+ff(v13)*s3f+ff(v14)*s3+ff(v15)*s3t);        
  } 
  else {return ("fail");}
  }


// lambda returns the value lambda_{a_0}(v)
          
proc lambda (vector v)
{
 if ( v[12]==0 and v[14]==0 and v[15]==0 and v[17]==0) 
  { poly p= lm4*v[1]+lm3*v[2]+lm2*v[3]+lm*v[4]+v[5]+lm*v[6]+lm2*v[7]+lm3*v[8]+lm4*v[9]+(lm-bt-bt*lm)*v[10]+(lm2-bt-bt*lm2)*v[11]+(lm3-bt-bt*lm3)*v[13]+(lm4-bt-bt*lm4)*v[16];
  return(p); }
  else {return ("fail")}
}  

 
vector v1= (s1+(bt-lm)*a0+(bt/2)*(a1+am1))/al; // al-eigenvalue for a0
vector u1=a1-lm*a0-v1-(a1-am1)/2;              // 0-eigenvalue for a0

vector uu=Times(u1, u1)[1];  // contains  1/(16*bt^2)*[10][10]
vector uv=Times(u1, v1)[1];  // contains  -1/(16*bt^2)*[10][10]
vector vv=Times(v1, v1)[1];  // contains  1/(16*bt^2)*[10][10]

// (uu+uv) is a 0-eigenvector for a0 and contains s2f. Since lambda(uu+uv)=0, we deduce lambda(s2f)

vector w=(uu+uv)-(uu+uv)[12]*s2f;
 
poly ls2f=-lambda(w)/((uu+uv)[12]);

// update the function lambda

proc lambda (vector v)
{
if ( v[14]==0 and v[15]==0 and v[17]==0) 
   {
   poly p= lm4*v[1]+lm3*v[2]+lm2*v[3]+lm*v[4]+v[5]+lm*v[6]+lm2*v[7]+lm3*v[8]+lm4*v[9]+(lm-bt-bt*lm)*v[10]+(lm2-bt-bt*lm2)*v[11]+(lm3-bt-bt*lm3)*v[13]+ls2f*v[12]+(lm4-bt-bt*lm4)*v[16];
  return(p); 
  } 
 else
 {return("fail")}
 }

  
// u1*u1 is a 0-eigenvector for a0, hence it has 0-projection on a0, and it contains 1/(al^2)*s1*s1 

poly ls1s1=-(al^2)*lambda(uu);; // lambda(s1*s1)


//
// Express a0*s2f
//
vector w=uu-vv -lambda(uu-vv)*a0; // this is a 0-eigenvector for a0

// Times(a0, w)  this vector must be zero
// it involves a0*s2f with coefficient (al-2*bt)/(2*al), hence we can express a0*s2f

Mult[5][12]=-2*al*Times(a0, w)[1]/(al-2*bt);
Mult[12][5]=Mult[5][12];

Mult[6][11]=f(Mult[5][12]); Mult[11][6]=Mult[6][11];
Mult[7][12]=tau1(Mult[5][12]); Mult[12][7]=Mult[7][12];
Mult[4][11]=tau0(f(Mult[5][12])); Mult[11][4]=Mult[4][11];


//
// Express s1*s1
//

// a0*(u1*u1+u1*v1)=al*u1*v1, where the right side has -s1*s1/al

vector s1s1=-al*(Times(a0,uu+uv)-al*uv);
Mult[10][10]=s1s1;


// THE FOLLOWING COMMANDS GIVE THE FORMULAS IN LEMMAS 7.1, 7.2, 7.3, 7.5


vector rel=s1s1-f(s1s1); // relation in Lemma 7.1

vector S2f=(rel[12]*s2f-rel)/rel[12];  // formula for s_{1,2} in Lemma 7.2

vector Rel=rel-tau0(rel); // first equation in Lemma 7.3
vector Relf=f(Rel);       // second equation in Lemma 7.3

vector rel1=Times(a0, S2f)-Mult[5][12];  // equation in Lemma 7.5


//
// function update substitutes S2f in place of s2f
//
proc update(vector v)
{
return(v-v[12]*s2f+v[12]*S2f);
}



// THE FOLLOWING COMMANDS GIVE THE FORMULA IN LEMMA 7.6(v)

vector v2= (s2+(bt-lm2)*a0+(bt/2)*(a2+am2))/al; // al-eigenvector for a0
vector u2=a2-lm2*a0-v2-(a2-am2)/2;              // 0-eigenvector for a0

Times(u1,u2);

vector u1u2=Times(u1, u2)[1];  // contains  1/(16*bt^2)*[10][11]
vector u1v2=Times(u1, v2)[1];  // contains  -1/(16*bt^2)*[10][11]
vector v1v2=Times(v1, v2)[1];  // contains  1/(16*bt^2)*[10][11]

// (u1u2+u1v2) is a 0-al-eigenvector for a0 and contains 3/16*(s3f+s3t): we compute lambda(s3f+s3t)
// s3f-s3t is a bt-eigenvector for a0, hence lambda(s3f)=lambda(s3t)

vector w=(u1u2+u1v2)-3/16*(s3f+s3t);
 
poly ls3f=-lambda(w)*8/3;


// update lambda

proc lambda (vector v)
{
if ( v[17]==0) 
   {
   poly p= lm4*v[1]+lm3*v[2]+lm2*v[3]+lm*v[4]+v[5]+lm*v[6]+lm2*v[7]+lm3*v[8]+lm4*v[9]+(lm-bt-bt*lm)*v[10]+(lm2-bt-bt*lm2)*v[11]+(lm3-bt-bt*lm3)*v[13]+ls2f*v[12]+ls3f*v[14]+ls3f*v[15]+(lm4-bt-bt*lm4)*v[16];
  return(p); 
  } 
 else
 {return("fail")}
 }



// u1*u2 is a 0-eigenvector for a0, hence it has 0-projection on a0, and it contains 1/(16*bt^2)*s1*s2 

poly ls1s2=-(16*bt^2)*lambda(u1u2);; // lambda(s1*s2)

Times(a0,u1u2+u1v2); // it contains a0*(s3f+s3t)

//
// Express a0*(s3f+s3t)
//

vector w=u1u2-v1v2 -lambda(u1u2-v1v2)*a0; // this is a 0-eigenvector for a0

// Times(a0, w)  this vector must be zero
// it involves a0*(s3f+s3t) with coefficient 1/8, hence we can express a0*(s3f+s3t)

vector a0s3ft=-8*Times(a0, w)[1];

Mult[5][14]=1/2*(bt*(s3f-s3t)+a0s3ft); 
Mult[14][5]=Mult[5][14];

Mult[5][15]=1/2*(-bt*(s3f-s3t)+a0s3ft);
Mult[15][5]=Mult[5][15];



//
// Express s1*s2
//

// a0*(u1*u2+u1*v2)=al*u1*v2, where the product al*u1*v2 contains -1/al*s1*s2

vector s1s2=-al*(Times(a0,u1u2+u1v2)-al*u1v2);// formula in Lemma 7.6(v)
Mult[10][11]=s1s2;
Mult[11][10]=Mult[10][11];

//
// Express s2*s2
//

// a0*(u2*u2+u2*v2)=al*u2*v2, where the product al*u2*v2 contains -1/al*s2*s2

vector u2u2=Times(u2, u2)[1];; // contains [11][11] with coefficient 1/al^2
vector v2u2=Times(v2, u2)[1];; // contains [11][11] with coefficient -1/al^2
vector v2v2=Times(v2, v2)[1];; // contains [11][11] with coefficient 1/al^2

vector p22=u2u2+v2u2;;

//  p22 is the sum of a 0 and an al-eigenvector for a0, hence lambda(p22)=0, p22 contains 3/8*s42

poly lms42=-8/3*lambda(p22-3/8*s42); 

// update lambda

proc lambda (vector v)
   {
   poly p= lm4*v[1]+lm3*v[2]+lm2*v[3]+lm*v[4]+v[5]+lm*v[6]+lm2*v[7]+lm3*v[8]+lm4*v[9]+(lm-bt-bt*lm)*v[10]+(lm2-bt-bt*lm2)*v[11]+(lm3-bt-bt*lm3)*v[13]+ls2f*v[12]+ls3f*v[14]+ls3f*v[15]+(lm4-bt-bt*lm4)*v[16]+lms42*v[17];
  return(p); 
  } 


vector w2=u2u2-v2v2+lambda(u2u2-v2v2)*a0; // this is a 0-eigenvector for a0

// Times(a0, w2) contains 1/4*a0*s42

Mult[5][17]=-4*Times(a0, w2)[1];
Mult[17][5]=Mult[5][17];


Mult[11][11]=-al*(Times(a0,p22)-al*v2u2);


 // function lambdb returns the value of lambda_{a_1}(v)

proc lambdb(vector v)
{
poly l=lambda(f(v));
 return(ff(l)); 
}



//THE FOLLOWING COMMANDS GIVE THE FORMULAS FOR LEMMAS 7.8 AND 7.9

vector rel2=s1s2-tau1(s1s2); // equation in Lemma 7.8
//contains (15*bt^2)/4*s3t

vector rel3=Times(s1, S2f)-f(Times(s1,s2)); // equation in Lemma 7.9


// Assume the characteristic is not 5

// express s3f

vector S3f=-4/(15*bt^2)*(rel3-(15*bt^2)/4*s3f);
vector S3t=tau0(S3f);


proc update(vector v)
{ 
return( v-v[12]*s2f+v[12]*S2f-v[14]*s3f+v[14]*S3f-v[15]*s3t+v[15]*S3f);
}



//THE FOLLOWING COMMAND GIVES THE FORMULA FOR LEMMAS 7.11

vector rel23=rel2-tau0(rel3); // equation in Lemma 7.11

Mult[6][13]=f(Mult[5][14]);  // a1*s3
Mult[13][6]=f(Mult[5][14]);
Mult[6][14]=f(Mult[5][13]);


// for every i in {1, ..., 11}, the function Formula returns the formulas in Lemma 7.i. For i=4, it gives the products a0*s2f and a1*s2; for i=7, it gives the expression of s2*s2

proc Formula(int n)
{
list k=Rel,Relf;
list k1=Times(a0,S2f), Times(a1,s2);
list k2=S3f, S3t;
list h= rel,S2f,k, k1, rel1,s1s2, Times(s2,s2), rel2,rel3, k2,rel23;
 
return(h[n]);
}



//CHECK OF CLAIM 2 IN THE PROOF OF THEOREM 7.12

map val=F, (18*bt-1)/8, lmf, lm2, lm2f, lm3, lm3f, lm4, lm4f;

poly pRelf5=Relf[5];
val(pRelf5); // (40*bt^2-14*bt+1)/4*lmf+(-3*bt^2)*lm2+(-100*bt^2+20*bt-1)/64



//CHECK OF CLAIM 3 IN THE PROOF OF THEOREM 7.12: CASE lm2f=P(lm)

map val=F, lm, (18*bt-1)/8, lm2, (40*bt^2-14*bt+1)/(12*bt^2)*lm-(100*bt^2-20*bt+1)/(192*bt^2), lm3, lm3f, lm4, lm4f;

proc Val(vector v)
{ 
poly v1=v[1];
poly v2=v[2];
poly v3=v[3];
poly v4=v[4];
poly v5=v[5];
poly v6=v[6];
poly v7=v[7];
poly v8=v[8];
poly v9=v[9];
poly v10=v[10];
poly v11=v[11];
poly v12=v[12];
poly v13=v[13];
poly v14=v[14];
poly v15=v[15];
poly v16=v[16];
poly v17=v[17];

return( val(v1)*am4+val(v2)*am3+val(v3)*am2+val(v4)*am1+val(v5)*a0+val(v6)*a1+val(v7)*a2+val(v8)*a3+val(v9)*a4+val(v10)*s1+val(v11)*s2+val(v12)*s2f+val(v13)*s3+val(v14)*s3f+val(v15)*s3t+val(v16)*s4+val(v17)*s42);
}


Val(rel1); //  Equation (62)


// if lm<>1/16, then Val(rel1)[6]<>0 and we get dim<4
// Assume lm==1/16


map val1=F, 1/16, (18*bt-1)/8, lm2, (40*bt^2-14*bt+1)/(12*bt^2)*1/16-(100*bt^2-20*bt+1)/(192*bt^2), lm3, lm3f, lm4, lm4f;

proc Val1(vector v)
{ 
poly v1=v[1];
poly v2=v[2];
poly v3=v[3];
poly v4=v[4];
poly v5=v[5];
poly v6=v[6];
poly v7=v[7];
poly v8=v[8];
poly v9=v[9];
poly v10=v[10];
poly v11=v[11];
poly v12=v[12];
poly v13=v[13];
poly v14=v[14];
poly v15=v[15];
poly v16=v[16];
poly v17=v[17];

return( val1(v1)*am4+val1(v2)*am3+val1(v3)*am2+val1(v4)*am1+val1(v5)*a0+val1(v6)*a1+val1(v7)*a2+val1(v8)*a3+val1(v9)*a4+val1(v10)*s1+val1(v11)*s2+val1(v12)*s2f+val1(v13)*s3+val1(v14)*s3f+val1(v15)*s3t+val1(v16)*s4+val1(v17)*s42);
}

Val1(Val1(rel1)); // [(-6*bt+1)/3*lm2+(1056*bt^3-344*bt^2+33*bt-1)/(384*bt^2)]*a0 

poly Lm2=(1056*bt^3-344*bt^2+33*bt-1)/(384*bt^2)*3/(6*bt-1); // Equation (64)

// using this value of lm2 we get the following


map val2=F, 1/16, (18*bt-1)/8, (1056*bt^3-344*bt^2+33*bt-1)/(384*bt^2)*3/(6*bt-1), (40*bt^2-14*bt+1)/(12*bt^2)*1/16-(100*bt^2-20*bt+1)/(192*bt^2), lm3, lm3f, lm4, lm4f;

proc Val2(vector v)
{ 
poly v1=v[1];
poly v2=v[2];
poly v3=v[3];
poly v4=v[4];
poly v5=v[5];
poly v6=v[6];
poly v7=v[7];
poly v8=v[8];
poly v9=v[9];
poly v10=v[10];
poly v11=v[11];
poly v12=v[12];
poly v13=v[13];
poly v14=v[14];
poly v15=v[15];
poly v16=v[16];
poly v17=v[17];

return( val2(v1)*am4+val2(v2)*am3+val2(v3)*am2+val2(v4)*am1+val2(v5)*a0+val2(v6)*a1+val2(v7)*a2+val2(v8)*a3+val2(v9)*a4+val2(v10)*s1+val2(v11)*s2+val2(v12)*s2f+val2(v13)*s3+val2(v14)*s3f+val2(v15)*s3t+val2(v16)*s4+val2(v17)*s42);
}

Val2(rel2)[9]; // (96*bt^3-3*bt^2)/16 =3*bt^2*(32*bt-1)/16

Val2(rel3)[3];  // (-24*bt^3-3*bt^2)/16 = -3*bt^2*(8*bt+1)/16 = -15/65536 if bt=1/32



//CHECK OF CLAIM 5 IN THE PROOF OF THEOREM 7.12
 

map val3=F, (18*bt-1)/8, (18*bt-1)/8, lm2, lm2f, lm3, lm3f, lm4, lm4f;

proc Val3(vector v)
{ 
poly v1=v[1];
poly v2=v[2];
poly v3=v[3];
poly v4=v[4];
poly v5=v[5];
poly v6=v[6];
poly v7=v[7];
poly v8=v[8];
poly v9=v[9];
poly v10=v[10];
poly v11=v[11];
poly v12=v[12];
poly v13=v[13];
poly v14=v[14];
poly v15=v[15];
poly v16=v[16];
poly v17=v[17];
 return( val3(v1)*am4+val3(v2)*am3+val3(v3)*am2+val3(v4)*am1+val3(v5)*a0+val3(v6)*a1+val3(v7)*a2+val3(v8)*a3+val3(v9)*a4+val3(v10)*s1+val3(v11)*s2+val3(v12)*s2f+val3(v13)*s3+val3(v14)*s3f+val3(v15)*s3t+val3(v16)*s4+val3(v17)*s42);
}

poly P=(40*bt^2-14*bt+1)/(12*bt^2)*lm-(100*bt^2-20*bt+1)/(192*bt^2);
poly om=val3(P);


map val4=F, (18*bt-1)/8, (18*bt-1)/8, om, om, lm3, lm3f, lm4, lm4f;

proc Val4(vector v)
{ 
poly v1=v[1];
poly v2=v[2];
poly v3=v[3];
poly v4=v[4];
poly v5=v[5];
poly v6=v[6];
poly v7=v[7];
poly v8=v[8];
poly v9=v[9];
poly v10=v[10];
poly v11=v[11];
poly v12=v[12];
poly v13=v[13];
poly v14=v[14];
poly v15=v[15];
poly v16=v[16];
poly v17=v[17];
 return( val4(v1)*am4+val4(v2)*am3+val4(v3)*am2+val4(v4)*am1+val4(v5)*a0+val4(v6)*a1+val4(v7)*a2+val4(v8)*a3+val4(v9)*a4+val4(v10)*s1+val4(v11)*s2+val4(v12)*s2f+val4(v13)*s3+val4(v14)*s3f+val4(v15)*s3t+val4(v16)*s4+val4(v17)*s42);
}

Val4(rel1);  

Val4(f(rel1));


// assume now bt=1/12

poly bt=1/12;


Val4(S2f); // S2f=s2

poly L=lambda(rel3);
val4(L);  // -1/768*lm3+1/768*lm3f => lm3=lm3f


map val5=F, (18*bt-1)/8, (18*bt-1)/8, om, om, lm3, lm3, lm4, lm4f;

proc Val5(vector v)
{ 
poly v1=v[1];
poly v2=v[2];
poly v3=v[3];
poly v4=v[4];
poly v5=v[5];
poly v6=v[6];
poly v7=v[7];
poly v8=v[8];
poly v9=v[9];
poly v10=v[10];
poly v11=v[11];
poly v12=v[12];
poly v13=v[13];
poly v14=v[14];
poly v15=v[15];
poly v16=v[16];
poly v17=v[17];
 return( val5(v1)*am4+val5(v2)*am3+val5(v3)*am2+val5(v4)*am1+val5(v5)*a0+val5(v6)*a1+val5(v7)*a2+val5(v8)*a3+val5(v9)*a4+val5(v10)*s1+val5(v11)*s2+val5(v12)*s2f+val5(v13)*s3+val5(v14)*s3f+val5(v15)*s3t+val5(v16)*s4+val5(v17)*s42);
}

2304/5*Val5(update(rel23));

vector A4=8/5*(6*lm3-1)*(a2-am1)-a3+am2+am3;
vector Am4=tau0(A4);
vector A5=am2+am1-8/5*(6*lm3-1)*(a0-a3)-A4;

Val5(S3f); // Equation (67)

Val5(S3t); // Equation (68)

vector relnew=update(Times(a1, S3f)-Times(a1,s3f)); // 

Val5(update(Val5(relnew))); // Equation (69)

12*Val5(update(Val5(relnew-tau0(relnew))));  // we can express am3 as a liner combination of am2, am1, a1, a2,a3


// for every i in {62, 66, 67, 68}, the function Eq returns the second member of Equation (i) in the proof of Theorem 7.12. 

proc Eq(int n)
{
 if (n==62)
{ return(Val(rel1));}
 if (n==66)
{ return(A4);}
 if (n==67)
{ return(Val5(S3f));}
 if (n==68)
{ return(Val5(S3t));}
else
{ return("see the code")}
}









