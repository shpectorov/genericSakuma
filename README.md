# 2-Generated Axial Algebras of Monster Type

This is a Singular code for performing the computations needed for the paper 

[C.Franchi, M. Mainardis, and S. Shpectorov, 2-generated axial algebras of Monster type, arXiv:2101.10379](https://doi.org/10.48550/arXiv.2101.10379)

# Content

- genericsakuma.s constructing the universal 2-generated primitive axial algebra of Monster type $(\alpha, \beta)$ when $\alpha\not \in$ { $2\beta, 4\beta$ } and producing the four polynomials poly1, poly2, poly3, poly4 (corresponding to the polynomials $p_1,p_2,p_3, p_4$ in the paper). Use this to obtain the formulas in Section 6 of the paper and produce the four polynomials poly1, poly2, poly3, poly4 (corresponding to the polynomials $p_1,p_2,p_3, p_4$ in Theorem 1.2 and Section 8).
- universalalgebra-al=4bt.s computing the formulas in Section 7 of the paper.
- checksSection9.s performing the computations needed in Section 9.
- majorana.s performing the computations needed in Section 10.

# Code-paper dictionary

The following correspondence between the symbols used in the code and those used in the paper hold:
- al corresponds to $\alpha$ and bt corresponds to $\beta$
- lm, lmf, lm2, lm2f, lm3, lm3f, lm4, lm4f correspond to $\lambda_1$, $\lambda_1^f$, $\lambda_2$, $\lambda_2^f$, $\lambda_3$, $\lambda_3^f$, $\lambda_4$, $\lambda_4^f$,  respectively;
- am4, am3, am2, am1, a0, a1, a2, a3, a4 correspond to $a_{-4}$, $a_{-3}$, $a_{-2}$, $a_{-1}$, $a_0$, $a_1$, $a_2$, $a_3$, $a_4$, respectively;
- s1, s2, s3, s4 correspond to $s_{0,1}$, $s_{0,2}$, $s_{0,3}$, $s_{0,4}$, respectively;
- s2f, s3f, s3t, s42 correspond to $s_{1,2}$, $s_{1,3}$, $s_{2,3}$, $s_{2,4}$, respectively;
- poly1, poly2, poly3, poly4 correspond to the polynomials $p_1$, $p_2$, $p_3$, $p_4$,  respectively.
 

# Examples

## genericsakuma.s

To produce the polynomial $p_i$ (resp $q_i$)  run the file genericsakuma.s in Singular and type

    P(i); (resp Q(i);)

Note that Singular returns gen(i) for the ith entry of ( am2, am1, a0, a1, a2, s1, s2, s2f ).   

## universalalgebra-al=4bt.s

To produce the displayed formulas in Lemma 7.i ($i \in$ { $1,..., 11$ }) run the file universalalgebra-al=4bt.s in Singular  and  type

    Formula(i);    

Note that Singular returns gen(i) for the ith entry of ( am4, am3, am2, am1, a0, a1, a2, a3, a4, s1, s2, s2f, s3, s3f, s3t, s4, s42 ). 
For other products use the function Times: for example, to compute $a_0s_{1,3}$ in Lemma 7.6(iii) type

    Times(a0,s3f);

To produce the value of the second member of Equation (i) ($i \in$ { $62, 66, 67, 68$ }) run the file universalalgebra-al=4bt.s in Singular  and  type

    Eq(i);    

The other formulas in Section 7 can be easily computed or found in the code.  

# checksSection9.s

To check the result o Lemma 9.1, run the file genericsakuma.s and type
 
    factorize(resultant(q1,q2,lm2)); 
    
to obtain the irreducible factors of the resultant of $q_1$ and $q_2$ with respect to $\lambda_2$. Considering the linear factors, we get the five values $1$, $0$, $\beta/2$, $\alpha/2$, $(3\alpha^2+3\alpha\beta-\alpha-2\beta)/(8\alpha-4)$ for $\lambda$ in the field.
For each of the previous values we compute the corresponding value of $\lambda_2$ by typing

    map val1=F, 1,1, lm2,lm2;
    factorize(val1(q1)); 
    factorize(val1(q2)); 


    map val2=F, 0,0, lm2,lm2; 
    factorize(val2(q1));  
    factorize(val2(q2));  
