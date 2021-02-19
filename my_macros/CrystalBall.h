double crystalBallLow(double* x, double* par)
{
  //[0] = N
  //[1] = mean
  //[2] = sigma
  //[3] = alpha
  //[4] = n
  
  double xx = x[0];
  double mean = par[1];
  double sigma = par[2];
  double alpha = par[3];
  double n = par[4];
  
  if( (xx-mean)/sigma <= -1.*fabs(alpha) )  
  {
    double A = pow(n/fabs(alpha), n) * exp(-0.5 * alpha*alpha);
    double B = n/fabs(alpha) - fabs(alpha);
    
    return par[0] * A * pow(B - (xx-mean)/sigma, -1.*n);
  }
  
  else
  {
    return par[0] * exp(-1. * (xx-mean)*(xx-mean) / (2*sigma*sigma) );
  } 
  
}


double crystalBallHigh(double* x, double* par)
{
  //[0] = N
  //[1] = mean
  //[2] = sigma
  //[3] = alpha
  //[4] = n
  
  double xx = x[0];
  double mean = par[1];
  double sigma = par[2];
  double alpha = par[3];
  double n = par[4];
  
  if( (xx-mean)/sigma >= 1.*fabs(alpha) )  
  {
    double A = pow(n/fabs(alpha), n) * exp(-0.5 * alpha*alpha);
    double B = n/fabs(alpha) - fabs(alpha);
    
    return par[0] * A * pow(B + (xx-mean)/sigma, -1.*n);
  }
  
  else
  {
    return par[0] * exp(-1. * (xx-mean)*(xx-mean) / (2*sigma*sigma) );
  } 
  
}



double crystalBallDouble(double* x, double* par)
{
  //[0] = N
  //[1] = mean
  //[2] = sigma
  //[3] = alpha
  //[4] = n
  
  double xx = x[0];
  double mean = par[1];
  double sigma = par[2];
  double alpha1 = par[3];
  double alpha2 = par[4];
  double n1 = par[5];
  double n2 = par[6];
  

  if( (xx-mean)/sigma <= -1.*fabs(alpha1) )  
  {
    double A = pow(n1/fabs(alpha1), n1) * exp(-0.5 * alpha1*alpha1);
    double B = n1/fabs(alpha1) - fabs(alpha1);
    
    return par[0] * A * pow(B - (xx-mean)/sigma, -1.*n1);
  }

  else if( (xx-mean)/sigma >= 1.*fabs(alpha2) )  
  {
    double A = pow(n2/fabs(alpha2), n1) * exp(-0.5 * alpha2*alpha2);
    double B = n2/fabs(alpha2) - fabs(alpha2);
    
    return par[0] * A * pow(B + (xx-mean)/sigma, -1.*n2);
  }
  
  else
  {
    return par[0] * exp(-1. * (xx-mean)*(xx-mean) / (2*sigma*sigma) );
  } 
  
}


double fullBananaPlusLin(double* x, double* par)
{
  

  double xx = x[0];
  
  double slopeOff 	= par[0];
  double slopeCoeff 	= par[1];
  double slopeQuadCoeff = par[2];
  
  double low_limit  = par[3];
  double up_limit   = par[4];
  
  double linCoeff   = par[5];
  double linOffset  = 0;
  double risingSlope  = par[6];
  double risingSlopeQuad  = par[7];
  double risingSlopeCube  = 0;
  /*
  if (xx <= low_limit || xx >= up_limit )
  {
       TF1::RejectPoint();
       return 0;
  }*/
 
  if (xx <= low_limit)
  {
    return (xx*linCoeff + linOffset);
  }

  else if (xx >= up_limit)
  {
    return (xx*linCoeff + linOffset);
  }
  else
  {
//     return (par[0]/xx + par[1]/xx/xx);
//     return (par[0]);// + par[1] / x);
    return (slopeQuadCoeff/xx/xx + slopeCoeff/xx + slopeOff + risingSlopeCube*xx*xx*xx + risingSlopeQuad*xx*xx + risingSlope*xx);
  }
  
}

