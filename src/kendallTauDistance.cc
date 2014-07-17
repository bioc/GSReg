#include <R.h>
#include "math.h"

#define max(x,y)  ((x)<(y))? (y):(x)
// dist = kendalltaudist(vect);
// dist(i,j)=kendall-tau-dist(vect[,i],vect[,j])
//
extern "C"
{
void kendalltaudist( double *vect, int *np,int *mp, double *dist)
{
    //I/O variables
    //n number of elements in vect and m the number of vectors
    //vect to be processed
    // comp    should be of size dimN^2*dimM
    int n = *np;
    int m = *mp;
    bool zk,zl;
    //Internal variables
    int i,j,k,l;//indices
    for(k=0;k<m;k++)
      for(l=k+1;l<m;l++)
      {
        for(i=0;i<n;i++)
            for(j=i+1;j<n;j++)
            {  
                zk = vect[i+k*n] < vect[j+k*n];
                zl = vect[i+l*n] < vect[j+l*n];
                if( zk != zl)
                {
                    dist[k+l*m]++;
                    dist[l+k*m]++;
                    
                }

            }   
      }
}
}
