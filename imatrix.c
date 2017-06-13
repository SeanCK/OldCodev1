/* > c.imatrix */ 


#include <stdlib.h>



int **imatrix(int nrh,int nch) { 

   int i; 

   int **m; 

 

   m=(int **) malloc((unsigned) (nrh+1)*sizeof(int*)); 



 

   for(i=0;i<=nrh;i++) { 

      m[i]=(int *) malloc((unsigned) (nch+1)*sizeof(int)); 


   } 

   return m; 

} 





 

void free_imatrix(int **m,int nrh,int nch) { 

   int i; 

 

   for(i=nrh;i>=0;i--) free((char*) (m[i])); 

   free((char*) (m)); 

} 
