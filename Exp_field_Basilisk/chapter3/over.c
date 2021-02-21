#include<stdio.h>

int main(void){
   int initt;
   initt = 2147483648;
   printf("This is the max int value %d, after 1 more value larger %d, if 2 more %d \n", initt, initt+1, initt+2);

   float overflow;
   float overflow1;
   overflow = 0.1234567;
   overflow1 = overflow*100e100;
   printf("The value of overflow is %f, however the value overflow to %f\n", overflow, overflow1);

   return 0;

}
