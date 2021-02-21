#include<stdio.h>

int main(void){
   float salary;
   printf("enter your desired monthly salary");
   printf("$________\b\b\b\b\b\b\b");

   scanf("%f", &salary);
   printf("\n \t$ %.2f a month is $%.2f one year.", salary, salary*12.0);

   printf("\r Gee!\n");
   /* the \b is used to back space*/
   return 0;


}
