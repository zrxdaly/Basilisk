#include<stdio.h>

int main(void){
   int age;
   double sec = 3.156e7;
   printf("please enter your age:\n");
   scanf("%d",&age);
   printf("You are %d years old, the overall second your life is %e\n", age, age*sec);

   return 0;

}
