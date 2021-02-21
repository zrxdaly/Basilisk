#include<stdio.h>

int main(){
   int value;
   printf("please enter a value:\n");
   scanf("%d", &value);

   printf("The value you entered is %d, if it is transfer to char %c.\n", value, value);

   return 0;

}
