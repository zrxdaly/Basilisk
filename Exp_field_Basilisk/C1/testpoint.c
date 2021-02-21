#include<stdio.h>

int main(){
    
    int x = 1, y = 2, z[10];

    int *p;

    p = &x;

    printf("%d %d %d %d\n", p, y = *p, &p, *p);

}