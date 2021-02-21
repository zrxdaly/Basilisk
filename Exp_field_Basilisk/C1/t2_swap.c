#include<stdio.h>

void swap(int *px, int *py);


int main(){
    int a = 4,b = 6;
    swap(&a, &b);
    printf("%d %d\n", a, b);
}

void swap(int *px, int *py){
    int temp;

    temp = *px;
    *px = *py;
    *py = temp;

}