#include<stdio.h>

int fahr(int fahr);

int main(){
    int fah, temp;
    for (fah=0;fah<301;fah=fah+30){
        temp = fahr(fah);
        printf("%d %d\n", fah, temp);
    }
}

int fahr(int f){
    int temp;
    temp = 5*(f-32)/9;
    return temp;
}