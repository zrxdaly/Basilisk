#include<stdio.h>

int lower(int b);

int main(){
    int a, d;
    float b;
    char c;
    a = 37;
    b = 1.11;
    c = 'G';
    d = 21;
    printf("%d\n", a & d);
    printf("%d %f %f\n", a, b, a+b);
    printf("%c\n", lower(c));
    c = a;
    a = c;
    printf("%c %c\n", a, c);
}

int lower(int c){
    // if (c>= 'A' && c<='Z'){
    //     return c+'a'-'A';
    // }
    // else{
    //     return c;
    // }
    return (c>'A'&& c<='Z')? (c+'a'-'A'):(c);
}