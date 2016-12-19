%Compare 3 titres model M1(Titre.A) M1.5(Titre.A with different pre-existing titre) Fix(Titre.D)

X = [1 2 3 4; 1 2 3 4; 1 2 3 4];

A = [35.5 19.6 18.3 7.1;
     37.4 20.4 18.5 5.2;
     12.9 17.9 16.2 7.9
];

B = [2 2 2 4;
     2 2 2 2;
     2 2 2 4
];

C = A - B;

L = [ 9.8  4.9 5.4 1.6;
     10.2  5.3 5.2 1.4;
      2.9  3.8 4.4 2.1
];

U = [ 7.2 4.3 4.6 1.4;
      6.7 3.9 4.4 1.2;
      3.0 4.1 4.5 1.5
];
errorbar(X',C',L',U')
UB