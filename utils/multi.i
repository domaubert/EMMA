WARNING freeing builtin function is_scalar
WARNING freeing builtin function is_vector
WARNING freeing builtin function is_matrix
WARNING freeing builtin function is_real
WARNING freeing builtin function is_complex
WARNING freeing builtin function is_integer
WARNING freeing builtin function is_numerical
 Copyright (c) 2005.  The Regents of the University of California.
 All rights reserved.  Yorick 2.2.01 ready.  For help type 'help'
> #include "utils/readamr.i"
> info,oct2cube
 func oct2cube(fname,lvl,field,&time,ncpu=,execut=)
> i=90;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);
Casting Rays on 64x64x64 cube from data/grid.00090.p00000
size= 2097160
tsim=3.791667e+01
0oct=0x7fe863f80010
done with 5401 octs
dumping data/grid.00090.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00090 6 701 data/grid.00090.f701 1 \
0"
> ws1
> pli,e(,,32)
> i=90;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 64x64x64 cube from data/grid.00090.p00000
size= 2097160
tsim=3.791667e+01
0oct=0x7fe863f80010
done with 5401 octs
dumping data/grid.00090.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00090 6 701 data/grid.00090.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00090.p00000
size= 2097160
tsim=3.791667e+01
0oct=0x7fe863f80010
done with 5401 octs
dumping data/grid.00090.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00090 6 705 data/grid.00090.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00090.p00000
size= 2097160
tsim=3.791667e+01
0oct=0x7fe863f80010
done with 5401 octs
dumping data/grid.00090.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00090 6 0 data/grid.00090.f0 1 0"
> pli,x(,,32)
> ws
> i=90;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 64x64x64 cube from data/grid.00090.p00000
size= 2097160
tsim=3.791667e+01
0oct=0x7fe863f80010
done with 5401 octs
dumping data/grid.00090.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00090 6 701 data/grid.00090.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00090.p00000
size= 2097160
tsim=3.791667e+01
0oct=0x7fe863f80010
done with 5401 octs
dumping data/grid.00090.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00090 6 706 data/grid.00090.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00090.p00000
size= 2097160
tsim=3.791667e+01
0oct=0x7fe863f80010
done with 5401 octs
dumping data/grid.00090.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00090 6 0 data/grid.00090.f0 1 0"
> pli,x(,,32)
> plotamr(l(,,32))
6
5
4
3
2
1
[]
> palette,"idl-33.gp"
> i=9;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 64x64x64 cube from data/grid.00009.p00000
size= 2097160
tsim=3.754167e+02
0oct=0x7fde5ff53010
done with 8241 octs
dumping data/grid.00009.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00009 6 701 data/grid.00009.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00009.p00000
size= 2097160
tsim=3.754167e+02
0oct=0x7fde5ff53010
done with 8241 octs
dumping data/grid.00009.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00009 6 706 data/grid.00009.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00009.p00000
size= 2097160
tsim=3.754167e+02
0oct=0x7fde5ff53010
done with 8241 octs
dumping data/grid.00009.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00009 6 0 data/grid.00009.f0 1 0"
> plotamr(l(,,32))
6
5
4
3
2
1
[]
> ws
> pli,x(,,32)
> plotamr(l(,,32))
6
5
4
3
2
1
[]
> i=12;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 64x64x64 cube from data/grid.00012.p00000
size= 2097160
tsim=5.004167e+02
0oct=0x7fde5ff53010
done with 9153 octs
dumping data/grid.00012.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00012 6 701 data/grid.00012.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00012.p00000
size= 2097160
tsim=5.004167e+02
0oct=0x7fde5ff53010
done with 9153 octs
dumping data/grid.00012.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00012 6 706 data/grid.00012.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00012.p00000
size= 2097160
tsim=5.004167e+02
0oct=0x7fde5ff53010
done with 9153 octs
dumping data/grid.00012.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00012 6 0 data/grid.00012.f0 1 0"
> ws1
> palette,"idl-33.gp"
> pli,x(,,32)
> plotamr(l(,,32))
6
5
4
3
2
1
[]
> pli,e(,,32)
> ws
> pli,log10(e(,,32))
> plotamr(l(,,32))
6
5
4
3
2
1
[]
> i=13;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 64x64x64 cube from data/grid.00013.p00000
size= 2097160
tsim=5.420833e+02
0oct=0x7fde5ff53010
done with 9249 octs
dumping data/grid.00013.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00013 6 701 data/grid.00013.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00013.p00000
size= 2097160
tsim=5.420833e+02
0oct=0x7fde5ff53010
done with 9249 octs
dumping data/grid.00013.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00013 6 706 data/grid.00013.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00013.p00000
size= 2097160
tsim=5.420833e+02
0oct=0x7fde5ff53010
done with 9249 octs
dumping data/grid.00013.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00013 6 0 data/grid.00013.f0 1 0"
> ws
> pli,log10(e(,,32))
> plotamr(l(,,32))
6
5
4
3
2
1
[]
> ws
> pli,x(,,32)
> plotamr(l(,,32))
6
5
4
3
2
1
[]
> i=15;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 64x64x64 cube from data/grid.00015.p00000
size= 2097160
tsim=6.254167e+02
0oct=0x7fde5ff53010
done with 9897 octs
dumping data/grid.00015.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00015 6 701 data/grid.00015.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00015.p00000
size= 2097160
tsim=6.254167e+02
0oct=0x7fde5ff53010
done with 9897 octs
dumping data/grid.00015.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00015 6 706 data/grid.00015.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00015.p00000
size= 2097160
tsim=6.254167e+02
0oct=0x7fde5ff53010
done with 9897 octs
dumping data/grid.00015.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00015 6 0 data/grid.00015.f0 1 0"
> ws
> pli,x(,,32)
> plotamr(l(,,32))
6
5
4
3
2
1
[]
> i=59;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 64x64x64 cube from data/grid.00059.p00000
size= 2097160
tsim=2.458750e+03
0oct=0x7fde5ff53010
done with 17137 octs
dumping data/grid.00059.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00059 6 701 data/grid.00059.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00059.p00000
size= 2097160
tsim=2.458750e+03
0oct=0x7fde5ff53010
done with 17137 octs
dumping data/grid.00059.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00059 6 706 data/grid.00059.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00059.p00000
size= 2097160
tsim=2.458750e+03
0oct=0x7fde5ff53010
done with 17137 octs
dumping data/grid.00059.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00059 6 0 data/grid.00059.f0 1 0"
> ws
> pli,x(,,32)
> plotamr(l(,,32))
6
5
4
3
2
1
[]
> ws
> pli,e(,,32)
> plotamr(l(,,32))
6
5
4
3
2
1
[]
> ws
> plotamr(l(,,32))
6
5
4
3
2
1
[]
> pli,e(,,32)
> ws
> pli,log10(e(,,32))
> plotamr(l(,,32))
6
5
4
3
2
1
[]
> ws
> plotamr(l(,,32))
6
5
4
3
2
1
[]
> ws
> plg,x(32,,32)
> logxy,0,1
> ws
> plg,x(32,,32)
> logxy,0,1
> PL,x(32,,32)
> i=5;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=2.500000e+00
0oct=0x7f39ec26e010
done with 4841 octs
dumping data/grid.00005.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 701 data/grid.00005.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=2.500000e+00
0oct=0x7f39ec26e010
done with 4841 octs
dumping data/grid.00005.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 706 data/grid.00005.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=2.500000e+00
0oct=0x7f39ec26e010
done with 4841 octs
dumping data/grid.00005.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 0 data/grid.00005.f0 1 0"
> ws
> plotamr(l(,,32))
6
5
4
3
2
1
[]
> i=20;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 64x64x64 cube from data/grid.00020.p00000
size= 2097160
tsim=8.750000e+00
0oct=0x7fa0816b3010
done with 7265 octs
dumping data/grid.00020.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00020 6 701 data/grid.00020.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00020.p00000
size= 2097160
tsim=8.750000e+00
0oct=0x7fa0816b3010
done with 7265 octs
dumping data/grid.00020.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00020 6 706 data/grid.00020.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00020.p00000
size= 2097160
tsim=8.750000e+00
0oct=0x7fa0816b3010
done with 7265 octs
dumping data/grid.00020.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00020 6 0 data/grid.00020.f0 1 0"
> ws
> plotamr(l(,,32))
6
5
4
3
2
1
[]
> i=20;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 64x64x64 cube from data/grid.00020.p00000
size= 2097160
tsim=8.750000e+00
0oct=0x7f73b36df010
done with 5401 octs
dumping data/grid.00020.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00020 6 701 data/grid.00020.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00020.p00000
size= 2097160
tsim=8.750000e+00
0oct=0x7f73b36df010
done with 5401 octs
dumping data/grid.00020.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00020 6 706 data/grid.00020.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00020.p00000
size= 2097160
tsim=8.750000e+00
0oct=0x7f73b36df010
done with 5401 octs
dumping data/grid.00020.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00020 6 0 data/grid.00020.f0 1 0"
> ws
> plotamr(l(,,32))
6
5
4
3
2
1
[]
> i=20;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 64x64x64 cube from data/grid.00020.p00000
size= 2097160
tsim=8.750000e+00
0oct=0x7f935f942010
done with 5257 octs
dumping data/grid.00020.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00020 6 701 data/grid.00020.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00020.p00000
size= 2097160
tsim=8.750000e+00
0oct=0x7f935f942010
done with 5257 octs
dumping data/grid.00020.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00020 6 706 data/grid.00020.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00020.p00000
size= 2097160
tsim=8.750000e+00
0oct=0x7f935f942010
done with 5257 octs
dumping data/grid.00020.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00020 6 0 data/grid.00020.f0 1 0"
> ws
> plg,x(32,,32)
> ws
> plotamr(l(,,32))
6
5
4
3
2
1
[]
> i=100;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 64x64x64 cube from data/grid.00100.p00000
size= 2097160
tsim=4.208333e+01
0oct=0x7f1063359010
done with 6097 octs
dumping data/grid.00100.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00100 6 701 data/grid.00100.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00100.p00000
size= 2097160
tsim=4.208333e+01
0oct=0x7f1063359010
done with 6097 octs
dumping data/grid.00100.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00100 6 706 data/grid.00100.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00100.p00000
size= 2097160
tsim=4.208333e+01
0oct=0x7f1063359010
done with 6097 octs
dumping data/grid.00100.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00100 6 0 data/grid.00100.f0 1 0"
> ws
> plotamr(l(,,32))
6
5
4
3
2
1
[]
> i=160;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 64x64x64 cube from data/grid.00160.p00000
size= 2097160
tsim=6.708333e+01
0oct=0x7f1063359010
done with 6481 octs
dumping data/grid.00160.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00160 6 701 data/grid.00160.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00160.p00000
size= 2097160
tsim=6.708333e+01
0oct=0x7f1063359010
done with 6481 octs
dumping data/grid.00160.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00160 6 706 data/grid.00160.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00160.p00000
size= 2097160
tsim=6.708333e+01
0oct=0x7f1063359010
done with 6481 octs
dumping data/grid.00160.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00160 6 0 data/grid.00160.f0 1 0"
> ws
> plotamr(l(,,32))
6
5
4
3
2
1
[]
> i=200;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 64x64x64 cube from data/grid.00200.p00000
size= 2097160
tsim=8.375000e+01
0oct=0x7f1063359010
done with 6729 octs
dumping data/grid.00200.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00200 6 701 data/grid.00200.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00200.p00000
size= 2097160
tsim=8.375000e+01
0oct=0x7f1063359010
done with 6729 octs
dumping data/grid.00200.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00200 6 706 data/grid.00200.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00200.p00000
size= 2097160
tsim=8.375000e+01
0oct=0x7f1063359010
done with 6729 octs
dumping data/grid.00200.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00200 6 0 data/grid.00200.f0 1 0"
> ws
> plotamr(l(,,32))
6
5
4
3
2
1
[]
> i=160;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 64x64x64 cube from data/grid.00160.p00000
size= 2097160
tsim=6.708333e+01
0oct=0x7f1063359010
done with 6481 octs
dumping data/grid.00160.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00160 6 701 data/grid.00160.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00160.p00000
size= 2097160
tsim=6.708333e+01
0oct=0x7f1063359010
done with 6481 octs
dumping data/grid.00160.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00160 6 706 data/grid.00160.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00160.p00000
size= 2097160
tsim=6.708333e+01
0oct=0x7f1063359010
done with 6481 octs
dumping data/grid.00160.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00160 6 0 data/grid.00160.f0 1 0"
> ws
> plotamr(l(,,32))
6
5
4
3
2
1
[]
> i=200;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 64x64x64 cube from data/grid.00200.p00000
size= 2097160
tsim=8.375000e+01
0oct=0x7f1063359010
done with 6729 octs
dumping data/grid.00200.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00200 6 701 data/grid.00200.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00200.p00000
size= 2097160
tsim=8.375000e+01
0oct=0x7f1063359010
done with 6729 octs
dumping data/grid.00200.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00200 6 706 data/grid.00200.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00200.p00000
size= 2097160
tsim=8.375000e+01
0oct=0x7f1063359010
done with 6729 octs
dumping data/grid.00200.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00200 6 0 data/grid.00200.f0 1 0"
> ws
> plotamr(l(,,32))
6
5
4
3
2
1
[]
> pli,x(,,32)
> plotamr(l(,,32))
6
5
4
3
2
1
[]
> ws
> i=200;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 64x64x64 cube from data/grid.00200.p00000
size= 2097160
tsim=8.375000e+01
0oct=0x7f2510212010
done with 6729 octs
dumping data/grid.00200.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00200 6 701 data/grid.00200.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00200.p00000
size= 2097160
tsim=8.375000e+01
0oct=0x7f2510212010
done with 6729 octs
dumping data/grid.00200.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00200 6 706 data/grid.00200.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00200.p00000
size= 2097160
tsim=8.375000e+01
0oct=0x7f2510212010
done with 6729 octs
dumping data/grid.00200.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00200 6 0 data/grid.00200.f0 1 0"
> ws
> plotamr(l(,,32))
6
5
4
3
2
1
[]
> i=200;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 64x64x64 cube from data/grid.00200.p00000
size= 2097160
tsim=8.375000e+01
0oct=0x7f8adadce010
done with 6577 octs
dumping data/grid.00200.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00200 6 701 data/grid.00200.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00200.p00000
size= 2097160
tsim=8.375000e+01
0oct=0x7f8adadce010
done with 6577 octs
dumping data/grid.00200.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00200 6 706 data/grid.00200.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00200.p00000
size= 2097160
tsim=8.375000e+01
0oct=0x7f8adadce010
done with 6577 octs
dumping data/grid.00200.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00200 6 0 data/grid.00200.f0 1 0"
> ws
> plotamr(l(,,32))
6
5
4
3
2
1
[]
> pli,x(,,32)
> plotamr(l(,,32))
6
5
4
3
2
1
[]
> ws
> plg,x(32,,32)
> pPL,x(32,,32)
ERROR (*main*) attempt to call non-function or index non-array
WARNING source code unavailable (try dbdis function)
now at pc= 1 (of 18), failed at pc= 12
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> PL,x(32,,32)
> i=9;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 128x128x128 cube from data/grid.00009.p00000
size= 16777224
tsim=4.166667e+00
0oct=0x7f246b293010
done with 6681 octs
dumping data/grid.00009.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00009 7 701 data/grid.00009.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00009.p00000
size= 16777224
tsim=4.166667e+00
0oct=0x7f246b293010
done with 6681 octs
dumping data/grid.00009.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00009 7 706 data/grid.00009.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00009.p00000
size= 16777224
tsim=4.166667e+00
0oct=0x7f246b293010
done with 6681 octs
dumping data/grid.00009.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00009 7 0 data/grid.00009.f0 1 0"
> ws
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> pli,x(,,64)
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> i=30;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 128x128x128 cube from data/grid.00030.p00000
size= 16777224
tsim=1.291667e+01
0oct=0x7f19e46ca010
done with 8185 octs
dumping data/grid.00030.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00030 7 701 data/grid.00030.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00030.p00000
size= 16777224
tsim=1.291667e+01
0oct=0x7f19e46ca010
done with 8185 octs
dumping data/grid.00030.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00030 7 706 data/grid.00030.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00030.p00000
size= 16777224
tsim=1.291667e+01
0oct=0x7f19e46ca010
done with 8185 octs
dumping data/grid.00030.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00030 7 0 data/grid.00030.f0 1 0"
> ws
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> pli,x(,,64)
> pli,e(,,32)
> pli,e(,,64)
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> i=77;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 128x128x128 cube from data/grid.00077.p00000
size= 16777224
tsim=3.250000e+01
0oct=0x7f19e46ca010
done with 10753 octs
dumping data/grid.00077.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00077 7 701 data/grid.00077.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00077.p00000
size= 16777224
tsim=3.250000e+01
0oct=0x7f19e46ca010
done with 10753 octs
dumping data/grid.00077.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00077 7 706 data/grid.00077.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00077.p00000
size= 16777224
tsim=3.250000e+01
0oct=0x7f19e46ca010
done with 10753 octs
dumping data/grid.00077.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00077 7 0 data/grid.00077.f0 1 0"
> ws
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> pli,e(,,64)
> pli,log10(e(,,64))
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> ws
> plg,x(64,,64)
> PL,x(64,,64)
> logxy,0,1
> ws
> plg,l(64,,64)
> i=146;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 128x128x128 cube from data/grid.00146.p00000
size= 16777224
tsim=6.125000e+01
0oct=0x7f19e46ca010
done with 13049 octs
dumping data/grid.00146.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00146 7 701 data/grid.00146.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00146.p00000
size= 16777224
tsim=6.125000e+01
0oct=0x7f19e46ca010
done with 13049 octs
dumping data/grid.00146.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00146 7 706 data/grid.00146.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00146.p00000
size= 16777224
tsim=6.125000e+01
0oct=0x7f19e46ca010
done with 13049 octs
dumping data/grid.00146.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00146 7 0 data/grid.00146.f0 1 0"
> w
[]
> ws
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> pli,x(,,64)
> i=180;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 128x128x128 cube from data/grid.00180.p00000
size= 16777224
tsim=7.541667e+01
0oct=0x7f19e46ca010
done with 13849 octs
dumping data/grid.00180.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00180 7 701 data/grid.00180.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00180.p00000
size= 16777224
tsim=7.541667e+01
0oct=0x7f19e46ca010
done with 13849 octs
dumping data/grid.00180.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00180 7 706 data/grid.00180.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00180.p00000
size= 16777224
tsim=7.541667e+01
0oct=0x7f19e46ca010
done with 13849 octs
dumping data/grid.00180.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00180 7 0 data/grid.00180.f0 1 0"
> ws
> pli,x(,,64)
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> pli,log10(e(,,64))
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> ws
> i=200;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 128x128x128 cube from data/grid.00200.p00000
size= 16777224
tsim=8.375000e+01
0oct=0x7f19e46ca010
done with 14017 octs
dumping data/grid.00200.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00200 7 701 data/grid.00200.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00200.p00000
size= 16777224
tsim=8.375000e+01
0oct=0x7f19e46ca010
done with 14017 octs
dumping data/grid.00200.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00200 7 706 data/grid.00200.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00200.p00000
size= 16777224
tsim=8.375000e+01
0oct=0x7f19e46ca010
done with 14017 octs
dumping data/grid.00200.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00200 7 0 data/grid.00200.f0 1 0"
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> i=3;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 128x128x128 cube from data/grid.00003.p00000
size= 16777224
tsim=1.254167e+02
0oct=0x7fb661254010
done with 16705 octs
dumping data/grid.00003.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00003 7 701 data/grid.00003.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00003.p00000
size= 16777224
tsim=1.254167e+02
0oct=0x7fb661254010
done with 16705 octs
dumping data/grid.00003.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00003 7 706 data/grid.00003.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00003.p00000
size= 16777224
tsim=1.254167e+02
0oct=0x7fb661254010
done with 16705 octs
dumping data/grid.00003.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00003 7 0 data/grid.00003.f0 1 0"
> ws
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> pli,log10(e(,,64))
> pli,x(,,64)
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> pli,x(,,64)
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> pli,x(,,64)
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> pli,x(,,64)
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> pli,x(,,64)
> 
> pli,log10(1.-x(,,64))
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> rvp
> i=4;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=1.670833e+02
0oct=0x7fb661254010
done with 18697 octs
dumping data/grid.00004.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 701 data/grid.00004.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=1.670833e+02
0oct=0x7fb661254010
done with 18697 octs
dumping data/grid.00004.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 706 data/grid.00004.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=1.670833e+02
0oct=0x7fb661254010
done with 18697 octs
dumping data/grid.00004.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 0 data/grid.00004.f0 1 0"
> ws
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> pli,log10(1.-x(,,64))
> pli,x(,,64)
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> ws
> pli,log10(e(,,64))
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> i=11;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 128x128x128 cube from data/grid.00011.p00000
size= 16777224
tsim=4.587500e+02
0oct=0x7fb661254010
done with 28761 octs
dumping data/grid.00011.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00011 7 701 data/grid.00011.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00011.p00000
size= 16777224
tsim=4.587500e+02
0oct=0x7fb661254010
done with 28761 octs
dumping data/grid.00011.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00011 7 706 data/grid.00011.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00011.p00000
size= 16777224
tsim=4.587500e+02
0oct=0x7fb661254010
done with 28761 octs
dumping data/grid.00011.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00011 7 0 data/grid.00011.f0 1 0"
> ws
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> pli,log10(e(,,64))
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> ws
> pli,x(,,64)
> pltitle,"H ionized fraction"
> ws,1
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> pltitle,"AMR Mesh"
> lim
[]
> limits
> pli,log10(e(,,64))
> i=17;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 128x128x128 cube from data/grid.00017.p00000
size= 16777224
tsim=7.087500e+02
0oct=0x7fb661254010
done with 34817 octs
dumping data/grid.00017.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00017 7 701 data/grid.00017.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00017.p00000
size= 16777224
tsim=7.087500e+02
0oct=0x7fb661254010
done with 34817 octs
dumping data/grid.00017.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00017 7 706 data/grid.00017.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00017.p00000
size= 16777224
tsim=7.087500e+02
0oct=0x7fb661254010
done with 34817 octs
dumping data/grid.00017.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00017 7 0 data/grid.00017.f0 1 0"
> ws
> wkl
> ws1
> pli,x(,,64)
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> ws
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> pli,x(,,64)
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> palette,"idl-33.gp"
> ws
> pli,x(,,64)
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> ws
> pli,e(,,64)
> pli,log10(e(,,64))
> i=25;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 128x128x128 cube from data/grid.00025.p00000
size= 16777224
tsim=1.042083e+03
0oct=0x7fb661254010
done with 41353 octs
dumping data/grid.00025.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00025 7 701 data/grid.00025.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00025.p00000
size= 16777224
tsim=1.042083e+03
0oct=0x7fb661254010
done with 41353 octs
dumping data/grid.00025.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00025 7 706 data/grid.00025.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00025.p00000
size= 16777224
tsim=1.042083e+03
0oct=0x7fb661254010
done with 41353 octs
dumping data/grid.00025.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00025 7 0 data/grid.00025.f0 1 0"
> ws
> pli,log10(e(,,64))
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> ws
> pli,x(,,64)
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> i=25;lvl=8;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 256x256x256 cube from data/grid.00025.p00000
size= 134217736
tsim=1.042083e+03
0oct=0x7fb661254010
done with 41353 octs
dumping data/grid.00025.f701.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00025 8 701 data/grid.00025.f701 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00025.p00000
size= 134217736
tsim=1.042083e+03
0oct=0x7fb661254010
done with 41353 octs
dumping data/grid.00025.f706.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00025 8 706 data/grid.00025.f706 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00025.p00000
size= 134217736
tsim=1.042083e+03
0oct=0x7fb661254010
  C-c"~/Project/Quartz/utils/oct2grid data/grid.00025 8 0 data/grid.00025.f0 1 0"
> i=2;lvl=8;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 256x256x256 cube from data/grid.00002.p00000
size= 134217736
tsim=1.250000e+00
0oct=0x7fffee978010
done with 9089 octs
dumping data/grid.00002.f701.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 8 701 data/grid.00002.f701 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00002.p00000
size= 134217736
tsim=1.250000e+00
0oct=0x7fffee978010
done with 9089 octs
dumping data/grid.00002.f706.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 8 706 data/grid.00002.f706 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00002.p00000
size= 134217736
tsim=1.250000e+00
0oct=0x7fffee978010
done with 9089 octs
dumping data/grid.00002.f0.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 8 0 data/grid.00002.f0 1 0"
> ws
> plotamr(l(,,128))
8
7
6
5
4
3
2
1
[]
> i=2;lvl=8;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 256x256x256 cube from data/grid.00002.p00000
size= 134217736
tsim=1.250000e+00
0oct=0x7fffee978010
done with 9305 octs
dumping data/grid.00002.f701.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 8 701 data/grid.00002.f701 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00002.p00000
size= 134217736
tsim=1.250000e+00
0oct=0x7fffee978010
done with 9305 octs
dumping data/grid.00002.f706.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 8 706 data/grid.00002.f706 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00002.p00000
size= 134217736
tsim=1.250000e+00
0oct=0x7fffee978010
done with 9305 octs
dumping data/grid.00002.f0.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 8 0 data/grid.00002.f0 1 0"
> ws
> plotamr(l(,,128))
8
7
6
5
4
3
2
1
[]
> i=2;lvl=8;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 256x256x256 cube from data/grid.00002.p00000
size= 134217736
tsim=1.250000e+00
0oct=0x7fffee978010
done with 9305 octs
dumping data/grid.00002.f701.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 8 701 data/grid.00002.f701 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00002.p00000
size= 134217736
tsim=1.250000e+00
0oct=0x7fffee978010
done with 9305 octs
dumping data/grid.00002.f706.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 8 706 data/grid.00002.f706 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00002.p00000
size= 134217736
tsim=1.250000e+00
0oct=0x7fffee978010
done with 9305 octs
dumping data/grid.00002.f0.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 8 0 data/grid.00002.f0 1 0"
ws
> pli,x(,,64)
> pli,x(,,128)
> ws
> pli,x(,,128)
> plotamr(l(,,128))
8
7
6
5
4
3
2
1
[]
> i=6;lvl=8;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 256x256x256 cube from data/grid.00006.p00000
size= 134217736
tsim=2.504167e+02
0oct=0x7fb661254010
done with 21913 octs
dumping data/grid.00006.f701.p00000 with nmap=16777216
  C-c"~/Project/Quartz/utils/oct2grid data/grid.00006 8 701 data/grid.00006.f701 1 \
0"
ERROR (readcube) encountered end-of-file before read completed
  LINE: 31  FILE: /home/daubert/Project/Quartz/utils/readamr.i
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> i=6;lvl=8;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 256x256x256 cube from data/grid.00006.p00000
size= 134217736
tsim=2.541667e+01
0oct=0x7fffee978010
done with 34689 octs
dumping data/grid.00006.f701.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 8 701 data/grid.00006.f701 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00006.p00000
size= 134217736
tsim=2.541667e+01
0oct=0x7fffee978010
done with 34689 octs
dumping data/grid.00006.f706.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 8 706 data/grid.00006.f706 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00006.p00000
size= 134217736
tsim=2.541667e+01
0oct=0x7fffee978010
done with 34689 octs
dumping data/grid.00006.f0.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 8 0 data/grid.00006.f0 1 0"
> ws1
> pli,x(,,128)
> plotamr(l(,,128))
8
7
6
5
4
3
2
1
[]
> ws
> plotamr(l(,,128))
8
7
6
5
4
3
2
1
[]
> ws
> pli,x(,,128)
> palette,"idl-33.gp"
> plotamr(l(,,128))
8
7
6
5
4
3
2
1
[]
> i=20;lvl=8;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 256x256x256 cube from data/grid.00020.p00000
size= 134217736
tsim=8.375000e+01
0oct=0x7fffee978010
ws
done with 68305 octs
dumping data/grid.00020.f701.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00020 8 701 data/grid.00020.f701 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00020.p00000
size= 134217736
tsim=8.375000e+01
0oct=0x7fffee978010
done with 68305 octs
dumping data/grid.00020.f706.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00020 8 706 data/grid.00020.f706 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00020.p00000
size= 134217736
tsim=8.375000e+01
0oct=0x7fffee978010
done with 68305 octs
dumping data/grid.00020.f0.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00020 8 0 data/grid.00020.f0 1 0"
> pli,x(,,128)
> plotamr(l(,,128))
8
7
6
5
4
3
2
1
[]
> ws
> pli,x(,,128)
> pli,log10(e(,,64))
> plotamr(l(,,128))
8
7
6
5
4
3
2
1
[]
> ws
> pli,log10(e(,,128))
> plotamr(l(,,128))
8
7
6
5
4
3
2
1
[]
> i=20;lvl=8;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707);
Casting Rays on 256x256x256 cube from data/grid.00020.p00000
size= 134217736
tsim=8.375000e+01
0oct=0x7f0f53efa010
done with 4681 octs
dumping data/grid.00020.f701.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00020 8 701 data/grid.00020.f701 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00020.p00000
size= 134217736
tsim=8.375000e+01
0oct=0x7f0f53efa010
  C-c"~/Project/Quartz/utils/oct2grid data/grid.00020 8 706 data/grid.00020.f706 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00020.p00000
size= 134217736
tsim=8.375000e+01
0oct=0x7f0f53efa010
  C-c"~/Project/Quartz/utils/oct2grid data/grid.00020 8 0 data/grid.00020.f0 1 0"
Casting Rays on 256x256x256 cube from data/grid.00020.p00000
size= 134217736
tsim=8.375000e+01
0oct=0x7f0f53efa010
  C-c"~/Project/Quartz/utils/oct2grid data/grid.00020 8 707 data/grid.00020.f707 1 \
0"
> i=20;lvl=5;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707);
Casting Rays on 32x32x32 cube from data/grid.00020.p00000
size= 262152
tsim=8.375000e+01
0oct=0x7f0f53efa010
done with 4681 octs
dumping data/grid.00020.f701.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00020 5 701 data/grid.00020.f701 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00020.p00000
size= 262152
tsim=8.375000e+01
0oct=0x7f0f53efa010
done with 4681 octs
dumping data/grid.00020.f706.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00020 5 706 data/grid.00020.f706 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00020.p00000
size= 262152
tsim=8.375000e+01
0oct=0x7f0f53efa010
done with 4681 octs
dumping data/grid.00020.f0.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00020 5 0 data/grid.00020.f0 1 0"
Casting Rays on 32x32x32 cube from data/grid.00020.p00000
size= 262152
tsim=8.375000e+01
0oct=0x7f0f53efa010
done with 4681 octs
dumping data/grid.00020.f707.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00020 5 707 data/grid.00020.f707 1 \
0"
> ws
> info,x
 array(double,32,32,32)
> pli,x(,,16)
> pli,t(,,16)
> ws
> plg,x(16,,16)
> plg,t(16,,16)
> logxy,0,1
> i=5;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707);
Casting Rays on 128x128x128 cube from data/grid.00005.p00000
size= 16777224
tsim=2.125000e+01
0oct=0x7f21c36c0010
done with 17185 octs
dumping data/grid.00005.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 7 701 data/grid.00005.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00005.p00000
size= 16777224
tsim=2.125000e+01
0oct=0x7f21c36c0010
done with 17185 octs
dumping data/grid.00005.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 7 706 data/grid.00005.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00005.p00000
size= 16777224
tsim=2.125000e+01
0oct=0x7f21c36c0010
done with 17185 octs
dumping data/grid.00005.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 7 0 data/grid.00005.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00005.p00000
size= 16777224
tsim=2.125000e+01
0oct=0x7f21c36c0010
done with 17185 octs
dumping data/grid.00005.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 7 707 data/grid.00005.f707 1 \
0"
> ws
>
> ws1
> info,t
 array(double,128,128,128)
> pli,t(,,64)
> plg,t(64,,64)
> ws
> plg,t(64,,64)
> logxy,0,1
> ws
> pli,t(,,64)
> palette,"idl-33.gp"
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> ws
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> i=5;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707);
Casting Rays on 128x128x128 cube from data/grid.00005.p00000
size= 16777224
tsim=2.125000e+01
0oct=0x7fba4627f010
done with 17185 octs
dumping data/grid.00005.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 7 701 data/grid.00005.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00005.p00000
size= 16777224
tsim=2.125000e+01
0oct=0x7fba4627f010
done with 17185 octs
dumping data/grid.00005.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 7 706 data/grid.00005.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00005.p00000
size= 16777224
tsim=2.125000e+01
0oct=0x7fba4627f010
done with 17185 octs
dumping data/grid.00005.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 7 0 data/grid.00005.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00005.p00000
size= 16777224
tsim=2.125000e+01
0oct=0x7fba4627f010
done with 17185 octs
dumping data/grid.00005.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 7 707 data/grid.00005.f707 1 \
0"
> ws
> ws1
> palette,"idl-33.gp"
> info,x
 array(double,128,128,128)
> pli,x(,,64)
> pli,e(,,64)
> pli,log10(e(,,64))
> pli,log10(e(,,70))
> pli,log10(e(,,75))
> pli,x(,,75)
> pli,log10(x(,,75))
> i=7;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707);
Casting Rays on 128x128x128 cube from data/grid.00007.p00000
size= 16777224
tsim=2.958333e+01
0oct=0x7fba4627f010
done with 19284 octs
dumping data/grid.00007.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00007 7 701 data/grid.00007.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00007.p00000
size= 16777224
tsim=2.958333e+01
0oct=0x7fba4627f010
done with 19284 octs
dumping data/grid.00007.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00007 7 706 data/grid.00007.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00007.p00000
size= 16777224
tsim=2.958333e+01
0oct=0x7fba4627f010
done with 19284 octs
dumping data/grid.00007.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00007 7 0 data/grid.00007.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00007.p00000
size= 16777224
tsim=2.958333e+01
0oct=0x7fba4627f010
done with 19284 octs
dumping data/grid.00007.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00007 7 707 data/grid.00007.f707 1 \
0"
> ws1
> pli,log10(x(,,75))
> pli,log10(x(,,70))
> pli,log10(x(,,80))
> rvp
> rvp
> palette,"idl-33.gp"
> pli,log10(e(,,80))
> plotamr(l(,,80))
7
6
5
4
3
2
1
[]
> ws
> plk,log10(e(,,80))
> plotamr(l(,,80))
7
6
5
4
3
2
1
[]
> i=9;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707);
Casting Rays on 128x128x128 cube from data/grid.00009.p00000
size= 16777224
tsim=3.791667e+01
0oct=0x7fba4627f010
done with 21598 octs
dumping data/grid.00009.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00009 7 701 data/grid.00009.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00009.p00000
size= 16777224
tsim=3.791667e+01
0oct=0x7fba4627f010
done with 21598 octs
dumping data/grid.00009.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00009 7 706 data/grid.00009.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00009.p00000
size= 16777224
tsim=3.791667e+01
0oct=0x7fba4627f010
done with 21598 octs
dumping data/grid.00009.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00009 7 0 data/grid.00009.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00009.p00000
size= 16777224
tsim=3.791667e+01
0oct=0x7fba4627f010
done with 21598 octs
dumping data/grid.00009.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00009 7 707 data/grid.00009.f707 1 \
0"
> wS1
[]
> ws1
> pli,log10(e(,,80))
> plotamr(l(,,80))
7
6
5
4
3
2
1
[]
> pli,x(,,80)
> pli,log10(x(,,80))
> plotamr(l(,,80))
7
6
5
4
3
2
1
[]
> ws
> pli,log10(e(,,80))
> pli,e(,,80)
> i=10;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707);
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=4.208333e+01
0oct=0x7fba4627f010
done with 22777 octs
dumping data/grid.00010.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 701 data/grid.00010.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=4.208333e+01
0oct=0x7fba4627f010
done with 22777 octs
dumping data/grid.00010.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 706 data/grid.00010.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=4.208333e+01
0oct=0x7fba4627f010
done with 22777 octs
dumping data/grid.00010.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 0 data/grid.00010.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=4.208333e+01
0oct=0x7fba4627f010
done with 22777 octs
dumping data/grid.00010.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 707 data/grid.00010.f707 1 \
0"
> pli,x(,,80)
> pli,log10(x(,,80))
> plotamr(l(,,80))
7
6
5
4
3
2
1
[]
> ws
> pli,log10(e(,,80))
> plotamr(l(,,80))
7
6
5
4
3
2
1
[]
> i=12;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707);
Casting Rays on 128x128x128 cube from data/grid.00012.p00000
size= 16777224
tsim=5.041667e+01
0oct=0x7fba4627f010
done with 24495 octs
dumping data/grid.00012.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00012 7 701 data/grid.00012.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00012.p00000
size= 16777224
tsim=5.041667e+01
0oct=0x7fba4627f010
done with 24495 octs
dumping data/grid.00012.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00012 7 706 data/grid.00012.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00012.p00000
size= 16777224
tsim=5.041667e+01
0oct=0x7fba4627f010
done with 24495 octs
dumping data/grid.00012.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00012 7 0 data/grid.00012.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00012.p00000
size= 16777224
tsim=5.041667e+01
0oct=0x7fba4627f010
done with 24495 octs
dumping data/grid.00012.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00012 7 707 data/grid.00012.f707 1 \
0"
> pli,x(,,80)
> pli,x(,,64)
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> ws
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> i=14;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707);
Casting Rays on 128x128x128 cube from data/grid.00014.p00000
size= 16777224
tsim=5.875000e+01
0oct=0x7fba4627f010
done with 26145 octs
dumping data/grid.00014.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00014 7 701 data/grid.00014.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00014.p00000
size= 16777224
tsim=5.875000e+01
0oct=0x7fba4627f010
done with 26145 octs
dumping data/grid.00014.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00014 7 706 data/grid.00014.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00014.p00000
size= 16777224
tsim=5.875000e+01
0oct=0x7fba4627f010
done with 26145 octs
dumping data/grid.00014.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00014 7 0 data/grid.00014.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00014.p00000
size= 16777224
tsim=5.875000e+01
0oct=0x7fba4627f010
done with 26145 octs
dumping data/grid.00014.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00014 7 707 data/grid.00014.f707 1 \
0"
> pli,x(,,64)
> pli,x(,,90)
> plotamr(l(,,90))
7
6
5
4
3
2
1
[]
> palette,"idl-33.gp"
> rvp
> rvp
> pli,t(,,64)
> pli,t(,,90)
> pli,tlog10((,,90))
SYNTAX: syntax error near ,,90))
SYNTAX: syntax error near <EOF>
> pli,log10(t(,,90))
> ws
> pli,log10(t(,,90))
> pli,log10(e(,,90))
> pli,log10(e(,,64))
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> i=4;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707);
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=1.708333e+01
0oct=0x7ff2a3ee8010
done with 9981 octs
dumping data/grid.00004.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 701 data/grid.00004.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=1.708333e+01
0oct=0x7ff2a3ee8010
done with 9981 octs
dumping data/grid.00004.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 706 data/grid.00004.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=1.708333e+01
0oct=0x7ff2a3ee8010
done with 9981 octs
dumping data/grid.00004.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 0 data/grid.00004.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=1.708333e+01
0oct=0x7ff2a3ee8010
done with 9981 octs
dumping data/grid.00004.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 707 data/grid.00004.f707 1 \
0"
> ws
> ws1
> pli,np(,,64)
ERROR (*main*) attempt to call non-function or index non-array
WARNING source code unavailable (try dbdis function)
now at pc= 1 (of 19), failed at pc= 9
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> i=4;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707);n=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=1.708333e+01
0oct=0x7ff2a3ee8010
done with 9981 octs
dumping data/grid.00004.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 701 data/grid.00004.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=1.708333e+01
0oct=0x7ff2a3ee8010
done with 9981 octs
dumping data/grid.00004.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 706 data/grid.00004.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=1.708333e+01
0oct=0x7ff2a3ee8010
done with 9981 octs
dumping data/grid.00004.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 0 data/grid.00004.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=1.708333e+01
0oct=0x7ff2a3ee8010
done with 9981 octs
dumping data/grid.00004.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 707 data/grid.00004.f707 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=1.708333e+01
0oct=0x7ff2a3ee8010
done with 9981 octs
dumping data/grid.00004.f708.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 708 data/grid.00004.f708 1 \
0"
> info,nh
 []
> info,np
 []
> info,n
 array(double,128,128,128)
> pli,n(,,64)
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> ws
> pli,x(,,64)
> pli,log10(e(,,64))
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> ws
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> i=10;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707);n=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=4.208333e+01
0oct=0x7ff2a3ee8010
done with 13665 octs
dumping data/grid.00010.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 701 data/grid.00010.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=4.208333e+01
0oct=0x7ff2a3ee8010
done with 13665 octs
dumping data/grid.00010.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 706 data/grid.00010.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=4.208333e+01
0oct=0x7ff2a3ee8010
done with 13665 octs
dumping data/grid.00010.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 0 data/grid.00010.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=4.208333e+01
0oct=0x7ff2a3ee8010
done with 13665 octs
dumping data/grid.00010.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 707 data/grid.00010.f707 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=4.208333e+01
0oct=0x7ff2a3ee8010
done with 13665 octs
dumping data/grid.00010.f708.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 708 data/grid.00010.f708 1 \
0"
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> ws
> pli,log10(e(,,64))
> pli,x(,,64)
> palette,"idl-33.gp"
> i=12;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707);n=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 128x128x128 cube from data/grid.00012.p00000
size= 16777224
tsim=5.041667e+01
0oct=0x7ff2a3ee8010
done with 14329 octs
dumping data/grid.00012.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00012 7 701 data/grid.00012.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00012.p00000
size= 16777224
tsim=5.041667e+01
0oct=0x7ff2a3ee8010
done with 14329 octs
dumping data/grid.00012.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00012 7 706 data/grid.00012.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00012.p00000
size= 16777224
tsim=5.041667e+01
0oct=0x7ff2a3ee8010
done with 14329 octs
dumping data/grid.00012.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00012 7 0 data/grid.00012.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00012.p00000
size= 16777224
tsim=5.041667e+01
0oct=0x7ff2a3ee8010
done with 14329 octs
dumping data/grid.00012.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00012 7 707 data/grid.00012.f707 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00012.p00000
size= 16777224
tsim=5.041667e+01
0oct=0x7ff2a3ee8010
done with 14329 octs
dumping data/grid.00012.f708.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00012 7 708 data/grid.00012.f708 1 \
0"
> pli,x(,,64)
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> ws
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> ws
> pli,x(,,64)
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> i=14;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707);n=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 128x128x128 cube from data/grid.00014.p00000
size= 16777224
tsim=5.875000e+01
0oct=0x7ff2a3ee8010
done with 15649 octs
dumping data/grid.00014.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00014 7 701 data/grid.00014.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00014.p00000
size= 16777224
tsim=5.875000e+01
0oct=0x7ff2a3ee8010
done with 15649 octs
dumping data/grid.00014.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00014 7 706 data/grid.00014.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00014.p00000
size= 16777224
tsim=5.875000e+01
0oct=0x7ff2a3ee8010
done with 15649 octs
dumping data/grid.00014.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00014 7 0 data/grid.00014.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00014.p00000
size= 16777224
tsim=5.875000e+01
0oct=0x7ff2a3ee8010
done with 15649 octs
dumping data/grid.00014.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00014 7 707 data/grid.00014.f707 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00014.p00000
size= 16777224
tsim=5.875000e+01
0oct=0x7ff2a3ee8010
done with 15649 octs
dumping data/grid.00014.f708.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00014 7 708 data/grid.00014.f708 1 \
0"
> ws
> pli,x(,,64)
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> pli,log10(e(,,64))
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> pli,log10(e(,,64))
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> i=16;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707);n=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 128x128x128 cube from data/grid.00016.p00000
size= 16777224
tsim=6.708333e+01
0oct=0x7ff2a3ee8010
done with 16521 octs
dumping data/grid.00016.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00016 7 701 data/grid.00016.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00016.p00000
size= 16777224
tsim=6.708333e+01
0oct=0x7ff2a3ee8010
done with 16521 octs
dumping data/grid.00016.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00016 7 706 data/grid.00016.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00016.p00000
size= 16777224
tsim=6.708333e+01
0oct=0x7ff2a3ee8010
done with 16521 octs
dumping data/grid.00016.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00016 7 0 data/grid.00016.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00016.p00000
size= 16777224
tsim=6.708333e+01
0oct=0x7ff2a3ee8010
done with 16521 octs
dumping data/grid.00016.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00016 7 707 data/grid.00016.f707 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00016.p00000
size= 16777224
tsim=6.708333e+01
0oct=0x7ff2a3ee8010
done with 16521 octs
dumping data/grid.00016.f708.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00016 7 708 data/grid.00016.f708 1 \
0"
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> ws
> pli,log10(e(,,64))
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> ws
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> i=2;lvl=8;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707);n=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 256x256x256 cube from data/grid.00002.p00000
size= 134217736
tsim=2.100000e+01
0oct=0x7f4214566010
done with 44417 octs
dumping data/grid.00002.f701.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 8 701 data/grid.00002.f701 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00002.p00000
size= 134217736
tsim=2.100000e+01
0oct=0x7f4214566010
done with 44417 octs
dumping data/grid.00002.f706.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 8 706 data/grid.00002.f706 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00002.p00000
size= 134217736
tsim=2.100000e+01
0oct=0x7f4214566010
done with 44417 octs
dumping data/grid.00002.f0.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 8 0 data/grid.00002.f0 1 0"
Casting Rays on 256x256x256 cube from data/grid.00002.p00000
size= 134217736
tsim=2.100000e+01
0oct=0x7f4214566010
done with 44417 octs
dumping data/grid.00002.f707.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 8 707 data/grid.00002.f707 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00002.p00000
size= 134217736
tsim=2.100000e+01
0oct=0x7f4214566010
done with 44417 octs
dumping data/grid.00002.f708.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 8 708 data/grid.00002.f708 1 \
0"
> pli,log10(e(,,64))
> limits
> pli,log10(e(,,128))
> pli,x(,,128)
> plotamr(l(,,128))
8
7
6
5
4
3
2
1
[]
> i=2;lvl=9;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707);n=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 512x512x512 cube from data/grid.00002.p00000
size= 1073741832
tsim=2.100000e+01
0oct=0x7f4214566010
done with 56053 octs
dumping data/grid.00002.f701.p00000 with nmap=134217728
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 9 701 data/grid.00002.f701 1 \
0"
Casting Rays on 512x512x512 cube from data/grid.00002.p00000
size= 1073741832
tsim=2.100000e+01
0oct=0x7f4214566010
done with 56053 octs
dumping data/grid.00002.f706.p00000 with nmap=134217728
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 9 706 data/grid.00002.f706 1 \
0"
Casting Rays on 512x512x512 cube from data/grid.00002.p00000
size= 1073741832
tsim=2.100000e+01
0oct=0x7f4214566010
done with 56053 octs
dumping data/grid.00002.f0.p00000 with nmap=134217728
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 9 0 data/grid.00002.f0 1 0"
Casting Rays on 512x512x512 cube from data/grid.00002.p00000
size= 1073741832
tsim=2.100000e+01
0oct=0x7f4214566010
done with 56053 octs
dumping data/grid.00002.f707.p00000 with nmap=134217728
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 9 707 data/grid.00002.f707 1 \
0"
Casting Rays on 512x512x512 cube from data/grid.00002.p00000
size= 1073741832
tsim=2.100000e+01
0oct=0x7f4214566010
done with 56053 octs
dumping data/grid.00002.f708.p00000 with nmap=134217728
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 9 708 data/grid.00002.f708 1 \
0"
> ws
> pli,x(,,256)
> plotamr(l(,,256))
9
8
7
6
5
4
3
2
1
[]
> t=[]
> n=[]
> i=6;lvl=9;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 512x512x512 cube from data/grid.00006.p00000
size= 1073741832
tsim=6.100000e+01
0oct=0x7f4214566010
done with 92869 octs
dumping data/grid.00006.f701.p00000 with nmap=134217728
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 9 701 data/grid.00006.f701 1 \
0"
Casting Rays on 512x512x512 cube from data/grid.00006.p00000
size= 1073741832
tsim=6.100000e+01
0oct=0x7f4214566010
done with 92869 octs
dumping data/grid.00006.f706.p00000 with nmap=134217728
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 9 706 data/grid.00006.f706 1 \
0"
Casting Rays on 512x512x512 cube from data/grid.00006.p00000
size= 1073741832
tsim=6.100000e+01
0oct=0x7f4214566010
done with 92869 octs
dumping data/grid.00006.f0.p00000 with nmap=134217728
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 9 0 data/grid.00006.f0 1 0"
> ws
> pli,x(,,256)
> plotamr(l(,,128))
6
5
4
3
2
1
[]
> ws
> plotamr(l(,,256))
9
8
7
6
5
4
3
2
1
[]
> pli,x(,,256)
> pli,x(,,256)
> plotamr(l(,,256))
9
8
7
6
5
4
3
2
1
[]
> i=9;lvl=9;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 512x512x512 cube from data/grid.00009.p00000
size= 1073741832
tsim=9.100000e+01
0oct=0x7f4214566010
done with 125321 octs
dumping data/grid.00009.f701.p00000 with nmap=134217728
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00009 9 701 data/grid.00009.f701 1 \
0"
Casting Rays on 512x512x512 cube from data/grid.00009.p00000
size= 1073741832
tsim=9.100000e+01
0oct=0x7f4214566010
done with 125321 octs
dumping data/grid.00009.f706.p00000 with nmap=134217728
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00009 9 706 data/grid.00009.f706 1 \
0"
Casting Rays on 512x512x512 cube from data/grid.00009.p00000
size= 1073741832
tsim=9.100000e+01
0oct=0x7f4214566010
done with 125321 octs
dumping data/grid.00009.f0.p00000 with nmap=134217728
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00009 9 0 data/grid.00009.f0 1 0"
> ws
> pli,x(,,256)
> plotamr(l(,,256))
9
8
7
6
5
4
3
2
1
[]
> ws
> pli,log10(e(,,256))
> i=11;lvl=9;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 512x512x512 cube from data/grid.00011.p00000
size= 1073741832
tsim=1.110000e+02
0oct=0x7f4214566010
done with 150237 octs
dumping data/grid.00011.f701.p00000 with nmap=134217728
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00011 9 701 data/grid.00011.f701 1 \
0"
Casting Rays on 512x512x512 cube from data/grid.00011.p00000
size= 1073741832
tsim=1.110000e+02
0oct=0x7f4214566010
done with 150237 octs
dumping data/grid.00011.f706.p00000 with nmap=134217728
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00011 9 706 data/grid.00011.f706 1 \
0"
Casting Rays on 512x512x512 cube from data/grid.00011.p00000
size= 1073741832
tsim=1.110000e+02
0oct=0x7f4214566010
done with 150237 octs
dumping data/grid.00011.f0.p00000 with nmap=134217728
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00011 9 0 data/grid.00011.f0 1 0"
> pli,log10(e(,,256))
> plotamr(l(,,256))
9
8
7
6
5
4
3
2
1
[]
> pli,log10(x(,,256))
> plotamr(l(,,256))
9
8
7
6
5
4
3
2
1
[]
> pli,log10(1.-x(,,256))
> ws
> pli,log10(1.-x(,,256))
> plotamr(l(,,256))
9
8
7
6
5
4
3
2
1
[]
> i=14;lvl=9;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 512x512x512 cube from data/grid.00014.p00000
size= 1073741832
tsim=1.410000e+02
0oct=0x7f4214566010
done with 187657 octs
dumping data/grid.00014.f701.p00000 with nmap=134217728
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00014 9 701 data/grid.00014.f701 1 \
0"
Casting Rays on 512x512x512 cube from data/grid.00014.p00000
size= 1073741832
tsim=1.410000e+02
0oct=0x7f4214566010
done with 187657 octs
dumping data/grid.00014.f706.p00000 with nmap=134217728
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00014 9 706 data/grid.00014.f706 1 \
0"
Casting Rays on 512x512x512 cube from data/grid.00014.p00000
size= 1073741832
tsim=1.410000e+02
0oct=0x7f4214566010
done with 187657 octs
dumping data/grid.00014.f0.p00000 with nmap=134217728
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00014 9 0 data/grid.00014.f0 1 0"
> ws
> pli,log10(1.-x(,,256))
> plotamr(l(,,256))
9
8
7
6
5
4
3
2
1
[]
> ws
> pli,log10(e(,,256))
> plotamr(l(,,256))
9
8
7
6
5
4
3
2
1
[]
> i=15;lvl=9;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 512x512x512 cube from data/grid.00015.p00000
size= 1073741832
tsim=1.510000e+02
0oct=0x7f4214566010
done with 200537 octs
dumping data/grid.00015.f701.p00000 with nmap=134217728
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00015 9 701 data/grid.00015.f701 1 \
0"
Casting Rays on 512x512x512 cube from data/grid.00015.p00000
size= 1073741832
tsim=1.510000e+02
0oct=0x7f4214566010
done with 200537 octs
dumping data/grid.00015.f706.p00000 with nmap=134217728
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00015 9 706 data/grid.00015.f706 1 \
0"
Casting Rays on 512x512x512 cube from data/grid.00015.p00000
size= 1073741832
tsim=1.510000e+02
0oct=0x7f4214566010
done with 200537 octs
dumping data/grid.00015.f0.p00000 with nmap=134217728
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00015 9 0 data/grid.00015.f0 1 0"
> ws1
> pli,log10(e(,,256))
> plotamr(l(,,256))
9
8
7
6
5
4
3
2
1
[]
> palette,"rainbow.gp"
> ws
> pli,log10(1.-x(,,256))
> plotamr(l(,,256))
9
8
7
6
5
4
3
2
1
[]
> 3.*(100e3/3.085677e22)^2/(8.*pi*6.67384e-11)
1.87847e-26
> 3.*(100e3/3.085677e22)^2/(8.*pi*6.67384e-11)*0.7^2
9.20451e-27
> 3.*(100e3/3.085677e22)^2/(8.*pi*6.67384e-11)*0.74^2
1.02865e-26
> 3.*(100e3/3.085677e22)^2/(8.*pi*6.67384e-11)*0.74^2*0.3
3.08595e-27
> 3.*(100e3/3.085677e22)^2/(8.*pi*6.67384e-11)*0.74^2*0.045
4.62893e-28
> 3.*(100e3/3.085677e22)^2/(8.*pi*6.67384e-11)*0.74^2*0.045/1.67e-27
0.277181
> 3.*(100e3/3.085677e22)^2/(8.*pi*6.67384e-11)*0.7^2*0.045/1.67e-27
0.248026
> 3.*(100e3/3.085677e22)^2/(8.*pi*6.67384e-11)*0.74^2*0.045/1.67e-27
0.277181
> 3.*(100e3/3.085677e22)^2/(8.*pi*6.67384e-11)*0.74^2*0.045/1.67e-27*0.75
0.207886
> 3.*(100e3/3.085677e22)^2/(8.*pi*6.67384e-11)*0.74^2*0.045
4.62893e-28
> 3.*(100e3/3.085677e22)^2/(8.*pi*6.67384e-11)*0.74^2*0.045*(100*3.085677e22)^3
1.35998e+46
> 1.67e-23/3.*(100e3/3.085677e22)^2/(8.*pi*6.67384e-11)*0.74^2*0.045*(100*3.085677e22)^3
2.52352e+22
> 1.67e-23/(3.*(100e3/3.085677e22)^2/(8.*pi*6.67384e-11)*0.74^2*0.045*(100*3.085677e22)^3)
1.22796e-69
> 3.*(100e3/3.085677e22)^2/(8.*pi*6.67384e-11)
1.87847e-26
> i=0;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f5da806b010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 701 data/grid.00000.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f5da806b010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 706 data/grid.00000.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f5da806b010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 0 data/grid.00000.f0 1 0"
> i=0;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f5da806b010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 701 data/grid.00000.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f5da806b010
  C-c"~/Project/Quartz/utils/oct2grid data/grid.00000 6 706 data/grid.00000.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f5da806b010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 0 data/grid.00000.f0 1 0"
> i=0;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f5da806b010
done with 37449 octs
dumping data/grid.00000.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 701 data/grid.00000.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f5da806b010
done with 37449 octs
dumping data/grid.00000.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 706 data/grid.00000.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f5da806b010
done with 37449 octs
dumping data/grid.00000.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 0 data/grid.00000.f0 1 0"
> ws1
> info,e
 array(double,64,64,64)
> i=15;lvl=9;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 512x512x512 cube from data/grid.00015.p00000
size= 1073741832
tsim=1.510000e+02
0oct=0x7f4214566010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00015 9 701 data/grid.00015.f701 1 \
0"
Casting Rays on 512x512x512 cube from data/grid.00015.p00000
size= 1073741832
tsim=1.510000e+02
0oct=0x7f4214566010
  C-c"~/Project/Quartz/utils/oct2grid data/grid.00015 9 706 data/grid.00015.f706 1 \
0"
  C-cERROR (reform) Keyboard interrupt received (SIGINT)
  LINE: 272  FILE: /usr/lib/yorick/i/Eric/utils.i
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> i=15;lvl=9;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 512x512x512 cube from data/grid.00015.p00000
size= 1073741832
tsim=1.510000e+02
0oct=0x7f4214566010
done with 200537 octs
dumping data/grid.00015.f701.p00000 with nmap=134217728
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00015 9 701 data/grid.00015.f701 1 \
0"
Casting Rays on 512x512x512 cube from data/grid.00015.p00000
size= 1073741832
tsim=1.510000e+02
0oct=0x7f4214566010
done with 200537 octs
dumping data/grid.00015.f706.p00000 with nmap=134217728
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00015 9 706 data/grid.00015.f706 1 \
0"
Casting Rays on 512x512x512 cube from data/grid.00015.p00000
size= 1073741832
tsim=1.510000e+02
0oct=0x7f4214566010
done with 200537 octs
dumping data/grid.00015.f0.p00000 with nmap=134217728
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00015 9 0 data/grid.00015.f0 1 0"
> pli,log10(1.-x(,,256))
> ws1
> pli,log10(1.-x(,,256))
> 
> palette,"rainbow.gp"
> ws1,1
> plotamr(l(,,256))
9
8
7
6
5
4
3
2
1
[]
> ws1,2
> pli,log10(e(,,256))
> palette,"rainbow.gp"
> ws1,3
> pli,log10(1.-x(,,256))
> rvp
> plotamr(l(,,256))
9
8
7
6
5
4
3
2
1
[]
> palette,"rainbow.gp"
> wkl
> ws1
> wkl
> ws
> i=0;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0);
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f5da806b010
done with 37449 octs
dumping data/grid.00000.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 701 data/grid.00000.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f5da806b010
done with 37449 octs
dumping data/grid.00000.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 706 data/grid.00000.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f5da806b010
done with 37449 octs
dumping data/grid.00000.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 0 data/grid.00000.f0 1 0"
> INFO,x
 1: array(double,64,64,64) min=0 max=0 avg=0 std=0
> info,oct2grid
 []
> i=6;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);
Casting Rays on 64x64x64 cube from data/grid.00006.p00000
size= 2097160
tsim=1.488135e-01
0oct=0x7f9afe896010
done with 37449 octs
dumping data/grid.00006.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 6 701 data/grid.00006.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00006.p00000
size= 2097160
tsim=1.488135e-01
0oct=0x7f9afe896010
done with 37449 octs
dumping data/grid.00006.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 6 706 data/grid.00006.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00006.p00000
size= 2097160
tsim=1.488135e-01
0oct=0x7f9afe896010
done with 37449 octs
dumping data/grid.00006.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 6 101 data/grid.00006.f101 1 \
0"
> ws1
> info,d
 array(double,64,64,64)
> INFO,d
 1: array(double,64,64,64) min=0.00368437 max=1.48746 avg=0.045 std=0.0332946
> pli,d(,,32)
> #include "utils/readpart.i"
> info,readpart
 func readpart(fname,&time)
> i=6;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);p=readpart(swrite(format="data/part.%05d",i));
Casting Rays on 64x64x64 cube from data/grid.00006.p00000
size= 2097160
tsim=1.488135e-01
0oct=0x7f9afe896010
done with 37449 octs
dumping data/grid.00006.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 6 701 data/grid.00006.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00006.p00000
size= 2097160
tsim=1.488135e-01
0oct=0x7f9afe896010
done with 37449 octs
dumping data/grid.00006.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 6 706 data/grid.00006.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00006.p00000
size= 2097160
tsim=1.488135e-01
0oct=0x7f9afe896010
done with 37449 octs
dumping data/grid.00006.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 6 101 data/grid.00006.f101 1 \
0"
ERROR (readpart) cannot open file data/part.00006 (mode rb)
  LINE: 6  FILE: /home/daubert/Project/Quartz/utils/readpart.i
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> i=6;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);p=readpart(swrite(format="data/part.%05d.p00000",i));
Casting Rays on 64x64x64 cube from data/grid.00006.p00000
size= 2097160
tsim=1.488135e-01
0oct=0x7f9afe896010
done with 37449 octs
dumping data/grid.00006.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 6 701 data/grid.00006.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00006.p00000
size= 2097160
tsim=1.488135e-01
0oct=0x7f9afe896010
done with 37449 octs
dumping data/grid.00006.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 6 706 data/grid.00006.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00006.p00000
size= 2097160
tsim=1.488135e-01
0oct=0x7f9afe896010
done with 37449 octs
dumping data/grid.00006.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 6 101 data/grid.00006.f101 1 \
0"
 found 262144 particles
> pl,p(2,),p(1,)
> ws
> pl,p(2,),p(1,)
> ws
> pl,p(2,),p(1,)
> pli,d(,,32)
> pl,p(2,),p(1,)
> ws
> pli,d(,,avg)
> pl,p(2,),p(1,)
> ws
> pli,d(,,avg),0,1,0,1
> pli,d(,,avg),0,0,1,1
> pl,p(2,),p(1,)
> INFO,d
 1: array(double,64,64,64) min=0.00368437 max=1.48746 avg=0.045 std=0.0332946
> INFO,x
 1: array(double,64,64,64) min=0.0012 max=0.0012 avg=0.0012 std=6.79795e-16
> i=13;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);p=readpart(swrite(format="data/part.%05d.p00000",i));
Casting Rays on 64x64x64 cube from data/grid.00013.p00000
size= 2097160
tsim=2.516389e-01
0oct=0x7f9afe896010
done with 37449 octs
dumping data/grid.00013.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00013 6 701 data/grid.00013.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00013.p00000
size= 2097160
tsim=2.516389e-01
0oct=0x7f9afe896010
done with 37449 octs
dumping data/grid.00013.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00013 6 706 data/grid.00013.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00013.p00000
size= 2097160
tsim=2.516389e-01
0oct=0x7f9afe896010
done with 37449 octs
dumping data/grid.00013.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00013 6 101 data/grid.00013.f101 1 \
0"
 found 262144 particles
> INFO,d
 1: array(double,64,64,64) min=0.00225843 max=4.38328 avg=0.045 std=0.0696811
> i=13;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);p=readpart(swrite(format="data/part.%05d.p00000",i));
Casting Rays on 64x64x64 cube from data/grid.00013.p00000
size= 2097160
tsim=2.516389e-01
0oct=0x7f9afe896010
done with 37449 octs
dumping data/grid.00013.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00013 6 701 data/grid.00013.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00013.p00000
size= 2097160
tsim=2.516389e-01
0oct=0x7f9afe896010
done with 37449 octs
dumping data/grid.00013.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00013 6 706 data/grid.00013.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00013.p00000
size= 2097160
tsim=2.516389e-01
0oct=0x7f9afe896010
done with 37449 octs
dumping data/grid.00013.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00013 6 101 data/grid.00013.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00013.p00000
size= 2097160
tsim=2.516389e-01
0oct=0x7f9afe896010
done with 37449 octs
dumping data/grid.00013.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00013 6 1 data/grid.00013.f1 1 0"
 found 262144 particles
> INFO,dt
 1: array(double,64,64,64) min=-0.997038 max=114.391 avg=-1.48978e-08 std=1.53649
> i=6;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);p=readpart(swrite(format="data/part.%05d.p00000",i));
Casting Rays on 64x64x64 cube from data/grid.00006.p00000
size= 2097160
tsim=1.488135e-01
0oct=0x7fb37f09f010
done with 37449 octs
dumping data/grid.00006.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 6 701 data/grid.00006.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00006.p00000
size= 2097160
tsim=1.488135e-01
0oct=0x7fb37f09f010
done with 37449 octs
dumping data/grid.00006.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 6 706 data/grid.00006.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00006.p00000
size= 2097160
tsim=1.488135e-01
0oct=0x7fb37f09f010
done with 37449 octs
dumping data/grid.00006.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 6 101 data/grid.00006.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00006.p00000
size= 2097160
tsim=1.488135e-01
0oct=0x7fb37f09f010
done with 37449 octs
dumping data/grid.00006.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 6 1 data/grid.00006.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00006.p00000
size= 2097160
tsim=1.488135e-01
0oct=0x7fb37f09f010
done with 37449 octs
dumping data/grid.00006.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 6 705 data/grid.00006.f705 1 \
0"
 found 262144 particles
> INFO,src
 1: array(double,64,64,64) min=0 max=0 avg=0 std=0
> INFO,e
 1: array(double,64,64,64) min=4.21112e+62 max=1.39107e+65 avg=5.56079e+63 std=3.60729e+63
> i=6;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);p=readpart(swrite(format="data/part.%05d.p00000",i));
Casting Rays on 64x64x64 cube from data/grid.00006.p00000
size= 2097160
tsim=1.488135e-01
0oct=0x7fb37f09f010
done with 37449 octs
dumping data/grid.00006.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 6 701 data/grid.00006.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00006.p00000
size= 2097160
tsim=1.488135e-01
0oct=0x7fb37f09f010
done with 37449 octs
dumping data/grid.00006.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 6 706 data/grid.00006.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00006.p00000
size= 2097160
tsim=1.488135e-01
0oct=0x7fb37f09f010
done with 37449 octs
dumping data/grid.00006.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 6 101 data/grid.00006.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00006.p00000
size= 2097160
tsim=1.488135e-01
0oct=0x7fb37f09f010
done with 37449 octs
dumping data/grid.00006.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 6 1 data/grid.00006.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00006.p00000
size= 2097160
tsim=1.488135e-01
0oct=0x7fb37f09f010
done with 37449 octs
dumping data/grid.00006.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 6 705 data/grid.00006.f705 1 \
0"
 found 262144 particles
> ws1
> INFO,src
 1: array(double,64,64,64) min=0 max=0 avg=0 std=0
> pli,e(,,80)
ERROR (*main*) index overreach beyond array bounds
WARNING source code unavailable (try dbdis function)
now at pc= 1 (of 19), failed at pc= 9
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> INFO,e
 1: array(double,64,64,64) min=4.21112e+62 max=1.39107e+65 avg=5.56079e+63 std=3.60729e+63
> pli,e(,,max)
> i=7;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);p=readpart(swrite(format="data/part.%05d.p00000",i));
Casting Rays on 64x64x64 cube from data/grid.00007.p00000
size= 2097160
tsim=1.656215e-01
0oct=0x7fb37f09f010
done with 37449 octs
dumping data/grid.00007.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00007 6 701 data/grid.00007.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00007.p00000
size= 2097160
tsim=1.656215e-01
0oct=0x7fb37f09f010
done with 37449 octs
dumping data/grid.00007.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00007 6 706 data/grid.00007.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00007.p00000
size= 2097160
tsim=1.656215e-01
0oct=0x7fb37f09f010
done with 37449 octs
dumping data/grid.00007.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00007 6 101 data/grid.00007.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00007.p00000
size= 2097160
tsim=1.656215e-01
0oct=0x7fb37f09f010
done with 37449 octs
dumping data/grid.00007.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00007 6 1 data/grid.00007.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00007.p00000
size= 2097160
tsim=1.656215e-01
0oct=0x7fb37f09f010
done with 37449 octs
dumping data/grid.00007.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00007 6 705 data/grid.00007.f705 1 \
0"
 found 262144 particles
> INFO,e
 1: array(double,64,64,64) min=4.07065e+62 max=1.49421e+65 avg=5.56707e+63 std=3.6637e+63
> pli,e(,,max)
> ws
> pli,e(,,avg)
> i=0;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);p=readpart(swrite(format="data/part.%05d.p00000",i));
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fb37f09f010
done with 37449 octs
dumping data/grid.00000.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 701 data/grid.00000.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fb37f09f010
done with 37449 octs
dumping data/grid.00000.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 706 data/grid.00000.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fb37f09f010
done with 37449 octs
dumping data/grid.00000.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 101 data/grid.00000.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fb37f09f010
done with 37449 octs
dumping data/grid.00000.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 1 data/grid.00000.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fb37f09f010
done with 37449 octs
dumping data/grid.00000.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 705 data/grid.00000.f705 1 \
0"
 found 262144 particles
> INFO,e
 1: array(double,64,64,64) min=0 max=0 avg=0 std=0
> INFO,src
 1: array(double,64,64,64) min=0 max=0 avg=0 std=0
> i=0;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);p=readpart(swrite(format="data/part.%05d.p00000",i));
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fb37f09f010
done with 37449 octs
dumping data/grid.00000.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 701 data/grid.00000.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fb37f09f010
done with 37449 octs
dumping data/grid.00000.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 706 data/grid.00000.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fb37f09f010
done with 37449 octs
dumping data/grid.00000.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 101 data/grid.00000.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fb37f09f010
done with 37449 octs
dumping data/grid.00000.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 1 data/grid.00000.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fb37f09f010
done with 37449 octs
dumping data/grid.00000.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 705 data/grid.00000.f705 1 \
0"
 found 262144 particles
> INFO,d
 1: array(double,64,64,64) min=0.0120465 max=0.0854171 avg=0.045 std=0.00822227
> i=6;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);p=readpart(swrite(format="data/part.%05d.p00000",i));
Casting Rays on 64x64x64 cube from data/grid.00006.p00000
size= 2097160
tsim=1.488135e-01
0oct=0x7fb37f09f010
done with 37449 octs
dumping data/grid.00006.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 6 701 data/grid.00006.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00006.p00000
size= 2097160
tsim=1.488135e-01
0oct=0x7fb37f09f010
done with 37449 octs
dumping data/grid.00006.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 6 706 data/grid.00006.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00006.p00000
size= 2097160
tsim=1.488135e-01
0oct=0x7fb37f09f010
done with 37449 octs
dumping data/grid.00006.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 6 101 data/grid.00006.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00006.p00000
size= 2097160
tsim=1.488135e-01
0oct=0x7fb37f09f010
done with 37449 octs
dumping data/grid.00006.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 6 1 data/grid.00006.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00006.p00000
size= 2097160
tsim=1.488135e-01
0oct=0x7fb37f09f010
done with 37449 octs
dumping data/grid.00006.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 6 705 data/grid.00006.f705 1 \
0"
 found 262144 particles
> INFO,src
 1: array(double,64,64,64) min=0 max=0 avg=0 std=0
> INFO,e
 1: array(double,64,64,64) min=4.21112e+62 max=1.39107e+65 avg=5.56079e+63 std=3.60729e+63
> i=0;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f9f8ea0e010
done with 37449 octs
dumping data/grid.00000.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 701 data/grid.00000.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f9f8ea0e010
done with 37449 octs
dumping data/grid.00000.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 706 data/grid.00000.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f9f8ea0e010
done with 37449 octs
dumping data/grid.00000.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 101 data/grid.00000.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f9f8ea0e010
done with 37449 octs
dumping data/grid.00000.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 1 data/grid.00000.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f9f8ea0e010
done with 37449 octs
dumping data/grid.00000.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 705 data/grid.00000.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f9f8ea0e010
done with 37449 octs
dumping data/grid.00000.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 708 data/grid.00000.f708 1 \
0"
> INFO,nh
 1: array(double,64,64,64) min=0 max=0 avg=0 std=0
> i=0;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 701 data/grid.00000.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 706 data/grid.00000.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 101 data/grid.00000.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 1 data/grid.00000.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 705 data/grid.00000.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 708 data/grid.00000.f708 1 \
0"
> INFO,src
 1: array(double,64,64,64) min=0 max=0 avg=0 std=0
> i=1;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00001.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 701 data/grid.00001.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00001.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 706 data/grid.00001.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00001.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 101 data/grid.00001.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00001.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 1 data/grid.00001.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00001.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 705 data/grid.00001.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00001.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 708 data/grid.00001.f708 1 \
0"
> INFO,src
 1: array(double,64,64,64) min=0 max=0 avg=0 std=0
> i=1;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00001.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 701 data/grid.00001.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00001.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 706 data/grid.00001.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00001.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 101 data/grid.00001.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00001.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 1 data/grid.00001.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00001.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 705 data/grid.00001.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00001.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 708 data/grid.00001.f708 1 \
0"
> INFO,src
 1: array(double,64,64,64) min=0 max=0 avg=0 std=0
> i=0;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 701 data/grid.00000.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 706 data/grid.00000.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 101 data/grid.00000.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 1 data/grid.00000.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 705 data/grid.00000.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 708 data/grid.00000.f708 1 \
0"
> INFO,src
 1: array(double,64,64,64) min=0 max=0 avg=0 std=0
> i=0;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 701 data/grid.00000.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 706 data/grid.00000.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 101 data/grid.00000.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 1 data/grid.00000.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 705 data/grid.00000.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 708 data/grid.00000.f708 1 \
0"
> INFO,src
 1: array(double,64,64,64) min=0 max=2.31112e+77 avg=7.05298e+72 std=1.27671e+75
> ws1
> pli,src(,,max)
> pli,nh(,,max)
> pli,src(,,max)
> i=0;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 701 data/grid.00000.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 706 data/grid.00000.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 101 data/grid.00000.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 1 data/grid.00000.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 705 data/grid.00000.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 708 data/grid.00000.f708 1 \
0"
> INFO,src
 1: array(double,64,64,64) min=0 max=2.31112e+77 avg=7.05298e+72 std=1.27671e+75
> i=0;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 701 data/grid.00000.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 706 data/grid.00000.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 101 data/grid.00000.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 1 data/grid.00000.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 705 data/grid.00000.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 708 data/grid.00000.f708 1 \
0"
> INFO,src
 1: array(double,64,64,64) min=0 max=2.31112e+77 avg=7.05298e+72 std=1.27671e+75
> i=0;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 701 data/grid.00000.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 706 data/grid.00000.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 101 data/grid.00000.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 1 data/grid.00000.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 705 data/grid.00000.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 708 data/grid.00000.f708 1 \
0"
> INFO,src
 1: array(double,64,64,64) min=0 max=2.31112e+77 avg=7.05298e+72 std=1.27671e+75
> ws
> pli,e(,,avg)
> INFO,e
 1: array(double,64,64,64) min=0 max=0 avg=0 std=0
> INFO,x
 1: array(double,64,64,64) ERROR (_stat_worker) Floating point interrupt (SIGFPE)
  LINE: 916  FILE: /usr/lib/yorick/i/Eric/utils.i
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> INFO,nh
 1: array(double,64,64,64) min=5.3943e+71 max=3.82489e+72 avg=2.01505e+72 std=3.68185e+71
> i=0;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 701 data/grid.00000.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 706 data/grid.00000.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 101 data/grid.00000.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 1 data/grid.00000.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 705 data/grid.00000.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 708 data/grid.00000.f708 1 \
0"
> INFO,e
 1: array(double,64,64,64) min=0 max=0 avg=0 std=0
> INFO,xion
 1: void, []
> INFO,x
 1: array(double,64,64,64) min=0 max=0 avg=0 std=0
> i=0;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffef2a0010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 701 data/grid.00000.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffef2a0010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 706 data/grid.00000.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffef2a0010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 101 data/grid.00000.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffef2a0010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 1 data/grid.00000.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffef2a0010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 705 data/grid.00000.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffef2a0010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 708 data/grid.00000.f708 1 \
0"
> size= 2097160
> i=0;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffef2a0010
done with 37449 octs
dumping data/grid.00000.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 701 data/grid.00000.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffef2a0010
done with 37449 octs
dumping data/grid.00000.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 706 data/grid.00000.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffef2a0010
done with 37449 octs
dumping data/grid.00000.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 101 data/grid.00000.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffef2a0010
done with 37449 octs
dumping data/grid.00000.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 1 data/grid.00000.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffef2a0010
done with 37449 octs
dumping data/grid.00000.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 705 data/grid.00000.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffef2a0010
done with 37449 octs
dumping data/grid.00000.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 708 data/grid.00000.f708 1 \
0"
> ws
> INFO,src
 1: array(double,64,64,64) min=0 max=2.31112e+77 avg=7.05298e+72 std=1.27671e+75
> INFO,e
 1: array(double,64,64,64) min=0 max=3.08726e+76 avg=9.42158e+71 std=1.70546e+74
> pli,e(,,avg)
> i=1;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7fffef2a0010
done with 37449 octs
dumping data/grid.00001.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 701 data/grid.00001.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7fffef2a0010
done with 37449 octs
dumping data/grid.00001.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 706 data/grid.00001.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7fffef2a0010
done with 37449 octs
dumping data/grid.00001.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 101 data/grid.00001.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7fffef2a0010
done with 37449 octs
dumping data/grid.00001.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 1 data/grid.00001.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7fffef2a0010
done with 37449 octs
dumping data/grid.00001.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 705 data/grid.00001.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7fffef2a0010
done with 37449 octs
dumping data/grid.00001.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 708 data/grid.00001.f708 1 \
0"
> pli,e(,,avg)
> pli,log(e(,,max))
ERROR (*main*) Floating point interrupt (SIGFPE)
WARNING source code unavailable (try dbdis function)
now at pc= 1 (of 23), failed at pc= 13
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> INFO,e
 1: array(double,64,64,64) min=0 max=5.91094e+76 avg=4.55857e+72 std=4.3022e+74
> i=5;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=4.513571e-02
0oct=0x7fffef2a0010
done with 37449 octs
dumping data/grid.00005.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 701 data/grid.00005.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=4.513571e-02
0oct=0x7fffef2a0010
done with 37449 octs
dumping data/grid.00005.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 706 data/grid.00005.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=4.513571e-02
0oct=0x7fffef2a0010
done with 37449 octs
dumping data/grid.00005.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 101 data/grid.00005.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=4.513571e-02
0oct=0x7fffef2a0010
done with 37449 octs
dumping data/grid.00005.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 1 data/grid.00005.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=4.513571e-02
0oct=0x7fffef2a0010
done with 37449 octs
dumping data/grid.00005.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 705 data/grid.00005.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=4.513571e-02
0oct=0x7fffef2a0010
done with 37449 octs
dumping data/grid.00005.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 708 data/grid.00005.f708 1 \
0"
> pli,log(e(,,max))
ERROR (*main*) Floating point interrupt (SIGFPE)
WARNING source code unavailable (try dbdis function)
now at pc= 1 (of 23), failed at pc= 13
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> pli,e(,,avg)
> pli,e(,,avg)
> pli,e(,,max)
> i=10;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=5.462054e-02
0oct=0x7fffef2a0010
done with 37449 octs
dumping data/grid.00010.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 701 data/grid.00010.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=5.462054e-02
0oct=0x7fffef2a0010
done with 37449 octs
dumping data/grid.00010.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 706 data/grid.00010.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=5.462054e-02
0oct=0x7fffef2a0010
done with 37449 octs
dumping data/grid.00010.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 101 data/grid.00010.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=5.462054e-02
0oct=0x7fffef2a0010
done with 37449 octs
dumping data/grid.00010.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 1 data/grid.00010.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=5.462054e-02
0oct=0x7fffef2a0010
done with 37449 octs
dumping data/grid.00010.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 705 data/grid.00010.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=5.462054e-02
0oct=0x7fffef2a0010
done with 37449 octs
dumping data/grid.00010.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 708 data/grid.00010.f708 1 \
0"
> pli,e(,,max)
> info,x
 array(double,64,64,64)
> pli,x(,,256)
ERROR (*main*) index overreach beyond array bounds
WARNING source code unavailable (try dbdis function)
now at pc= 1 (of 19), failed at pc= 9
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> pli,x(,,max)
> pli,e(,,max)
> pli,log(e(,,max)+1e-10)
> pli,log(e(,,max)+1e-20)
> palette,"idl-33.gp"
> rvp
> ws
> pli,e(,,max)
> ws
> i=10;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=5.462054e-02
0oct=0x7fffef2a0010
done with 37449 octs
dumping data/grid.00010.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 701 data/grid.00010.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=5.462054e-02
0oct=0x7fffef2a0010
done with 37449 octs
dumping data/grid.00010.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 706 data/grid.00010.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=5.462054e-02
0oct=0x7fffef2a0010
done with 37449 octs
dumping data/grid.00010.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 101 data/grid.00010.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=5.462054e-02
0oct=0x7fffef2a0010
done with 37449 octs
dumping data/grid.00010.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 1 data/grid.00010.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=5.462054e-02
0oct=0x7fffef2a0010
done with 37449 octs
dumping data/grid.00010.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 705 data/grid.00010.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=5.462054e-02
0oct=0x7fffef2a0010
done with 37449 octs
dumping data/grid.00010.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 708 data/grid.00010.f708 1 \
0"
> pli,e(,,max)
> i=0;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 701 data/grid.00000.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
  C-c"~/Project/Quartz/utils/oct2grid data/grid.00000 6 706 data/grid.00000.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 101 data/grid.00000.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 1 data/grid.00000.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 705 data/grid.00000.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 708 data/grid.00000.f708 1 \
0"
> i=0;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 701 data/grid.00000.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 706 data/grid.00000.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 101 data/grid.00000.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 1 data/grid.00000.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 705 data/grid.00000.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 708 data/grid.00000.f708 1 \
0"
> ws
> INFO,x
 1: array(double,64,64,64) min=0 max=0 avg=0 std=0
> INFO,e
 1: array(double,64,64,64) min=0 max=0 avg=0 std=0
> i=0;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 701 data/grid.00000.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 706 data/grid.00000.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 101 data/grid.00000.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 1 data/grid.00000.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 705 data/grid.00000.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00000.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 708 data/grid.00000.f708 1 \
0"
> ws
> info,x
 array(double,64,64,64)
> INFO,x
 1: array(double,64,64,64) min=0.00120006 max=0.00169802 avg=0.00120025 std=2.70566e-06
> pli,x(,,max)
> i=1;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00001.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 701 data/grid.00001.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00001.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 706 data/grid.00001.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00001.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 101 data/grid.00001.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00001.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 1 data/grid.00001.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00001.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 705 data/grid.00001.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00001.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 708 data/grid.00001.f708 1 \
0"
> ws
> pli,x(,,max)
> plk,x(,,max)
> INFO?x
SYNTAX: syntax error near x
SYNTAX: syntax error near <EOF>
> INFO,x
 1: array(double,64,64,64) min=0.00120012 max=0.00227348 avg=0.00120052 std=8.20972e-06
> ws
> i=4;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00004.p00000
size= 2097160
tsim=4.324466e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00004.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 6 701 data/grid.00004.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00004.p00000
size= 2097160
tsim=4.324466e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00004.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 6 706 data/grid.00004.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00004.p00000
size= 2097160
tsim=4.324466e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00004.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 6 101 data/grid.00004.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00004.p00000
size= 2097160
tsim=4.324466e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00004.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 6 1 data/grid.00004.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00004.p00000
size= 2097160
tsim=4.324466e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00004.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 6 705 data/grid.00004.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00004.p00000
size= 2097160
tsim=4.324466e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00004.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 6 708 data/grid.00004.f708 1 \
0"
> pli,x(,,max)
> i=5;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=4.513571e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00005.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 701 data/grid.00005.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=4.513571e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00005.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 706 data/grid.00005.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=4.513571e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00005.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 101 data/grid.00005.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=4.513571e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00005.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 1 data/grid.00005.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=4.513571e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00005.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 705 data/grid.00005.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=4.513571e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00005.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 708 data/grid.00005.f708 1 \
0"
> pli,x(,,max)
> INFO,x
 1: array(double,64,64,64) min=0.00120027 max=0.00547878 avg=0.00120485 std=9.76767e-05
> i=1;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00001.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 701 data/grid.00001.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00001.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 706 data/grid.00001.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00001.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 101 data/grid.00001.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00001.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 1 data/grid.00001.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00001.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 705 data/grid.00001.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00001.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 708 data/grid.00001.f708 1 \
0"
> INFO,x
 1: array(double,64,64,64) min=0.00120012 max=0.999848 avg=0.00131243 std=0.0100052
> ws
> pli,x(,,max)
> i=2;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00002.p00000
size= 2097160
tsim=3.947236e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00002.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 6 701 data/grid.00002.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00002.p00000
size= 2097160
tsim=3.947236e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00002.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 6 706 data/grid.00002.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00002.p00000
size= 2097160
tsim=3.947236e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00002.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 6 101 data/grid.00002.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00002.p00000
size= 2097160
tsim=3.947236e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00002.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 6 1 data/grid.00002.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00002.p00000
size= 2097160
tsim=3.947236e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00002.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 6 705 data/grid.00002.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00002.p00000
size= 2097160
tsim=3.947236e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00002.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 6 708 data/grid.00002.f708 1 \
0"
> pli,x(,,max)
> i=5;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=4.513571e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00005.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 701 data/grid.00005.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=4.513571e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00005.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 706 data/grid.00005.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=4.513571e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00005.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 101 data/grid.00005.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=4.513571e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00005.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 1 data/grid.00005.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=4.513571e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00005.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 705 data/grid.00005.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=4.513571e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00005.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 708 data/grid.00005.f708 1 \
0"
> pli,x(,,max)
> pli,x(,,avg)
> pli,x(,,max)
> info,nh
 array(double,64,64,64)
> pli,nh(,,max)
> pli,(nh*(1.-x))(,,max)
> pli,(nh*x))(,,max)
SYNTAX: syntax error near )(,,max)
SYNTAX: syntax error near <EOF>
> 
> pli,(nh*x)(,,max)
> ws
> pli,x(,,max)
> i=10;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=5.462054e-02
0oct=0x7fffef2a0010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 701 data/grid.00010.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=5.462054e-02
0oct=0x7fffef2a0010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 706 data/grid.00010.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=5.462054e-02
0oct=0x7fffef2a0010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 101 data/grid.00010.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=5.462054e-02
0oct=0x7fffef2a0010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 1 data/grid.00010.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=5.462054e-02
0oct=0x7fffef2a0010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 705 data/grid.00010.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=5.462054e-02
0oct=0x7fffef2a0010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 708 data/grid.00010.f708 1 \
0"
> i=9;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00009.p00000
size= 2097160
tsim=5.271985e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00009.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00009 6 701 data/grid.00009.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00009.p00000
size= 2097160
tsim=5.271985e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00009.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00009 6 706 data/grid.00009.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00009.p00000
size= 2097160I
tsim=5.271985e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00009.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00009 6 101 data/grid.00009.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00009.p00000
size= 2097160
tsim=5.271985e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00009.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00009 6 1 data/grid.00009.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00009.p00000
size= 2097160
tsim=5.271985e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00009.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00009 6 705 data/grid.00009.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00009.p00000
size= 2097160
tsim=5.271985e-02
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00009.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00009 6 708 data/grid.00009.f708 1 \
0"
> pli,x(,,max)
> pli,x(,,avg)
> INFO?d
SYNTAX: syntax error near d
SYNTAX: syntax error near <EOF>
> INFO,d
 1: array(double,64,64,64) min=0.00916048 max=0.14405 avg=0.045 std=0.0114803
> i=15;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00015.p00000
size= 2097160
tsim=2.768507e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00015.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00015 6 701 data/grid.00015.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00015.p00000
size= 2097160
tsim=2.768507e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00015.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00015 6 706 data/grid.00015.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00015.p00000
size= 2097160
tsim=2.768507e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00015.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00015 6 101 data/grid.00015.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00015.p00000
size= 2097160
tsim=2.768507e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00015.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00015 6 1 data/grid.00015.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00015.p00000
size= 2097160
tsim=2.768507e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00015.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00015 6 705 data/grid.00015.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00015.p00000
size= 2097160
tsim=2.768507e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00015.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00015 6 708 data/grid.00015.f708 1 \
0"
> pli,x(,,avg)
> i=20;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00020.p00000
size= 2097160
tsim=3.375011e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00020.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00020 6 701 data/grid.00020.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00020.p00000
size= 2097160
tsim=3.375011e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00020.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00020 6 706 data/grid.00020.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00020.p00000
size= 2097160
tsim=3.375011e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00020.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00020 6 101 data/grid.00020.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00020.p00000
size= 2097160
tsim=3.375011e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00020.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00020 6 1 data/grid.00020.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00020.p00000
size= 2097160
tsim=3.375011e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00020.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00020 6 705 data/grid.00020.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00020.p00000
size= 2097160
tsim=3.375011e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00020.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00020 6 708 data/grid.00020.f708 1 \
0"
> pli,x(,,avg)
> i=40;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00040.p00000
size= 2097160
tsim=5.327715e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00040.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00040 6 701 data/grid.00040.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00040.p00000
size= 2097160
tsim=5.327715e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00040.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00040 6 706 data/grid.00040.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00040.p00000
size= 2097160
tsim=5.327715e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00040.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00040 6 101 data/grid.00040.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00040.p00000
size= 2097160
tsim=5.327715e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00040.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00040 6 1 data/grid.00040.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00040.p00000
size= 2097160
tsim=5.327715e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00040.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00040 6 705 data/grid.00040.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00040.p00000
size= 2097160
tsim=5.327715e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00040.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00040 6 708 data/grid.00040.f708 1 \
0"
> pli,x(,,avg)
> i=80;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00080.p00000
size= 2097160
tsim=9.943900e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00080.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00080 6 701 data/grid.00080.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00080.p00000
size= 2097160
tsim=9.943900e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00080.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00080 6 706 data/grid.00080.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00080.p00000
size= 2097160
tsim=9.943900e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00080.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00080 6 101 data/grid.00080.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00080.p00000
size= 2097160
tsim=9.943900e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00080.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00080 6 1 data/grid.00080.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00080.p00000
size= 2097160
tsim=9.943900e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00080.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00080 6 705 data/grid.00080.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00080.p00000
size= 2097160
tsim=9.943900e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00080.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00080 6 708 data/grid.00080.f708 1 \
0"
> pli,x(,,avg)
> ws
> INFO,x
 1: array(double,64,64,64) min=0.00120005 max=1 avg=0.221354 std=0.402377
> pli,x(,,max)
> pli,x(,,avg)
> pli,log(1.-x(,,avg))
> plk,log(1.-x(,,avg))
> ws
> pli,x(,,avg)
> pli,x(,,64)
> pli,x(,,avg)
> pli,x(,,64)
> pli,nh(,,64)
> pli,log10(nh(,,64))
> pli,x(,,64)
> pli,log10(nh(,,64))
> pli,x(,,64)
> pli,x(,,64)
> INFO,x
 1: array(double,64,64,64) min=0.00120005 max=1 avg=0.221354 std=0.402377
> pli,x(,,32)
> pli,nh(,,32)
> pli,log10(nh(,,32))
> pli,log10(nh(,,32)*(1.-x)(,,32))
> ws1
> pli,log10(nh(,,32)*(1.-x)(,,32))
> plk,log10(nh(,,32)*(1.-x)(,,32))
> ws
> pli,log10(nh(,,32))
> palette,"idl-33.gp"
> plk,log10(nh(,,32)*(1.-x)(,,32))
> pli,nh(,,32)
> pli,nh(,,32)*x(,,32)
> pli,nh(,,32)
> pli,nh(,,32)*x(,,32)
> pli,nh(,,32)
> pli,nh(,,32)*x(,,32)
> pli,nh(,,32)
> pli,nh(,,32)*x(,,32)
> pli,log10(nh(,,32))
> ws1,0
> pli,x(,,32)
> ws1,1
> pli,log10(nh(,,32))
> pli,src(,,32)
> ws1
> pli,src(,,32)
> ws
> pli,src(,,32)
> ws
> pli,src(,,12)
> ws
> ws1
> pli,x(,,32)
> palette,"idl-33.gp"
> pli,log10(nh(,,32))
> plcont,x(,,32)
> help,plc
 builtin plc()
 /* DOCUMENT plc, z, y, x, levs=z_values
          or plc, z, y, x, ireg, levs=z_values
          or plc, z, levs=z_values
      plots a contours of Z on the mesh Y versus X.  Y, X, and IREG are
      as for plm.  The Z array must have the same shape as Y and X.
      The function being contoured takes the value Z at each point
      (X,Y) -- that is, the Z array is presumed to be point-centered.
      The Y, X, and IREG arguments may all be omitted to default to the
      mesh set by the most recent plmesh call.
      The LEVS keyword is a list of the values of Z at which you want
      contour curves.  The default is eight contours spanning the
      range of Z.
      See plfc if you want to color the regions between contours.
      The following keywords are legal (each has a separate help entry):
    KEYWORDS: legend, hide
              type, width, color, smooth
              marks, marker, mspace, mphase
              smooth, triangle, region
    SEE ALSO: plg, plm, plc, plv, plf, pli, plt, pldj, plfp, plmesh, plfc
              contour, spann, limits, logxy, range, fma, hcp
  */
 defined at:  LINE: 732  FILE: /usr/lib/yorick/i0/graph.i
> ws
> pli,log10(nh(,,32))
> plc,x(,,32)
ERROR (*main*) no default mesh exists to define y and x -- use plmesh
WARNING source code unavailable (try dbdis function)
now at pc= 1 (of 19), failed at pc= 13
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> plc,x(,,32),span(0,1,64)(-:1:64,)
> 
> ws1
> pli,log10(nh(,,30))
> palette,"idl-33.gp"
> plc,x(,,30),span(0,64,64)(-:1:64,),span(0,64,64)(,-:1:64)
> i=70;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00070.p00000
size= 2097160
tsim=8.731168e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00070.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00070 6 701 data/grid.00070.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00070.p00000
size= 2097160
tsim=8.731168e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00070.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00070 6 706 data/grid.00070.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00070.p00000
size= 2097160
tsim=8.731168e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00070.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00070 6 101 data/grid.00070.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00070.p00000
size= 2097160
tsim=8.731168e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00070.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00070 6 1 data/grid.00070.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00070.p00000
size= 2097160
tsim=8.731168e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00070.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00070 6 705 data/grid.00070.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00070.p00000
size= 2097160
tsim=8.731168e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00070.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00070 6 708 data/grid.00070.f708 1 \
0"
> ws
> pli,x(,,32)
> INFO,x
 1: array(double,64,64,64) min=0.00120002 max=1 avg=0.160699 std=0.354003
> ws
> plotamr(l(,,256))
ERROR (*main*) index overreach beyond array bounds
WARNING source code unavailable (try dbdis function)
now at pc= 1 (of 18), failed at pc= 9
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> pli,nh(,,32)*x(,,32)
> pli,log10(nh(,,32)*x(,,32))
> pli,log10(nh(,,32)*(1.-x)(,,32))
> pli,log10(nh(,,10)*(1.-x)(,,10))
> plk,log10(nh(,,10)*(1.-x)(,,10))
> ws1
> plk,log10(nh(,,10)*(1.-x)(,,10))
> plk,log10(nh(,,)*(1.-x)(,,))(,,avg)
> plk,log10(nh(,,)*(1.-x)(,,))(,,max)
> plk,log10(nh(,,)*(1.-x)(,,))(,,avg)
> ws
> plk,log10(nh(,,10)*(1.-x)(,,10))
> palette,"idl-33.gp"
> rvp
> rvp
> rvp
> rvp
> rvp
> i=74;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00074.p00000
size= 2097160
tsim=9.210266e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00074.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00074 6 701 data/grid.00074.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00074.p00000
size= 2097160
tsim=9.210266e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00074.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00074 6 706 data/grid.00074.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00074.p00000
size= 2097160
tsim=9.210266e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00074.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00074 6 101 data/grid.00074.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00074.p00000
size= 2097160
tsim=9.210266e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00074.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00074 6 1 data/grid.00074.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00074.p00000
size= 2097160
tsim=9.210266e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00074.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00074 6 705 data/grid.00074.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00074.p00000
size= 2097160
tsim=9.210266e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00074.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00074 6 708 data/grid.00074.f708 1 \
0"
> plk,log10(nh(,,10)*(1.-x)(,,10))
> plk,log10(nh(,,64)*(1.-x)(,,64))
> plk,log10(nh(,12,)*(1.-x)(,,64))
> plk,log10(nh(,12,)*(1.-x)(,12,))
> i=13;plk,log10(nh(i,,)*(1.-x)(i,,))
> i=12;plk,log10(nh(i,,)*(1.-x)(i,,))
> i=40;plk,log10(nh(i,,)*(1.-x)(i,,))
> i=50;plk,log10(nh(i,,)*(1.-x)(i,,))
> i=2;plk,log10(nh(i,,)*(1.-x)(i,,))
> i=2;plk,log10(nh(i,,)*(1.-x)(avg,,))
> i=2;plk,log10(nh(i,,)*(1.-x)(,avg,))
> i=10;plk,log10(nh(i,,)*(1.-x)(,avg,))
> i=20;plk,log10(nh(i,,)*(1.-x)(,avg,))
> i=30;plk,log10(nh(i,,)*(1.-x)(i,,))
> i=30;plk,log10(nh(,,)*(1.-x)(,,))(avg,,)
> i=75;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00075.p00000
size= 2097160
tsim=9.331907e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00075.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00075 6 701 data/grid.00075.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00075.p00000
size= 2097160
tsim=9.331907e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00075.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00075 6 706 data/grid.00075.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00075.p00000
size= 20971601./
tsim=9.331907e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00075.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00075 6 101 data/grid.00075.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00075.p00000
size= 2097160
tsim=9.331907e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00075.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00075 6 1 data/grid.00075.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00075.p00000
size= 2097160
tsim=9.331907e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00075.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00075 6 705 data/grid.00075.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00075.p00000
size= 2097160
tsim=9.331907e-01
0oct=0x7fffee050010
done with 37449 octs
dumping data/grid.00075.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00075 6 708 data/grid.00075.f708 1 \
0"
> plk,log10(nh(,12,)*(1.-x)(,12,))
> plk,smooth(log10(nh(,12,)*(1.-x)(,12,)))
> size= 20971601./
cont> 
cont> qsd
ERROR (*main*) bad data type(s) in binary /
WARNING source code unavailable (try dbdis function)
now at pc= 1 (of 12), failed at pc= 5
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> 
dbug> debxit
[]
dbug> dbexit
> 10./64^3
3.8147e-05
> d=oct2cube("data/grid.00008",6,1)
Casting Rays on 64x64x64 cube from data/grid.00008.p00000
size= 2097160
tsim=5.082127e-02
0oct=0x7fa0f0e22010
done with 37449 octs
dumping data/grid.00008.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00008 6 1 data/grid.00008.f1 1 0"
> INFO,s
 1: void, []
> INFO,d
 1: array(double,64,64,64) min=-0.62958 max=2.23098 avg=-1.48978e-08 std=0.232425
> d=oct2cube("data/grid.00006",7,101)
Casting Rays on 128x128x128 cube from data/grid.00006.p00000
size= 16777224
tsim=1.479818e-01
0oct=0x7f77fd316010
done with 38102 octs
dumping data/grid.00006.f101.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 7 101 data/grid.00006.f101 1 \
0"
> ws
> pli,d(,,max)
> pli,log(d(,,max))
> l=oct2cube("data/grid.00006",7,0)
Casting Rays on 128x128x128 cube from data/grid.00006.p00000
size= 16777224
tsim=1.479818e-01
0oct=0x7f77fd316010
done with 38102 octs
dumping data/grid.00006.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 7 0 data/grid.00006.f0 1 0"
> plotamr(l(,,max))
7
6
5
4
3
2
1
[]
> ws
> pli,log(d(,,avg))
> plotamr(l(,,max))
7
6
5
4
3
2
1
[]
> ws1
> pli,log(d(,,))(,,avg)
> rvp
> pli,log(d(,,))(,,max)
> plotamr(l(,,max))
7
6
5
4
3
2
1
[]
> d=oct2cube("data/grid.00010",7,101)
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=1.999441e-01
0oct=0x7f77fd316010
done with 42988 octs
dumping data/grid.00010.f101.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 101 data/grid.00010.f101 1 \
0"
> ws
> pli,d(,,max)
> pli,log(d(,,))(,,max)
> l=oct2cube("data/grid.00010",7,0)
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=1.999441e-01
0oct=0x7f77fd316010
done with 42988 octs
dumping data/grid.00010.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 0 data/grid.00010.f0 1 0"
> plotamr(l(,,max))
7
6
5
4
3
2
1
[]
> rvp
> palette,"gray.gp"
> rvp
> ws
> pli,log(d(,,))(,,max)
> www=where2(l==7)
> info,www
 array(long,3,44312)
> INFO,www(3,)
 1: array(long,44312) min=1 max=128 avg=68.3218 std=31.9992
> h=histo1d(www(3,),span(0,128,128))
> h(mxx)
67
> ws
> pli,d(,,67)
> plotamr(l(,,67))
7
6
5
4
3
2
1
[]
> ws
> pli,log10(d(,,67))
> plotamr(l(,,67))
7
6
5
4
3
2
1
[]
> i=5;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=1.316201e-01
0oct=0x7f0e87c6c010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 701 data/grid.00005.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=1.316201e-01
0oct=0x7f0e87c6c010
  C-c"~/Project/Quartz/utils/oct2grid data/grid.00005 6 706 data/grid.00005.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=1.316201e-01
0oct=0x7f0e87c6c010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 101 data/grid.00005.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=1.316201e-01
0oct=0x7f0e87c6c010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 1 data/grid.00005.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=1.316201e-01
0oct=0x7f0e87c6c010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 705 data/grid.00005.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=1.316201e-01
0oct=0x7f0e87c6c010
  C-c"~/Project/Quartz/utils/oct2grid data/grid.00005 6 708 data/grid.00005.f708 1 \
0"
> i=5;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=1.316201e-01
0oct=0x7f0e87c6c010
done with 37449 octs
dumping data/grid.00005.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 701 data/grid.00005.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=1.316201e-01
0oct=0x7f0e87c6c010
done with 37449 octs
dumping data/grid.00005.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 706 data/grid.00005.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=1.316201e-01
0oct=0x7f0e87c6c010
done with 37449 octs
dumping data/grid.00005.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 101 data/grid.00005.f101 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=1.316201e-01
0oct=0x7f0e87c6c010
done with 37449 octs
dumping data/grid.00005.f1.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 1 data/grid.00005.f1 1 0"
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=1.316201e-01
0oct=0x7f0e87c6c010
done with 37449 octs
dumping data/grid.00005.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 705 data/grid.00005.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=1.316201e-01
0oct=0x7f0e87c6c010
done with 37449 octs
dumping data/grid.00005.f708.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 708 data/grid.00005.f708 1 \
0"
> ws1
> pli,nh(,,64)
> pli,nh(,,max)
> www=where2(l==7)
> INFO,www
 1: array(long,3,44312) min=1 max=128 avg=64.8952 std=35.0117
> h=histo1d(www(3,),span(0,128,128))
> h(mxx)
67
> ws
> pli,nh(,,67)
ERROR (*main*) index overreach beyond array bounds
WARNING source code unavailable (try dbdis function)
now at pc= 1 (of 19), failed at pc= 9
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> i=5;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 128x128x128 cube from data/grid.00005.p00000
size= 16777224
tsim=1.316201e-01
0oct=0x7f0e87c6c010
done with 37602 octs
dumping data/grid.00005.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 7 701 data/grid.00005.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00005.p00000
size= 16777224
tsim=1.316201e-01
0oct=0x7f0e87c6c010
done with 37602 octs
dumping data/grid.00005.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 7 706 data/grid.00005.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00005.p00000
size= 16777224
tsim=1.316201e-01
0oct=0x7f0e87c6c010
done with 37602 octs
dumping data/grid.00005.f101.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 7 101 data/grid.00005.f101 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00005.p00000
size= 16777224
tsim=1.316201e-01
0oct=0x7f0e87c6c010
done with 37602 octs
dumping data/grid.00005.f1.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 7 1 data/grid.00005.f1 1 0"
Casting Rays on 128x128x128 cube from data/grid.00005.p00000
size= 16777224
tsim=1.316201e-01
0oct=0x7f0e87c6c010
done with 37602 octs
dumping data/grid.00005.f705.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 7 705 data/grid.00005.f705 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00005.p00000
size= 16777224
tsim=1.316201e-01
0oct=0x7f0e87c6c010
done with 37602 octs
dumping data/grid.00005.f708.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 7 708 data/grid.00005.f708 1 \
0"
> www=where2(l==7)
> h=histo1d(www(3,),span(0,128,128))
> h(mxx)
67
> pli,nh(,,67)
> plotamr(l(,,67))
7
6
5
4
3
2
1
[]
> pli,x(,,67)
> pli,e(,,67)
> pli,log10(e(,,67))
> plotamr(l(,,67))
7
6
5
4
3
2
1
[]
> ws
> pli,x(,,67)
> pli,log(x(,,67))
> pli,log(1.-x(,,67))
> plotamr(l(,,67))
7
6
5
4
3
2
1
[]
> ws
> pli,log(1.-x(,,67))
> palette,"gray.gp"
> palette,"idl-33.gp"
> pli,x(,,max)
> INFO,x
 1: array(double,128,128,128) min=0.00120003 max=1 avg=0.0102601 std=0.0921594
> i=4;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=1.130159e-01
0oct=0x7fffcbaff010
done with 37476 octs
dumping data/grid.00004.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 701 data/grid.00004.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=1.130159e-01
0oct=0x7fffcbaff010
done with 37476 octs
dumping data/grid.00004.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 706 data/grid.00004.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=1.130159e-01
0oct=0x7fffcbaff010
done with 37476 octs
dumping data/grid.00004.f101.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 101 data/grid.00004.f101 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=1.130159e-01
0oct=0x7fffcbaff010
done with 37476 octs
dumping data/grid.00004.f1.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 1 data/grid.00004.f1 1 0"
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=1.130159e-01
0oct=0x7fffcbaff010
done with 37476 octs
dumping data/grid.00004.f705.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 705 data/grid.00004.f705 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=1.130159e-01
0oct=0x7fffcbaff010
ws
done with 37476 octs
dumping data/grid.00004.f708.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 708 data/grid.00004.f708 1 \
0"
> www=where2(l==7)
> h=histo1d(www(3,),span(0,128,128))
> h(mxx)
67
> pli,log(1.-x(,,67))
> plotamr(l(,,67))
7
6
5
4
3
2
1
[]
> info,nh
 array(double,128,128,128)
> pli,nh(,,67)
> pli,nh(,,67)
> plc,x(,,67),span(0,128,128)(-:1:128,),span(0,128,128)(,-:1:128)
> INFO,nh
 1: array(double,128,128,128) min=2.14168e+71 max=3.15379e+73 avg=2.01505e+72 std=1.06828e+72
> i=5;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 128x128x128 cube from data/grid.00005.p00000
size= 16777224
tsim=1.316201e-01
0oct=0x7fffcbaff010
done with 37602 octs
dumping data/grid.00005.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 7 701 data/grid.00005.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00005.p00000
size= 16777224
tsim=1.316201e-01
0oct=0x7fffcbaff010
done with 37602 octs
dumping data/grid.00005.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 7 706 data/grid.00005.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00005.p00000
size= 16777224
tsim=1.316201e-01
0oct=0x7fffcbaff010
done with 37602 octs
dumping data/grid.00005.f101.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 7 101 data/grid.00005.f101 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00005.p00000
size= 16777224
tsim=1.316201e-01
0oct=0x7fffcbaff010
done with 37602 octs
dumping data/grid.00005.f1.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 7 1 data/grid.00005.f1 1 0"
Casting Rays on 128x128x128 cube from data/grid.00005.p00000
size= 16777224
tsim=1.316201e-01
0oct=0x7fffcbaff010
done with 37602 octs
dumping data/grid.00005.f705.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 7 705 data/grid.00005.f705 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00005.p00000
size= 16777224
tsim=1.316201e-01
0oct=0x7fffcbaff010
done with 37602 octs
dumping data/grid.00005.f708.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 7 708 data/grid.00005.f708 1 \
0"
> ws1
> www=where2(l==7)
> h=histo1d(www(3,),span(0,128,128))
> h(mxx)
67
> pli,nh(,,67)
> pli,x(,,max)
> ws
> pli,log(1.-x(,,67))
> plotamr(l(,,67))
7
6
5
4
3
2
1
[]
> rvp
> rvp
> rvp
> rvp
> plotamr(l(,,67))
7
6
5
4
3
2
1
[]
> ws
> plotamr(l(,,67))
7
6
5
4
3
2
1
[]
> pli,x(,,67)
> plotamr(l(,,67))
7
6
5
4
3
2
1
[]
> palette,"idl-33.gp"
> pli,log(1.-x(,,67))
> plotamr(l(,,67))
7
6
5
4
3
2
1
[]
> i=0;lvl=5;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 32x32x32 cube from data/grid.00000.p00000
size= 262152
tsim=9.090909e-02
0oct=0x7fffe84c3010
done with 4681 octs
dumping data/grid.00000.f701.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 5 701 data/grid.00000.f701 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00000.p00000
size= 262152
tsim=9.090909e-02
0oct=0x7fffe84c3010
done with 4681 octs
dumping data/grid.00000.f706.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 5 706 data/grid.00000.f706 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00000.p00000
size= 262152
tsim=9.090909e-02
0oct=0x7fffe84c3010
done with 4681 octs
dumping data/grid.00000.f101.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 5 101 data/grid.00000.f101 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00000.p00000
size= 262152
tsim=9.090909e-02
0oct=0x7fffe84c3010
done with 4681 octs
dumping data/grid.00000.f1.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 5 1 data/grid.00000.f1 1 0"
Casting Rays on 32x32x32 cube from data/grid.00000.p00000
size= 262152
tsim=9.090909e-02
0oct=0x7fffe84c3010
done with 4681 octs
dumping data/grid.00000.f705.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 5 705 data/grid.00000.f705 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00000.p00000
size= 262152
tsim=9.090909e-02
0oct=0x7fffe84c3010
done with 4681 octs
dumping data/grid.00000.f708.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 5 708 data/grid.00000.f708 1 \
0"
> info,src
 array(double,32,32,32)
> ws1
> pli,src(,,16)
> pli,e(,,16)
> i=0;lvl=5;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 32x32x32 cube from data/grid.00000.p00000
size= 262152
tsim=9.090909e-02
0oct=0x7f25193f8010
done with 4681 octs
dumping data/grid.00000.f701.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 5 701 data/grid.00000.f701 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00000.p00000
size= 262152
tsim=9.090909e-02
0oct=0x7f25193f8010
done with 4681 octs
dumping data/grid.00000.f706.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 5 706 data/grid.00000.f706 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00000.p00000
size= 262152
tsim=9.090909e-02
0oct=0x7f25193f8010
done with 4681 octs
dumping data/grid.00000.f101.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 5 101 data/grid.00000.f101 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00000.p00000
size= 262152
tsim=9.090909e-02
0oct=0x7f25193f8010
done with 4681 octs
dumping data/grid.00000.f1.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 5 1 data/grid.00000.f1 1 0"
Casting Rays on 32x32x32 cube from data/grid.00000.p00000
size= 262152
tsim=9.090909e-02
0oct=0x7f25193f8010
done with 4681 octs
dumping data/grid.00000.f705.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 5 705 data/grid.00000.f705 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00000.p00000
size= 262152
tsim=9.090909e-02
0oct=0x7f25193f8010
done with 4681 octs
dumping data/grid.00000.f708.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 5 708 data/grid.00000.f708 1 \
0"
> ws
> pli,e(,,16)
> i=1;lvl=5;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 32x32x32 cube from data/grid.00001.p00000
size= 262152
tsim=9.508573e-02
0oct=0x7f25193f8010
done with 4681 octs
dumping data/grid.00001.f701.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 5 701 data/grid.00001.f701 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00001.p00000
size= 262152
tsim=9.508573e-02
0oct=0x7f25193f8010
done with 4681 octs
dumping data/grid.00001.f706.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 5 706 data/grid.00001.f706 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00001.p00000
size= 262152
tsim=9.508573e-02
0oct=0x7f25193f8010
done with 4681 octs
dumping data/grid.00001.f101.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 5 101 data/grid.00001.f101 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00001.p00000
size= 262152
tsim=9.508573e-02
0oct=0x7f25193f8010
done with 4681 octs
dumping data/grid.00001.f1.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 5 1 data/grid.00001.f1 1 0"
Casting Rays on 32x32x32 cube from data/grid.00001.p00000
size= 262152
tsim=9.508573e-02
0oct=0x7f25193f8010
done with 4681 octs
dumping data/grid.00001.f705.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 5 705 data/grid.00001.f705 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00001.p00000
size= 262152
tsim=9.508573e-02
0oct=0x7f25193f8010
done with 4681 octs
dumping data/grid.00001.f708.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 5 708 data/grid.00001.f708 1 \
0"
> pli,e(,,16)
> ws
> pli,e(,,16)
> i=10;lvl=5;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101);dt=oct2cube(swrite(format="data/grid.%05d",i),lvl,1);src=oct2cube(swrite(format="data/grid.%05d",i),lvl,705);nh=oct2cube(swrite(format="data/grid.%05d",i),lvl,708);
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=1.367939e-01
0oct=0x7fc6ffd1e010
done with 4681 octs
dumping data/grid.00010.f701.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 701 data/grid.00010.f701 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=1.367939e-01
0oct=0x7fc6ffd1e010
done with 4681 octs
dumping data/grid.00010.f706.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 706 data/grid.00010.f706 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=1.367939e-01
0oct=0x7fc6ffd1e010
done with 4681 octs
dumping data/grid.00010.f101.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 101 data/grid.00010.f101 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=1.367939e-01
0oct=0x7fc6ffd1e010
done with 4681 octs
dumping data/grid.00010.f1.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 1 data/grid.00010.f1 1 0"
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=1.367939e-01
0oct=0x7fc6ffd1e010
done with 4681 octs
dumping data/grid.00010.f705.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 705 data/grid.00010.f705 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=1.367939e-01
0oct=0x7fc6ffd1e010
done with 4681 octs
dumping data/grid.00010.f708.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 708 data/grid.00010.f708 1 \
0"
> ws
> pli,e(,,16)
> pli,f(,,16)
ERROR (*main*) attempt to call non-function or index non-array
WARNING source code unavailable (try dbdis function)
now at pc= 1 (of 19), failed at pc= 9
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> i=10;lvl=5;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702);
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=1.367939e-01
0oct=0x7fc6ffd1e010
done with 4681 octs
dumping data/grid.00010.f701.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 701 data/grid.00010.f701 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=1.367939e-01
0oct=0x7fc6ffd1e010
done with 4681 octs
dumping data/grid.00010.f702.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 702 data/grid.00010.f702 1 \
0"
> ws
> pli,fx(,,16)
> ws
> plg,e(,16,16)
> plg,e(,16,16)/fx(,16,16)
ERROR (*main*) Floating point interrupt (SIGFPE)
WARNING source code unavailable (try dbdis function)
now at pc= 1 (of 30), failed at pc= 21
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> ws
> plg,e(,16,16)/(fx(,16,16)+1e-3)
> i=1;lvl=5;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702);
Casting Rays on 32x32x32 cube from data/grid.00001.p00000
size= 262152
tsim=9.508573e-02
0oct=0x7fc6ffd1e010
done with 4681 octs
dumping data/grid.00001.f701.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 5 701 data/grid.00001.f701 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00001.p00000
size= 262152
tsim=9.508573e-02
0oct=0x7fc6ffd1e010
done with 4681 octs
dumping data/grid.00001.f702.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 5 702 data/grid.00001.f702 1 \
0"
> ws
> plg,fx(,16,16)
> plg,e(,16,16)
> logxy,0,1
> plg,e(,16,16),color="red"
> i=1;lvl=5;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);
Casting Rays on 32x32x32 cube from data/grid.00001.p00000
size= 262152
tsim=9.508573e-02
0oct=0x7fc6ffd1e010
done with 4681 octs
dumping data/grid.00001.f701.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 5 701 data/grid.00001.f701 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00001.p00000
size= 262152
tsim=9.508573e-02
0oct=0x7fc6ffd1e010
done with 4681 octs
dumping data/grid.00001.f702.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 5 702 data/grid.00001.f702 1 \
0"
> a
0.0950857
> ws
> i=10;lvl=5;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702);
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=1.367939e-01
0oct=0x7f2337bca010
done with 4681 octs
dumping data/grid.00010.f701.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 701 data/grid.00010.f701 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=1.367939e-01
0oct=0x7f2337bca010
done with 4681 octs
dumping data/grid.00010.f702.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 702 data/grid.00010.f702 1 \
0"
> i=10;lvl=5;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=1.367939e-01
0oct=0x7f2337bca010
done with 4681 octs
dumping data/grid.00010.f701.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 701 data/grid.00010.f701 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=1.367939e-01
0oct=0x7f2337bca010
done with 4681 octs
dumping data/grid.00010.f702.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 702 data/grid.00010.f702 1 \
0"
> a
0.136794
> plg,e(,16,16)/a^3,color="red"
> ws
> plg,e(,16,16)*0.08511784,color="red"
> ws
> plg,e(,16,16)*0.08511784,color="red"
> plg,fx(,16,16)
> logxy,0,1
> i=50;lvl=5;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);
Casting Rays on 32x32x32 cube from data/grid.00050.p00000
size= 262152
tsim=4.114215e-01
0oct=0x7f0c95ee7010
done with 4681 octs
dumping data/grid.00050.f701.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 5 701 data/grid.00050.f701 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00050.p00000
size= 262152
tsim=4.114215e-01
0oct=0x7f0c95ee7010
done with 4681 octs
dumping data/grid.00050.f702.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 5 702 data/grid.00050.f702 1 \
0"
> ws
> plg,e(,16,16)*0.08511784,color="red"
> plg,fx(,16,16)
> a
0.411422
> i=100;lvl=5;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);
Casting Rays on 32x32x32 cube from data/grid.00100.p00000
size= 262152
tsim=9.585990e-01
0oct=0x7f0c95ee7010
done with 4681 octs
dumping data/grid.00100.f701.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00100 5 701 data/grid.00100.f701 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00100.p00000
size= 262152
tsim=9.585990e-01
0oct=0x7f0c95ee7010
done with 4681 octs
dumping data/grid.00100.f702.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00100 5 702 data/grid.00100.f702 1 \
0"
> plg,e(,16,16)*0.5964,color="red"
> ws
> plg,e(,16,16)*0.5964,color="red"
> plg,fx(,16,16)
> logxy,0,1
> i=100;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);
Casting Rays on 64x64x64 cube from data/grid.00100.p00000
size= 2097160
tsim=4.081296e-01
0oct=0x7fb749981010
done with 37449 octs
dumping data/grid.00100.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00100 6 701 data/grid.00100.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00100.p00000
size= 2097160
tsim=4.081296e-01
0oct=0x7fb749981010
done with 37449 octs
dumping data/grid.00100.f702.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00100 6 702 data/grid.00100.f702 1 \
0"
> ws
> plg,fx(,16,16)
> plg,e(,16,16)*0.2539,color="red"
> logxy,0,1
> i=100;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);
Casting Rays on 64x64x64 cube from data/grid.00100.p00000
size= 2097160
tsim=4.081296e-01
0oct=0x7fb749981010
done with 12483 octs
dumping data/grid.00100.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00100 6 701 data/grid.00100.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00100.p00000
size= 2097160
tsim=4.081296e-01
0oct=0x7fb749981010
done with 12483 octs
dumping data/grid.00100.f702.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00100 6 702 data/grid.00100.f702 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00100.p00000
size= 2097160
tsim=4.081296e-01
0oct=0x7fb749981010
done with 12483 octs
dumping data/grid.00100.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00100 6 706 data/grid.00100.f706 1 \
0"
> ws
> i=100;lvl=5;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);
Casting Rays on 32x32x32 cube from data/grid.00100.p00000
size= 262152
tsim=4.081296e-01
0oct=0x7fb749981010
done with 1560 octs
dumping data/grid.00100.f701.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00100 5 701 data/grid.00100.f701 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00100.p00000
size= 262152
tsim=4.081296e-01
0oct=0x7fb749981010
done with 1560 octs
dumping data/grid.00100.f702.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00100 5 702 data/grid.00100.f702 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00100.p00000
size= 262152
tsim=4.081296e-01
0oct=0x7fb749981010
done with 1560 octs
dumping data/grid.00100.f706.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00100 5 706 data/grid.00100.f706 1 \
0"
> pli,x(,,16)
> INFO?x
SYNTAX: syntax error near x
SYNTAX: syntax error near <EOF>
> INFO,x
 1: array(double,32,32,32) min=0 max=0 avg=0 std=0
> pli,e(,,16)
> ws
> pli,e(,,16)
> INFO?e
SYNTAX: syntax error near e
SYNTAX: syntax error near <EOF>
> i=10;lvl=5;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=9.585990e-01
0oct=0x7f721dcaf010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 701 data/grid.00010.f701 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=9.585990e-01
0oct=0x7f721dcaf010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 702 data/grid.00010.f702 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=9.585990e-01
0oct=0x7f721dcaf010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 706 data/grid.00010.f706 1 \
0"
> i=10;lvl=5;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=9.585990e-01
0oct=0x7f721dcaf010
done with 4681 octs
dumping data/grid.00010.f701.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 701 data/grid.00010.f701 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=9.585990e-01
0oct=0x7f721dcaf010
done with 4681 octs
dumping data/grid.00010.f702.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 702 data/grid.00010.f702 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=9.585990e-01
0oct=0x7f721dcaf010
done with 4681 octs
dumping data/grid.00010.f706.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 706 data/grid.00010.f706 1 \
0"
> ws1
> pli,e(,,16)
> pli,fx(,,16)
> pli,x(,,16)
> ws
> INFO,x
 1: array(double,32,32,32) min=0.00239258 max=0.00239276 avg=0.00239258 std=2.86917e-09
> i=10;lvl=5;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=9.585990e-01
0oct=0x7fcf1b795010
done with 4681 octs
dumping data/grid.00010.f701.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 701 data/grid.00010.f701 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=9.585990e-01
0oct=0x7fcf1b795010
done with 4681 octs
dumping data/grid.00010.f702.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 702 data/grid.00010.f702 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=9.585990e-01
0oct=0x7fcf1b795010
done with 4681 octs
dumping data/grid.00010.f706.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 706 data/grid.00010.f706 1 \
0"
> INFO,x
 1: array(double,32,32,32) min=0.00239258 max=0.00257578 avg=0.00239262 std=2.86218e-06
> ws
> plg,e(,16,16)*0.5964,color="red"
> plg,fx(,16,16)
> ws
> plg,e(,16,16)*0.5964,color="red"
> plg,fx(,16,16)
> i=10;lvl=5;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=9.585990e-01
0oct=0x7f3106726010
done with 4681 octs
dumping data/grid.00010.f701.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 701 data/grid.00010.f701 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=9.585990e-01
0oct=0x7f3106726010
done with 4681 octs
dumping data/grid.00010.f702.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 702 data/grid.00010.f702 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=9.585990e-01
0oct=0x7f3106726010
done with 4681 octs
dumping data/grid.00010.f706.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 706 data/grid.00010.f706 1 \
0"
> ws
> INFO,x
 1: array(double,32,32,32) min=0.00239258 max=1 avg=0.048449 std=0.203814
> ws
> pli,x(,,16)
> ws
> plg,e(,16,16)*0.5964,color="red"
> plg,fx(,16,16)
> logxy,0,1
> i=10;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=9.575022e-01
0oct=0x7f351a513010
done with 8905 octs
dumping data/grid.00010.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 701 data/grid.00010.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=9.575022e-01
0oct=0x7f351a513010
done with 8905 octs
dumping data/grid.00010.f702.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 702 data/grid.00010.f702 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=9.575022e-01
0oct=0x7f351a513010
done with 8905 octs
dumping data/grid.00010.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 706 data/grid.00010.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=9.575022e-01
0oct=0x7f351a513010
done with 8905 octs
dumping data/grid.00010.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 0 data/grid.00010.f0 1 0"
> ws
> plotamr(l(,,32))
6
5
4
3
2
1
[]
> pli,x(,,16)
> ws
> pli,x(,,32)
> plotamr(l(,,32))
6
5
4
3
2
1
[]
> palette,"idl-33.gp"
> ws
> plg,e(,16,16)*0.5917065,color="red"
> ws
> plg,e(,32,32)*0.5917065,color="red"
> 

> plg,fx(,32,32)
> logxy,0,1
> i=10;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=9.585990e-01
0oct=0x7ffd661ec010
done with 4681 octs
dumping data/grid.00010.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 701 data/grid.00010.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=9.585990e-01
0oct=0x7ffd661ec010
done with 4681 octs
dumping data/grid.00010.f702.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 702 data/grid.00010.f702 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=9.585990e-01
0oct=0x7ffd661ec010
done with 4681 octs
dumping data/grid.00010.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 706 data/grid.00010.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=9.585990e-01
0oct=0x7ffd661ec010
done with 4681 octs
dumping data/grid.00010.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 0 data/grid.00010.f0 1 0"
> ws
> pli,e(,,16)
> info,e
 array(double,64,64,64)
> pli,e(,,32)
> pli,x(,,32)
> INFO,x
 1: array(double,64,64,64) min=0.00120074 max=0.00156154 avg=0.00120088 std=5.75862e-06
> i=10;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7f73dbb70010
done with 47065 octs
dumping data/grid.00010.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 701 data/grid.00010.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7f73dbb70010
done with 47065 octs
dumping data/grid.00010.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 702 data/grid.00010.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7f73dbb70010
done with 47065 octs
dumping data/grid.00010.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 706 data/grid.00010.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7f73dbb70010
done with 47065 octs
dumping data/grid.00010.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 0 data/grid.00010.f0 1 0"
> ws
> pli,x(,,64)
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> ws
> plg,e(,64,64)*0.594,color="red"
> ws
> plg,e(,64,64)*0.594,color="red"
> plg,fx(64,64)
> logxy,0,1
> ws
> plg,fx(64,64)
> info,fx
 array(double,128,128,128)
> ws
> plg,fx(,64,64)
> plg,e(,64,64)*0.594,color="red"
> logxy,0,1
> ws
> pli,e(,,64)
> pli,log10(e(,,64))
> i=5;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);
Casting Rays on 128x128x128 cube from data/grid.00005.p00000
size= 16777224
tsim=4.096090e-01
0oct=0x7f73dbb70010
done with 12521 octs
dumping data/grid.00005.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 7 701 data/grid.00005.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00005.p00000
size= 16777224
tsim=4.096090e-01
0oct=0x7f73dbb70010
done with 12521 octs
dumping data/grid.00005.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 7 702 data/grid.00005.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00005.p00000
size= 16777224
tsim=4.096090e-01
0oct=0x7f73dbb70010
done with 12521 octs
dumping data/grid.00005.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 7 706 data/grid.00005.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00005.p00000
size= 16777224
tsim=4.096090e-01
0oct=0x7f73dbb70010
done with 12521 octs
dumping data/grid.00005.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 7 0 data/grid.00005.f0 1 0"
> pli,log10(e(,,64))
> i=8;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);
Casting Rays on 128x128x128 cube from data/grid.00008.p00000
size= 16777224
tsim=7.082656e-01
0oct=0x7f73dbb70010
done with 28249 octs
dumping data/grid.00008.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00008 7 701 data/grid.00008.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00008.p00000
size= 16777224
tsim=7.082656e-01
0oct=0x7f73dbb70010
done with 28249 octs
dumping data/grid.00008.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00008 7 702 data/grid.00008.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00008.p00000
size= 16777224
tsim=7.082656e-01
0oct=0x7f73dbb70010
done with 28249 octs
dumping data/grid.00008.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00008 7 706 data/grid.00008.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00008.p00000
size= 16777224
tsim=7.082656e-01
0oct=0x7f73dbb70010
done with 28249 octs
dumping data/grid.00008.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00008 7 0 data/grid.00008.f0 1 0"
> pli,log10(e(,,64))
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> ws
> pli,x(,,64)
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> ws
> pli,e(,,64)
> pli,log10(e(,,64))
> ws
> plg,e(,64,64)*0.594,color="red"
> logxy,0,1
> pli,fx(,,16)
> plg,fx(,64,64)
> done with 28249 octs
SYNTAX: syntax error near with 28249 octs
SYNTAX: syntax error near <EOF>
> i=8;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);
Casting Rays on 128x128x128 cube from data/grid.00008.p00000
size= 16777224
tsim=7.082656e-01
0oct=0x7f829b285010
done with 29025 octs
dumping data/grid.00008.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00008 7 701 data/grid.00008.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00008.p00000
size= 16777224
tsim=7.082656e-01
0oct=0x7f829b285010
done with 29025 octs
dumping data/grid.00008.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00008 7 702 data/grid.00008.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00008.p00000
size= 16777224
tsim=7.082656e-01
0oct=0x7f829b285010
done with 29025 octs
dumping data/grid.00008.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00008 7 706 data/grid.00008.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00008.p00000
size= 16777224
tsim=7.082656e-01
0oct=0x7f829b285010
done with 29025 octs
dumping data/grid.00008.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00008 7 0 data/grid.00008.f0 1 0"
> ws
> pli,log10(e(,,64))
> i=6;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);
Casting Rays on 128x128x128 cube from data/grid.00006.p00000
size= 16777224
tsim=5.001621e-01
0oct=0x7f829b285010
done with 16513 octs
dumping data/grid.00006.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 7 701 data/grid.00006.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00006.p00000
size= 16777224
tsim=5.001621e-01
0oct=0x7f829b285010
done with 16513 octs
dumping data/grid.00006.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 7 702 data/grid.00006.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00006.p00000
size= 16777224
tsim=5.001621e-01
0oct=0x7f829b285010
done with 16513 octs
dumping data/grid.00006.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 7 706 data/grid.00006.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00006.p00000
size= 16777224
tsim=5.001621e-01
0oct=0x7f829b285010
done with 16513 octs
dumping data/grid.00006.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 7 0 data/grid.00006.f0 1 0"
> pli,log10(e(,,64))
> ws
> plg,e(,64,64)*0.594,color="red"
> logxy,0,1
> ws1
> plg,e(,64,64)*0.594,color="red"
> logxy,0,1
> i=2;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);
Casting Rays on 128x128x128 cube from data/grid.00002.p00000
size= 16777224
tsim=1.919873e-01
0oct=0x7f829b285010
done with 6681 octs
dumping data/grid.00002.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 7 701 data/grid.00002.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00002.p00000
size= 16777224
tsim=1.919873e-01
0oct=0x7f829b285010
done with 6681 octs
dumping data/grid.00002.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 7 702 data/grid.00002.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00002.p00000
size= 16777224
tsim=1.919873e-01
0oct=0x7f829b285010
done with 6681 octs
dumping data/grid.00002.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 7 706 data/grid.00002.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00002.p00000
size= 16777224
tsim=1.919873e-01
0oct=0x7f829b285010
done with 6681 octs
dumping data/grid.00002.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 7 0 data/grid.00002.f0 1 0"
> plg,e(,64,64)*0.594,color="red"
> i=4;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=3.280592e-01
0oct=0x7f829b285010
done with 9689 octs
dumping data/grid.00004.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 701 data/grid.00004.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=3.280592e-01
0oct=0x7f829b285010
done with 9689 octs
dumping data/grid.00004.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 702 data/grid.00004.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=3.280592e-01
0oct=0x7f829b285010
done with 9689 octs
dumping data/grid.00004.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 706 data/grid.00004.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=3.280592e-01
0oct=0x7f829b285010
done with 9689 octs
dumping data/grid.00004.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 0 data/grid.00004.f0 1 0"
> plg,e(,64,64)*0.594,color="red"
> i=5;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);
Casting Rays on 128x128x128 cube from data/grid.00005.p00000
size= 16777224
tsim=4.096090e-01
0oct=0x7f829b285010
done with 12593 octs
dumping data/grid.00005.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 7 701 data/grid.00005.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00005.p00000
size= 16777224
tsim=4.096090e-01
0oct=0x7f829b285010
done with 12593 octs
dumping data/grid.00005.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 7 702 data/grid.00005.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00005.p00000
size= 16777224
tsim=4.096090e-01
0oct=0x7f829b285010
done with 12593 octs
dumping data/grid.00005.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 7 706 data/grid.00005.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00005.p00000
size= 16777224
tsim=4.096090e-01
0oct=0x7f829b285010
done with 12593 octs
dumping data/grid.00005.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 7 0 data/grid.00005.f0 1 0"
> plg,e(,64,64)*0.594,color="red"
> i=8;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);
Casting Rays on 128x128x128 cube from data/grid.00008.p00000
size= 16777224
tsim=7.082656e-01
0oct=0x7f829b285010
done with 29025 octs
dumping data/grid.00008.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00008 7 701 data/grid.00008.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00008.p00000
size= 16777224
tsim=7.082656e-01
0oct=0x7f829b285010
done with 29025 octs
dumping data/grid.00008.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00008 7 702 data/grid.00008.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00008.p00000
size= 16777224
tsim=7.082656e-01
0oct=0x7f829b285010
done with 29025 octs
dumping data/grid.00008.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00008 7 706 data/grid.00008.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00008.p00000
size= 16777224
tsim=7.082656e-01
0oct=0x7f829b285010
done with 29025 octs
dumping data/grid.00008.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00008 7 0 data/grid.00008.f0 1 0"
> plg,e(,64,64)*0.594,color="red"
> i=9;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);
Casting Rays on 128x128x128 cube from data/grid.00009.p00000
size= 16777224
tsim=8.258107e-01
0oct=0x7f829b285010
done with 38265 octs
dumping data/grid.00009.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00009 7 701 data/grid.00009.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00009.p00000
size= 16777224
tsim=8.258107e-01
0oct=0x7f829b285010
done with 38265 octs
dumping data/grid.00009.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00009 7 702 data/grid.00009.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00009.p00000
size= 16777224
tsim=8.258107e-01
0oct=0x7f829b285010
done with 38265 octs
dumping data/grid.00009.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00009 7 706 data/grid.00009.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00009.p00000
size= 16777224
tsim=8.258107e-01
0oct=0x7f829b285010
done with 38265 octs
dumping data/grid.00009.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00009 7 0 data/grid.00009.f0 1 0"
> plg,e(,64,64)*0.594,color="red"
> plg,e(,64,64)*0.594,color="blue"
> i=10;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7f829b285010
done with 49241 octs
dumping data/grid.00010.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 701 data/grid.00010.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7f829b285010
done with 49241 octs
dumping data/grid.00010.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 702 data/grid.00010.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7f829b285010
done with 49241 octs
dumping data/grid.00010.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 706 data/grid.00010.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7f829b285010
done with 49241 octs
dumping data/grid.00010.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 0 data/grid.00010.f0 1 0"
> plg,e(,64,64)*0.594,color="blue"
> ws
> pli,log10(e(,,64))
> i=10;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fec54827010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 701 data/grid.00010.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fec54827010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 702 data/grid.00010.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fec54827010
ws
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 706 data/grid.00010.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fec54827010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 0 data/grid.00010.f0 1 0"
> i=10;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fec54827010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 701 data/grid.00010.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fec54827010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 702 data/grid.00010.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fec54827010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 706 data/grid.00010.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fec54827010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 0 data/grid.00010.f0 1 0"
> i=10;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fec54827010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 701 data/grid.00010.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fec54827010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 702 data/grid.00010.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fec54827010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 706 data/grid.00010.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fec54827010
Segmentation fault (core dumped)
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 0 data/grid.00010.f0 1 0"
> i=10;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fec54827010
done with 70457 octs
dumping data/grid.00010.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 701 data/grid.00010.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fec54827010
done with 70457 octs
dumping data/grid.00010.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 702 data/grid.00010.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fec54827010
done with 70457 octs
dumping data/grid.00010.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 706 data/grid.00010.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fec54827010
done with 70457 octs
dumping data/grid.00010.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 0 data/grid.00010.f0 1 0"
> pli,log10(e(,,64))
> ws
> pli,x(,,64)
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> rvp
> rvp
> ws
> pli,log10(e(,,64))
> palette,"idl-33.gp"
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> i=10;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fec54827010
done with 70457 octs
dumping data/grid.00010.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 701 data/grid.00010.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fec54827010
done with 70457 octs
dumping data/grid.00010.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 702 data/grid.00010.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fec54827010
done with 70457 octs
dumping data/grid.00010.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 706 data/grid.00010.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fec54827010
done with 70457 octs
dumping data/grid.00010.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 0 data/grid.00010.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fec54827010
done with 70457 octs
dumping data/grid.00010.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 707 data/grid.00010.f707 1 \
0"
> ws
> pli,t(,,64)
> ws
> plg,t(64,,64)
> i=8;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);
Casting Rays on 128x128x128 cube from data/grid.00008.p00000
size= 16777224
tsim=7.082656e-01
0oct=0x7fec54827010
done with 40849 octs
dumping data/grid.00008.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00008 7 701 data/grid.00008.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00008.p00000
size= 16777224
tsim=7.082656e-01
0oct=0x7fec54827010
done with 40849 octs
dumping data/grid.00008.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00008 7 702 data/grid.00008.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00008.p00000
size= 16777224
tsim=7.082656e-01
0oct=0x7fec54827010
done with 40849 octs
dumping data/grid.00008.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00008 7 706 data/grid.00008.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00008.p00000
size= 16777224
tsim=7.082656e-01
0oct=0x7fec54827010
done with 40849 octs
dumping data/grid.00008.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00008 7 0 data/grid.00008.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00008.p00000
size= 16777224
tsim=7.082656e-01
0oct=0x7fec54827010
done with 40849 octs
dumping data/grid.00008.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00008 7 707 data/grid.00008.f707 1 \
0"
> pli,t(,,64)
> ws*
cont> dbexit
ERROR (*main*) bad data type(s) in binary *
WARNING source code unavailable (try dbdis function)
now at pc= 1 (of 11), failed at pc= 5
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> plg,t(64,,64)
> i=1;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);
Casting Rays on 128x128x128 cube from data/grid.00001.p00000
size= 16777224
tsim=1.374774e-01
0oct=0x7fec54827010
done with 6681 octs
dumping data/grid.00001.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 7 701 data/grid.00001.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00001.p00000
size= 16777224
tsim=1.374774e-01
0oct=0x7fec54827010
done with 6681 octs
dumping data/grid.00001.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 7 702 data/grid.00001.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00001.p00000
size= 16777224
tsim=1.374774e-01
0oct=0x7fec54827010
done with 6681 octs
dumping data/grid.00001.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 7 706 data/grid.00001.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00001.p00000
size= 16777224
tsim=1.374774e-01
0oct=0x7fec54827010
done with 6681 octs
dumping data/grid.00001.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 7 0 data/grid.00001.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00001.p00000
size= 16777224
tsim=1.374774e-01
0oct=0x7fec54827010
done with 6681 octs
dumping data/grid.00001.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 7 707 data/grid.00001.f707 1 \
0"
> plg,t(64,,64)
> i=1;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);
Casting Rays on 128x128x128 cube from data/grid.00001.p00000
size= 16777224
tsim=1.374774e-01
0oct=0x7fec54827010
done with 6681 octs
dumping data/grid.00001.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 7 701 data/grid.00001.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00001.p00000
size= 16777224
tsim=1.374774e-01
0oct=0x7fec54827010
done with 6681 octs
dumping data/grid.00001.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 7 702 data/grid.00001.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00001.p00000
size= 16777224
tsim=1.374774e-01
0oct=0x7fec54827010
done with 6681 octs
dumping data/grid.00001.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 7 706 data/grid.00001.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00001.p00000
size= 16777224
tsim=1.374774e-01
0oct=0x7fec54827010
done with 6681 octs
dumping data/grid.00001.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 7 0 data/grid.00001.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00001.p00000
size= 16777224
tsim=1.374774e-01
0oct=0x7fec54827010
done with 6681 octs
dumping data/grid.00001.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 7 707 data/grid.00001.f707 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00001.p00000
size= 16777224
tsim=1.374774e-01
0oct=0x7fec54827010
done with 6681 octs
dumping data/grid.00001.f705.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 7 705 data/grid.00001.f705 1 \
0"
> INFO,s
 1: array(double,128,128,128) min=0 max=1.43088e+75 avg=3.49335e+71 std=2.23547e+73
> i=10;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fec54827010
done with 70457 octs
dumping data/grid.00010.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 701 data/grid.00010.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fec54827010
done with 70457 octs
dumping data/grid.00010.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 702 data/grid.00010.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fec54827010
done with 70457 octs
dumping data/grid.00010.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 706 data/grid.00010.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fec54827010
done with 70457 octs
dumping data/grid.00010.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 0 data/grid.00010.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fec54827010
done with 70457 octs
dumping data/grid.00010.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 707 data/grid.00010.f707 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fec54827010
done with 70457 octs
dumping data/grid.00010.f705.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 705 data/grid.00010.f705 1 \
0"
> INFO,s
 1: array(double,128,128,128) min=0 max=1.43088e+75 avg=3.49335e+71 std=2.23547e+73
> i=10;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fe6fe06b010
done with 113977 octs
dumping data/grid.00010.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 701 data/grid.00010.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fe6fe06b010
done with 113977 octs
dumping data/grid.00010.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 702 data/grid.00010.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fe6fe06b010
done with 113977 octs
dumping data/grid.00010.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 706 data/grid.00010.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fe6fe06b010
done with 113977 octs
dumping data/grid.00010.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 0 data/grid.00010.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fe6fe06b010
done with 113977 octs
dumping data/grid.00010.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 707 data/grid.00010.f707 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7fe6fe06b010
done with 113977 octs
dumping data/grid.00010.f705.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 705 data/grid.00010.f705 1 \
0"
> ws1
> pli,x(,,64)
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> ws
> pli,t(,,64)
> ws
> plg,t(64,,64)
> plg,e(,64,64)*0.594,color="blue"
> ws
> plg,e(,64,64)*0.594,color="blue"
> logxy,0,1
> plg,fx(,64,64)
> ws
> pli,t(,,64)
> i=10;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.585990e-01
0oct=0x7f6b39fd8010
done with 4681 octs
dumping data/grid.00010.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 701 data/grid.00010.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.585990e-01
0oct=0x7f6b39fd8010
done with 4681 octs
dumping data/grid.00010.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 702 data/grid.00010.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.585990e-01
0oct=0x7f6b39fd8010
done with 4681 octs
dumping data/grid.00010.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 706 data/grid.00010.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.585990e-01
0oct=0x7f6b39fd8010
done with 4681 octs
dumping data/grid.00010.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 0 data/grid.00010.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.585990e-01
0oct=0x7f6b39fd8010
done with 4681 octs
dumping data/grid.00010.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 707 data/grid.00010.f707 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.585990e-01
0oct=0x7f6b39fd8010
done with 4681 octs
dumping data/grid.00010.f705.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 705 data/grid.00010.f705 1 \
0"
> plg,t(64,,64)
> ws
> plg,t(64,,64)
> ws
> plg,e(,64,64),color="blue"
> logxy,0,1
> ws
> plg,t(64,,64)
> logxy,0,1
> i=10;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.575022e-01
0oct=0x7fe19c0c7010
done with 20985 octs
dumping data/grid.00010.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 701 data/grid.00010.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.575022e-01
0oct=0x7fe19c0c7010
done with 20985 octs
dumping data/grid.00010.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 702 data/grid.00010.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.575022e-01
0oct=0x7fe19c0c7010
done with 20985 octs
dumping data/grid.00010.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 706 data/grid.00010.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.575022e-01
0oct=0x7fe19c0c7010
done with 20985 octs
dumping data/grid.00010.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 0 data/grid.00010.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.575022e-01
0oct=0x7fe19c0c7010
done with 20985 octs
dumping data/grid.00010.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 707 data/grid.00010.f707 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.575022e-01
0oct=0x7fe19c0c7010
done with 20985 octs
dumping data/grid.00010.f705.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 705 data/grid.00010.f705 1 \
0"
> plg,t(64,,64),color="red"
> i=10;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7f064d4e8010
done with 113977 octs
dumping data/grid.00010.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 701 data/grid.00010.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7f064d4e8010
done with 113977 octs
dumping data/grid.00010.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 702 data/grid.00010.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7f064d4e8010
done with 113977 octs
dumping data/grid.00010.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 706 data/grid.00010.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7f064d4e8010
done with 113977 octs
dumping data/grid.00010.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 0 data/grid.00010.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7f064d4e8010
done with 113977 octs
dumping data/grid.00010.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 707 data/grid.00010.f707 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7f064d4e8010
done with 113977 octs
dumping data/grid.00010.f705.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 705 data/grid.00010.f705 1 \
0"
> plg,t(64,,64),color="blue"
> i=10;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7f3b2e5fe010
done with 113977 octs
dumping data/grid.00010.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 701 data/grid.00010.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7f3b2e5fe010
done with 113977 octs
dumping data/grid.00010.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 702 data/grid.00010.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7f3b2e5fe010
done with 113977 octs
dumping data/grid.00010.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 706 data/grid.00010.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7f3b2e5fe010
done with 113977 octs
dumping data/grid.00010.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 0 data/grid.00010.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7f3b2e5fe010
done with 113977 octs
dumping data/grid.00010.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 707 data/grid.00010.f707 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7f3b2e5fe010
done with 113977 octs
dumping data/grid.00010.f705.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 705 data/grid.00010.f705 1 \
0"
> plg,t(64,,64),color="black"
> plg,t(64,,64),color="green"
> i=6;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);
Casting Rays on 128x128x128 cube from data/grid.00006.p00000
size= 16777224
tsim=5.001621e-01
0oct=0x7f2d1f5d5010
done with 53513 octs
dumping data/grid.00006.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 7 701 data/grid.00006.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00006.p00000
size= 16777224
tsim=5.001621e-01
0oct=0x7f2d1f5d5010
done with 53513 octs
dumping data/grid.00006.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 7 702 data/grid.00006.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00006.p00000
size= 16777224
tsim=5.001621e-01
0oct=0x7f2d1f5d5010
done with 53513 octs
dumping data/grid.00006.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 7 706 data/grid.00006.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00006.p00000
size= 16777224
tsim=5.001621e-01
0oct=0x7f2d1f5d5010
done with 53513 octs
dumping data/grid.00006.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 7 0 data/grid.00006.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00006.p00000
size= 16777224
tsim=5.001621e-01
0oct=0x7f2d1f5d5010
done with 53513 octs
dumping data/grid.00006.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 7 707 data/grid.00006.f707 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00006.p00000
size= 16777224
tsim=5.001621e-01
0oct=0x7f2d1f5d5010
done with 53513 octs
dumping data/grid.00006.f705.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00006 7 705 data/grid.00006.f705 1 \
0"
> plg,t(64,,64),color="green"
> i=7;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);
Casting Rays on 128x128x128 cube from data/grid.00007.p00000
size= 16777224
tsim=5.997152e-01
0oct=0x7f2d1f5d5010
done with 72857 octs
dumping data/grid.00007.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00007 7 701 data/grid.00007.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00007.p00000
size= 16777224
tsim=5.997152e-01
0oct=0x7f2d1f5d5010
done with 72857 octs
dumping data/grid.00007.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00007 7 702 data/grid.00007.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00007.p00000
size= 16777224
tsim=5.997152e-01
0oct=0x7f2d1f5d5010
done with 72857 octs
dumping data/grid.00007.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00007 7 706 data/grid.00007.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00007.p00000
size= 16777224
tsim=5.997152e-01
0oct=0x7f2d1f5d5010
done with 72857 octs
dumping data/grid.00007.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00007 7 0 data/grid.00007.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00007.p00000
size= 16777224
tsim=5.997152e-01
0oct=0x7f2d1f5d5010
done with 72857 octs
dumping data/grid.00007.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00007 7 707 data/grid.00007.f707 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00007.p00000
size= 16777224
tsim=5.997152e-01
0oct=0x7f2d1f5d5010
done with 72857 octs
dumping data/grid.00007.f705.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00007 7 705 data/grid.00007.f705 1 \
0"
> plg,t(64,,64),color="black"
> i=10;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7f2d1f5d5010
done with 122017 octs
dumping data/grid.00010.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 701 data/grid.00010.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7f2d1f5d5010
done with 122017 octs
dumping data/grid.00010.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 702 data/grid.00010.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7f2d1f5d5010
done with 122017 octs
dumping data/grid.00010.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 706 data/grid.00010.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7f2d1f5d5010
done with 122017 octs
dumping data/grid.00010.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 0 data/grid.00010.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7f2d1f5d5010
done with 122017 octs
dumping data/grid.00010.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 707 data/grid.00010.f707 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=9.523488e-01
0oct=0x7f2d1f5d5010
done with 122017 octs
dumping data/grid.00010.f705.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 705 data/grid.00010.f705 1 \
0"
> plg,t(64,,64),color="red"
> ws
> plg,t(64,,64),color="red"
> i=9;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);
Casting Rays on 128x128x128 cube from data/grid.00009.p00000
size= 16777224
tsim=8.258107e-01
0oct=0x7f2d1f5d5010
done with 115481 octs
dumping data/grid.00009.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00009 7 701 data/grid.00009.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00009.p00000
size= 16777224
tsim=8.258107e-01
0oct=0x7f2d1f5d5010
done with 115481 octs
dumping data/grid.00009.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00009 7 702 data/grid.00009.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00009.p00000
size= 16777224
tsim=8.258107e-01
0oct=0x7f2d1f5d5010
done with 115481 octs
dumping data/grid.00009.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00009 7 706 data/grid.00009.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00009.p00000
size= 16777224
tsim=8.258107e-01
0oct=0x7f2d1f5d5010
done with 115481 octs
dumping data/grid.00009.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00009 7 0 data/grid.00009.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00009.p00000
size= 16777224
tsim=8.258107e-01
0oct=0x7f2d1f5d5010
done with 115481 octs
dumping data/grid.00009.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00009 7 707 data/grid.00009.f707 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00009.p00000
size= 16777224
tsim=8.258107e-01
0oct=0x7f2d1f5d5010
done with 115481 octs
dumping data/grid.00009.f705.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00009 7 705 data/grid.00009.f705 1 \
0"
> plg,t(64,,64),color="blue"
> logxy,0,1
> i=8;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);
Casting Rays on 128x128x128 cube from data/grid.00008.p00000
size= 16777224
tsim=7.082656e-01
0oct=0x7f2d1f5d5010
done with 94193 octs
dumping data/grid.00008.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00008 7 701 data/grid.00008.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00008.p00000
size= 16777224
tsim=7.082656e-01
0oct=0x7f2d1f5d5010
done with 94193 octs
dumping data/grid.00008.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00008 7 702 data/grid.00008.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00008.p00000
size= 16777224
tsim=7.082656e-01
0oct=0x7f2d1f5d5010
done with 94193 octs
dumping data/grid.00008.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00008 7 706 data/grid.00008.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00008.p00000
size= 16777224
tsim=7.082656e-01
0oct=0x7f2d1f5d5010
done with 94193 octs
dumping data/grid.00008.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00008 7 0 data/grid.00008.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00008.p00000
size= 16777224
tsim=7.082656e-01
0oct=0x7f2d1f5d5010
done with 94193 octs
dumping data/grid.00008.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00008 7 707 data/grid.00008.f707 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00008.p00000
size= 16777224
tsim=7.082656e-01
0oct=0x7f2d1f5d5010
done with 94193 octs
dumping data/grid.00008.f705.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00008 7 705 data/grid.00008.f705 1 \
0"
> plg,t(64,,64),color="black"
> i=4;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=3.280592e-01
0oct=0x7f2d1f5d5010
done with 26473 octs
dumping data/grid.00004.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 701 data/grid.00004.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=3.280592e-01
0oct=0x7f2d1f5d5010
done with 26473 octs
dumping data/grid.00004.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 702 data/grid.00004.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=3.280592e-01
0oct=0x7f2d1f5d5010
done with 26473 octs
dumping data/grid.00004.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 706 data/grid.00004.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=3.280592e-01
0oct=0x7f2d1f5d5010
done with 26473 octs
dumping data/grid.00004.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 0 data/grid.00004.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=3.280592e-01
0oct=0x7f2d1f5d5010
done with 26473 octs
dumping data/grid.00004.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 707 data/grid.00004.f707 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=3.280592e-01
0oct=0x7f2d1f5d5010
done with 26473 octs
dumping data/grid.00004.f705.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 705 data/grid.00004.f705 1 \
0"
> plg,t(64,,64),color="black"
> i=5;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=1.307416e-01
0oct=0x7f42694ba010
done with 37449 octs
dumping data/grid.00005.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 701 data/grid.00005.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=1.307416e-01
0oct=0x7f42694ba010
done with 37449 octs
dumping data/grid.00005.f702.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 702 data/grid.00005.f702 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=1.307416e-01
0oct=0x7f42694ba010
done with 37449 octs
dumping data/grid.00005.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 706 data/grid.00005.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=1.307416e-01
0oct=0x7f42694ba010
done with 37449 octs
dumping data/grid.00005.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 0 data/grid.00005.f0 1 0"
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=1.307416e-01
0oct=0x7f42694ba010
done with 37449 octs
dumping data/grid.00005.f707.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 707 data/grid.00005.f707 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00005.p00000
size= 2097160
tsim=1.307416e-01
0oct=0x7f42694ba010
done with 37449 octs
dumping data/grid.00005.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00005 6 705 data/grid.00005.f705 1 \
0"
> ws
> pli,x(,,32)
ERROR (*main*) Floating point interrupt (SIGFPE)
WARNING source code unavailable (try dbdis function)
now at pc= 1 (of 19), failed at pc= 13
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> INFO,x
 1: array(double,64,64,64) ERROR (_stat_worker) Floating point interrupt (SIGFPE)
  LINE: 916  FILE: /usr/lib/yorick/i/Eric/utils.i
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> i=1;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=5.462054e-02
0oct=0x7f2505f45010
done with 37449 octs
dumping data/grid.00001.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 701 data/grid.00001.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=5.462054e-02
0oct=0x7f2505f45010
done with 37449 octs
dumping data/grid.00001.f702.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 702 data/grid.00001.f702 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=5.462054e-02
0oct=0x7f2505f45010
done with 37449 octs
dumping data/grid.00001.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 706 data/grid.00001.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=5.462054e-02
0oct=0x7f2505f45010
done with 37449 octs
dumping data/grid.00001.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 0 data/grid.00001.f0 1 0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=5.462054e-02
0oct=0x7f2505f45010
done with 37449 octs
dumping data/grid.00001.f707.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 707 data/grid.00001.f707 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=5.462054e-02
0oct=0x7f2505f45010
done with 37449 octs
dumping data/grid.00001.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 705 data/grid.00001.f705 1 \
0"
> INFO,x
 1: array(double,64,64,64) ERROR (_stat_worker) Floating point interrupt (SIGFPE)
  LINE: 916  FILE: /usr/lib/yorick/i/Eric/utils.i
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> i=1;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7f5fa8a69010
done with 37449 octs
dumping data/grid.00001.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 701 data/grid.00001.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7f5fa8a69010
done with 37449 octs
dumping data/grid.00001.f702.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 702 data/grid.00001.f702 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7f5fa8a69010
done with 37449 octs
dumping data/grid.00001.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 706 data/grid.00001.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7f5fa8a69010
done with 37449 octs
dumping data/grid.00001.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 0 data/grid.00001.f0 1 0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7f5fa8a69010
done with 37449 octs
dumping data/grid.00001.f707.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 707 data/grid.00001.f707 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=3.759265e-02
0oct=0x7f5fa8a69010
done with 37449 octs
dumping data/grid.00001.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 705 data/grid.00001.f705 1 \
0"
> INFO,x
 1: array(double,64,64,64) ERROR (_stat_worker) Floating point interrupt (SIGFPE)
  LINE: 916  FILE: /usr/lib/yorick/i/Eric/utils.i
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> i=0;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f5fa8a69010
done with 37449 octs
dumping data/grid.00000.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 701 data/grid.00000.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f5fa8a69010
done with 37449 octs
dumping data/grid.00000.f702.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 702 data/grid.00000.f702 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f5fa8a69010
done with 37449 octs
dumping data/grid.00000.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 706 data/grid.00000.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f5fa8a69010
done with 37449 octs
dumping data/grid.00000.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 0 data/grid.00000.f0 1 0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f5fa8a69010
done with 37449 octs
dumping data/grid.00000.f707.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 707 data/grid.00000.f707 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f5fa8a69010
done with 37449 octs
dumping data/grid.00000.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 705 data/grid.00000.f705 1 \
0"
> INFO,x
 1: array(double,64,64,64) ERROR (_stat_worker) Floating point interrupt (SIGFPE)
  LINE: 916  FILE: /usr/lib/yorick/i/Eric/utils.i
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> INFO,e
 1: array(double,64,64,64) min=0 max=0 avg=0 std=0
> INFO,src
 1: array(double,32,32,32) min=0 max=1.43088e+68 avg=3.49335e+64 std=2.23547e+66
> INFO?s
SYNTAX: syntax error near s
SYNTAX: syntax error near <EOF>
> INFO,s
 1: array(double,64,64,64) min=0 max=2.31112e+78 avg=7.05298e+73 std=1.27671e+76
> i=0;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f430c6b9010
done with 37449 octs
dumping data/grid.00000.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 701 data/grid.00000.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f430c6b9010
done with 37449 octs
dumping data/grid.00000.f702.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 702 data/grid.00000.f702 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f430c6b9010
done with 37449 octs
dumping data/grid.00000.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 706 data/grid.00000.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f430c6b9010
done with 37449 octs
dumping data/grid.00000.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 0 data/grid.00000.f0 1 0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f430c6b9010
done with 37449 octs
dumping data/grid.00000.f707.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 707 data/grid.00000.f707 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f430c6b9010
done with 37449 octs
dumping data/grid.00000.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 705 data/grid.00000.f705 1 \
0"
> INFO,e
 1: array(double,64,64,64) min=0 max=0 avg=0 std=0
> INFO,s
 1: array(double,64,64,64) min=0 max=2.31112e+78 avg=7.05298e+73 std=1.27671e+76
> INFO,x
 1: array(double,64,64,64) ERROR (_stat_worker) Floating point interrupt (SIGFPE)
  LINE: 916  FILE: /usr/lib/yorick/i/Eric/utils.i
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> i=0;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f2b6f656010
done with 37449 octs
dumping data/grid.00000.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 701 data/grid.00000.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f2b6f656010
done with 37449 octs
dumping data/grid.00000.f702.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 702 data/grid.00000.f702 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f2b6f656010
done with 37449 octs
dumping data/grid.00000.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 706 data/grid.00000.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f2b6f656010
done with 37449 octs
dumping data/grid.00000.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 0 data/grid.00000.f0 1 0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f2b6f656010
done with 37449 octs
dumping data/grid.00000.f707.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 707 data/grid.00000.f707 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7f2b6f656010
done with 37449 octs
dumping data/grid.00000.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 705 data/grid.00000.f705 1 \
0"
> INFO,x
 1: array(double,64,64,64) min=0 max=0 avg=0 std=0
> INFO,e
 1: array(double,64,64,64) min=0 max=3.08726e+77 avg=9.42158e+72 std=1.70546e+75
> pli,e(,,max)
> i=3;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);
Casting Rays on 64x64x64 cube from data/grid.00003.p00000
size= 2097160
tsim=9.289027e-02
0oct=0x7f3eefbc1010
done with 37449 octs
dumping data/grid.00003.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00003 6 701 data/grid.00003.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00003.p00000
size= 2097160
tsim=9.289027e-02
0oct=0x7f3eefbc1010
done with 37449 octs
dumping data/grid.00003.f702.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00003 6 702 data/grid.00003.f702 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00003.p00000
size= 2097160
tsim=9.289027e-02
0oct=0x7f3eefbc1010
done with 37449 octs
dumping data/grid.00003.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00003 6 706 data/grid.00003.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00003.p00000
size= 2097160
tsim=9.289027e-02
0oct=0x7f3eefbc1010
done with 37449 octs
dumping data/grid.00003.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00003 6 0 data/grid.00003.f0 1 0"
Casting Rays on 64x64x64 cube from data/grid.00003.p00000
size= 2097160
tsim=9.289027e-02
0oct=0x7f3eefbc1010
done with 37449 octs
dumping data/grid.00003.f707.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00003 6 707 data/grid.00003.f707 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00003.p00000
size= 2097160
tsim=9.289027e-02
0oct=0x7f3eefbc1010
done with 37449 octs
dumping data/grid.00003.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00003 6 705 data/grid.00003.f705 1 \
0"
> INFO,x
 1: array(double,64,64,64) min=0.0012 max=1 avg=0.00120764 std=0.00275881
> ws
> pli,x(,,max)
> pli,e(,,max)
> w
[]
> ws
> pli,e(,,max)
> i=4;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);
Casting Rays on 64x64x64 cube from data/grid.00004.p00000
size= 2097160
tsim=1.119884e-01
0oct=0x7f3eefbc1010
done with 37449 octs
dumping data/grid.00004.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 6 701 data/grid.00004.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00004.p00000
size= 2097160
tsim=1.119884e-01
0oct=0x7f3eefbc1010
done with 37449 octs
dumping data/grid.00004.f702.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 6 702 data/grid.00004.f702 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00004.p00000
size= 2097160
tsim=1.119884e-01
0oct=0x7f3eefbc1010
done with 37449 octs
dumping data/grid.00004.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 6 706 data/grid.00004.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00004.p00000
size= 2097160
tsim=1.119884e-01
0oct=0x7f3eefbc1010
done with 37449 octs
dumping data/grid.00004.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 6 0 data/grid.00004.f0 1 0"
Casting Rays on 64x64x64 cube from data/grid.00004.p00000
size= 2097160
tsim=1.119884e-01
0oct=0x7f3eefbc1010
done with 37449 octs
dumping data/grid.00004.f707.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 6 707 data/grid.00004.f707 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00004.p00000
size= 2097160
tsim=1.119884e-01
0oct=0x7f3eefbc1010
done with 37449 octs
dumping data/grid.00004.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 6 705 data/grid.00004.f705 1 \
0"
> ws1
> pli,x(,,max)
> pli,e(,,max)
> ws
> pli,log(e(,,max))
> i=10;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=2.114747e-01
0oct=0x7f3eefbc1010
done with 37449 octs
dumping data/grid.00010.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 701 data/grid.00010.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=2.114747e-01
0oct=0x7f3eefbc1010
done with 37449 octs
dumping data/grid.00010.f702.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 702 data/grid.00010.f702 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=2.114747e-01
0oct=0x7f3eefbc1010
done with 37449 octs
dumping data/grid.00010.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 706 data/grid.00010.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=2.114747e-01
0oct=0x7f3eefbc1010
done with 37449 octs
dumping data/grid.00010.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 0 data/grid.00010.f0 1 0"
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=2.114747e-01
0oct=0x7f3eefbc1010
done with 37449 octs
dumping data/grid.00010.f707.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 707 data/grid.00010.f707 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=2.114747e-01
0oct=0x7f3eefbc1010
done with 37449 octs
dumping data/grid.00010.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 705 data/grid.00010.f705 1 \
0"
> ws
> pli,x(,,avg)
> pli,x(,,max)
> ws
> INFO,x
 1: array(double,64,64,64) min=0.00120006 max=1 avg=0.54798 std=0.480817
> pli,x(,,max)
> pli,x(,,avg)
> pli,x(,,32)
> pli,log(e(,,max))
> ws
> pli,log(e(,,max))
> pli,log(e(,,avg))
> wkl
> ws1
> wkl
> ws
> pli,log(e(,,avg))
> 10./64^3
3.8147e-05
> i=4;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=1.130031e-01
0oct=0x7f90f836b010
done with 37557 octs
dumping data/grid.00004.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 701 data/grid.00004.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=1.130031e-01
0oct=0x7f90f836b010
done with 37557 octs
dumping data/grid.00004.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 702 data/grid.00004.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=1.130031e-01
0oct=0x7f90f836b010
done with 37557 octs
dumping data/grid.00004.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 706 data/grid.00004.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=1.130031e-01
0oct=0x7f90f836b010
done with 37557 octs
dumping data/grid.00004.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 0 data/grid.00004.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=1.130031e-01
0oct=0x7f90f836b010
done with 37557 octs
dumping data/grid.00004.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 707 data/grid.00004.f707 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00004.p00000
size= 16777224
tsim=1.130031e-01
0oct=0x7f90f836b010
done with 37557 octs
dumping data/grid.00004.f705.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 7 705 data/grid.00004.f705 1 \
0"
> INFO,x
 1: array(double,128,128,128) min=0.00120002 max=1 avg=0.00275634 std=0.0390931
> ws1
> pli,x(,,32)
> ws
> pli,x(,,avg)
> plotamr(l(,,64))
6
5
4
3
2
1
[]
> ws
> www=where2(l==7)
> h=histo1d(www(3,),span(0,128,128))
>
> 
> h(mxx)
115
> pli,x(,,115)
> plotamr(l(,,115))
7
6
5
4
3
2
1
[]
> ws
> plotamr(l(,,115))
7
6
5
4
3
2
1
[]
> pli,log(e(,,115))
> palette,"idl-33.gp"
> plotamr(l(,,115))
7
6
5
4
3
2
1
[]
> INFO,t
 1: array(double,128,128,128) min=10000 max=10000 avg=10000 std=0
> ws
> pli,t(,,64)
> INFO,x
 1: array(double,128,128,128) min=0.00120002 max=1 avg=0.00275634 std=0.0390931
> ws1
> i=10;lvl=5;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=4.474665e-01
0oct=0x7fffd168d010
done with 4681 octs
dumping data/grid.00010.f701.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 701 data/grid.00010.f701 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=4.474665e-01
0oct=0x7fffd168d010
done with 4681 octs
dumping data/grid.00010.f702.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 702 data/grid.00010.f702 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=4.474665e-01
0oct=0x7fffd168d010
done with 4681 octs
dumping data/grid.00010.f706.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 706 data/grid.00010.f706 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=4.474665e-01
0oct=0x7fffd168d010
done with 4681 octs
dumping data/grid.00010.f0.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 0 data/grid.00010.f0 1 0"
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=4.474665e-01
0oct=0x7fffd168d010
done with 4681 octs
dumping data/grid.00010.f707.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 707 data/grid.00010.f707 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=4.474665e-01
0oct=0x7fffd168d010
done with 4681 octs
dumping data/grid.00010.f705.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 705 data/grid.00010.f705 1 \
0"
> ws1
> pli,x(,,16)
> INFO,x
 1: array(double,32,32,32) min=0 max=0 avg=0 std=0
> pli,e(,,max)
> pli,e(,,16)
> pli,log(e(,,16))
> pli,e(,,16)
> 10./32^3
0.000305176
> i=10;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=4.100093e-01
0oct=0x7fffd168d010
done with 7730 octs
dumping data/grid.00010.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 701 data/grid.00010.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=4.100093e-01
0oct=0x7fffd168d010
done with 7730 octs
dumping data/grid.00010.f702.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 702 data/grid.00010.f702 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=4.100093e-01
0oct=0x7fffd168d010
done with 7730 octs
dumping data/grid.00010.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 706 data/grid.00010.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=4.100093e-01
0oct=0x7fffd168d010
done with 7730 octs
dumping data/grid.00010.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 0 data/grid.00010.f0 1 0"
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=4.100093e-01
0oct=0x7fffd168d010
done with 7730 octs
dumping data/grid.00010.f707.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 707 data/grid.00010.f707 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00010.p00000
size= 2097160
tsim=4.100093e-01
0oct=0x7fffd168d010
done with 7730 octs
dumping data/grid.00010.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 6 705 data/grid.00010.f705 1 \
0"
> ws
> pli,e(,,32)
> plotamr(l(,,32))
6
5
4
3
2
1
[]
> i=10;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=3.824121e-01
0oct=0x7fffd168d010
done with 9251 octs
dumping data/grid.00010.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 701 data/grid.00010.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=3.824121e-01
0oct=0x7fffd168d010
done with 9251 octs
dumping data/grid.00010.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 702 data/grid.00010.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=3.824121e-01
0oct=0x7fffd168d010
done with 9251 octs
dumping data/grid.00010.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 706 data/grid.00010.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=3.824121e-01
0oct=0x7fffd168d010
done with 9251 octs
dumping data/grid.00010.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 0 data/grid.00010.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=3.824121e-01
0oct=0x7fffd168d010
done with 9251 octs
dumping data/grid.00010.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 707 data/grid.00010.f707 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=3.824121e-01
0oct=0x7fffd168d010
done with 9251 octs
dumping data/grid.00010.f705.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 705 data/grid.00010.f705 1 \
0"
> pli,e(,,64)
> ws
> pli,e(,,64)
> plotamr(l(,,64))
6
5
4
3
2
1
[]
> www=where2(l==7)
> h=histo1d(www(3,),span(0,128,128))
> h(mxx)
103
> ws
> pli,e(,,103)
> plotamr(l(,,103))
7
6
5
4
3
2
1
[]
> info,n
 []
> info,nh
 array(double,32,32,32)
>
> 
> 
> i=10;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101,a);
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=3.824121e-01
0oct=0x7fffd168d010
done with 9251 octs
dumping data/grid.00010.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 701 data/grid.00010.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=3.824121e-01
0oct=0x7fffd168d010
done with 9251 octs
dumping data/grid.00010.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 702 data/grid.00010.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=3.824121e-01
0oct=0x7fffd168d010
done with 9251 octs
dumping data/grid.00010.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 706 data/grid.00010.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=3.824121e-01
0oct=0x7fffd168d010
done with 9251 octs
dumping data/grid.00010.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 0 data/grid.00010.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=3.824121e-01
0oct=0x7fffd168d010
done with 9251 octs
dumping data/grid.00010.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 707 data/grid.00010.f707 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=3.824121e-01
0oct=0x7fffd168d010
done with 9251 octs
dumping data/grid.00010.f705.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 705 data/grid.00010.f705 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00010.p00000
size= 16777224
tsim=3.824121e-01
0oct=0x7fffd168d010
done with 9251 octs
dumping data/grid.00010.f101.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 7 101 data/grid.00010.f101 1 \
0"
> ws
> pli,d(,,103)
> pli,log(d(,,103))
> plotamr(l(,,103))
7
6
5
4
3
2
1
[]
> pli,s(,,103)
> sqrt(3.3^2+5.54^2+1.056^2)
6.53427
> sqrt(3.3^2+5.54^2+10.56^2)
12.3732
> i=41;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101,a);
Casting Rays on 128x128x128 cube from data/grid.00041.p00000
size= 16777224
tsim=2.355476e-01
0oct=0x7fffd168d010
done with 5259 octs
dumping data/grid.00041.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00041 7 701 data/grid.00041.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00041.p00000
size= 16777224
tsim=2.355476e-01
0oct=0x7fffd168d010
done with 5259 octs
dumping data/grid.00041.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00041 7 702 data/grid.00041.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00041.p00000
size= 16777224
tsim=2.355476e-01
0oct=0x7fffd168d010
done with 5259 octs
dumping data/grid.00041.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00041 7 706 data/grid.00041.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00041.p00000
size= 16777224
tsim=2.355476e-01
0oct=0x7fffd168d010
done with 5259 octs
dumping data/grid.00041.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00041 7 0 data/grid.00041.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00041.p00000
size= 16777224
tsim=2.355476e-01
0oct=0x7fffd168d010
done with 5259 octs
dumping data/grid.00041.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00041 7 707 data/grid.00041.f707 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00041.p00000
size= 16777224
tsim=2.355476e-01
0oct=0x7fffd168d010
done with 5259 octs
dumping data/grid.00041.f705.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00041 7 705 data/grid.00041.f705 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00041.p00000
size= 16777224
tsim=2.355476e-01
0oct=0x7fffd168d010
done with 5259 octs
dumping data/grid.00041.f101.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00041 7 101 data/grid.00041.f101 1 \
0"
> ws
> INFO,x
 1: array(double,128,128,128) min=0 max=0 avg=0 std=0
> INFO,e
 1: array(double,128,128,128) min=0 max=3.08074e+84 avg=1.01154e+81 std=3.92241e+82
> ws1
> pli,e(,,103)
> pli,log(e(,,103))
ERROR (*main*) Floating point interrupt (SIGFPE)
WARNING source code unavailable (try dbdis function)
now at pc= 1 (of 23), failed at pc= 13
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> pli,log10(e(,,103)+1e-5)
> plotamr(l(,,103))
6
5
4
3
2
1
[]
> ws
> www=where2(l==7)
> h=histo1d(www(3,),span(0,128,128))
> h(mxx)
37
> pli,log10(e(,,37)+1e-5)
> plotamr(l(,,37))
7
6
5
4
3
2
1
[]
> ws
> plotamr(l(,,37))
7
6
5
4
3
2
1
[]
> pli,log10(e(,,37)+1e-5)
> palette,"idl-33.gp"
> pli,log10(e(,,37)+1e-5)
> plotamr(l(,,37))
7
6
5
4
3
2
1
[]
> ws
> pli,log10(e(,,37)+1e-5)
> plg,e(1,,37)
> ws
> plg,e(1,,37)
> plg,log(e(1,,37))
ERROR (*main*) Floating point interrupt (SIGFPE)
WARNING source code unavailable (try dbdis function)
now at pc= 1 (of 24), failed at pc= 14
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> ws
> plotamr(l(,,37))
7
6
5
4
3
2
1
[]
> plg,e(8,,37)
> ws
> plg,e(8,,37)
> plg,e(,8,37)
> ws
> pli,fx(,,37)
> ws
> pli,e(,,37)
> plg,e(,8,37)
> ws
> plg,e(,8,37)
> plg,log10(e(,8,37)+1e-10⁾
> )
> ws
> plg,log10(e(,8,37)+1e-10)
> pli,fx(,,37)ws
SYNTAX: syntax error near ws
SYNTAX: syntax error near <EOF>
> ws
> pli,fx(,,37)ws
SYNTAX: syntax error near ws
SYNTAX: syntax error near <EOF>
> pli,fx(,,37)
> plotamr(l(,,37))
7
6
5
4
3
2
1
[]
> ws
> pli,log10(e(,,37)+1e-5)
> pli,log10(e(,,37)+1e-10)
> plotamr(l(,,37))
7
6
5
4
3
2
1
[]
> i=24;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101,a);
Casting Rays on 128x128x128 cube from data/grid.00024.p00000
size= 16777224
tsim=1.629220e-01
0oct=0x7fffd168d010
done with 4681 octs
dumping data/grid.00024.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00024 7 701 data/grid.00024.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00024.p00000
size= 16777224
tsim=1.629220e-01
0oct=0x7fffd168d010
done with 4681 octs
dumping data/grid.00024.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00024 7 702 data/grid.00024.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00024.p00000
size= 16777224
tsim=1.629220e-01
0oct=0x7fffd168d010
done with 4681 octs
dumping data/grid.00024.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00024 7 706 data/grid.00024.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00024.p00000
size= 16777224
tsim=1.629220e-01
0oct=0x7fffd168d010
done with 4681 octs
dumping data/grid.00024.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00024 7 0 data/grid.00024.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00024.p00000
size= 16777224
tsim=1.629220e-01
0oct=0x7fffd168d010
done with 4681 octs
dumping data/grid.00024.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00024 7 707 data/grid.00024.f707 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00024.p00000
size= 16777224
tsim=1.629220e-01
0oct=0x7fffd168d010
done with 4681 octs
dumping data/grid.00024.f705.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00024 7 705 data/grid.00024.f705 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00024.p00000
size= 16777224
tsim=1.629220e-01
0oct=0x7fffd168d010
done with 4681 octs
dumping data/grid.00024.f101.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00024 7 101 data/grid.00024.f101 1 \
0"
> ws
> pli,log10(e(,,37)+1e-10)
> www=where2(l==7)
> h(mxx)
37
> h=histo1d(www(3,),span(0,128,128))
ERROR (*main*) attempt to call non-function or index non-array
WARNING source code unavailable (try dbdis function)
now at pc= 1 (of 28), failed at pc= 8
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> www=where2(l==6)
> h=histo1d(www(3,),span(0,128,128))
ERROR (*main*) attempt to call non-function or index non-array
WARNING source code unavailable (try dbdis function)
now at pc= 1 (of 28), failed at pc= 8
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> i=41;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101,a);
Casting Rays on 128x128x128 cube from data/grid.00041.p00000
size= 16777224
tsim=2.355476e-01
0oct=0x7fffd168d010
done with 5259 octs
dumping data/grid.00041.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00041 7 701 data/grid.00041.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00041.p00000
size= 16777224
tsim=2.355476e-01
0oct=0x7fffd168d010
done with 5259 octs
dumping data/grid.00041.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00041 7 702 data/grid.00041.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00041.p00000
size= 16777224
tsim=2.355476e-01
0oct=0x7fffd168d010
done with 5259 octs
dumping data/grid.00041.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00041 7 706 data/grid.00041.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00041.p00000
size= 16777224
tsim=2.355476e-01
0oct=0x7fffd168d010
done with 5259 octs
dumping data/grid.00041.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00041 7 0 data/grid.00041.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00041.p00000
size= 16777224
tsim=2.355476e-01
0oct=0x7fffd168d010
done with 5259 octs
dumping data/grid.00041.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00041 7 707 data/grid.00041.f707 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00041.p00000
size= 16777224
tsim=2.355476e-01
0oct=0x7fffd168d010
done with 5259 octs
dumping data/grid.00041.f705.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00041 7 705 data/grid.00041.f705 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00041.p00000
size= 16777224
tsim=2.355476e-01
0oct=0x7fffd168d010
done with 5259 octs
dumping data/grid.00041.f101.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00041 7 101 data/grid.00041.f101 1 \
0"
> ws
> pli,log10(e(,,37)+1e-10)
> plotamr(l(,,37))
7
6
5
4
3
2
1
[]
> i=10;lvl=5;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101,a);
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=1.367939e-01
0oct=0x7fffe84c3010
done with 4681 octs
dumping data/grid.00010.f701.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 701 data/grid.00010.f701 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=1.367939e-01
0oct=0x7fffe84c3010
done with 4681 octs
dumping data/grid.00010.f702.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 702 data/grid.00010.f702 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=1.367939e-01
0oct=0x7fffe84c3010
done with 4681 octs
dumping data/grid.00010.f706.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 706 data/grid.00010.f706 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=1.367939e-01
0oct=0x7fffe84c3010
done with 4681 octs
dumping data/grid.00010.f0.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 0 data/grid.00010.f0 1 0"
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=1.367939e-01
0oct=0x7fffe84c3010
done with 4681 octs
dumping data/grid.00010.f707.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 707 data/grid.00010.f707 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=1.367939e-01
0oct=0x7fffe84c3010
done with 4681 octs
dumping data/grid.00010.f705.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 705 data/grid.00010.f705 1 \
0"
Casting Rays on 32x32x32 cube from data/grid.00010.p00000
size= 262152
tsim=1.367939e-01
0oct=0x7fffe84c3010
done with 4681 octs
dumping data/grid.00010.f101.p00000 with nmap=32768
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00010 5 101 data/grid.00010.f101 1 \
0"
> ws1
> pli,e(,,16)
> pli,log(e(,,16))
ERROR (*main*) Floating point interrupt (SIGFPE)
WARNING source code unavailable (try dbdis function)
now at pc= 1 (of 23), failed at pc= 13
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> pli,log(e(,,16)+1e-10)
> ws
> plg,e(,16,16))
SYNTAX: syntax error near )
SYNTAX: syntax error near <EOF>
> plg,e(,16,16)
> i=50;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101,a);
Casting Rays on 64x64x64 cube from data/grid.00050.p00000
size= 2097160
tsim=4.125402e-01
0oct=0x7fffe84c3010
done with 4745 octs
dumping data/grid.00050.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 6 701 data/grid.00050.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00050.p00000
size= 2097160
tsim=4.125402e-01
0oct=0x7fffe84c3010
done with 4745 octs
dumping data/grid.00050.f702.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 6 702 data/grid.00050.f702 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00050.p00000
size= 2097160
tsim=4.125402e-01
0oct=0x7fffe84c3010
done with 4745 octs
dumping data/grid.00050.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 6 706 data/grid.00050.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00050.p00000
size= 2097160
tsim=4.125402e-01
0oct=0x7fffe84c3010
done with 4745 octs
dumping data/grid.00050.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 6 0 data/grid.00050.f0 1 0"
Casting Rays on 64x64x64 cube from data/grid.00050.p00000
size= 2097160
tsim=4.125402e-01
0oct=0x7fffe84c3010
done with 4745 octs
dumping data/grid.00050.f707.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 6 707 data/grid.00050.f707 1 \
0"i
Casting Rays on 64x64x64 cube from data/grid.00050.p00000
size= 2097160
tsim=4.125402e-01
0oct=0x7fffe84c3010
done with 4745 octs
dumping data/grid.00050.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 6 705 data/grid.00050.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00050.p00000
size= 2097160
tsim=4.125402e-01
0oct=0x7fffe84c3010
done with 4745 octs
dumping data/grid.00050.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 6 101 data/grid.00050.f101 1 \
0"
> ws1
> pli,log(e(,,32)+1e-10)
> rvp
> palette,"idl-33.gp"
> plotamr(l(,,32))
6
5
4
3
2
1
[]
> ws
> plg,e(,32,32)
> logxy,0,1
> 0"i
SYNTAX: missing close "
SYNTAX: syntax error near <EOF>
> i=50;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101,a);
Casting Rays on 64x64x64 cube from data/grid.00050.p00000
size= 2097160
tsim=4.125402e-01
0oct=0x7fffe84c3010
done with 5769 octs
dumping data/grid.00050.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 6 701 data/grid.00050.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00050.p00000
size= 2097160
tsim=4.125402e-01
0oct=0x7fffe84c3010
done with 5769 octs
dumping data/grid.00050.f702.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 6 702 data/grid.00050.f702 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00050.p00000
size= 2097160
tsim=4.125402e-01
0oct=0x7fffe84c3010
done with 5769 octs
dumping data/grid.00050.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 6 706 data/grid.00050.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00050.p00000
size= 2097160
tsim=4.125402e-01
0oct=0x7fffe84c3010
done with 5769 octs
dumping data/grid.00050.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 6 0 data/grid.00050.f0 1 0"
Casting Rays on 64x64x64 cube from data/grid.00050.p00000
size= 2097160
tsim=4.125402e-01
0oct=0x7fffe84c3010
done with 5769 octs
dumping data/grid.00050.f707.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 6 707 data/grid.00050.f707 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00050.p00000
size= 2097160
tsim=4.125402e-01
0oct=0x7fffe84c3010
done with 5769 octs
dumping data/grid.00050.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 6 705 data/grid.00050.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00050.p00000
size= 2097160
tsim=4.125402e-01
0oct=0x7fffe84c3010
done with 5769 octs
dumping data/grid.00050.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 6 101 data/grid.00050.f101 1 \
0"
> ws1
> pli,log(e(,,32)+1e-10)
> plotamr(l(,,32))
6
5
4
3
2
1
[]
> rvp
> rvp
> ws
> pli,e(,,16)
> pli,e(,,32)
> pli,e(,,32)
> pli,log(e(,,32)+1e-10)
> pli,e(,,32)
> ws
> plg,e(,32,32)
> logxy,0,1
> e(,32,32)
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5.88438e+72,8.82749e+72,1.20305e+73,
1.57927e+73,2.03616e+73,2.6136e+73,3.38974e+73,4.52888e+73,6.40864e+73,
1.00338e+74,1.87766e+74,2.15869e+74,2.15869e+74,1.87766e+74,1.00338e+74,
6.40864e+73,4.52888e+73,3.38974e+73,2.6136e+73,2.03616e+73,1.57927e+73,
1.20305e+73,8.82749e+72,5.88438e+72,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
> i=50;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101,a);
Casting Rays on 128x128x128 cube from data/grid.00050.p00000
size= 16777224
tsim=4.096090e-01
0oct=0x7fffe84c3010
done with 14313 octs
dumping data/grid.00050.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 7 701 data/grid.00050.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00050.p00000
size= 16777224
tsim=4.096090e-01
0oct=0x7fffe84c3010
done with 14313 octs
dumping data/grid.00050.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 7 702 data/grid.00050.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00050.p00000
size= 16777224
tsim=4.096090e-01
0oct=0x7fffe84c3010
done with 14313 octs
dumping data/grid.00050.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 7 706 data/grid.00050.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00050.p00000
size= 16777224
tsim=4.096090e-01
0oct=0x7fffe84c3010
done with 14313 octs
dumping data/grid.00050.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 7 0 data/grid.00050.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00050.p00000
size= 16777224
tsim=4.096090e-01
0oct=0x7fffe84c3010
done with 14313 octs
dumping data/grid.00050.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 7 707 data/grid.00050.f707 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00050.p00000
size= 16777224
tsim=4.096090e-01
0oct=0x7fffe84c3010
done with 14313 octs
dumping data/grid.00050.f705.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 7 705 data/grid.00050.f705 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00050.p00000
size= 16777224
tsim=4.096090e-01
0oct=0x7fffe84c3010
done with 14313 octs
dumping data/grid.00050.f101.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 7 101 data/grid.00050.f101 1 \
0"
> ws1
> palette,"idl-33.gp"
> pli,e(,,32)
> pli,e(,,64)
> plotamr(l(,,32))
5
4
3
2
1
[]
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> ws
> plg,e(,64,64)
> logxy,0,1
> PL,e(,64,64)
> plg,l(,64,64)*1e70
> i=50;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101,a);
Casting Rays on 64x64x64 cube from data/grid.00050.p00000
size= 2097160
tsim=4.096090e-01
0oct=0x7fffe84c3010
done with 7777 octs
dumping data/grid.00050.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 6 701 data/grid.00050.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00050.p00000
size= 2097160
tsim=4.096090e-01
0oct=0x7fffe84c3010
done with 7777 octs
dumping data/grid.00050.f702.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 6 702 data/grid.00050.f702 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00050.p00000
size= 2097160
tsim=4.096090e-01
0oct=0x7fffe84c3010
done with 7777 octs
dumping data/grid.00050.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 6 706 data/grid.00050.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00050.p00000
size= 2097160
tsim=4.096090e-01
0oct=0x7fffe84c3010
done with 7777 octs
dumping data/grid.00050.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 6 0 data/grid.00050.f0 1 0"
Casting Rays on 64x64x64 cube from data/grid.00050.p00000
size= 2097160
tsim=4.096090e-01
0oct=0x7fffe84c3010
done with 7777 octs
dumping data/grid.00050.f707.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 6 707 data/grid.00050.f707 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00050.p00000
size= 2097160
tsim=4.096090e-01
0oct=0x7fffe84c3010
done with 7777 octs
dumping data/grid.00050.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 6 705 data/grid.00050.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00050.p00000
size= 2097160
tsim=4.096090e-01
0oct=0x7fffe84c3010
done with 7777 octs
dumping data/grid.00050.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 6 101 data/grid.00050.f101 1 \
0"
> ws
> pli,e(,,32)
> plg,e(,32,32)
> ws
> plg,e(,32,32)
> logxy,0,1
> plg,l(,32,32)*1e70
> i=50;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101,a);
Casting Rays on 128x128x128 cube from data/grid.00050.p00000
size= 16777224
tsim=4.096090e-01
0oct=0x7fffe84c3010
done with 14625 octs
dumping data/grid.00050.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 7 701 data/grid.00050.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00050.p00000
size= 16777224
tsim=4.096090e-01
0oct=0x7fffe84c3010
done with 14625 octs
dumping data/grid.00050.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 7 702 data/grid.00050.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00050.p00000
size= 16777224
tsim=4.096090e-01
0oct=0x7fffe84c3010
done with 14625 octs
dumping data/grid.00050.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 7 706 data/grid.00050.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00050.p00000
size= 16777224
tsim=4.096090e-01
0oct=0x7fffe84c3010
done with 14625 octs
dumping data/grid.00050.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 7 0 data/grid.00050.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00050.p00000
size= 16777224
tsim=4.096090e-01
0oct=0x7fffe84c3010
done with 14625 octs
dumping data/grid.00050.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 7 707 data/grid.00050.f707 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00050.p00000
size= 16777224
tsim=4.096090e-01
0oct=0x7fffe84c3010
done with 14625 octs
dumping data/grid.00050.f705.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 7 705 data/grid.00050.f705 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00050.p00000
size= 16777224
tsim=4.096090e-01
0oct=0x7fffe84c3010
done with 14625 octs
dumping data/grid.00050.f101.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 7 101 data/grid.00050.f101 1 \
0"
> 
> ws
> pli,e(,,64)
> pli,log(e(,,64)+1e-10)
> pli,log(e(,,32)+1e-10)
> ws
> pli,log(e(,,32)+1e-10)
> plotamr(l(,,64))
7
6
5
4
3
2
1
[]
> ws
> plg,e(,64,64)
> PL,e(,64,64)
> logxy,0,1
> ws
> 10./32^3
0.000305176
> i=50;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101,a);
Casting Rays on 64x64x64 cube from data/grid.00050.p00000
size= 2097160
tsim=2.681333e-01
0oct=0x7fffd168d010
done with 5488 octs
dumping data/grid.00050.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 6 701 data/grid.00050.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00050.p00000
size= 2097160
tsim=2.681333e-01
0oct=0x7fffd168d010
done with 5488 octs
dumping data/grid.00050.f702.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 6 702 data/grid.00050.f702 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00050.p00000
size= 2097160
tsim=2.681333e-01
0oct=0x7fffd168d010
done with 5488 octs
dumping data/grid.00050.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 6 706 data/grid.00050.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00050.p00000
size= 2097160
tsim=2.681333e-01
0oct=0x7fffd168d010
done with 5488 octs
dumping data/grid.00050.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 6 0 data/grid.00050.f0 1 0"
Casting Rays on 64x64x64 cube from data/grid.00050.p00000
size= 2097160
tsim=2.681333e-01
0oct=0x7fffd168d010
done with 5488 octs
dumping data/grid.00050.f707.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 6 707 data/grid.00050.f707 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00050.p00000
size= 2097160
tsim=2.681333e-01
0oct=0x7fffd168d010
done with 5488 octs
dumping data/grid.00050.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 6 705 data/grid.00050.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00050.p00000
size= 2097160
tsim=2.681333e-01
0oct=0x7fffd168d010
done with 5488 octs
dumping data/grid.00050.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 6 101 data/grid.00050.f101 1 \
0"
> ws1
> www=where2(l==6)
> h=histo1d(www(3,),span(0,128,128))
> h(mxx)
53
> pli,e(,,53)
> plotamr(l(,,53))
6
5
4
3
2
1
[]
> i=50;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101,a);
Casting Rays on 128x128x128 cube from data/grid.00050.p00000
size= 16777224
tsim=2.655791e-01
0oct=0x7fffd168d010
done with 5794 octs
dumping data/grid.00050.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 7 701 data/grid.00050.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00050.p00000
size= 16777224
tsim=2.655791e-01
0oct=0x7fffd168d010
done with 5794 octs
dumping data/grid.00050.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 7 702 data/grid.00050.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00050.p00000
size= 16777224
tsim=2.655791e-01
0oct=0x7fffd168d010
done with 5794 octs
dumping data/grid.00050.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 7 706 data/grid.00050.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00050.p00000
size= 16777224
tsim=2.655791e-01
0oct=0x7fffd168d010
done with 5794 octs
dumping data/grid.00050.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 7 0 data/grid.00050.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00050.p00000
size= 16777224
tsim=2.655791e-01
0oct=0x7fffd168d010
done with 5794 octs
dumping data/grid.00050.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 7 707 data/grid.00050.f707 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00050.p00000
size= 16777224
tsim=2.655791e-01
0oct=0x7fffd168d010
done with 5794 octs
dumping data/grid.00050.f705.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 7 705 data/grid.00050.f705 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00050.p00000
size= 16777224
tsim=2.655791e-01
0oct=0x7fffd168d010
done with 5794 octs
dumping data/grid.00050.f101.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00050 7 101 data/grid.00050.f101 1 \
0"
> www=where2(l==7)
> h=histo1d(www(3,),span(0,128,128))
> h(mxx)
39
> ws
> pli,e(,,39)
> plotamr(l(,,39))
7
6
5
4
3
2
1
[]
> pli,log(e(,,39))
> plk,log(e(,,39))
> rvp
> palette,"idl-33.gp"
> i=100;lvl=7;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101,a);
Casting Rays on 128x128x128 cube from data/grid.00100.p00000
size= 16777224
tsim=3.824121e-01
0oct=0x7fffd168d010
done with 9251 octs
dumping data/grid.00100.f701.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00100 7 701 data/grid.00100.f701 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00100.p00000
size= 16777224
tsim=3.824121e-01
0oct=0x7fffd168d010
done with 9251 octs
dumping data/grid.00100.f702.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00100 7 702 data/grid.00100.f702 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00100.p00000
size= 16777224
tsim=3.824121e-01
0oct=0x7fffd168d010
done with 9251 octs
dumping data/grid.00100.f706.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00100 7 706 data/grid.00100.f706 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00100.p00000
size= 16777224
tsim=3.824121e-01
0oct=0x7fffd168d010
done with 9251 octs
dumping data/grid.00100.f0.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00100 7 0 data/grid.00100.f0 1 0"
Casting Rays on 128x128x128 cube from data/grid.00100.p00000
size= 16777224
tsim=3.824121e-01
0oct=0x7fffd168d010
done with 9251 octs
dumping data/grid.00100.f707.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00100 7 707 data/grid.00100.f707 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00100.p00000
size= 16777224
tsim=3.824121e-01
0oct=0x7fffd168d010
done with 9251 octs
dumping data/grid.00100.f705.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00100 7 705 data/grid.00100.f705 1 \
0"
Casting Rays on 128x128x128 cube from data/grid.00100.p00000
size= 16777224
tsim=3.824121e-01
0oct=0x7fffd168d010
done with 9251 octs
dumping data/grid.00100.f101.p00000 with nmap=2097152
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00100 7 101 data/grid.00100.f101 1 \
0"
> ws
> pli,log(e(,,39))
> plotamr(l(,,39))
7
6
5
4
3
2
1
[]
> i=2;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101,a);
Casting Rays on 64x64x64 cube from data/grid.00002.p00000
size= 2097160
tsim=7.372642e-02
0oct=0x7fffd168d010
done with 37449 octs
dumping data/grid.00002.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 6 701 data/grid.00002.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00002.p00000
size= 2097160
tsim=7.372642e-02
0oct=0x7fffd168d010
done with 37449 octs
dumping data/grid.00002.f702.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 6 702 data/grid.00002.f702 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00002.p00000
size= 2097160
tsim=7.372642e-02
0oct=0x7fffd168d010
done with 37449 octs
dumping data/grid.00002.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 6 706 data/grid.00002.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00002.p00000
size= 2097160
tsim=7.372642e-02
0oct=0x7fffd168d010
done with 37449 octs
dumping data/grid.00002.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 6 0 data/grid.00002.f0 1 0"
Casting Rays on 64x64x64 cube from data/grid.00002.p00000
size= 2097160
tsim=7.372642e-02
0oct=0x7fffd168d010
done with 37449 octs
dumping data/grid.00002.f707.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 6 707 data/grid.00002.f707 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00002.p00000
size= 2097160
tsim=7.372642e-02
0oct=0x7fffd168d010
done with 37449 octs
dumping data/grid.00002.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 6 705 data/grid.00002.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00002.p00000
size= 2097160
tsim=7.372642e-02
0oct=0x7fffd168d010
done with 37449 octs
dumping data/grid.00002.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00002 6 101 data/grid.00002.f101 1 \
0"
> ws
> pli,e(,,32)
> pli,log(e(,,32))
ERROR (*main*) Floating point interrupt (SIGFPE)
WARNING source code unavailable (try dbdis function)
now at pc= 1 (of 23), failed at pc= 13
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> pli,log(e(,,32)+1e-10)
> ws
> pli,e(,,32)
> pli,e(,,max)
> wsws
[]
> i=4;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101,a);
Casting Rays on 64x64x64 cube from data/grid.00004.p00000
size= 2097160
tsim=1.119884e-01
0oct=0x7fffd168d010
done with 37449 octs
dumping data/grid.00004.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 6 701 data/grid.00004.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00004.p00000
size= 2097160
tsim=1.119884e-01
0oct=0x7fffd168d010
done with 37449 octs
dumping data/grid.00004.f702.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 6 702 data/grid.00004.f702 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00004.p00000
size= 2097160
tsim=1.119884e-01
0oct=0x7fffd168d010
done with 37449 octs
dumping data/grid.00004.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 6 706 data/grid.00004.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00004.p00000
size= 2097160
tsim=1.119884e-01
0oct=0x7fffd168d010
done with 37449 octs
dumping data/grid.00004.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 6 0 data/grid.00004.f0 1 0"
Casting Rays on 64x64x64 cube from data/grid.00004.p00000
size= 2097160
tsim=1.119884e-01
0oct=0x7fffd168d010
done with 37449 octs
dumping data/grid.00004.f707.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 6 707 data/grid.00004.f707 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00004.p00000
size= 2097160
tsim=1.119884e-01
0oct=0x7fffd168d010
done with 37449 octs
dumping data/grid.00004.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 6 705 data/grid.00004.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00004.p00000
size= 2097160
tsim=1.119884e-01
0oct=0x7fffd168d010
done with 37449 octs
dumping data/grid.00004.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00004 6 101 data/grid.00004.f101 1 \
0"
> ws
> pli,e(,,32)
> pli,log(e(,,32)+1e-10)
> i=1;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101,a);
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=5.462054e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00001.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 701 data/grid.00001.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=5.462054e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00001.f702.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 702 data/grid.00001.f702 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=5.462054e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00001.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 706 data/grid.00001.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=5.462054e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00001.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 0 data/grid.00001.f0 1 0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=5.462054e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00001.f707.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 707 data/grid.00001.f707 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=5.462054e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00001.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 705 data/grid.00001.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00001.p00000
size= 2097160
tsim=5.462054e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00001.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00001 6 101 data/grid.00001.f101 1 \
0"
> ws
> pli,e(,,32)
> pli,e(,,max)
> INFO,e
 1: array(double,64,64,64) min=0 max=0 avg=0 std=0
> INFO?x
SYNTAX: syntax error near x
SYNTAX: syntax error near <EOF>
> INFO?x
SYNTAX: syntax error near x
SYNTAX: syntax error near <EOF>
> INFO,x
 1: array(double,64,64,64) ERROR (_stat_worker) Floating point interrupt (SIGFPE)
  LINE: 916  FILE: /usr/lib/yorick/i/Eric/utils.i
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> INFO?s
SYNTAX: syntax error near s
SYNTAX: syntax error near <EOF>
> INFO,s
 1: array(double,64,64,64) min=0 max=6.89501e+84 avg=5.6366e+82 std=6.2086e+83
> i=0;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101,a);
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00000.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 701 data/grid.00000.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00000.f702.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 702 data/grid.00000.f702 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00000.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 706 data/grid.00000.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00000.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 0 data/grid.00000.f0 1 0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00000.f707.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 707 data/grid.00000.f707 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00000.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 705 data/grid.00000.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00000.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 101 data/grid.00000.f101 1 \
0"
> INFO,e
 1: array(double,64,64,64) min=0 max=0 avg=0 std=0
> INFO,x
 1: array(double,64,64,64) ERROR (_stat_worker) Floating point interrupt (SIGFPE)
  LINE: 916  FILE: /usr/lib/yorick/i/Eric/utils.i
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> i=0;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101,a);
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00000.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 701 data/grid.00000.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00000.f702.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 702 data/grid.00000.f702 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00000.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 706 data/grid.00000.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00000.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 0 data/grid.00000.f0 1 0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00000.f707.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 707 data/grid.00000.f707 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00000.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 705 data/grid.00000.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00000.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 101 data/grid.00000.f101 1 \
0"
> INFO,x
 1: array(double,64,64,64) ERROR (_stat_worker) Floating point interrupt (SIGFPE)
  LINE: 916  FILE: /usr/lib/yorick/i/Eric/utils.i
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> INFO,e
 1: array(double,64,64,64) min=0 max=0 avg=0 std=0
> INFO,t
 1: array(double,64,64,64) min=0 max=0 avg=0 std=0
> INFO,s
 1: array(double,64,64,64) min=0 max=2.94855e+84 avg=8.99826e+79 std=1.62883e+82
> i=0;lvl=6;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101,a);
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00000.f701.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 701 data/grid.00000.f701 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00000.f702.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 702 data/grid.00000.f702 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00000.f706.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 706 data/grid.00000.f706 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00000.f0.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 0 data/grid.00000.f0 1 0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00000.f707.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 707 data/grid.00000.f707 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00000.f705.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 705 data/grid.00000.f705 1 \
0"
Casting Rays on 64x64x64 cube from data/grid.00000.p00000
size= 2097160
tsim=3.571846e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00000.f101.p00000 with nmap=262144
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00000 6 101 data/grid.00000.f101 1 \
0"
> INFO,x
 1: array(double,64,64,64) ERROR (_stat_worker) Floating point interrupt (SIGFPE)
  LINE: 916  FILE: /usr/lib/yorick/i/Eric/utils.i
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> INFO,temp
 1: void, []
> INFO,t
 1: array(double,64,64,64) min=0 max=0 avg=0 std=0
> i=20;lvl=8;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101,a);
Casting Rays on 256x256x256 cube from data/grid.00020.p00000
size= 134217736
tsim=6.583581e-02
0oct=0x7fffcbaff010
done with 111032 octs
dumping data/grid.00020.f701.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00020 8 701 data/grid.00020.f701 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00020.p00000
size= 134217736
tsim=6.583581e-02
0oct=0x7fffcbaff010
done with 111032 octs
dumping data/grid.00020.f702.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00020 8 702 data/grid.00020.f702 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00020.p00000
size= 134217736
tsim=6.583581e-02
0oct=0x7fffcbaff010
done with 111032 octs
dumping data/grid.00020.f706.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00020 8 706 data/grid.00020.f706 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00020.p00000
size= 134217736
tsim=6.583581e-02
0oct=0x7fffcbaff010
done with 111032 octs
dumping data/grid.00020.f0.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00020 8 0 data/grid.00020.f0 1 0"
Casting Rays on 256x256x256 cube from data/grid.00020.p00000
size= 134217736
tsim=6.583581e-02
0oct=0x7fffcbaff010
done with 111032 octs
dumping data/grid.00020.f707.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00020 8 707 data/grid.00020.f707 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00020.p00000
size= 134217736
tsim=6.583581e-02
0oct=0x7fffcbaff010
done with 111032 octs
dumping data/grid.00020.f705.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00020 8 705 data/grid.00020.f705 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00020.p00000
size= 134217736
tsim=6.583581e-02
0oct=0x7fffcbaff010
done with 111032 octs
dumping data/grid.00020.f101.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00020 8 101 data/grid.00020.f101 1 \
0"
> www=where2(l==8)
> h=histo1d(www(3,),span(0,256,256))
> h(mxx)
141
> ws
> ws1
> pli,log(d(,,141))
> palette,"idl-33.gp"
> plotamr(l(,,141))
8
7
6
5
4
3
2
1
[]
> ws
> pli,x(,,141)
> pli,log(1.-x(,,141))
> plk,log(1.-x(,,141))
> plotamr(l(,,141))
8
7
6
5
4
3
2
1
[]
> ws
> pli,log(1.-x(,,141))
> ws
> pli,log(1.-x(,,141))
> plk,log(1.-x(,,141))
> ws
> pli,s(,,141)
> i=21;lvl=8;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101,a);
Casting Rays on 256x256x256 cube from data/grid.00021.p00000
size= 134217736
tsim=6.756658e-02
0oct=0x7fffcbaff010
done with 86931 octs
dumping data/grid.00021.f701.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00021 8 701 data/grid.00021.f701 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00021.p00000
size= 134217736
tsim=6.756658e-02
0oct=0x7fffcbaff010
done with 86931 octs
dumping data/grid.00021.f702.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00021 8 702 data/grid.00021.f702 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00021.p00000
size= 134217736
tsim=6.756658e-02
0oct=0x7fffcbaff010
done with 86931 octs
dumping data/grid.00021.f706.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00021 8 706 data/grid.00021.f706 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00021.p00000
size= 134217736
tsim=6.756658e-02
0oct=0x7fffcbaff010
done with 86931 octs
dumping data/grid.00021.f0.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00021 8 0 data/grid.00021.f0 1 0"
Casting Rays on 256x256x256 cube from data/grid.00021.p00000
size= 134217736
tsim=6.756658e-02
0oct=0x7fffcbaff010
done with 86931 octs
dumping data/grid.00021.f707.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00021 8 707 data/grid.00021.f707 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00021.p00000
size= 134217736
tsim=6.756658e-02
0oct=0x7fffcbaff010
done with 86931 octs
dumping data/grid.00021.f705.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00021 8 705 data/grid.00021.f705 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00021.p00000
size= 134217736
tsim=6.756658e-02
0oct=0x7fffcbaff010
done with 86931 octs
dumping data/grid.00021.f101.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00021 8 101 data/grid.00021.f101 1 \
0"
> ws1
> pli,log(1.-x(,,141))
> palette,"idl-33.gp"
> plcont,log(1.-x(,,141))
> ws
> pli,log(1.-x(,,141))
> plct,log(1.-x(,,141))
> ws
> pli,log(1.-x(,,141))
> a
0.0675666
> 1./a-1
13.8002
> ws
> pli,d(,,103)
> pli,d(,,141)
> pli,log(d(,,141))
> plotamr(l(,,141),lmin=7)
8
7
[]
> www=where2(l==8)
> h=histo1d(www(3,),span(0,256,256))
> h(mxx)
133
> ws
> pli,d(,,133)
> pli,log(d(,,133))
> plotamr(l(,,133),lmin=7)
8
7
[]
> ws
> pli,log(1.-x(,,133))
> i=22;lvl=8;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101,a);
Casting Rays on 256x256x256 cube from data/grid.00022.p00000
size= 134217736
tsim=6.926227e-02
0oct=0x7fffcbaff010
done with 123679 octs
dumping data/grid.00022.f701.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00022 8 701 data/grid.00022.f701 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00022.p00000
size= 134217736
tsim=6.926227e-02
0oct=0x7fffcbaff010
done with 123679 octs
dumping data/grid.00022.f702.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00022 8 702 data/grid.00022.f702 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00022.p00000
size= 134217736
tsim=6.926227e-02
0oct=0x7fffcbaff010
done with 123679 octs
dumping data/grid.00022.f706.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00022 8 706 data/grid.00022.f706 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00022.p00000
size= 134217736
tsim=6.926227e-02
0oct=0x7fffcbaff010
done with 123679 octs
dumping data/grid.00022.f0.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00022 8 0 data/grid.00022.f0 1 0"
Casting Rays on 256x256x256 cube from data/grid.00022.p00000
size= 134217736
tsim=6.926227e-02
0oct=0x7fffcbaff010
done with 123679 octs
dumping data/grid.00022.f707.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00022 8 707 data/grid.00022.f707 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00022.p00000
size= 134217736
tsim=6.926227e-02
0oct=0x7fffcbaff010
done with 123679 octs
dumping data/grid.00022.f705.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00022 8 705 data/grid.00022.f705 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00022.p00000
size= 134217736
tsim=6.926227e-02
0oct=0x7fffcbaff010
done with 123679 octs
dumping data/grid.00022.f101.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00022 8 101 data/grid.00022.f101 1 \
0"
> pli,log(1.-x(,,133))
> pli,d(,,133)
> pli,log10(d(,,133))
> pli,log10(d(,,1))
> pli,log10(d(1,,))
> i=23;lvl=8;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101,a);
Casting Rays on 256x256x256 cube from data/grid.00023.p00000
size= 134217736
tsim=7.947448e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00023.f701.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00023 8 701 data/grid.00023.f701 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00023.p00000
size= 134217736
tsim=7.947448e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00023.f702.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00023 8 702 data/grid.00023.f702 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00023.p00000
size= 134217736
tsim=7.947448e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00023.f706.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00023 8 706 data/grid.00023.f706 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00023.p00000
size= 134217736
tsim=7.947448e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00023.f0.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00023 8 0 data/grid.00023.f0 1 0"
Casting Rays on 256x256x256 cube from data/grid.00023.p00000
size= 134217736
tsim=7.947448e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00023.f707.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00023 8 707 data/grid.00023.f707 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00023.p00000
size= 134217736
tsim=7.947448e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00023.f705.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00023 8 705 data/grid.00023.f705 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00023.p00000
size= 134217736
tsim=7.947448e-02
0oct=0x7fffcbaff010
done with 37449 octs
dumping data/grid.00023.f101.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00023 8 101 data/grid.00023.f101 1 \
0"
> pli,log10(d(1,,))
ERROR (*main*) Floating point interrupt (SIGFPE)
WARNING source code unavailable (try dbdis function)
now at pc= 1 (of 23), failed at pc= 13
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> i=22;lvl=8;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101,a);
Casting Rays on 256x256x256 cube from data/grid.00022.p00000
size= 134217736
tsim=6.926227e-02
0oct=0x7fffcbaff010
done with 123679 octs
dumping data/grid.00022.f701.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00022 8 701 data/grid.00022.f701 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00022.p00000
size= 134217736
tsim=6.926227e-02
0oct=0x7fffcbaff010
done with 123679 octs
dumping data/grid.00022.f702.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00022 8 702 data/grid.00022.f702 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00022.p00000
size= 134217736
tsim=6.926227e-02
0oct=0x7fffcbaff010
done with 123679 octs
dumping data/grid.00022.f706.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00022 8 706 data/grid.00022.f706 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00022.p00000
size= 134217736
tsim=6.926227e-02
0oct=0x7fffcbaff010
done with 123679 octs
dumping data/grid.00022.f0.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00022 8 0 data/grid.00022.f0 1 0"
Casting Rays on 256x256x256 cube from data/grid.00022.p00000
size= 134217736
tsim=6.926227e-02
0oct=0x7fffcbaff010
done with 123679 octs
dumping data/grid.00022.f707.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00022 8 707 data/grid.00022.f707 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00022.p00000
size= 134217736
tsim=6.926227e-02
0oct=0x7fffcbaff010
done with 123679 octs
dumping data/grid.00022.f705.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00022 8 705 data/grid.00022.f705 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00022.p00000
size= 134217736
tsim=6.926227e-02
0oct=0x7fffcbaff010
done with 123679 octs
dumping data/grid.00022.f101.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00022 8 101 data/grid.00022.f101 1 \
0"
> pli,log10(d(1,,))
ERROR (*main*) Floating point interrupt (SIGFPE)
WARNING source code unavailable (try dbdis function)
now at pc= 1 (of 23), failed at pc= 13
 To enter debug mode, type <RETURN> now (then dbexit to get out)
> ws
> i=23;lvl=8;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101,a);
Casting Rays on 256x256x256 cube from data/grid.00023.p00000
size= 134217736
tsim=7.093570e-02
0oct=0x7fffcbaff010
done with 128400 octs
dumping data/grid.00023.f701.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00023 8 701 data/grid.00023.f701 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00023.p00000
size= 134217736
tsim=7.093570e-02
0oct=0x7fffcbaff010
done with 128400 octs
dumping data/grid.00023.f702.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00023 8 702 data/grid.00023.f702 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00023.p00000
size= 134217736
tsim=7.093570e-02
0oct=0x7fffcbaff010
done with 128400 octs
dumping data/grid.00023.f706.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00023 8 706 data/grid.00023.f706 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00023.p00000
size= 134217736
tsim=7.093570e-02
0oct=0x7fffcbaff010
done with 128400 octs
dumping data/grid.00023.f0.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00023 8 0 data/grid.00023.f0 1 0"
Casting Rays on 256x256x256 cube from data/grid.00023.p00000
size= 134217736
tsim=7.093570e-02
0oct=0x7fffcbaff010
done with 128400 octs
dumping data/grid.00023.f707.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00023 8 707 data/grid.00023.f707 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00023.p00000
size= 134217736
tsim=7.093570e-02
0oct=0x7fffcbaff010
done with 128400 octs
dumping data/grid.00023.f705.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00023 8 705 data/grid.00023.f705 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00023.p00000
size= 134217736
tsim=7.093570e-02
0oct=0x7fffcbaff010
done with 128400 octs
dumping data/grid.00023.f101.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00023 8 101 data/grid.00023.f101 1 \
0"
> pli,log10(d(1,,))
> pli,log10(d(10,,))
> pli,log10(d(0,,))
> pli,log10(d(10,,))
> pli,log10(d(,0,))
> plotamr(l(,0,),lmin=7)
8
7
[]
> plct,log(1.-x(,0,))
> ws
> pli,log10(d(,0,))
> plct,log(1.-x(,0,))
> ws
> pli,log(1.-x(,0,))
> ws1
> pli,log(1.-x(,0,))
> plotamr(l(,0,),lmin=7)
8
7
[]
> ws
> pli,e(,0,)
> pli,log10(e(,0,))
> ws
> pli,d(,,133)
> pli,d(,,1)
> ws
> pli,log10(d(,1,))
> pli,log10(d(1,,))
> ws
> pli,log10(d(1,,))
> pli,log10(d(,,1))
> pli,log10(d(,32,))
> plotamr(l(,32,),lmin=7)
8
7
[]
> ws
> i=27;lvl=8;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101,a);
Casting Rays on 256x256x256 cube from data/grid.00027.p00000
size= 134217736
tsim=7.713703e-02
0oct=0x7fffcbaff010
done with 100821 octs
dumping data/grid.00027.f701.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00027 8 701 data/grid.00027.f701 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00027.p00000
size= 134217736
tsim=7.713703e-02
0oct=0x7fffcbaff010
done with 100821 octs
dumping data/grid.00027.f702.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00027 8 702 data/grid.00027.f702 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00027.p00000
size= 134217736
tsim=7.713703e-02
0oct=0x7fffcbaff010
done with 100821 octs
dumping data/grid.00027.f706.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00027 8 706 data/grid.00027.f706 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00027.p00000
size= 134217736
tsim=7.713703e-02
0oct=0x7fffcbaff010
done with 100821 octs
dumping data/grid.00027.f0.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00027 8 0 data/grid.00027.f0 1 0"
Casting Rays on 256x256x256 cube from data/grid.00027.p00000
size= 134217736
tsim=7.713703e-02
0oct=0x7fffcbaff010
done with 100821 octs
dumping data/grid.00027.f707.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00027 8 707 data/grid.00027.f707 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00027.p00000
size= 134217736
tsim=7.713703e-02
0oct=0x7fffcbaff010
done with 100821 octs
dumping data/grid.00027.f705.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00027 8 705 data/grid.00027.f705 1 \
0"
Casting Rays on 256x256x256 cube from data/grid.00027.p00000
size= 134217736
tsim=7.713703e-02
0oct=0x7fffcbaff010
done with 100821 octs
dumping data/grid.00027.f101.p00000 with nmap=16777216
status=0
"~/Project/Quartz/utils/oct2grid data/grid.00027 8 101 data/grid.00027.f101 1 \
0"
> pli,log10(d(,32,))
> palette,"idl-33.gp"
> pli,log10(d(1,,))
> pli,log10(d(,,1))
> pli,log10(d(,,133))
> plotamr(l(,,133),lmin=7)
8
7
[]
> ws
> pli,log(1.-x(,,133))
> plotamr(l(,,133),lmin=7)
8
7
[]
> ws
> plotamr(l(,,133),lmin=7)
8
7
[]
> plotamr(l(,,133),lmin=7)
8
7
[]
> ws
> pli,log(1.-x(,,133))
> pli,log(1.-x(,,1))
> pli,log(1.-x(,144,))
> ws
> pli,log10(d(,,133))
> plct,log(1.-x(,,133))
> ws
> pli,log(1.-x(,144,))
> plct,log(d(,,133))
> rvp
> rvp
> palette,"idl-02.gp"
> ws
> i=27;lvl=8;e=oct2cube(swrite(format="data/grid.%05d",i),lvl,701);fx=oct2cube(swrite(format="data/grid.%05d",i),lvl,702,a);x=oct2cube(swrite(format="data/grid.%05d",i),lvl,706,a);l=oct2cube(swrite(format="data/grid.%05d",i),lvl,0,a);t=oct2cube(swrite(format="data/grid.%05d",i),lvl,707,a);s=oct2cube(swrite(format="data/grid.%05d",i),lvl,705,a);d=oct2cube(swrite(format="data/grid.%05d",i),lvl,101,a);
