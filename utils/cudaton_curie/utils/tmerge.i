

x=array(float,512,512,512);
s=readsnap("out/snap.00009.p00000",nbnd=16);
x(:256,:256,:256)=(*s.xion);
s=readsnap("out/snap.00009.p00001",nbnd=16);
x(257:,:256,:256)=(*s.xion);
s=readsnap("out/snap.00009.p00002",nbnd=16);
x(:256,257:,:256)=(*s.xion);
s=readsnap("out/snap.00009.p00003",nbnd=16);
x(257:,257:,:256)=(*s.xion);

s=readsnap("out/snap.00009.p00004",nbnd=16);
x(:256,:256,257:)=(*s.xion);
s=readsnap("out/snap.00009.p00005",nbnd=16);
x(257:,:256,257:)=(*s.xion);
s=readsnap("out/snap.00009.p00006",nbnd=16);
x(:256,257:,257:)=(*s.xion);
s=readsnap("out/snap.00009.p00007",nbnd=16);
x(257:,257:,257:)=(*s.xion);


