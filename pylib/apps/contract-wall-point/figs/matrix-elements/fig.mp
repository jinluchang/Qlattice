prologues := 3;
outputtemplate := "%j-%c.eps";

verbatimtex
\documentclass{article}
\usepackage[intlimits]{amsmath}
\usepackage{amsfonts}
\begin{document}
etex

input feynmp;

def drawbase_meson =
begingroup;
pickup pencircle scaled 1pt;
u := 5.0pt;
%
x1=x2;
x3=x4;
%
y1=y3;
y2=y4;
%
y1-y5=0.5(y5-y6)=y6-y2;
x3-x6=0.5(x6-x5)=x5-x1;
y5=y5p;
y6=y6p;
x5p=x6;
x6p=x5;
%
z1=(0,0);
x3-x1=20u;
y1-y2=10u;
%
%label.top(btex $0$ etex, z5+(2,0));
%label.bot(btex $x$ etex, z6+(-2,0));
%
%label.lft(btex $d$ etex, z1);
%label.lft(btex $\bar{u}$ etex, z2);
label.lft(btex $K/\pi$ etex, 0.5(z1+z2));
%label.rt(btex $u$ etex, z3);
%label.rt(btex $\bar{d}$ etex, z4);
label.rt(btex $K/\pi$ etex, 0.5(z3+z4));
%
endgroup;
enddef;

def drawbase_decay =
begingroup;
pickup pencircle scaled 1pt;
u := 5.0pt;
%
x1=x2;
x3=x4;
%
y1=y3;
y2=y4;
%
y1-y5=0.5(y5-y6)=y6-y2;
x3-x6=0.5(x6-x5)=x5-x1;
y5=y5p;
y6=y6p;
x5p=x6;
x6p=x5;
%
z1=(0,0);
x3-x1=20u;
y1-y2=10u;
%
%label.top(btex $0$ etex, z5+(2,0));
%label.bot(btex $x$ etex, z6+(-2,0));
%
%label.lft(btex $d$ etex, z1);
%label.lft(btex $\bar{u}$ etex, z2);
%label.lft(btex $K/\pi$ etex, 0.5(z1+z2));
%label.rt(btex $u$ etex, z3);
%label.rt(btex $\bar{d}$ etex, z4);
label.rt(btex $K/\pi$ etex, 0.5(z3+z4));
%
endgroup;
enddef;

beginfig(1);
pickup pencircle scaled 0.5bp;
%
drawbase_meson;
%
dotlabel.top(btex $J_{\mu}(x)$ etex, z5);
dotlabel.bot(btex $J_{\nu}(y)$ etex, z6);
draw_fermion z1..z5..z5p..z3;
draw_fermion z4..z6..z6p..z2;
%
endfig;

beginfig(2);
%
drawbase_meson;
%
dotlabel.rt(btex $J_{\mu}(x)$ etex, z5);
dotlabel.lft(btex $J_{\nu}(y)$ etex, z6);
draw_fermion z1..z5..z6p..z2;
draw_fermion z4..z6..z5p..z3;
%
endfig;

beginfig(3);
%
drawbase_meson;
%
dotlabel.top(btex $J_{\mu}(x)$ etex, z5);
dotlabel.top(btex $J_{\nu}(y)$ etex, z5p);
draw_fermion z1..z5..z5p..z3;
draw_fermion z4..z6..z6p..z2;
%
endfig;

beginfig(4);
%
drawbase_meson;
%
z1l=(0.5(z5+z5p)+(0,5u));
lradius:=3u;
draw_fermion z1l..(z1l+(lradius,lradius))..(z1l+(2lradius,0lradius))..(z1l+(lradius,-lradius))..cycle;
dotlabel.lft(btex $J_{\mu}(x)$ etex, z1l);
dotlabel.bot(btex $J_{\nu}(y)$ etex, z5);
draw_fermion z1..z5..z5p..z3;
draw_fermion z4..z6..z6p..z2;
%
endfig;

beginfig(5);
%
drawbase_meson;
%
lradius:=3u;
z1l=(0.5(z5+z5p)+(0,5u) + (-lradius,0));
draw_fermion (z1l+(lradius,-0.8lradius))..z1l..(z1l+(lradius,0.8lradius))..(z1l+(2lradius,0lradius))..cycle;
dotlabel.lft(btex $J_{\mu}(x)$ etex, z1l);
dotlabel.rt(btex $J_{\nu}(y)$ etex, (z1l+(2lradius,0lradius)));
draw_fermion z1..z5..z5p..z3;
draw_fermion z4..z6..z6p..z2;
%
endfig;

beginfig(6);
%
drawbase_decay;
%
dotlabel.ulft(btex $J_{\mu}(x)$ etex, z5);
dotlabel.llft(btex $J_{\nu}(y)$ etex, z6p);
%
draw_fermion z4..z6..z6p..z5..z5p..z3;
%
endfig;



end
