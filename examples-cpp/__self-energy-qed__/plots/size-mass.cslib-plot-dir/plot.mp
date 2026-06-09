%GNUPLOT Metapost output: Sat Sep 28 22:53:49 2019

prologues:=3;

%% Add \documentclass and \begin{dcoument} for latex
%% NB you should set the environment variable TEX to the name of your
%% latex executable (normally latex) inorder for metapost to work
%% or run
%% mpost --tex=latex ...

% BEGPRE
verbatimtex
\documentclass{article}
\usepackage[intlimits]{amsmath}
\usepackage{amsfonts}
\begin{document}
etex
% ENDPRE

warningcheck:=0;
defaultmpt:=mpt:=4;
th:=.6;
%% Have nice sharp joins on our lines
linecap:=butt;
linejoin:=mitered;

def scalepen expr n = pickup pencircle scaled (n*th) enddef;
def ptsize expr n = mpt:=n*defaultmpt enddef;


textmag:= 1.000;
vardef makepic(expr str) =
  if picture str : str scaled textmag
  % otherwise a string
  else: str infont defaultfont scaled (defaultscale*textmag)
  fi
enddef;

def infontsize(expr str, size) =
  infont str scaled (size / fontsize str)
enddef;

color currentcolor; currentcolor:=black;
picture currentdash; currentdash:=dashpattern(on 1);
color fillcolor;
boolean colorlines,dashedlines;
colorlines:=true;
dashedlines:=true;

def _wc = withpen currentpen withcolor currentcolor dashed currentdash enddef;
def _ac = addto currentpicture enddef;
def _sms = scaled mpt shifted enddef;
% drawing point-types
def gpdraw (expr n, x, y) =
  if n<0: _ac contour fullcircle _sms (x,y)
  elseif (n=1) or (n=3):
    _ac doublepath ptpath[n] _sms (x,y) _wc;
    _ac doublepath ptpath[n] rotated 90 _sms (x,y) _wc
  elseif n<6: _ac doublepath ptpath[n] _sms (x,y) _wc
  else: _ac contour ptpath[n] _sms (x,y) _wc
  fi
enddef;

% the point shapes
path ptpath[];
%diamond
ptpath0 = ptpath6 = (-1/2,0)--(0,-1/2)--(1/2,0)--(0,1/2)--cycle;
% plus sign
ptpath1 = (-1/2,0)--(1/2,0);
% square
ptpath2 = ptpath7 = (-1/2,-1/2)--(1/2,-1/2)--(1/2,1/2)--(-1/2,1/2)--cycle;
% cross
ptpath3 := (-1/2,-1/2)--(1/2,1/2);
% circle:
ptpath4 = ptpath8:= fullcircle;
% triangle
ptpath5 = ptpath9 := (0,1/2)--(-1/2,-1/2)--(1/2,-1/2)--cycle;

def linetype expr n =
  currentcolor:= if colorlines : col[n] else: black fi;
  currentdash:= if dashedlines : lt[n] else: dashpattern(on 1) fi;
  if n = -1 :
      drawoptions(withcolor currentcolor withpen (currentpen scaled .5));
  else :
    drawoptions(_wc);
  fi
enddef;

% dash patterns
picture lt[];
lt[-2]:=lt[-1]:=lt0:=dashpattern(on 1);
lt1=dashpattern(on 2 off 2); % dashes
lt2=dashpattern(on 2 off 2 on 0.2 off 2); %dash-dot
lt3=lt1 scaled 1.414;
lt4=lt2 scaled 1.414;
lt5=lt1 scaled 2;
lt6:=lt2 scaled 2;
lt7=dashpattern(on 0.2 off 2); %dots

color col[],cyan, magenta, yellow;
cyan=blue+green; magenta=red+blue;yellow=green+red;
col[-2]:=col[-1]:=col0:=black;
col1:=red;
col2:=(.2,.2,1); %blue
col3:=(1,.66,0); %orange
col4:=.85*green;
col5:=.9*magenta;
col6:=0.85*cyan;
col7:=.85*yellow;

%placing text
picture GPtext;
def put_text(expr pic, x, y, r, j) =
  GPtext:=makepic(pic);
  GPtext:=GPtext shifted
    if j = 1: (-(ulcorner GPtext + llcorner GPtext)/2)
    elseif j = 2: (-center GPtext)
    else: (-(urcorner GPtext + lrcorner GPtext)/2)
    fi
    rotated r;
  addto currentpicture also GPtext shifted (x,y)
enddef;

beginfig(0);
w:=5.000in;h:=3.000in;
a:=w/1200.0;b:=h/720.0;
scalepen 1; ptsize 1.000;linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (169.3a,117.1b)--(185.9a,117.1b);
draw (666.8a,117.1b)--(650.2a,117.1b);
put_text( btex $0.001$ etex, 151.6a, 117.1b, 0, 3);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (169.3a,148.6b)--(177.6a,148.6b);
draw (666.8a,148.6b)--(658.5a,148.6b);
draw (169.3a,167.0b)--(177.6a,167.0b);
draw (666.8a,167.0b)--(658.5a,167.0b);
draw (169.3a,180.0b)--(177.6a,180.0b);
draw (666.8a,180.0b)--(658.5a,180.0b);
draw (169.3a,190.2b)--(177.6a,190.2b);
draw (666.8a,190.2b)--(658.5a,190.2b);
draw (169.3a,198.4b)--(177.6a,198.4b);
draw (666.8a,198.4b)--(658.5a,198.4b);
draw (169.3a,205.4b)--(177.6a,205.4b);
draw (666.8a,205.4b)--(658.5a,205.4b);
draw (169.3a,211.5b)--(177.6a,211.5b);
draw (666.8a,211.5b)--(658.5a,211.5b);
draw (169.3a,216.9b)--(177.6a,216.9b);
draw (666.8a,216.9b)--(658.5a,216.9b);
draw (169.3a,221.6b)--(185.9a,221.6b);
draw (666.8a,221.6b)--(650.2a,221.6b);
put_text( btex $0.01$ etex, 151.6a, 221.6b, 0, 3);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (169.3a,253.1b)--(177.6a,253.1b);
draw (666.8a,253.1b)--(658.5a,253.1b);
draw (169.3a,271.5b)--(177.6a,271.5b);
draw (666.8a,271.5b)--(658.5a,271.5b);
draw (169.3a,284.6b)--(177.6a,284.6b);
draw (666.8a,284.6b)--(658.5a,284.6b);
draw (169.3a,294.7b)--(177.6a,294.7b);
draw (666.8a,294.7b)--(658.5a,294.7b);
draw (169.3a,303.0b)--(177.6a,303.0b);
draw (666.8a,303.0b)--(658.5a,303.0b);
draw (169.3a,310.0b)--(177.6a,310.0b);
draw (666.8a,310.0b)--(658.5a,310.0b);
draw (169.3a,316.0b)--(177.6a,316.0b);
draw (666.8a,316.0b)--(658.5a,316.0b);
draw (169.3a,321.4b)--(177.6a,321.4b);
draw (666.8a,321.4b)--(658.5a,321.4b);
draw (169.3a,326.2b)--(185.9a,326.2b);
draw (666.8a,326.2b)--(650.2a,326.2b);
put_text( btex $0.1$ etex, 151.6a, 326.2b, 0, 3);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (169.3a,357.6b)--(177.6a,357.6b);
draw (666.8a,357.6b)--(658.5a,357.6b);
draw (169.3a,376.0b)--(177.6a,376.0b);
draw (666.8a,376.0b)--(658.5a,376.0b);
draw (169.3a,389.1b)--(177.6a,389.1b);
draw (666.8a,389.1b)--(658.5a,389.1b);
draw (169.3a,399.2b)--(177.6a,399.2b);
draw (666.8a,399.2b)--(658.5a,399.2b);
draw (169.3a,407.5b)--(177.6a,407.5b);
draw (666.8a,407.5b)--(658.5a,407.5b);
draw (169.3a,414.5b)--(177.6a,414.5b);
draw (666.8a,414.5b)--(658.5a,414.5b);
draw (169.3a,420.6b)--(177.6a,420.6b);
draw (666.8a,420.6b)--(658.5a,420.6b);
draw (169.3a,425.9b)--(177.6a,425.9b);
draw (666.8a,425.9b)--(658.5a,425.9b);
draw (169.3a,430.7b)--(185.9a,430.7b);
draw (666.8a,430.7b)--(650.2a,430.7b);
put_text( btex $1$ etex, 151.6a, 430.7b, 0, 3);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (169.3a,117.1b)--(169.3a,133.7b);
draw (169.3a,430.7b)--(169.3a,414.1b);
put_text( btex $1$ etex, 169.3a, 80.5b, 0, 2);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (224.6a,117.1b)--(224.6a,133.7b);
draw (224.6a,430.7b)--(224.6a,414.1b);
put_text( btex $2$ etex, 224.6a, 80.5b, 0, 2);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (279.9a,117.1b)--(279.9a,133.7b);
draw (279.9a,430.7b)--(279.9a,414.1b);
put_text( btex $3$ etex, 279.9a, 80.5b, 0, 2);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (335.1a,117.1b)--(335.1a,133.7b);
draw (335.1a,430.7b)--(335.1a,414.1b);
put_text( btex $4$ etex, 335.1a, 80.5b, 0, 2);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (390.4a,117.1b)--(390.4a,133.7b);
draw (390.4a,430.7b)--(390.4a,414.1b);
put_text( btex $5$ etex, 390.4a, 80.5b, 0, 2);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (445.7a,117.1b)--(445.7a,133.7b);
draw (445.7a,430.7b)--(445.7a,414.1b);
put_text( btex $6$ etex, 445.7a, 80.5b, 0, 2);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (501.0a,117.1b)--(501.0a,133.7b);
draw (501.0a,430.7b)--(501.0a,414.1b);
put_text( btex $7$ etex, 501.0a, 80.5b, 0, 2);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (556.2a,117.1b)--(556.2a,133.7b);
draw (556.2a,430.7b)--(556.2a,414.1b);
put_text( btex $8$ etex, 556.2a, 80.5b, 0, 2);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (611.5a,117.1b)--(611.5a,133.7b);
draw (611.5a,430.7b)--(611.5a,414.1b);
put_text( btex $9$ etex, 611.5a, 80.5b, 0, 2);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (666.8a,117.1b)--(666.8a,133.7b);
draw (666.8a,430.7b)--(666.8a,414.1b);
put_text( btex $10$ etex, 666.8a, 80.5b, 0, 2);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (169.3a,430.7b)--(169.3a,117.1b)--(666.8a,117.1b)--(666.8a,430.7b)--(169.3a,430.7b);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
put_text( btex FV correction to $\Delta m_\pi$ (MeV) etex, 26.5a, 273.9b, -270, 2);
linetype -2;
put_text( btex $L$ (fm) etex, 418.0a, 25.6b, 0, 2);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
currentcolor:=0.5804*red+0*green+0.8275*blue;
put_text( btex FV error with $t_s = L/2$ etex, 544.0a, 469.1b, 0, 3);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
currentcolor:=0.5804*red+0*green+0.8275*blue;
draw (561.7a,469.1b)--(649.1a,469.1b);
draw (169.3a,413.1b)--(179.5a,401.8b)--(201.3a,382.7b)--(223.1a,366.4b)--(244.9a,352.2b)
--(266.7a,339.4b)--(288.5a,327.6b)--(310.4a,316.6b)--(332.2a,306.3b)--(354.0a,296.4b)
--(375.8a,287.0b)--(397.6a,278.0b)--(419.4a,269.3b)--(441.3a,260.8b)--(463.1a,252.6b)
--(484.9a,244.6b)--(506.7a,236.9b)--(528.5a,229.2b)--(550.3a,221.8b)--(572.1a,214.5b)
--(594.0a,207.2b)--(615.8a,200.1b)--(637.6a,193.1b)--(659.4a,186.2b)--(666.8a,183.8b);
scalepen 2.000;
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
scalepen 1.000;
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (169.3a,430.7b)--(169.3a,117.1b)--(666.8a,117.1b)--(666.8a,430.7b)--(169.3a,430.7b);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
endfig;
% BEGPOST
verbatimtex
 \end{document}
etex
% ENDPOST
end.
