%GNUPLOT Metapost output: Sat Sep 28 22:53:50 2019

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
put_text( btex $-0.3$ etex, 151.6a, 117.1b, 0, 3);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (169.3a,160.6b)--(185.9a,160.6b);
draw (666.8a,160.6b)--(650.2a,160.6b);
put_text( btex $-0.25$ etex, 151.6a, 160.6b, 0, 3);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (169.3a,204.0b)--(185.9a,204.0b);
draw (666.8a,204.0b)--(650.2a,204.0b);
put_text( btex $-0.2$ etex, 151.6a, 204.0b, 0, 3);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (169.3a,247.5b)--(185.9a,247.5b);
draw (666.8a,247.5b)--(650.2a,247.5b);
put_text( btex $-0.15$ etex, 151.6a, 247.5b, 0, 3);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (169.3a,291.0b)--(185.9a,291.0b);
draw (666.8a,291.0b)--(650.2a,291.0b);
put_text( btex $-0.1$ etex, 151.6a, 291.0b, 0, 3);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (169.3a,334.5b)--(185.9a,334.5b);
draw (666.8a,334.5b)--(650.2a,334.5b);
put_text( btex $-0.05$ etex, 151.6a, 334.5b, 0, 3);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (169.3a,377.9b)--(185.9a,377.9b);
draw (666.8a,377.9b)--(650.2a,377.9b);
put_text( btex $0$ etex, 151.6a, 377.9b, 0, 3);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (169.3a,421.4b)--(185.9a,421.4b);
draw (666.8a,421.4b)--(650.2a,421.4b);
put_text( btex $0.05$ etex, 151.6a, 421.4b, 0, 3);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (169.3a,464.9b)--(185.9a,464.9b);
draw (666.8a,464.9b)--(650.2a,464.9b);
put_text( btex $0.1$ etex, 151.6a, 464.9b, 0, 3);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (169.3a,117.1b)--(169.3a,133.7b);
draw (169.3a,464.9b)--(169.3a,448.3b);
put_text( btex $1$ etex, 169.3a, 80.5b, 0, 2);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (231.5a,117.1b)--(231.5a,133.7b);
draw (231.5a,464.9b)--(231.5a,448.3b);
put_text( btex $1.5$ etex, 231.5a, 80.5b, 0, 2);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (293.7a,117.1b)--(293.7a,133.7b);
draw (293.7a,464.9b)--(293.7a,448.3b);
put_text( btex $2$ etex, 293.7a, 80.5b, 0, 2);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (355.9a,117.1b)--(355.9a,133.7b);
draw (355.9a,464.9b)--(355.9a,448.3b);
put_text( btex $2.5$ etex, 355.9a, 80.5b, 0, 2);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (418.1a,117.1b)--(418.1a,133.7b);
draw (418.1a,464.9b)--(418.1a,448.3b);
put_text( btex $3$ etex, 418.1a, 80.5b, 0, 2);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (480.2a,117.1b)--(480.2a,133.7b);
draw (480.2a,464.9b)--(480.2a,448.3b);
put_text( btex $3.5$ etex, 480.2a, 80.5b, 0, 2);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (542.4a,117.1b)--(542.4a,133.7b);
draw (542.4a,464.9b)--(542.4a,448.3b);
put_text( btex $4$ etex, 542.4a, 80.5b, 0, 2);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (604.6a,117.1b)--(604.6a,133.7b);
draw (604.6a,464.9b)--(604.6a,448.3b);
put_text( btex $4.5$ etex, 604.6a, 80.5b, 0, 2);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (666.8a,117.1b)--(666.8a,133.7b);
draw (666.8a,464.9b)--(666.8a,448.3b);
put_text( btex $5$ etex, 666.8a, 80.5b, 0, 2);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (169.3a,464.9b)--(169.3a,117.1b)--(666.8a,117.1b)--(666.8a,464.9b)--(169.3a,464.9b);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
put_text( btex FV correction to $\Delta m_\pi$ (MeV) etex, 26.5a, 291.0b, -270, 2);
linetype -2;
put_text( btex $t_s$ (fm) etex, 418.0a, 25.6b, 0, 2);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
currentcolor:=0.5804*red+0*green+0.8275*blue;
draw (169.3a,378.0b)--(174.3a,378.0b)--(179.4a,378.0b)--(184.4a,378.0b)--(189.4a,378.0b)
--(194.4a,378.0b)--(199.5a,378.0b)--(204.5a,378.0b)--(209.5a,378.0b)--(214.5a,378.0b)
--(219.6a,378.0b)--(224.6a,378.0b)--(229.6a,378.0b)--(234.6a,378.0b)--(239.7a,378.0b)
--(244.7a,378.0b)--(249.7a,378.0b)--(254.7a,378.0b)--(259.8a,378.0b)--(264.8a,378.0b)
--(269.8a,378.0b)--(274.8a,378.0b)--(279.9a,378.0b)--(284.9a,378.0b)--(289.9a,378.0b)
--(294.9a,378.0b)--(300.0a,378.0b)--(305.0a,378.0b)--(310.0a,378.0b)--(315.0a,378.0b)
--(320.1a,378.0b)--(325.1a,378.0b)--(330.1a,378.0b)--(335.1a,378.0b)--(340.2a,378.0b)
--(345.2a,378.0b)--(350.2a,378.0b)--(355.2a,378.0b)--(360.3a,378.0b)--(365.3a,378.0b)
--(370.3a,378.0b)--(375.3a,378.0b)--(380.4a,378.0b)--(385.4a,378.0b)--(390.4a,378.0b)
--(395.4a,378.0b)--(400.5a,378.0b)--(405.5a,378.0b)--(410.5a,378.0b)--(415.5a,378.0b)
--(420.6a,378.0b)--(425.6a,378.0b)--(430.6a,378.0b)--(435.6a,378.0b)--(440.7a,378.0b)
--(445.7a,378.0b)--(450.7a,378.0b)--(455.7a,378.0b)--(460.8a,378.0b)--(465.8a,378.0b)
--(470.8a,378.0b)--(475.8a,378.0b)--(480.9a,378.0b)--(485.9a,378.0b)--(490.9a,378.0b)
--(495.9a,378.0b)--(501.0a,378.0b)--(506.0a,378.0b)--(511.0a,378.0b)--(516.0a,378.0b)
--(521.1a,378.0b)--(526.1a,378.0b)--(531.1a,378.0b)--(536.1a,378.0b)--(541.2a,378.0b)
--(546.2a,378.0b)--(551.2a,378.0b)--(556.2a,378.0b)--(561.3a,378.0b)--(566.3a,378.0b)
--(571.3a,378.0b)--(576.3a,378.0b)--(581.4a,378.0b)--(586.4a,378.0b)--(591.4a,378.0b)
--(596.4a,378.0b)--(601.5a,378.0b)--(606.5a,378.0b)--(611.5a,378.0b)--(616.5a,378.0b)
--(621.6a,378.0b)--(626.6a,378.0b)--(631.6a,378.0b)--(636.6a,378.0b)--(641.7a,378.0b)
--(646.7a,378.0b)--(651.7a,378.0b)--(656.7a,378.0b)--(661.8a,378.0b)--(666.8a,378.0b);
linetype -2;
%MP_dashtype% type = 0
currentcolor:=0.3373*red+0.7059*green+0.9137*blue;
put_text( btex $\mathcal{I}^{\phantom{(s)}}$ $L=4.7$ fm etex, 544.0a, 613.1b, 0, 3);
linetype -2;
%MP_dashtype% type = 0
currentcolor:=0.3373*red+0.7059*green+0.9137*blue;
draw (561.7a,613.1b)--(649.1a,613.1b);
draw (169.3a,400.8b)--(192.2a,403.3b)--(216.7a,405.6b)--(241.3a,407.7b)--(265.8a,409.7b)
--(290.4a,411.5b)--(314.9a,413.2b)--(339.4a,414.7b)--(364.0a,416.0b)--(388.5a,417.3b)
--(413.1a,418.4b)--(437.6a,419.4b)--(462.1a,420.4b)--(486.7a,421.2b)--(511.2a,422.0b)
--(535.8a,422.7b)--(560.3a,423.3b)--(584.9a,423.9b)--(609.4a,424.5b)--(633.9a,425.0b)
--(658.5a,425.4b)--(666.8a,425.5b);
linetype -2;
%MP_dashtype% type = 1
linetype 1;
currentcolor:=0.3373*red+0.7059*green+0.9137*blue;
put_text( btex $\mathcal{I}^{(s)}$ $L=4.7$ fm etex, 544.0a, 576.5b, 0, 3);
linetype -2;
%MP_dashtype% type = 1
linetype 1;
currentcolor:=0.3373*red+0.7059*green+0.9137*blue;
draw (561.7a,576.5b)--(649.1a,576.5b);
draw (169.3a,197.2b)--(192.2a,219.7b)--(216.7a,238.7b)--(241.3a,254.2b)--(265.8a,267.3b)
--(290.4a,278.5b)--(314.9a,288.3b)--(339.4a,296.9b)--(364.0a,304.6b)--(388.5a,311.4b)
--(413.1a,317.6b)--(437.6a,323.3b)--(462.1a,328.4b)--(486.7a,333.1b)--(511.2a,337.4b)
--(535.8a,341.4b)--(560.3a,345.1b)--(584.9a,348.5b)--(609.4a,351.7b)--(633.9a,354.6b)
--(658.5a,357.4b)--(666.8a,358.2b);
linetype -2;
%MP_dashtype% type = 0
currentcolor:=0*red+0.6196*green+0.451*blue;
put_text( btex $\mathcal{I}^{\phantom{(s)}}$ $L=6.3$ fm etex, 544.0a, 539.9b, 0, 3);
linetype -2;
%MP_dashtype% type = 0
currentcolor:=0*red+0.6196*green+0.451*blue;
draw (561.7a,539.9b)--(649.1a,539.9b);
draw (169.3a,385.6b)--(192.2a,386.8b)--(216.7a,387.9b)--(241.3a,388.9b)--(265.8a,389.8b)
--(290.4a,390.7b)--(314.9a,391.6b)--(339.4a,392.4b)--(364.0a,393.1b)--(388.5a,393.9b)
--(413.1a,394.5b)--(437.6a,395.2b)--(462.1a,395.8b)--(486.7a,396.3b)--(511.2a,396.8b)
--(535.8a,397.3b)--(560.3a,397.8b)--(584.9a,398.2b)--(609.4a,398.6b)--(633.9a,399.0b)
--(658.5a,399.3b)--(666.8a,399.4b);
linetype -2;
%MP_dashtype% type = 1
linetype 1;
currentcolor:=0*red+0.6196*green+0.451*blue;
put_text( btex $\mathcal{I}^{(s)}$ $L=6.3$ fm etex, 544.0a, 503.3b, 0, 3);
linetype -2;
%MP_dashtype% type = 1
linetype 1;
currentcolor:=0*red+0.6196*green+0.451*blue;
draw (561.7a,503.3b)--(649.1a,503.3b);
draw (169.3a,189.4b)--(192.2a,210.5b)--(216.7a,228.0b)--(241.3a,242.1b)--(265.8a,253.9b)
--(290.4a,263.9b)--(314.9a,272.6b)--(339.4a,280.1b)--(364.0a,286.9b)--(388.5a,292.9b)
--(413.1a,298.3b)--(437.6a,303.2b)--(462.1a,307.7b)--(486.7a,311.9b)--(511.2a,315.7b)
--(535.8a,319.2b)--(560.3a,322.4b)--(584.9a,325.5b)--(609.4a,328.3b)--(633.9a,330.9b)
--(658.5a,333.4b)--(666.8a,334.1b);
scalepen 2.000;
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
scalepen 1.000;
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
draw (169.3a,464.9b)--(169.3a,117.1b)--(666.8a,117.1b)--(666.8a,464.9b)--(169.3a,464.9b);
linetype -2;
%MP_dashtype% DASHTYPE_SOLID
endfig;
% BEGPOST
verbatimtex
 \end{document}
etex
% ENDPOST
end.
