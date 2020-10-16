function [Lims,marks]=colorspy(A,flag,ch,markersize,colors,LimBounds);
% function [Lims,marks]=colorspy(A,[flag],[ch],[markersize],[colors],[LimBounds]);
% defaults: flag=0, ch='.', markersize=[], colors='rmbcgyw';
%           markersize & colors must have lengths 0, 1, or 7.
% author: Daniel Boley, version of Tue Jun  1 11:14:22 CDT 1999
% --- version 1.2 - re-arranged how 'marks' is filled in.
% --- version 1.2.1 - fixed bug in setting multiple markersizes.
% --- version 1.2.2 - used only 2 return values; plot handles no longer cell.
% --- version 1.2.3 - fixed bug: error if any subrange was empty.
% --- version 1.2.4 - fixed bug: error if any subrange was empty -- again.
% --- version 1.2.5 - automated checking for colordef black vs colordef white.
% --- version 1.2.6 - moved some subplot commands to make it work with octave
% SPY the matrix A (just like MATLAB's own function "spy"), except
% that the dots are colored "rmbcgyw" (see PLOT) according each entry's value.
% The 3nd arg (ch) specifies the character to use (see PLOT).
% The limits are set so each color has about nnz(A)/7 entries.
% flag: 0 = divide the range of min to max value into equal parts (default).
%       1 = divide the colors so each is equally represented.
%       4, 5 = like 0, 1 (resp.), but just plot the matrix.
%       6, 7 = like 0, 1 (resp.), but just plot the scale.
%       above + 8  = like above, but draw the most positive items last.
%                    (default: draw the extreme values last.)
%       above + 16 = like above, but reverse the order the items are drawn.
%       above + 32 = like above, but use a log scale on the matrix values.
%                    (for this, recommend also option 8, yielding flag=40.)
%       above + 64 = like above, but always put one color bdry at zero.
%       above +128 : Plot 3D data: A has 3 cols: does plot(A(:,1),A(:,2),'?'),
%                  : where '?' is varied according to value in A(:,3).
% (actually, the 2nd & 3rd arguments can be supplied in either order.)
%
% Returns the lower & upper limits used for each color.
%
% Examples (first try these examples on a small but sparse matrix A).
%   colorspy(A)        % use defaults.
%   colorspy(A,'x')    % use a char bigger than dots to display entries
%   colorspy(A,40,'x') % use log scale, drawing biggest items last.
%   colorspy(A,'x',40) % same as above.
%   colorspy(A,40)     % use log scale, but use dots (for large dimensionality)
% To get a black and white plot, try the following (adjust the sizes to taste):
%   colorspy(A,[],'x',(1:7)*2,'w');
%   colorspy(A,[],'x',(1:7),'w');hold on;colorspy(A,[],'o',(1:7),'w');hold off
% or even swap the roles of character & color ("help plot" for character key):
%   colorspy(A,'w',[],[],'.o+xpsd');
% 3D data: assume x,y,z are compatible vectors. Then  (as examples)
%   colorspy([x,y,z],128)      or    colorspy([x,y,z],128+4,'x')     etc.
% will plot points on the x-y plane, using the z value to set the color.


washeld = ishold;

if nargin<6, LimBounds=[];end;
if nargin<5, colors=[];end; 
if nargin<4, markersize=[];  end;
if nargin<3, ch=[];  end;   
if nargin<2, flag=[];end;   

args={flag,ch,markersize,colors,LimBounds,'','',[],[]};
isstrargs=[];for i=1:length(args),isstrargs(i)=ischar(args{i});end;
strs=find(isstrargs);
nums=find(isstrargs==0);

ch=args{strs(1)};
colors=args{strs(2)};
flag=args{nums(1)};
markersize=args{nums(2)};
LimBounds=args{nums(3)};

if (all(ismember(colors,'rmbcgywk')) & length(ch)==0),
   if length(colors)==1, ch='.o+xpsd';else,ch='.';end;
elseif (all(ismember(ch,'rmbcgywk')) & length(colors)==0),
   if length(ch)==1, colors='.o+xpsd';else,colors='.';end;
end;

defaultcolors='ygcrmbw';

if length(colors)==0,colors=defaultcolors;  end;
if length(ch)==0,   ch='.'; end;
if length(flag)==0, flag=0; end;
if length(LimBounds)==0, SetLimP=0; else, SetLimP=1; end;

if flag>=128,ThreeDimP=1;flag=flag-128;else, ThreeDimP=0;  end;
if flag>=64,ShowSignP=1;flag=flag-64;else, ShowSignP=0;  end;
if flag>=32,LogP=1;flag=flag-32;else, LogP=0;  end;
if flag>=16,ReverseP=1;flag=flag-16;else, ReverseP=0;  end;
if flag>=8,OrderedP=1;flag=flag-8;else, OrderedP=0;  end;
if flag>=4,PlotBothP=0;flag=flag-4;else, PlotBothP=1;  end;
if flag>=2,PlotSortedP=1;flag=flag-2;else, PlotSortedP=0;  end;
if flag>=1,DivideByRangeP=0;flag=flag-1;else, DivideByRangeP=1;  end;

if nargin<1, %% just print out default flags and return
   disp('usage: colorspy(A,flag,ch,colors,markersize,LimBounds);');
   disp('       assumes "colordef black"');
   disp('optional arguments and their default values:');
   flag
   ch
   colors
   markersize
   LimBounds
   disp('set colors=''w'' to get black & white plot with ch=''.o+xpsd'' .');
   return
end

if length(markersize)==1, markersize=[1,1,1,1,1,1,1]*markersize; end;
if length(colors)==1, o=colors;colors=[o,o,o,o,o,o,o];  end;
if length(ch)==1,ch=[ch,ch,ch,ch,ch,ch,ch];  end;
defaultaxescolor=get(gcf,'defaultaxescolor');
if isequal(defaultaxescolor,[0 0 0]),
   ch(find(ch=='k'))='w';
   colors(find(colors=='k'))='w';
   set(gcf,'color',[0 0 0]);
end;
if isequal(defaultaxescolor,[1 1 1]),
   ch(find(ch=='w'))='k';
   colors(find(colors=='w'))='k';
end;

seq=[4,5,3,6,2,7,1];
if OrderedP, seq=([1,2,3,4,5,6,7]); end;
if ReverseP, seq=seq([7,6,5,4,3,2,1]); end;

if ThreeDimP,
   IA=A(:,1);JA=A(:,2);Anz=A(:,3);
   n=ceil(max(IA)-min(IA));p=ceil(max(JA)-min(JA));
   %Imin=(min(IA));IA=IA-Imin;n=ceil(max(IA));
   %Jmin=(min(JA));JA=JA-Jmin;p=ceil(max(JA));
else,
   [n,p]=size(A);
   [IA,JA,Anz]=find(A);
end;
find(isfinite(Anz));IA=IA(ans);JA=JA(ans);Anz=Anz(ans);
nz=length(Anz);
if nz==0,
  disp([int2str(n),' by ',int2str(p),...
     ' input matrix all zero - nothing to plot']);
  return;
end;
if LogP,Anz=log10(abs(Anz));  end;
truenz=nz;
if ~ThreeDimP,nz=nz+1;Anz(nz)=max(Anz);IA(nz)=n+1;JA(nz)=-1;  end;

if ~DivideByRangeP || PlotSortedP || PlotBothP,
   [Asort,Isort]=sort(Anz);
   JAsort=JA(Isort);
   IAsort=IA(Isort);
end;
if ~DivideByRangeP,
   Iu=[ceil([0,1,2,3,4,5,6]*(nz-1)/7)+1,nz];
   Limits=[Asort(Iu(1:7)),Asort(Iu(2:8))]';
else,
   AnzMax=max(Anz);
   AnzMin=min(Anz);
   diff=(AnzMax-AnzMin)/7;
   Lims=[AnzMin*[1,1,1,1,1,1,1]+diff*[0,1,2,3,4,5,6],AnzMax];
   Limits=[Lims(1:7);Lims(2:8)];
end;
if ShowSignP && Limits(1,1)<0 && Limits(2,7)>0,
   [o,i]=min(abs(Limits(2,1:6)));
   Limits(2,i)=0;Limits(1,i+1)=0;
end;
if SetLimP,
   LimBounds=LimBounds(:)';
   if ~isequal(LimBounds,sort(LimBounds)),
      error('colorspy: limits out of order');
   end;
   if LogP,LimBounds=log10(LimBounds); end;
   if length(LimBounds)<8,
      LimBounds=[-inf,LimBounds,inf,inf,inf,inf,inf,inf,inf];
   end;
   Limits=[LimBounds(1:7);LimBounds(2:8)];
end;

if ~washeld, oldorient=orient;clf; end;
marks=zeros(7,2);

if PlotSortedP || PlotBothP,
 if PlotBothP, if washeld, hold on;end; subplot(1,2,2); end;
 for i=seq,
   if     i==1, ii=find((Asort <  Limits(2,1)));
   elseif i==7, ii=find((Asort >= Limits(1,7)));
   else,        ii=find((Asort >= Limits(1,i)) & (Asort < Limits(2,i)));
   end;
   if ~isempty(ii)
      if LogP,M=semilogy((ii),10.^Asort(ii),[colors(i),ch(i)]);
      else, M=plot((ii),Asort(ii),[colors(i),ch(i)]);
      end
      marks(i,2)=M;
      hold on
   end
 end;
 % ans=axis;axis([0 (nz+1) ans(3:4)]);
 xlabel(['nonzeros: ',int2str(truenz),' out of ',int2str(n*p)]);
 if ThreeDimP,xlabel([int2str(truenz),' points']);  end
end;
if ~PlotSortedP || PlotBothP,
 if PlotBothP,if washeld, hold on;end; subplot(1,2,1);  end;
 for i=seq,
   if     i==1, ii=find((Anz <  Limits(2,1)));
   elseif i==7, ii=find((Anz >= Limits(1,7)));
   else,        ii=find((Anz >= Limits(1,i)) & (Anz < Limits(2,i))); end;
   if ~isempty(ii),
      M=plot(JA(ii),IA(ii),[colors(i),ch(i)]);
      marks(i,1)=M;
      hold on; if ~ThreeDimP, axis('ij');  end;
   end;
 end;
 vz=(max(IA)-min(IA))*.05;hz=(max(JA)-min(JA))*.05;
 if ThreeDimP,axis([min(JA)-hz,max(JA)+hz,min(IA)-vz,max(IA)+vz]);
	 else,axis([0 (p+1) 0 (n+1)]);  end;
 xlabel(['Spy: ',int2str(n),' by ',int2str(p), ...
	 ', fill = ',num2str(100*(truenz)/(n*p)) , '%']);
 if ThreeDimP,xlabel([int2str(truenz),' points']);  end
end;

if length(markersize)==7,
   for i=1:7,
      if marks(i,1), set(marks(i,1),'MarkerSize',markersize(i)); end
      if marks(i,2), set(marks(i,2),'MarkerSize',markersize(i)); end
   end; 
end;

if ~washeld , hold off;  orient(oldorient); end;
if LogP,Limits=10.^(Limits);  end;
Lims=[Limits(1,1:7),Limits(2,7)];
if nargout==0,Lims=[(cellstr([ch;colors]'))';num2cell(Limits)];  end;
