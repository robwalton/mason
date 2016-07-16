function [Num,Den] = mason(NetFile,Start,Stop)
% mason.m
% This function takes a netfile describing a signal flow graph
% with symbolic coefficients and generates an equation representing
% the equivilent term between an independent input node, and dependent
% output node. Please see the *readme* file for a full description, and 
% an example netfile.
%
% Author :  Rob Walton
% Organisation : TRLabs and the University of Calgary, Alberta, Canada
% Date  :  January 25th 1999
% Revised : January 20th 2000 (Removed bug that caused the odd loop to be thrown away)
%
% Please email me at <walton@ieee.org> if you find any errors!
%
% USAGE:
%   [Numerator,Denominator] = mason(Netfile,StartNode,StopNode)
%
%   Netfile     - is a string with the name of the netfile to load
%   StartNode   - is the integer number describing the independent input node
%   StopNode    - is the integer number describing the dependent output node
%   Numerator   - is a string containing the equation for the Numerator
%   Denominator - is a string containing the equation for the Denominator 


%*** Load the the net description file into Net ***
% *** The first column in the file is the coefficient number (These must be in order starting
% at 1). The second and third column are start node and stop node numbers respectively. The last
% column is the name of the coefficient. ** See readme file **

% *** Net will be the first 3 columns of the file. Coeff will be the last column. ***
fid=fopen(NetFile);			% Open the file
if (fid==-1)
  fprintf('\n*** File, %s, not found ***\n\n',NetFile)
  return
end

Net=[];					% Initialize the Numerical Net matrix
line_number=0;				% Line number read from file
Coeff_Names={};				% Initialize cell array of strings

while 1 				% Loop Until end of file, read one line at a time
   line_number=line_number+1;           % Count which line we are on
   x=fscanf(fid,'%d',3);		% Read the three decimal numbers into x
   Coeff=fscanf(fid,'%s\n',1);		% Read the one coefficient name into coeff
   if isempty(x)			% Stop the loop if no data is left in file
     break
   end   
   Net(line_number,:)=transpose(x);     % Append one row to bottom of numerical matrix
   Coeff_Names{line_number}= Coeff;     % Append the coefficient name to the coefficient array
end
fclose(fid);				% Remember to close the file!
			

%*** Determine Number of Coefficients in net file ***
temp=size(Net);				% Determine the number of rows in the net matrix
Number_Coeff=temp(1);


%*** Find all the paths connecting the start and end nodes, along which ***
%*** no node is crossed more than once                                  ***

[PathCoeffList,PathNodeList]=findpaths(Start,Stop,[],[],Net);

% PathCoeffList and PathNodeList are matrixes where each row lists the number of the coefficients or the nodes visited respectively. Each row is padded with zeros to make them the same length.

%fprintf('\n- Path List -\n');
%print_paths(PathCoeffList);

%*** Determine all the first-Order loops ***
LoopCoeffList=[];			% Initialise list of coefficients in for each loop
LoopNodeList=[];			% Initialise list of nodes for each loop found

for index=1:Number_Coeff;      		% Find list of loops originating from each node
	% Get the matrix describing all the loops from at node #index
	[FoundLoops,FoundNodes]=findpaths(index,index,[],[],Net);
	LoopCoeffList=[LoopCoeffList;FoundLoops];	% Append Coefficients of loops
	LoopNodeList=[LoopNodeList;FoundNodes];		% Append nodes visited for each loop
end



% Remove duplicate loops
[LoopCoeffList,LoopNodeList]=RemoveDuplicateLoops(LoopCoeffList,LoopNodeList);

%fprintf('\n\n- 1st Order Loop List -\n');
%print_paths(LoopCoeffList);


%*** Convert to nomenclature used by Pozars RF Book  ***
% P{n} represents the nth path found connecting start to stop
% P{n}.Coeff is an array with a list of the coefficients passed through in order
% P{n}.Nodes is an array listing the nodes passed through in order. Including the end nodes
% NOTE: A cell array is used because the path lengths can be different resulting in different sized arrays



% *** Make a cell array of P the different length paths ***
temp=size(PathCoeffList);			% Determine number of paths
NumberPaths=temp(1);
if (NumberPaths==0);
  fprintf('\n*** There are no paths connecting those nodes ***\n')
  return
end

for index=1:NumberPaths				% Do each path separately
  Coeff=PathCoeffList(index,:);			% Read Coefficients for a path
  P{index}.Coeff=Coeff(1:sum(Coeff>0));		% Strip trailing zeros and put in struct
  Node=PathNodeList(index,:);			% Read node list for a path
  P{index}.Node=[Node(1:sum(Coeff>0)),Stop];    % Append trailing zeros and put in struct.
						% Append the Stop node onto the end of the node list
end


% *** Make a cell array of the first order paths, each row a different order ***
% *** The first column contains the number of paths of that order
% L{LoopOrder}.NumberLoops = number of loops of this order
% L{LoopOrder}.Coeff{n} = Coefficients of nth loop
% L{LoopOrder}.Node{n}  = Nodes of nth loop



temp=size(LoopCoeffList);
NumberLoops=temp(1);				% Determine number of first order paths
L{1}.NumberLoops=NumberLoops;			% Set number of loops in the L{1} struct
for index=1:NumberLoops				% Do each loop seperately
  Coeff=LoopCoeffList(index,:);			% Read coefficients for that loop
  L{1}.Coeff{index}=Coeff(1:sum(Coeff>0));	% Strip Trailing zeros and put in struct
  Node=LoopNodeList(index,:);			% Read Node list for loop
  L{1}.Node{index}=[Node(1:sum(Coeff>0)),Node(1)]; % Strip trailing zeros and put in struct
						% Append Stop node (same as first node in list
end




%*** Determine nth order loops ***
n=1;						% n is the order of loops we are finding

while 1  					% Loop until an order of loop is reached that is empty
  n=n+1;					% Count which order we are on
  L{n}.NumberLoops=0;				% Assume no loops of this order
 
  % compare each first order loop with each n-1th loop. If non touching add to the 
  % two combined to the nth loop.
  
  for first=1:L{1}.NumberLoops			% Take each first order loop
    for second=1:L{n-1}.NumberLoops		% Compare with each n-1th loop

      if not(AreTouchingLoops(L{1}.Node{first},L{n-1}.Node{second})) % Non Touching loops found
	% Determine if loop is a duplicate
	Duplicate=0;                       % A flag to indicate loop found is duplicate(0=not dup)
	for index=1:L{n}.NumberLoops       % Add this loop if it is not a duplicate entry
          if IsSameLoop([L{1}.Coeff{first}, L{n-1}.Coeff{second}],L{n}.Coeff{index}) %Duplicate found
            Duplicate=1;        		% Set the duplicate flag
          end
        end
        if (Duplicate==0)                       % Add the loop if not a duplicate
	  L{n}.NumberLoops=L{n}.NumberLoops+1;	% Increment the number of loops of that order
	  % For Node and Coeff structs. Append a new array describing the loop of order n found
          L{n}.Coeff{(L{n}.NumberLoops)}=[L{1}.Coeff{first}, L{n-1}.Coeff{second}];
          L{n}.Node{(L{n}.NumberLoops)}=[L{1}.Node{first}, L{n-1}.Node{second}];
        end
      end 
    end
  end

  if (L{n}.NumberLoops==0)			% If no loops of order n where found, then break
    break					% There will be no loops of order n+1
  end
end

% ***** Display File info *****
fprintf('\n-- Network Info --\n')
fprintf('Net File   : ');fprintf(NetFile);fprintf('\n');
fprintf('Start Node : %d\n',Start);
fprintf('Stop Node  : %d\n',Stop);


% ***** Display the paths found ******
fprintf('\n----- Paths -----\n')
for pathn=1:length(P)				% Look at each Path and display it's Coeff numbers
  fprintf('P%d : ',pathn);                       % on a different line
  fprintf('%d ',P{pathn}.Coeff);
  fprintf('\n');
end

% ***** Display all the loops found *****

for loop_order=1:length(L)-1            	% Look at each loop order (last order has no loops
  fprintf('\n- Order %d Loops -\n',loop_order)  % Print header describing loop order
  for loop_number=1:L{loop_order}.NumberLoops   % Look at each loop of that order
    fprintf('L%d%d : ',loop_order,loop_number)  % Display coefficients on a different line
    fprintf('%d ',L{loop_order}.Coeff{loop_number})
    fprintf('\n')
  end
end  


% *******************************************************
% ************ Generate the final equation **************
% *******************************************************


% For now the equations are written in terms of the coefficient number : c#
% the coefficient strings will be substituted later

% Determine Numerator
Num='';				% Assume Numerator is empty to start
for pathn=1:length(P)            % Add Each path and related 
  Num=sprintf('%s%s*(1', Num, CoeffToString(P{pathn}.Coeff));            %    Pn*(1 ...  
  for order=1:length(L)-1       % Append each order of sums of non-touching loops
    % if order is odd order append a minus, otherwise a plus
    if (rem(order,2)==1)
      Num=sprintf('%s-',Num);
    else
      Num=sprintf('%s+',Num);
    end
    % Append the sum of all the nontouching loops that  don't touch the current path
    Num=[Num,PrintSumsNotTouching(L,order,P{pathn}.Node)];
  end
  Num=sprintf('%s)+',Num);   	% Close the bracket around paths list of sums
end



Num=Num(1:length(Num)-1);       % Remove the extra plus sign on the end.NOTE using /b screws up the symbolic math later

% Determine Denominator
Den='1';			% Denominator always start with a zero
for order=1:length(L)-1		% Add order subtract the sum of each orders loops
  % if order is odd order append a minus, otherwise a plus, 
  if (rem(order,2)==1)
    Den=sprintf('%s-',Den);
  else
    Den=sprintf('%s+',Den);
  end
  %Add or subtract all the loops
  % KLUDGE: all the loops not touching the path with nodes 999999999 are added
  %         That definetly should be all of them! 
  Den=[Den,PrintSumsNotTouching(L,order,[9999999 999999])];  %Sums of all the loops of order order
end


fprintf('\nThe variables returned are strings describing\n')
fprintf('the numerator and Denominator of the transfer equation.\n')
fprintf('If you have the symbolic toolbox, use Denominator=sym(Denominator)\n');
fprintf('and Numerator=sym(Numerator) to make these symbolic equations.\n')
fprintf('You can now use simple(Numerator/Denominator) to boil the whole\n')
fprintf('thing down. You could also use simple(Numerator) to simplify the\n')
fprintf('Numerator on it'' own.\n\n')
% ****** Convert to Symbolic and do substitutions *******

for coeff_num=length(Coeff_Names):-1:1;	%Cycle through Coefficient array, substituting each one
  orig=sprintf('c%d',Net(coeff_num,1)); % for each line generate c[Coeff Number] to replace  
  Den=strrep(Den,orig,Coeff_Names{coeff_num});	% Replace all the c#s with the strings from net file
  Num=strrep(Num,orig,Coeff_Names{coeff_num});
end % This loop had to count down so there was no risk of C12 being replace by C1


%*************************************************************************************************
function Touching=AreTouchingLoops(Nodes1,Nodes2)
%*************************************************************************************************
% This function takes two arrays describing two sets of nodes visited(each padded with zeros).
% Return 1 if they are they are touching loops.

% Determine length of loop arrays with zeros removed
Loop1Length=sum(Nodes1>0);
Loop2Length=sum(Nodes2>0);

for first=1:Loop1Length
	for second=1:Loop2Length
		if (Nodes1(first)==Nodes2(second))
			Touching=1;
			return;
		end
	end
end

Touching=0;


%*************************************************************************************************
function StrMult=CoeffToString(Coefficients)
%*************************************************************************************************
% StrMult=CoeffToString(Coefficients)
% Coefficients is an array with coefficients c1,c2..cN

N=length(Coefficients);     			% Get length of string
StrMult=sprintf('c%d',Coefficients(1));		% Start with first coefficient
for n=2:N		% Append each coefficent in list with * before it
  StrMult=[StrMult, sprintf('*c'),sprintf('%d',Coefficients(n))];
end

%*************************************************************************************************
function [PathUp,NodesUp]=findpaths(StartNode,StopNode,Path,Nodes,Net)
%*************************************************************************************************
%[PathUp,NodesUp]=findpaths(StartNode,StopNode,Path,Nodes,Net)
%
%Iterative function to find path between StartNode and StopNode. Net is the array with the network
%list in it. Path is the single path to date for a given route through the tree. PathUp is a
%list of all paths terminated below that node that are sucesfull.
%Nodes is the list of nodes tvaersed so far on the way down

% Determine number of coefficients in net
temp=size(Net);
NumberCoeff=temp(1,1);

PathUp=[];
NodesUp=[];

% Terminate branch and return nothing if the Nodes to date contains repetitions.
for index=1:NumberCoeff
	if not(isempty(Nodes))  % Only compare if the the Path has data in it
		if (sum(Nodes==index)>1)
			PathUp=[];
		%	fprintf('Repeated Node : ');
		%	fprintf('%d ',Nodes);
		%	fprintf('\n');
			return
		end
	end
end

% Terminate branch and return path if start and stop nodes are the same
if ((StartNode==StopNode) & (length(Path>1)))
	PathUp=Path;
	NodesUp=Nodes;
	%fprintf('Sucessfull Path : ');
	%fprintf('%d ',Path);
	%fprintf('\n');
	return
end


% Check for all branches leaving StartNode, and call another iteration for them
for index=1:NumberCoeff
	if (StartNode==Net(index,2)) 
		% Iterate with appended coeff to path and new startnode
		[FoundPath,FoundNodes]=findpaths(Net(index,3),StopNode,[Path,Net(index,1)],[Nodes,StartNode],Net);
	        if not(isempty(FoundPath))  %Only append if not empty
			PathUp=[PathUp;[FoundPath,zeros(1,NumberCoeff+1-length(FoundPath))]];
                        NodesUp=[NodesUp;[FoundNodes,zeros(1,NumberCoeff+1-length(FoundPath))]];
		end
	end
end
%*************************************************************************************************
function Same=IsSameLoop(Loop1,Loop2)
%*************************************************************************************************
% This function takes two arrays describing two loops(Can be padded with zeros if desired).
% Return 1 if they are they describe the same circular loop.

% Determine length of loop arrays with zeros removed
Loop1Length=sum(Loop1>0);		% Count all the non zero terms
Loop2Length=sum(Loop2>0);

%Return 0 if different length since the loops can't be the same!
if (Loop1Length~=Loop2Length)
	Same=0;
	return
end

%They are the same length so see if the contain the same nodes, but in any order.

% sort the nodes and subtract the two vectors. The resulting vector componments will all be zero, only if the lists contains the same values.

if (sum(abs(sort(Loop1)-sort(Loop2)))==0)
	Same=1;				% Loops are the same
else
	Same=0;				% Loops are different
end


  
%*************************************************************************************************
function Str=PrintSumsNotTouching(L,order,Pnodes)
%*************************************************************************************************
% Str=PrintSumsNotTouching(L,path)
% L is the array of structs containing all the first and higher order loops.
% Pnodes is the array of nodes that decibing a path
%
% The function returns a string with the sum off all the loops of order order
% that do not touch the path. The sum is contained within a set of brackets

No_NonTouching=1;  			% Flag set so indacate no nontouching loops found
Str=('(');				% Open the first bracket
for n=1:L{order}.NumberLoops		% Look at each loop of thet order
  if not(AreTouchingLoops(Pnodes,L{order}.Node{n}))   		%The loop doesn't touch the path
    Str=sprintf('%s%s+',Str,CoeffToString(L{order}.Coeff{n}));% So add its coefficients
    No_NonTouching=0;			% Set flag to indacet a nontouching loop was found
  end
end
Str=Str(1:(length(Str)-1));		% Remove the extra plus sign (or open bracket
					% if not loops where found
Str=sprintf('%s)',Str);			% Append the closed bracket
%If the sum is zero return zero instead

if No_NonTouching==1			% If no loops foun then return '0' instead
  Str='0';
end

%*************************************************************************************************
function [LoopList,NodeList]=RemoveDuplicateLoops(LoopList,NodeList);
%*************************************************************************************************
% [LoopList,NodeList]=RemoveDuplicateLoops(LoopList,NodeList)
% Returns a subset of the LoopList matrix, where each duplicate row has been removed
% This function works on the initial loops description where each loop is a row in the
% the matrix padded with zeros to make them all the same length

temp=size(LoopList);		% Get number of loops
NumberLoops=temp(1);

% Compare each loop with all the loops after it. And remove its duplicates from after it.
first=1;			% The first loop
while (first<=NumberLoops)      % Loop until all the loops have been used as first loops
	second=first+1;		% Start the second loop to compare as being one ahead of first loop
	while (second<=NumberLoops)  % choose the second loop to compare as all the ones after first
		% Remove the extra loop found at the second loop index.
		% NOTE: A for loop was not used since the number of loops is not constant
		if (IsSameLoop(LoopList(first,:),LoopList(second,:))==1) %Remove row at second loop
 			LoopList=[LoopList(1:second-1,:);LoopList(second+1:NumberLoops,:)];
			NodeList=[NodeList(1:second-1,:);NodeList(second+1:NumberLoops,:)];
			NumberLoops=NumberLoops-1;       % Decrement the number o loops
		else
			second=second+1;	% Only increment second if no loop was removed
						% Otherwise a new loop is move to second
		end
	end
	first=first+1;				% increment the first loop pointer
end	
	






