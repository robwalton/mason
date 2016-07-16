Mason Rule Solver Program
~~~~~~~~~~~~~~~~~~~~~~~~~

Theory:
~~~~~~~

This program crunches down a signal flow graph to generate a symbolic 
equation for the equivalent phaser relating an output node and an input
node. For example:

   1     s21      2
    *------>------*------------.
    |             |            |
    |             |            | 
s11 v             ^ s22        v R2
    |             |            | 
    |             |            |  
    *------<------*------------'
   3     s12      4

There are four nodes 1,2,3 and 4, and there are 5 Coefficients S11,S21,
S12,S22 and R2. If we set node 1 as an independent input node, and choose
node 3 as a dependent node we get the following simplification:

           1
   a --->   *
            |
            |                      R2
            v     b/a= s11+s21* -------- *s12
            |                   1-s22*R2
            |           
   b <---   *
           3

If we set node 1 as the independent input node, and node 2 as the
dependent output node we get:

          1               2
   a --->  *------>------*  ---> b 

                 s21
         b/a = --------
               1-s22*R2

This program generates these equations for a given network and pair of nodes.


Specifying the Network
~~~~~~~~~~~~~~~~~~~~~~

A network description file is used to specify the topology of the
flow diagram. Each node is assigned a number, and each branch is
defined a coefficient number. Each branch is described as one line
in a file, defined as follows:

[Coefficient #]   [Start Node #]  [Stop Node #]   [Coefficient Name]

The Coefficient numbers must be in order (starting at 1 not 0). For example,
to describe the following diagram (where coefficient numbers are in brackets),

       1    (1) s21   2
        *------>------*------------.
        |             |            |
        |             |            | 
(4) s11 v             ^ (2) s22    v (5) R2
        |             |            | 
        |             |            |  
        *------<------*------------'
       3    (3) s12   4

we can use this file (where the white space are tabs, but any white 
space should suffice):

1	1	2	s21
2	4	2	s22
3	4	3	s12
4	1	3	s11
5	2	4	R2

This is saved as 'example.net' in this directory

The coefficient name can be any valid symbolic expression used in Maple. If
the expression has multiple terms then enclose them in brackets. For example:

   (-D*z^(-1)) or (1+B)


Using the Program
~~~~~~~~~~~~~~~~~

It is important that the lines in the net file be ordered so that the
coefficient numbers count from 1 up. Don't use 0 to number the coefficients
or nodes! Once you have made the net file, run 'mason.m' from Matlab, as 
described below:

USAGE:
   [Numerator,Denominator] = mason(Netfile,StartNode,StopNode)

   Netfile     - is a string with the name of the netfile to load
   StartNode   - is the integer number describing the independent input node
   StopNode    - is the integer number describing the dependent output node
   Numerator   - is a string containing the equation for the Numerator
   Denominator - is a string containing the equation for the Denominator 

Try out the example network! To recreate the above examples use:

   [Numerator,Denominator] = mason('example.net',1,3)
   [Numerator,Denominator] = mason('example.net',1,2)


File Details
~~~~~~~~~~~~

mason.m -- the mason's rule function. Type help mason for arguments.

example.net -- the example network described in the theory section.


About the program
~~~~~~~~~~~~~~~~~

I have tested this program with many different networks, and use it
a lot in my research. However, as with any simulator, don't assume it
is perfect! If you find any bugs, or if you find this program useful,
give me an email.

Rob Walton

walton@ieee.org
TRLabs, Calgary, Alberta, Canada.
