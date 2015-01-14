#Introduction to the MatLab programming environment
###Pestilli, Franco K310 Spring 2014 Indiana University Bloomington

The goal of the first week is to get the student familiarized with the MatLab programming environment. You can practice the content of the tutorials using MatLab via IUAnyware (please see the Syllabus for more details, see also https://kb.iu.edu/d/bbbu for instructions).

Mathworks the company making MatLab provides useful videos to help familiarize with the langue:
http://www.mathworks.com/videos/getting-started-with-matlab-68985.html 
http://www.mathworks.com/videos/working-with-arrays-in-matlab-69022.html 
http://www.mathworks.com/videos/writing-a-matlab-program-69023.html 
http://www.mathworks.com/videos/introducing-matlab-fundamental-classes-data-types-68991.html 

The MatLab documentation is accessible online fro free: http://www.mathworks.com/help/matlab/

Mathworks also hosts a portal called Matlabcentral, where files and questions can be submitted to a comunity of users: http://www.mathworks.com/matlabcentral/
It is foten helpful to find that many of the problems we might encounter have been previosuly solved by other individuals.

### The MatLab environment 
####Command line 
The MatLab interface contains several panels (see video [above](http://www.mathworks.com/videos/getting-started-with-matlab-68985.html)). The core panel of matlab is the command line. This is the subpanel where we can type commands and MatLab will execute. Command line is indicated by the symbol:
```
>>
```
####Variables 
They represent in the memory of the computer (RAM) numbers, strings, vectors, matrices, etc. Variable have names that be a combination of letters and numbers, but that cannot start with a number, for example:
```
  >> this_is_one_matlab_variable = 1; 
  >> thisIsOneMatLabVariable = 2;
  >> ThisIs1MatLabVariable = 3;
```
  
Variable exist only as long as MatLab is running. As soon as MatLab is turned off all variables are wiped out of memory. VAriables can be saved onto a file on your hard disk (the file has extension .mat) and the file can be loaded later. When the file is loaded all the variabes saved in the file will be accessible into the MatLab workspace. 

####Workspace 
MatLab holds all the variables created always in memory. The entirety of the variables kept in memory at every moment is called the workspace. 

####Functions and scripts 
MatLab can store a series of commands in a file on saved on your hard disk. MatLab files are saved to disk with the extension *.m There are two fundamentalfile types in MatLab.

- Scripts: These are the simples form of a written set of instructions. Scripts are simple text files saved with a .m exstension. Scripts access variables in the workspace directly, meaning the if a variable called 'VAR' exists in a script it will exist also in the workspace.
- Functions: These are files saved with the same extension as scripts (.m) but they contain slightly more advanced coding syntax. Functions have a predefined set of input and output arguments. The arguments of a function are the set of variables that the function recognizes and can operate upon. This is a fundamental difference with scripts. If a variable called VAR exist inside a function the variable does not exist in the workspace unless it is retuned as an output. The majority of the operations that we perform in matlab is performed by functions, more often functions that come with the MatLab package.
 
