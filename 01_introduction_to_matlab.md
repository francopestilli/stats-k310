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
- Command line: The MatLab interface contains several panels (see video [above](http://www.mathworks.com/videos/getting-started-with-matlab-68985.html)). The core panel of matlab is the command line. This is the subpanel where we can type commands and MatLab wil execute. Command line is indicated by the symbol >>.
- Variables: Represent in the memory of the computer (RAM) numbers, strings, vectors, matrices etc. Variable have names that be a combination of letters and numbers, but that cannot start with a numer, for example:
  this_is_one_variable, 
  thisIsOneVariable 
  ThisIs1Variable. 
Variable exist only as long a sMatLab is running. As soon as MatLab is turned off all variables are wiped out of memory. This is a fundamental difference between variables and scripts or functions, see below.
- Workspace: MatLab holds all the variables created always in memory. The entirety of the variables kept in memroy is called the workspace. 
- Functions and scripts: MatLab can store a series of commands in a file. MatLab files can be saved to disk with the extension *.m There are two fundamental 
- Functions accept input and return output. There are built-in functions (e.g. inv) that have no
readable source code, functions that come with MATLAB (e.g. mean) that have readable
source code, and user-written functions. Importantly, functions cannot access nor modify the
workspace except through inputs and outputs.
- Scripts are text files that consist of a series of MATLAB commands. Scripts are similar to
user-written functions (which are also text files consisting of a series of MATLAB commands).
However, the crucial difference is that scripts do not operate on inputs and outputs. Rather,
scripts have direct access to the workspace, so they can freely read variables from the
workspace and write variables to the workspace.
- Functions and scripts written by the user are saved as .m files. These are text files.
- MATLAB saves variables into .mat files. These are binary files that can be read by MATLAB.

