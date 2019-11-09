#!/bin/bash

##############################################################################################################
#                                                                                                            #
#     mgrid                                                                                                  #
#     Copyright 2010-2013 Matthew Denwood - matthewdenwood@mac.com - Version 4.01 - January 2013             #
#                                                                                                            #
#     This file is part of runjags                                                                           #
#                                                                                                            #
#     This script is a replacement for xgrid -job submit that provides some extra features, including        #
#     support for jobs with multiple tasks, use of environmental variables and easier targeting of jobs to   #
#     specific nodes or hardware requirements.  To run, this script requires XGRID_CONTROLLER_HOSTNAME and   #
#     XGRID_CONTROLLER_PASSWORD to be set as environmental variables.  You should also put this script in,   #
#     for example, /usr/local/bin/ by issuing the following commands in Terminal:                            #
#                                                                                                            #
#     [cd <path_to_download_folder>/runjags/inst/xgrid/]                                                     #
#     sudo cp mgrid.sh /usr/local/bin/mgrid                                                                  #
#                                                                                                            #
#     As well as being accessible to all functions in the runjags package, the script will then be invokable #
#     from the command line by using the command 'mgrid' (see the options by typing 'mgrid -?').  This isn't #
#     done as part of the default package install process, and will require an administrator's password.     #
#                                                                                                            #
#     This script also includes default ART script that ensures jobs are preferentially sent to agents with  #
#     more free processors (taking into account the number of jobs already running), some free memory, and   #
#     preferring server nodes (any node running Mac OS X server).  Additional information detailing the      #
#     amount of available RAM, processor speed, and 32/64 bit availability on the nodes is collected but     #
#     not used by default.  To change the way that tasks are distributed among the available nodes, you can  #
#     specify an alternative scoring script using the -z option.  An example scoring script is included in   #
#     the xgrid folder of the runjags package (note that this script uses some of the node specification     #
#     information which is passed to scoring scripts as arguments).                                          #
#                                                                                                            #
#     [Note: the Xgrid controller on Mac OS X 10.6 (Snow Leopard) server and earlier does not correctly      #
#     prioritise nodes with higher scores, due to a longstanding bug (although a score of 0 does prevent a   #
#     node being assigned a job).  I strongly recommend updating your Xgrid controller to OS X 10.7 (Lion)   #
#     server.]                                                                                               #
#                                                                                                            #
#     USAGE                                                                                                  #
#     mgrid [-s stdin] [-i indir] [-d jobid] [-e email-address] [-a art-path | -z art-path] [-b batchname]   #
#     [-c arch] [-g osvers] [-v "env_variable=env_value [...]" ] [-q] [-h "node_task1[:node_task2...]"] [-f] #
#     [-r ram_required] [-n name] [-t number_of_tasks [-y]] [-w [-o file] [-x]] cmd [arg1 [...] ['$task']]   #
#     mgrid -l [-m] [-u] [-j number_of_jobs]                                                                 #
#     mgrid -?                                                                                               #
#     (Use 'mgrid' with no arguments to see the manual page)                                                 #
#                                                                                                            #
#     OPTIONS                                                                                                #
#     The following options are equivaent to those given for xgrid -job submit:                              #
#     -s  Use supplied file as standard input (equivalent to -si).  If the specified file does not exist,    #
#         the argumend is used as a command string to be passed as stdin to the specified command (with a    #
#         warning).  If no command is specified, the input is directed to /bin/bash                          #
#     -i  Input supplied directory (equivalent to -in).  If -w option specified, the same directory is used  #
#         as an output directory.                                                                            #
#     -d  Wait for the dependant job ID to finish (equivalent to -dids, with the limitation that only 1      #
#         dependant job can be specified)                                                                    #
#     -e  Use supplied email address to report status changes (equivalent to -email).  You can also set a    #
#         XGRID_EMAIL environmental variable to be used as the default email address for all submissions.    #
#         This argument overrides the shell variable if set (or use -e "" to temporarily disable emailing).  #
#     -a  Use supplied file as an additional ART script in conjunction with the default inbuild ranking      #
#         script (somewhat equivalent to -art).  The final art score returned to xgrid is the product of the #
#         inbuilt script and the user specified script.  This may be useful for specifying additional        #
#         dependencies (such as software dependencies) while preserving the scoring system for all nodes     #
#         that do meet the dependency. To override the default inbuilt script with a single specified ART    #
#         script, use the -z option instead.  Only one of -z and -a options can be supplied.                 #
#                                                                                                            #
#     Other options for xgrid -job submit are not supported.  The following options are also available:      #
#                                                                                                            #
#     -z  Use supplied file as a ranking script in place of the inbuilt defualt script. Specifying "none"    #
#         or "" disables the inbuilt script. Only one of -z and -a options can be supplied.                  #
#     -b  Produce a batch file with the specified file name and stop (does not submit to xgrid).             #
#     -c  Ensure that jobs are run only on intel (32 or 64 bit), or ppc machines - arch should be either     #
#         'intel' (32 or 64 bit), 'i386' (32 bit intel), 'x86_64' (64 bit intel), or 'ppc' (32 or 64 bit).   #
#     -g  Minimum Mac OS (major release) version for host nodes (format: "10.<major.vers>" eg "10.6").       #
#     -v  A space seperated string of environmental variable(s) in the form of "var_one=value_one.           #
#         var_two=value_two" to be set locally before running the command.                                   #
#     -q  Wait for server nodes to become available rather than running jobs on desktop nodes.  Note:  the   #
#         ranking script will look for an 'mgridserver.sh' executable on the remote host, and if present in  #
#         any standard path will run this executable and convert the output to an integer (ie either 0 or    #
#         1) to determine if the node should be regarded as a server. If no such executable is found then    #
#         the node will be assumed to be a server if it is running OS X server.                              #
#     -h  Use schedule hinting to allocate tasks to the specified node name(s).  For multiple tasks,         #
#         separate node names with a colon (no space) - the number of hints specified must match the number  #
#         of tasks unless only one hint is provided (which is copied to all tasks), or the -f option is      #
#         specified (see below), in which case any matching hostname is allowed to run any task.  The entire #
#         string must be quoted if any node name contains a space.  Colons in node names will produce either #
#         an error or unexpected results.  If the node name provided does not match any of the available     #
#         nodes then the first available node will be chosen (unless the -f option is also specified in      #
#         which case the job will hang).                                                                     #
#     -f  Force the controller to assign the job to the one of the node name(s) specified by -h by using an  #
#         ART script.  Note that ALL tasks will be run on the same node if only 1 node name is provided.     #
#     -r  The minimum amount of free RAM in MB that is required for the job.  If not supplied, 10MB is used  #
#         as a default amount to deter machines that have free processors but no free RAM from accepting     #
#         jobs by assigning them a base score of 1 (applies when using inbuilt default ranking script only). #
#     -n  The name to give the job (appears on xgrid admin etc).  If none is supplied, the current user and  #
#         the command supplied are combined to generate a jobname.                                           #
#     -t  The number of tasks being run.  The arguments should include one containing '$task' which denotes  #
#         the task number (this MUST be enclosed in single quotes or '\' used to escape '$'), otherwise the  #
#         task number will be passed as the last argument to the command (see below for exception).          #
#         -y  Suppress passing the task number as the last argument to the executable (this means that there #
#             is no way to differentiate tasks from inside the executable, but is useful when the command    #
#             doesn't have an argument that can be passed as the task number).  Ignored if no value supplied #
#             for -t option.                                                                                 #
#     -w  After submitting, wait synchronously for the job to complete and then retrieve the job output.  If #
#         an input directory was specified, the same directory will be used as an output directory.  Also    #
#         enables the following options:                                                                     #
#     	  -o  Write job output (stdout and stderr) to specified file (implies -w)                            #
#         -x  Delete the job from xgrid after retrieving results (implies -w)                                #
#         -p  Display the node ranking scores obtained from each node (implies -w)                           #
#                                                                                                            #
#     -l  Display a list of jobs currently on xgrid and exit.  The following arguments can also be given:    #
#         -m  Include a more detailed status of each of the jobs.                                            #
#         -u  Include the username and hostame of whoever submitted the job (provided the jobs were          #
#             submitted using mgrid).                                                                        #
#         -j  The number of jobs to retrieve information for.  By default, information on the last up to 10  #
#             jobs to be submitted is shown, specifier a higher number to see more jobs or 0 to see the      #
#             status of all jobs on that controller.                                                         #
#         All other arguments are ignored.                                                                   #
#                                                                                                            #
#     -?  Print a help/usage message and exit.  All other arguments are ignored.                             #
#                                                                                                            #
#     ARGUMENTS                                                                                              #
#     'cmd' is the command to be run on xgrid, and the remaining arguments are passed as arguments to this   #
#     command.  The special argument '$task' is used to denote the task number, and is appended to (any)     #
#     other commands if ntasks is specified (even if it is 1) and '$task' is not found among the other       #
#     arguments.  The '$task' variable can be embedded in other text to form an argument that changes with   #
#     the task number, for example mgrid -t 2 /usr/bin/cal -y '200$task' would print 2001 for task 1 and     #
#     2002 for task 2.  *NB* If using the '$task' special variable in this way, ensure that the argument is  #
#     enclosed in single quotes (NOT double quotes), or use backslash to escape it (as in "\$"), as the      #
#     shell will otherwise evaluate "$task" to "" on passing the argument to mgrid.  If cmd is not supplied  #
#     but a stdin file is supplied, then the cmd is assumed to be "/bin/bash".  If cmd is in the current     #
#     working directory and no input directory is supplied, then the file specified by cmd is copied as part #
#     of the xgrid job and "./" is prepended to cmd if necessary (to prevent this behaviour prepend the      #
#     specified command with "//" or ".//").                                                                 #
#                                                                                                            #
#     EXAMPLES                                                                                               #
#     Simple xgrid job submission:   mgrid /usr/bin/cal                                                      #
#     Check the status of the last 25 submitted jobs with detailed output:  mgrid -lmu -j 25                 #
#     Run a job synchronously:   mgrid -wx /usr/bin/cal                                                      #
#     Run a series of bash commands:   mgrid -s "sleep 2; echo hello world"                                  #
#     Force a job to run on one of two named nodes, and show the ART scores:                                 #
#          mgrid -p -f -h "node1:node2" /usr/bin/cal                                                         #
#     Run a job with 5 tasks and see which nodes it is assigned to (and all the profile scores):             #
#          mgrid -pxy -t 5 -s '/usr/sbin/system_profiler SPSoftwareDataType | /usr/bin/grep "Computer Name"' #
#     Run a job that executes a local shell script (script is automatically copied):                         #
#          printf '#!/bin/bash \n echo hello world \n' > cmd.sh                                              #
#          mgrid ./cmd.sh                                                                                    #
#                                                                                                            #
#     NOTES                                                                                                  #
#     This script is distributed 'as is', both FREELY and WITHOUT CHARGE, under the GNU general public       #
#     license (see http://www.gnu.org/copyleft/gpl.html).  I am therefore not liable for any damage to your  #
#     computer, xgrid cluster, or sanity caused by using it.                                                 #
#                                                                                                            #
#     If you find this script useful, or find any bugs, then feel free to email me at matthewdenwood@mac.com #
#     Paypal donations to the same address are also gratefully received and may encourage further            #
#     development of the software....                                                                        #
#                                                                                                            #
##############################################################################################################

#echo "Need to update interactive helps"

subid='mgrid_v4.01'


##### Some exit codes from /usr/include/sysexits.h
EX_OK=0
EX_USAGE=64
EX_UNAVAILABLE=69
EX_SOFTWARE=70
EX_CANTCREATE=73
#####


##### Ensure OS is later than Tiger (required for -job wait amongst other things):
osvers=$[`sysctl -n kern.osrelease | awk -F '.' '{print $1}'`-4]
if [ $osvers -lt 5 ]; then
	echo "Mac OS X 10.5 (Leopard) or later is required to use this script" >&2
	exit $EX_UNAVAILABLE
fi
#####


##### Arguments and options setup etc:

queue=0
tasks=0
thejobname=""
manhint=0
manhints=""
force=0
artrank=0
artscore=0
list=0
listcomment=0
listuser=0
printhelp=0
printerror=0
rankingdisabled=0
requiredram=10
pastjobs=10
email=$XGRID_EMAIL
waitafter=0
outfile=""
deleteafter=0
profile=0
notasks=0
minos="0"

# Undocumented variables:
customartname=$mgridcustomartname
if [ "$customartname" == "" ]; then
	customartname="customscore"
fi

# Unused arguments:  k
# Reserve k for kerberos?
        
while getopts "s:i:d:e:a:z:c:g:v:b:qh:fr:n:t:ywo:xplmuj:?" flag; do
if [ $flag == "s" ]; then
	stdin=$OPTARG
elif [ $flag == "i" ]; then
	indir=$OPTARG
elif [ $flag == "d" ]; then
	# Can only have 1 dependant job using this
	depjobs=$OPTARG
elif [ $flag == "e" ]; then
	email=$OPTARG
elif [ $flag == "a" ]; then
	artpath=$OPTARG
	artscore=1
	if [ $artrank == 1 ]; then
		echo "Input error:  cannot specify both -a and -z arguments simultaneously" >&2
		exit $EX_USAGE	
	fi
	if [ ! -f $artpath ]; then
		echo "Input error:  file at path specified to -a doesn't exist" >&2
		exit $EX_USAGE			
	fi	
elif [ $flag == "z" ]; then
	artrank=1
	artpath=$OPTARG
	if [ "$artpath" == "" -o "$artpath" == "none" ]; then
		artrank=0
		artpath=""
	fi
	rankingdisabled=1
	if [ $artscore == 1 ]; then
		echo "Input error:  cannot specify both -a and -z arguments simultaneously" >&2
		exit $EX_USAGE	
	fi
	if [ ! -f $artpath ]; then
		echo "Input error:  file at path specified to -z doesn't exist" >&2
		exit $EX_USAGE			
	fi
elif [ $flag == "c" ]; then
	arch=$OPTARG
	if [ "$arch" == "intel" -o "$arch" == "Intel" -o "$arch" == "INTEL" ]; then
		arch='intel'
	elif [ "$arch" == "intel32" -o "$arch" == "Intel32" -o "$arch" == "INTEL32" -o "$arch" == "i386" -o "$arch" == "I386" ]; then
		arch='i386'
	elif [ "$arch" == "intel64" -o "$arch" == "Intel64" -o "$arch" == "INTEL64" -o "$arch" == "x86_64" -o "$arch" == "X86_64" ]; then
		arch='x86_64'
	elif [ "$arch" == 'ppc' -o "$arch" == 'PPC' -o "$arch" == 'Ppc' ]; then
		arch='ppc'
	else
		echo "Unsupported cpu type '"$arch"'.  Specify one of 'intel', 'i386', 'x64_86' or 'ppc'." >&2
		exit $EX_USAGE
	fi
elif [ $flag == "g" ]; then
	minos=$OPTARG
	if [ "`echo $minos | grep '\.'`" == "" ]; then minos=`echo 10.$minos`; fi	
	minos=`echo $minos | awk -F . '{print $2}'`		
elif [ $flag == "v" ]; then
	env_var=$OPTARG
elif [ $flag == "b" ]; then
	batch=$OPTARG
elif [ $flag == "q" ]; then
	queue=1
elif [ $flag == "h" ]; then
	manhint=1
	manhints=$OPTARG
elif [ $flag == "f" ]; then
	force=1
elif [ $flag == "r" ]; then
	requiredram=$OPTARG
elif [ $flag == "n" ]; then
	thejobname=$OPTARG
elif [ $flag == "t" ]; then
	tasks=$OPTARG
elif [ $flag == "y" ]; then
	notasks=1
elif [ $flag == "w" ]; then
	waitafter=1
elif [ $flag == "o" ]; then
	waitafter=1
	outfile=$OPTARG
elif [ $flag == "x" ]; then
	waitafter=1
	deleteafter=1	
elif [ $flag == "p" ]; then
	waitafter=1
	profile=1
elif [ $flag == "l" ]; then
	list=1
elif [ $flag == "m" ]; then
	listcomment=1
elif [ $flag == "u" ]; then
	listuser=1
elif [ $flag == "j" ]; then
	pastjobs=$OPTARG
elif [ $flag == "?" ]; then
	printerror=1
fi
done


shift $(( $OPTIND-1 ))

nenv_var=`echo $env_var | awk -v p=2 {'print NF'}`
if [ $nenv_var -gt 0 ];then
	for (( i=1; i<=$nenv_var; i++ )); do
		env_vars[$i]=`echo $env_var | awk -v j=$i {'print $j'}`
	done
fi


# For checking environmental variables:
#for (( i=1; i<=$nenv_var; i++ )); do
#	echo ${env_vars[$i]} | awk -F "=" {'print $1'}
#	echo ${env_vars[$i]} | awk -F "=" {'print $2'}
#done

if [ $printerror == 1 ]; then
	#     mgrid [-s stdin] [-i indir] [-d jobid] [-e email-address] [-a art-path | -z art-path] [-b batchname]   #
	#     [-c arch] [-g osvers] [-v "env_variable=env_value [...]" ] [-q] [-h "node_task1[:node_task2...]"] [-f] #
	#     [-r ram_required] [-n name] [-t number_of_tasks [-y]] [-w [-o file] [-x]] cmd [arg1 [...] ['$task']]   #
	#     mgrid -l [-m] [-u] [-j number_of_jobs]                                                                 #
	#     mgrid -?                                                                                               #
	#     (Use 'mgrid' with no arguments to see the manual page)                                                 #
	
	printf "\nmgrid -- version 4.01, January 2013\nby Matthew Denwood (matthewdenwood@mac.com)\n\nUsage:\nmgrid [-s stdin] [-i indir] [-d jobid] [-e email-address] \n      [-a art-path | -z art-path] [-b batchname] [-c arch] \n      [-g osvers] [-v \"env_variable=env_value [...]\" ] [-q] \n      [-h \"node_task1[:node_task2...]\"] [-f] [-r ram_required_(MB)] \n      [-n name] [-t number_of_tasks [-y]] [-w [-o file] [-x]] \n      cmd [arg1 [...]] ['\$task']]\nmgrid -l [-m] [-u] [-j number_of_jobs] \nmgrid -?\n(Use 'mgrid' with no arguments to see the manual page)\n\n" >&2
	exit $EX_USAGE	
fi

if [ $rankingdisabled == 1 -a $profile == 1 ]; then
	echo "Error:  nodescores job requested but node ranking script disabled!" >&2
	exit $EX_USAGE
fi

if [ $# == 0 -a $profile == 0 -a $list == 0 -a "$stdin" == "" ]; then
	printhelp=1
fi

if [ $printhelp == 1 ]; then
	# Horrible hack of a function to allow overstriking (bold) for less (tput doesn't work...):
	embold ()
	{
	str=`echo $1`
	out=""
	for (( i=0; i<${#str}; i++ )); do
	out="$out${str:$i:1}\b${str:$i:1}"
	done

	echo "$out"
	return 0
	}


	printf "
`embold "mgrid -- version 4.01, January 2013"`
by Matthew Denwood (matthewdenwood@mac.com)

`embold SYNOPSIS`
A bash script that provides a replacement for `embold 'xgrid -job submit'` using batch
file submission with support for multi-task jobs, scheduler hinting,
environmental variables, and built-in ART scripts.

`embold USAGE`
`embold mgrid` [`embold -s` stdin] [`embold -i` indir] [`embold -d` jobid] [-\b-e\be email-address] 
   [`embold -a` art-path | `embold -z` art-path] [`embold -b` batchname] [`embold -c` arch]
   [`embold -g` osvers] [`embold -v` \"env_variable=env_value [...]\" ] [`embold -q`] 
   [`embold -h` \"node_task1[:node_task2...]\"] [`embold -f`] [`embold -r` ram_required] 
   [`embold -n` name] [`embold -t` number_of_tasks [`embold -y`]] [`embold -w` [`embold -o` file] [`embold -x`]] 
   `embold cmd` [arg1 [...] ['\$task']]   
`embold 'mgrid -l'` [`embold -m`] [`embold -u`] [`embold -j` number_of_jobs]
`embold 'mgrid -?'`

`embold OPTIONS`
The following options are equivaent to those given for xgrid -job submit:
`embold -s`  Use supplied file as standard input (equivalent to `embold -si`).  If the 
    specified file does not exist, the argumend is used as a command string
    to be passed as stdin to the specified command (with a warning).  If no 
    command is specified, the input is directed to /bin/bash 
`embold -i`  Input supplied directory (equivalent to `embold -in`).  If the `embold -w` option is 
    specified, the same directory is used as an output directory.
`embold -d`  Wait for the dependant job ID to finish (equivalent to `embold -dids`, with the 
    limitation that only 1 dependant job can be specified)
-\b-e\be  Use supplied email address to report status changes (equivalent to 
    `embold -email`). You can also set an XGRID_EMAIL environmental variable to be 
    used as the default email address for all submissions.  This argument 
    overrides the shell variable if set (or use -\b-e\be `embold \\\"\\\"` to temporarily disable 
    emailing).
`embold -a`  Use supplied file as an additional ART script in conjunction with the 
    default inbuild ranking script (somewhat equivalent to `embold -art`).  The 
    final art score returned to xgrid is the product of the inbuilt script and 
    the user specified script.  This may be useful for specifying additional 
    dependencies (such as software dependencies) while preserving the scoring 
    system for all nodes that do meet the dependency. To override the default 
    inbuilt script with a single specified ART script, use the `embold -z` option 
    instead. Only one of `embold -z` and `embold -a` options can be supplied. 

Other options for `embold xgrid` `embold -job` `embold submit` are not supported.  The following options 
are also available:
`embold -z`  Use supplied file as a ranking script in place of the inbuilt defualt 
    script. Specifying \"none\" or \"\" disables the inbuilt script. Only one of 
    `embold -z` and `embold -a` options can be supplied. 
`embold -b`  Produce a batch file with the specified file name and stop (does not submit 
    to xgrid). 
`embold -c`  Ensure that jobs are run only on intel (32 or 64 bit), or ppc machines - 
    arch should be either 'intel' (32 or 64 bit), 'i386' (32 bit intel), 
    'x86_64' (64 bit intel), or 'ppc' (32 or 64 bit). 
`embold -g`  Minimum Mac OS (major release) version for host nodes (format: 
    \"10.<major.vers>\" eg \"10.6\"). 
`embold -v`  A space seperated string of environmental variable(s) in the form of 
    \"var_one=value_one. var_two=value_two\" to be set locally before running 
    the command. 
`embold -q`  Wait for server nodes to become available rather than running jobs on 
    desktop nodes.  Note:  the ranking script will look for an 'mgridserver.sh' 
    executable on the remote host, and if present in any standard path will 
    run this executable and convert the output to an integer (ie either 0 or 1) 
    to determine if the node should be regarded as a server. If no such 
    executable is found then the node will be assumed to be a server if it is 
    running OS X server. 
`embold -h`  Use schedule hinting to allocate tasks to the specified node name(s).  
    For multiple tasks, separate node names with a colon (no space) - the 
    number of hints specified must match the number of tasks unless only one 
    hint is provided (which is copied to all tasks), or the `embold -f` option is 
    specified (see below), in which case any matching hostname is allowed to 
    run any task.  The entire string must be quoted if any node name contains 
    a space.  Colons in node names will produce either an error or unexpected 
    results.  If the node name provided does not match any of the available 
    nodes then the first available node will be chosen (unless the `embold -f` option 
    is also specified in which case the job will hang). 
`embold -f`  Force the controller to assign the job to the one of the node name(s) 
    specified by `embold -h` by using an ART script.  Note that ALL tasks will be run 
    on the same node if only 1 node name is provided. 
`embold -r`  The minimum amount of free RAM in MB that is required for the job.  If 
    not supplied, 10MB is used as a default amount to deter machines that 
    have free processors but no free RAM from accepting jobs by assigning 
    them a base score of 1 (applies when using inbuilt default ranking 
    script only). 
-\b-n\bn  The name to give the job (appears on xgrid admin etc).  If none is 
    supplied, the current user and the command supplied are combined to 
    generate a jobname. 
`embold -t`  The number of tasks being run.  The arguments should include one 
    containing '\$task' which denotes the task number (this MUST be enclosed 
    in single quotes or '\' used to escape '$'), otherwise the task number 
    will be passed as the last argument to the command (but see below for an
    exception, and more details in the arguments section). 
    `embold -y`  Suppress passing the task number as the last argument to the  
        executable (this means that there is no way to differentiate tasks from 
        inside the executable, but is useful when the command doesn't have an  
        argument that can be passed as the task number).  Ignored if no value  
        is supplied for `embold -t` option. 
`embold -w`  After submitting, wait synchronously for the job to complete and then 
    retrieve the job output.  If an input directory was specified, the same 
    directory will be used as an output directory.  Also enables the 
    following options:
    `embold -o`  Write job output (stdout and stderr) to specified file 
        (implies `embold -w`) 
    `embold -x`  Delete the job from xgrid after retrieving results 
        (implies `embold -w`) 
    `embold -p`  Display the node ranking scores obtained from each node 
        (implies `embold -w`) 

`embold -l`  Display a list of jobs currently on xgrid and exit.  The following
    arguments can also be given: 
    `embold -m`  Include a more detailed status of each of the jobs. 
    `embold -u`  Include the username and hostame of whoever submitted the job 
        (provided the jobs were submitted using mgrid). 
    `embold -j`  The number of jobs to retrieve information for.  By default, 
        information on the last up to 10 jobs to be submitted is shown, 
        specifier a higher number to see more jobs or 0 to see the 
        status of all jobs on that controller. 
    All other arguments are ignored.  
	
`embold -?`  Print a help/usage message and exit.  All other arguments are ignored.  

`embold ARGUMENTS`
'cmd' is the command to be run on xgrid, and the remaining arguments are passed
as arguments to this command.  The special argument '\$task' is used to denote 
the task number, and is appended to (any) other commands if ntasks is specified 
(even if it is 1) and '\$task' is not found among the other arguments.  The 
'\$task' variable can be embedded in other text to form an argument that changes 
with the task number, for example mgrid -t 2 /usr/bin/cal -y '200\$task' would 
print 2001 for task 1 and 2002 for task 2.  *NB* If using the '\$task' special 
variable in this way, ensure that the argument is enclosed in single quotes 
(NOT double quotes), or use backslash to escape it (as in \"\$\"), as the shell 
will otherwise evaluate \"\$task\" to \"\" on passing the argument to mgrid.  If 
cmd is not supplied but a stdin file is supplied, then the cmd is assumed to be 
\"/bin/bash\".  If cmd is in the current working directory and no input directory 
is supplied, then the file specified by cmd is copied as part of the xgrid job 
and \"./\" is prepended to cmd if necessary (to prevent this behaviour prepend 
the specified command with \"//\" or \".//\").  

`embold EXAMPLES`
Simple xgrid job submission:   
    `embold "mgrid /usr/bin/cal"`
Check the status of the last 25 submitted jobs with detailed output:  
    `embold "mgrid -lmu -j 25"`
Run a job synchronously:   
    `embold "mgrid -wx /usr/bin/cal"`
Run a series of bash commands:   
    `embold "mgrid -s \\\"sleep 2; echo hello world\\\""`
Force a job to run on one of two named nodes, and show the ART scores: 
    `embold "mgrid -p -f -h \\\"node1:node2\\\" /usr/bin/cal"`
Run a job with 5 tasks and see which nodes it is assigned to (and scores): 
    `embold "mgrid -pxy -t 5 -s '/usr/sbin/system_profiler SPSoftwareDataType | "`
	`embold "/usr/bin/grep \\\"Computer Name\\\"'"`

`embold NOTES`
This script is distributed 'as is', both FREELY and WITHOUT CHARGE, under the 
GNU general public license (see http://www.gnu.org/copyleft/gpl.html).  I am 
therefore not liable for any damage to your computer, xgrid cluster, or sanity 
caused by using it.  

If you find this script useful, or find any bugs, then feel free to email me at 
matthewdenwood@mac.com Paypal donations to the same address are also gratefully 
received and may encourage further development of the software....

" |  fold -s | less

	exit $EX_OK
fi


if [ $artscore == 1 -a $artrank == 1 ];then
	echo "Cannot specify both -a and -r options" >&2
	exit $EX_USAGE
fi


# I use quite a few (probably more than necessary) temporary files, so they should at least be secure:
tmpdir=`mktemp -d -t temp` && success=1 || success=0

if [ $success == 0 ]; then
	echo "Unable to create temporary working directory" >&2
	exit $EX_CANTCREATE
fi

# Assure the file is removed at program termination or after we received a signal:
trap 'rm -rf "$tmpdir" >/dev/null 2>&1' 0
trap "exit 2" 1 2 3 15

# Set up a temporary file I use:
tfile1=${tmpdir}/temp1

# Check that hostname and password are set as environmental variables:
if [ "$XGRID_CONTROLLER_HOSTNAME" == "" ]; then
	echo "Error:  The environmental variables XGRID_CONTROLLER_HOSTNAME (and, if required, XGRID_CONTROLLER_PASSWORD) are not set.  Use the export command to do this (probably in your .profile).  For a local job run, use the hostname ':private:'" | fold -s >&2
	exit $EX_USAGE
fi


# Check to see that the hostname and password are correct:
expect -c "
#		exp_internal 1
		spawn xgrid -job list
		expect {
			Password:  { send \r\n; interact }
			eof { exit }
		}
		exit
		" > $tfile1

if [ "`cat $tfile1 | grep 'Password:'`" == "" ]; then
	passblank=0
else
	passblank=1
fi
if [ "`cat $tfile1 | grep 'Unable to connect'`" == "" ]; then
	hostwrong=0
else
	hostwrong=1
fi
if [ "`cat $tfile1 | grep 'Authentication failed'`" == "" ]; then
	passwrong=0
else
	passwrong=1
fi
rm $tfile1
	
if [ $hostwrong == 1 ]; then
	echo "Error:  The Xgrid controller '$XGRID_CONTROLLER_HOSTNAME' could not be found.  Make sure the environmental variables XGRID_CONTROLLER_HOSTNAME (and, if required, XGRID_CONTROLLER_PASSWORD) are set to the match the details for the controller you are trying to connect to."  | fold -s >&2
	exit $EX_USAGE
fi
if [ $passwrong == 1 ]; then
	if [ $passblank == 1 ]; then
		echo "Error:  The blank Xgrid controller password specified is incorrect.  Use the export command to set the environmental variable XGRID_CONTROLLER_PASSWORD to the match the details for the controller you are trying to connect to."  | fold -s >&2
	else
		echo "Error:  The Xgrid controller password '$XGRID_CONTROLLER_PASSWORD' is incorrect.  Use the export command to set the environmental variable XGRID_CONTROLLER_PASSWORD to the match the details for the controller you are trying to connect to."  | fold -s >&2
	fi
	exit $EX_USAGE
fi
if [ $passblank == 1 -a $passwrong == 0 ]; then
	echo "Warning:  The password on the controller is set as blank and Xgrid is continually prompting for a password, which is likely to cause problems with using this script.  If you encounter problems you could try updating to the latest version of OS X, or just set a password on your controller."  | fold -s >&2
fi

cmd=$1

if [ "$cmd" == "" -a "$stdin" != "" ]; then
	cmd="/bin/bash"
fi

if [ "$cmd" == "" -a $profile == 1 ]; then
	outfile="/dev/null"
	deleteafter=1
	cmd="/usr/bin/cal"
	thejobname="nodescores"
fi

# If we don't have an outfile use a temporary file so we can re-read the results after deleting job:
catoutfile=0;
if [ "$outfile" == "" ]; then
	catoutfile=1;
	outfile="${tmpdir}/outfile"
	touch $outfile
fi

shift
nargs=$#

if [ $nargs -ge 1 ]; then
	for (( i=1; i<=$nargs; i++ )); do
		args[$i]=$1
		shift
	done
fi

if [ $tasks == 0 ]; then
	tasks=1
else
	# If ntasks < 0 (ie -t specified) then look for '$task' in any arguemt; if it isn't found then append an argument which is just $task unless -y argument is also supplied
	dolltaskf=0
	for (( i=1; i<= $nargs; i++ )); do
		task=1
		grepr=`echo ${args[$i]} | grep '$task'`
		if [ ! "$grepr" == "" ]; then
			dolltaskf=1
		fi
	done
	if [ $dolltaskf == 0 -a $notasks == 0 ]; then
		nargs=$[$nargs+1]
		args[$nargs]='$task'
	fi
fi

for (( i=1; i<=$tasks; i++ )); do
	hints[$i]=""
done

if [ $manhint == 1 ]; then
	nmanhints=`echo $manhints | awk -v p=2 'BEGIN { FS = ":" } ; { print NF }'`
	if [ "$manhints" == "" ]; then
		nmanhints=1
	fi
	if [ $nmanhints == 1 ]; then
	
		for (( i=1; i<=$tasks; i++ )); do
			hints[$i]=`echo $manhints | awk 'BEGIN { FS = ":" } ; { print $1 }'`
		done
	
	else
		
		if [ $nmanhints != $tasks -a $force == 0 ]; then
			echo "The number of manual schedule hints provided does not match the number of tasks - ensure the node names are separated by a colon"  | fold -s >&2
			exit $EX_USAGE
		fi
	
		for (( i=1; i<=$nmanhints; i++ )); do
			hints[$i]=`echo $manhints | awk -v j=$i 'BEGIN { FS = ":" } ; { print $j }'`
		done
	
	fi
fi

if [ $force == 1 -a $manhint == 0 ]; then
	echo "You must supply a node name using the -h argument when using the -f argument">&2
	exit $EX_USAGE
fi
if [ $force == 1 -a "$manhints" == "" ]; then
	echo "You must supply a non-blank node name to the -h argument when using the -f argument">&2
	exit $EX_USAGE
fi

if [ ! "$indir"=="" ]; then
	if [ ! -d "$indir" ]; then
		echo "The directory specified as the input directory does not exist">&2
		exit $EX_USAGE
	fi
fi
if [ ! "$artpath"=="" ]; then
	if [ ! -f "$artpath" ]; then
		echo "The file specified as art-path does not exist">&2
		exit $EX_USAGE
	fi
fi

if [ "$thejobname" == "" ]; then
	thejobname=`whoami`-$cmd
fi


##### End rguments and options setup etc




##################   This section is for the -l options only

if [ $list == 1 ]; then
	printf "Retrieving current job list from xgrid..."
	xgrid -job list >& $tfile1 && success=1 || success=0
	if [ $success == 1 ]; then
	
		ORIGIFS=$IFS
		#save the current value of ifs

		IFS=$(echo -en "\n\b\r")
		#reset ifs to end of line stuff

		exec 3<&0
		#save current value of stdin

		exec 0<$tfile1
		#set stdin to read from the temporary file
		
		jobs=0
		string=""
		while read line
		do
			linenospace=`echo $line | awk '{print $1}' | sed -e "s/,//g"`
			if [ ! "$linenospace" == "{" -a ! "$linenospace" == "}" -a ! "$linenospace" == ");" -a ! "$linenospace" == "jobList" ]; then
				string=`echo $string $linenospace`
				jobs=$(( $jobs+1 ))
				pid[$jobs]=$linenospace
			fi
		done
		
		exec 0<&3 3<&-
		# restore stdin and free up fd#6 for other processes to use

		IFS=$ORIGIFS
		# restore $IFS which was used to determine what the field separators are

		rm $tfile1
		tempstatusfile=${tmpdir}/tempstatusfile
		
		if [ $jobs == 0 ]; then
			printf "\rThere are no jobs currently listed on xgrid\n"
		else
			
			if [ $listuser == 1 ]; then
				if [ $listcomment == 1 ]; then
					echo "Name;ID;Status;Comment;User;Email" > $tempstatusfile
				else
					echo "Name;ID;Status;User;Email" > $tempstatusfile
				fi
			else
				if [ $listcomment == 1 ]; then
					echo "Name;ID;Status;Comment" > $tempstatusfile
				else
					echo "Name;ID;Status" > $tempstatusfile
				fi
			fi
			
			startat=1
			if [ $pastjobs != 0 ]; then
				if [ $[$jobs-$pastjobs] -gt 0 ]; then
					startat=$[$jobs-$pastjobs+1]
				fi
			fi
			
			for (( i=$startat; i<=$jobs; i++ )); do
				
				xgrid -job attributes -id ${pid[$i]} >& $tfile1 && success=1 || success=0
				if [ $success == 1 ]; then
					
					grepstring=`cat $tfile1 | grep "jobStatus"`
					status=`echo $grepstring | awk '{print $3}'`
					
					grepstring=`cat $tfile1 | grep "applicationIdentifier"`
					usern=`echo $grepstring | awk '{print $3}' | sed 's/;//g' | sed 's/"//g' | sed "s/'//g"`
					
					usere=`echo $usern | awk 'BEGIN { FS = "#@#" } ; { print $2 }'`
					if [ "$usere" == "" ]; then
						usere='Not avaialble'
					fi
					if [ "$usere" == "no_email" ]; then
						usere='None supplied'
					fi
					
					usern=`echo $usern | awk 'BEGIN { FS = "#@#" } ; { print $1 }'`
					
					if [ "$usern" == "com.apple.xgrid.cli" ]; then
						usern="Unknown"
					fi
					if [ "$usern" == "igrid_MD" ]; then
						usern="Unknown"
					fi
					if [ "$usern" == "mgrid_MD" ]; then
						usern="Unknown"
					fi
					
					grepstring=`cat $tfile1 | grep "percentDone"`
					percent=`echo $grepstring | awk '{print $3}' | sed 's/;//g' | sed 's/"//g' | awk '{print int(($1*10)+0.5)/10}'`
					
					grepstring=`cat $tfile1 | grep "dateSubmitted"`
					subon=`echo $grepstring | sed 's/;//g' | sed 's/"//g' | awk '{print $3 " at " $4}'`

					grepstring=`cat $tfile1 | grep "dateStarted"`
					startedon=`echo $grepstring | sed 's/;//g' | sed 's/"//g' | awk '{print $3 " at " $4}'`

					grepstring=`cat $tfile1 | grep "dateStopped"`
					stoppedon=`echo $grepstring | sed 's/;//g' | sed 's/"//g' | awk '{print $3 " at " $4}'`
					
					grepstring=`cat $tfile1 | grep "name"`
					jobname=`echo $grepstring | awk '{print $3}' | sed -e "s/[;\"]//g"`
					
					string="$jobname;${pid[$i]}"
					
					if [ "$status" == "Finished;" ]; then
						if [ $listcomment == 1 ]; then
							string=`echo "$string;Complete;Finised on: $stoppedon"`
						else
							string=`echo "$string;Complete"`
						fi
					elif [ "$status" == "Pending;" ]; then
						grepstring=`cat $tfile1 | grep "suspended"`

						if [ "`echo $grepstring`" == "" ]; then
							if [ $listcomment == 1 ]; then
								string=`echo "$string;Pending;Submitted on: $subon"`
							else
								string=`echo "$string;Pending;"`
							fi
						else
							if [ $listcomment == 1 ]; then
								string=`echo "$string;PAUSED;Submitted on: $subon"`
							else
								string=`echo "$string;PAUSED"`
							fi
						fi
					elif [ "$status" == "Canceled;" ]; then
						if [ $listcomment == 1 ]; then
							string=`echo "$string;Canceled;Submitted on: $subon"`
						else
							string=`echo "$string;Canceled"`
						fi
					elif [ "$status" == "Running;" ]; then
						grepstring=`cat $tfile1 | grep "activeCPUPower"`
						cpupower=`echo $grepstring | awk '{print $3}' | sed -e "s/;//g"`
						
						grepstring=`cat $tfile1 | grep "suspended"`

						if [ "`echo $grepstring`" == "" ]; then
							if [ $listcomment == 1 ]; then
								string=`echo "$string;RUNNING;${percent:0:4}% complete, started on: $startedon"`
							else
								string=`echo "$string;RUNNING"`
							fi
						else
							if [ $listcomment == 1 ]; then
								string=`echo "$string;PAUSED;${percent:0:4}% complete, started on: $startedon"`
							else
								string=`echo "$string;PAUSED"`
							fi
						fi
						
						#`", active CPU power: $cpupower"`
					elif [ "$status" == "Failed;" ]; then
						grepstring=`cat $tfile1 | grep "error"`
						error=`echo $grepstring | awk ' {
						for (i=3; i<=NF; i++)
						printf("%s ", $i)
						} ' | sed -e "s/;//g"`
						if [ $listcomment == 1 ]; then
							string=`echo "$string;FAILED;Returned error $error"`
						else
							string=`echo "$string;FAILED"`
						fi
					else
						if [ $listcomment == 1 ]; then
							string=`echo "$string;UNKNOWN;Xgrid returned an unrecognised status type for this job"`
						else
							string=`echo "$string;UNKNOWN"`
						fi
					fi
					
					if [ $listuser == 1 ]; then
						string=`echo "$string;$usern;$usere"`
					fi
					echo $string >> $tempstatusfile
				else
					echo "There was an error retrieving job information for process id ${pid[$i]}" >> $tempstatusfile
				fi
				rm $tfile1
			done
			#echo ""
			if [ $jobs -le $pastjobs -o $pastjobs == 0 ]; then				
				printf "\rThe following $jobs jobs are currently listed on xgrid:\n"
			else
				printf "\rThere are $jobs jobs currently listed on xgrid; these are the $pastjobs most recent:\n"
			fi
			#echo $string # echo ${pid[@]} would do the same thing anyway
			
			cat < $tempstatusfile | column -t -s ";" 
			rm $tempstatusfile
			#echo ""
			#echo "Use xgrid.results along with the job name to retrive results for a completed job"
		fi
		
		exit $EX_OK
		
	else
		echo "An error occured while contacting xgrid.  The following was returned:" >&2
		cat < $tfile1 >&2
		rm $tfile1
		exit $EX_UNAVAILABLE
	fi
	
	exit $EX_OK
fi

##################   End list option section



echo " "
echo "Generating xgrid job"
	
if [ ! -f "$stdin" -a "$stdin" != "" ]; then
	echo "[Specified stdin file doesn't exist; interpreting argument as a command string]"
	stdinfile=${tmpdir}/stdinfile
	echo $stdin > $stdinfile
	stdin=$stdinfile
fi


# Determine if cmd is a text file (not currently using):
istext=`file "$cmd" | grep "text"`
#echo $istext
# Determine if cmd is in the working directory:
incwd=`echo "$cmd" | sed -e "s/\.\///g" | grep "/"`
if [ "$incwd" == "" ]; then
	incwd=1
else
	incwd=0
fi
copycmd=0

# If cmd is in the working directory (and a text file?) then automagically include it in the batch file:
if [ $incwd == 1 -a "$indir" == "" ]; then
	copycmd=1
	dotslash=`echo "$cmd" | grep "\.\/"`
	if [ "$dotslash" == "" ]; then
		cmd=`echo ./$cmd`
	fi
fi

### Put together the scoring ART script (combination of default scoring script and -a custom script, or -z specified script with default scoring script as wrapper using no calcualtions):

#  Look for a (shell or) ruby scoring script in the local directory, then the xgrid directory.  It must be run as a function/script taking arguments xjobs, physical cores, logical cores, server, ram, cpu speed, bits (64 or 32), (ART) cpu type, (ART) server, (ART) os vers, (ART) custom and echoing the score.  If none is there then use an internal script (half of the previous one).  The profiling script returns all this information.  The scores are then worked out by bash before every task is assigned and the top score given a task until all tasks done.

if [ $artrank == 1 ]; then

	cat > ${tmpdir}/scoring.rb <<-EOA
	#!/usr/bin/ruby -w

	# arguments xjobs, physical cores, logical cores, server, ram, cpu speed, bits (64 or 32), (ART) cpu type, (ART) server, (ART) os vers, (ART) custom

	puts ARGV[10].to_i
	exit 0
	EOA

else
	
	if [ -f scoring.rb ]; then
		if [ $uog == 1 -a $manhint == 0 ]; then
			echo "Using local ruby file for grid scoring"
		fi
		cp scoring.rb ${tmpdir}/scoring.rb
	elif [ -f scoring.sh ]; then
		if [ $uog == 1 -a $manhint == 0 ]; then
			echo "Using local shell script for grid scoring"
		fi
		cp scoring.sh ${tmpdir}/scoring.rb
	elif [ -f "/Library/Application Support/mgrid/scoring.rb" -a $manhint == 0 -a ! "`whoami`" == "nobody" ]; then
		cp "/Library/Application Support/mgrid/scoring.rb" ${tmpdir}/scoring.rb
	elif [ -f "/Library/Application Support/mgrid/scoring.sh" -a $manhint == 0 -a ! "`whoami`" == "nobody" ]; then
		cp "/Library/Application Support/mgrid/scoring.sh" ${tmpdir}/scoring.rb
	else
	
		# The ` and $ all need to be escape charactered:
	
		cat > ${tmpdir}/scoring.rb <<-EOA
		#!/usr/bin/ruby -w

		# arguments xjobs, physical cores, logical cores, server, ram, cpu speed, bits (64 or 32), (ART) cpu type, (ART) server,  (ART) os vers, (ART) custom

		# Current jobs:
		xjobs = ARGV[0].to_i
		# Physical cores:
		ncpu = ARGV[1].to_i
		# Logical cores:
		nlogical = ARGV[2].to_i
		# 1 or 0:
		server = ARGV[3].to_i
		# Free RAM in MB:
		ram = ARGV[4].to_i
		# Clock speed in MHz:
		cpuspeed = ARGV[5].to_i
		# bits (64 or 32):
		bits = ARGV[6].to_i
		# Result of CPU type script (1 or 0; or 1 if not used):
		artcpu = ARGV[7].to_i
		# Result of server/desktop script (1 or 0; or 1 if not used):
		artserver = ARGV[8].to_i		
		# Result of os version script (1 or 0; or 1 if not used):
		artos = ARGV[9].to_i		
		# Result of the custom art score (or 1 if not used):
		artcustom = ARGV[10].to_i

		serverbonus = 6

		# Start with a base score of 2 in case no free processors calculated here is a lie!
		base = 2
		if (xjobs < nlogical)
		  # If there's any space left on the processor, check whether there
		  # are any physical cores left (ignoring one desktop core for the user!)
		  if (xjobs < ncpu + server - 1)
		    # If there are then rank by number of cores, giving a 4 bonus to
		    # servers, so that 4 physical cores will be used on each server
		    # before any desktops are touched
		    base = base + (ncpu - xjobs + serverbonus * server) * 1000
		  else
		    # If not, then rank by number of logical cores, again with bonus
		    # for servers to push jobs to them
		    base = base + (nlogical - xjobs + serverbonus * server) * 100
		  end
		  # Randomise to vary order of usage
		  base = base + Kernel.rand(50)
		end
		
		# Set score to 1 (or 0?) if the minimum free RAM requirement is not met (the dollarrequiredram thing is automatically swapped out for the right value when writing the script):
		if (ram < $requiredram)
			base = 1
		end
		
		# Take into account the ART scores (server/desktop * intel/ppc * os * custom.art - these will just be 1 if not used):
		base = base * artcustom * artserver * artcpu * artos
	
		puts base.to_i
	
		EOA
	
	fi
fi

# scoring.rb needs to be executable so I can use it locally:
chmod 755 ${tmpdir}/scoring.rb
scoring=${tmpdir}/scoring.rb

# The getinfo file is not meant to be user modified:

cat > ${tmpdir}/getinfo.rb <<-EOA
#!/usr/bin/ruby -w

# To get all sysctl info:
# sysctl -a

# Can extract number of physical cpu cores from sysctl
ncpu = \`/usr/sbin/sysctl hw.physicalcpu | /usr/bin/sed 's/hw.physicalcpu: //'\`.strip.to_i

# And the number of logical cores
nlogical = \`/usr/sbin/sysctl hw.ncpu | /usr/bin/sed 's/hw.ncpu: //'\`.strip.to_i
os=\`/usr/sbin/system_profiler SPSoftwareDataType | /usr/bin/grep 'System Version'\`.strip


# OLD approximation using tasks directory - problem is crashed jobs stay there....
#xjobs = \`/bin/ls -1 /var/xgrid/agent/tasks | /usr/bin/wc -l\`.strip.to_i - 1
#xjobs = \`/usr/bin/find /var/xgrid/agent/tasks -newerat '1 month ago' -and -not -newerBt '10 seconds ago' -maxdepth 1 -type d -name '??????*' | /usr/bin/wc -l\`.strip.to_i
#xjobs = 0 if (xjobs < 0)

# NEW approximation using iostat - much better and multi-task aware
#user = \`/usr/sbin/iostat -dCI -n 1 -w 5 -c 3 | /usr/bin/sed -e 's/us//g' -e 's/sy//g' -e 's/id//g' | /usr/bin/awk 'BEGIN {FS=" "}
#min=="" {
#min=max=\$4
#}
#{
#if (\$4 > max) {max = \$4};
#if (\$4 < min) {min = \$4};
#total += \$4
#count += 1
#}
#END {
#print min;
#}'\`.to_f

# Which corresponds to how many processors (allowing for a small amount of extraneous junk but otherwise rounded to the whole processor up)?
#xjobs = ((user*nlogical+95)/100).to_i

# NEWER approximation using vm.loadavg - much faster than iostat
user = \`/usr/sbin/sysctl vm.loadavg | /usr/bin/perl -wane 'print \$F[-3]'\`.to_f
xjobs = (user+0.5).to_i

# Is it a server?
server = 0
if (os =~ /Server/)
  server = 1
end
# Has the XGRID_SERVER override variable been set?
shellserver = \`mgridserver\`.strip
if (shellserver!="")
	server = shellserver.to_i
	if (shellserver=="true")
		server = 1
	end
	if (shellserver=="TRUE")
		server = 1
	end	
	if (shellserver=="yes")
		server = 1
	end
	if (shellserver=="YES")
		server = 1
	end	
end

# FREE MEMORY IN PAGES:
ram = \`vm_stat | awk '{ print \$3 }' | awk -F '\\\\n' 'BEGIN { RS = "" } { print \$2 }'\`.strip.to_i
# 1 page = 4096 bytes, so FREE RAM IN MB:
ram = (ram*4096) / (1024**2)

# TOTAL RAM (in MB):
totalram = \`/usr/sbin/sysctl hw.physmem | /usr/bin/sed 's/hw.physmem: //'\`.strip.to_i
totalram = totalram / (1024**2)

# CPU speed (in MHz):
cpuspeed = \`/usr/sbin/sysctl hw.cpufrequency | /usr/bin/sed 's/hw.cpufrequency: //'\`.strip.to_i
cpuspeed = cpuspeed / (1000**2)

# bits (64 or 32):
bits = \`/usr/sbin/sysctl hw.optional 2> /dev/null |/usr/bin/awk -F': ' '/64/ {print $2}'\`.strip.to_i
bits = (bits*32) + 32

# This swaps the spaces in names for _ and leaves an _ at the end of the name, but I deal with this in the shell script:
comname=\`/usr/sbin/system_profiler SPSoftwareDataType | /usr/bin/grep 'Computer Name' | /usr/bin/awk '{
for (i=3; i<=NF; i++)
printf("%s ", \$i)
}' | /usr/bin/sed -e 's/ /_/g' -e "s/'//g" -e "s/’//g" -e 's/"//g'\`.strip

artcustom=1
artserver=1
artcpu=1
artos=1
# If the serveronly.sh and/or cpu type and/or os and/or art.score.sh files exist then execute them and get the result:
if File.exist?('customart.sh')
  artcustom = \`./customart.sh\`.strip.to_i
end
if File.exist?('server.rb')
  artserver = \`./server.rb\`.strip.to_i
end
if File.exist?('arch.sh')
  artcpu = \`./arch.sh\`.strip.to_i
end
if File.exist?('minos.sh')
  artos = \`./minos.sh\`.strip.to_i
end

EOA

# Copy getinfo to prepare it for profiling job and ART scoring job:
cp ${tmpdir}/getinfo.rb ${tmpdir}/infoscoring.rb
# Append stuff to profiling job:
echo 'puts "\"#{comname}\" #{xjobs} #{ncpu} #{nlogical} #{server} #{ram} #{cpuspeed} #{bits} #{artcpu} #{artserver} #{artos} #{artcustom}"' >> ${tmpdir}/getinfo.rb
# Hexdump in (possibly local) scoring script and put in code to create the file locally and then run it:
printf '`echo "' >> ${tmpdir}/infoscoring.rb
xxd $scoring >> ${tmpdir}/infoscoring.rb
echo '" | xxd -r - scoring.rb`' >> ${tmpdir}/infoscoring.rb
echo '`chmod 755 scoring.rb`' >> ${tmpdir}/infoscoring.rb
echo 'puts `./scoring.rb #{xjobs} #{ncpu} #{nlogical} #{server} #{ram} #{cpuspeed} #{bits} #{artcpu} #{artserver} #{artos} #{artcustom}`.strip.to_i' >> ${tmpdir}/infoscoring.rb
echo '`rm scoring.rb`' >> ${tmpdir}/infoscoring.rb


if [ $queue == 1 ]; then
	
	cat > ${tmpdir}/server.rb <<-EOA
	#!/usr/bin/ruby -w

	os=\`/usr/sbin/system_profiler SPSoftwareDataType | /usr/bin/grep 'System Version'\`.strip

	# Is it a server?
	server = 0
	if (os =~ /Server/)
	  server = 1
	end
	# Has the XGRID_SERVER override variable been set?
	shellserver = \`mgridserver\`.strip
	if (shellserver!="")
		server = shellserver.to_i
		if (shellserver=="true")
			server = 1
		end
		if (shellserver=="TRUE")
			server = 1
		end	
		if (shellserver=="yes")
			server = 1
		end
		if (shellserver=="YES")
			server = 1
		end	
	end
	
	puts server

	EOA
fi

if [ "$arch" != "" ]; then
	
	cat > ${tmpdir}/arch.sh <<-EOA
	#!/bin/bash
	
	if [ "$arch" == "intel" ]; then
		if [ \`arch\` == "i386" ]; then
			echo 1
		else
			echo 0
		fi
	elif [ "$arch" == "ppc" ]; then
		if [ \`arch\` == "i386" ]; then
			echo 0
		else
			echo 1
		fi
	elif [ "$arch" == "i386" ]; then
		if [ \`arch\` == "i386" ]; then
			bits=\`uname -a | grep "i386"\`
			if [ "\$bits" == "" ]; then
				echo 0
			else
				echo 1
			fi			
		else
			echo 0
		fi
	elif [ "$arch" == "x86_64" ]; then
		if [ \`arch\` == "i386" ]; then
			bits=\`uname -a | grep "x86_64"\`
			if [ "\$bits" == "" ]; then
				echo 0
			else
				echo 1
			fi			
		else
			echo 0
		fi
	fi
	
	exit 0

	EOA
fi

if [ "$minos" != "0" ]; then
	
	cat > ${tmpdir}/minos.sh <<-EOA
	#!/bin/bash
	
	osvers=\`/usr/sbin/system_profiler SPSoftwareDataType | /usr/bin/grep 'System Version' | sed "s/Server //" | /usr/bin/awk '{print \$6}' | /usr/bin/awk -F . '{print \$2}'\`
	if [ \$osvers -lt $minos ]; then
		echo 0
	else
		echo 1
	fi
	exit 0

	EOA
fi

######  Profiling scripts are now ready to be used






###### Produce batch file for actual job:

batchname=${tmpdir}/batch

echo "{" > $batchname 
echo "jobSpecification = {" >> $batchname
echo "submissionIdentifier = \"$subid\";" >> $batchname
if [ "$email" == "" ]; then
	echo "applicationIdentifier = \"`whoami`"@"$HOSTNAME#@#no_email\";" >> $batchname
else
	echo "applicationIdentifier = \"`whoami`"@"$HOSTNAME#@#$email\";" >> $batchname
fi
echo "artConditions = {" >> $batchname
if [ $rankingdisabled == 0  ]; then
	echo "\"basescore\" = {" >> $batchname
	echo "artMin = 1;" >> $batchname
	echo "};" >> $batchname
fi
if [ ! "$artpath" == "" ]; then
	echo "\"$customartname\" = {" >> $batchname
	echo "artMin = 1;" >> $batchname
	echo "};" >> $batchname
fi
if [ $queue == 1 ]; then
	echo "\"serverscore\" = {" >> $batchname
	echo "artMin = 1;" >> $batchname
	echo "};" >> $batchname
fi
if [ "$arch" != "" ]; then
	echo "\"archscore\" = {" >> $batchname
	echo "artMin = 1;" >> $batchname
	echo "};" >> $batchname
fi
if [ "$minos" != "0" ]; then
	echo "\"osscore\" = {" >> $batchname
	echo "artMin = 1;" >> $batchname
	echo "};" >> $batchname
fi
if [ $force == 1 ]; then
	echo "\"forcescore\" = {" >> $batchname
	echo "artMin = 1;" >> $batchname
	echo "};" >> $batchname
fi
			
echo "};" >> $batchname
echo "artSpecifications = {" >> $batchname
if [ $artrank == 0 -a $rankingdisabled == 0 ]; then
	echo "\"basescore\" = {" >> $batchname
	echo "artData = <" >> $batchname
	hexdump  -v -e ' "" 4/1 "%02x" " "'  ${tmpdir}/infoscoring.rb >> $batchname
	echo ">;" >> $batchname
	echo "isExecutable = YES;" >> $batchname
	echo "};" >> $batchname
fi
if [ ! "$artpath" == "" ]; then
	echo "\"$customartname\" = {" >> $batchname
	echo "artData = <" >> $batchname
	hexdump  -v -e ' "" 4/1 "%02x" " "'  "$artpath" >> $batchname
	echo ">;" >> $batchname
	echo "isExecutable = YES;" >> $batchname
	echo "};" >> $batchname
fi
if [ $queue == 1 ]; then
	echo "\"serverscore\" = {" >> $batchname
	echo "artData = <" >> $batchname
	hexdump  -v -e ' "" 4/1 "%02x" " "'  ${tmpdir}/server.rb >> $batchname
	echo ">;" >> $batchname
	echo "isExecutable = YES;" >> $batchname
	echo "};" >> $batchname
fi
if [ "$arch" != "" ]; then
	echo "\"archscore\" = {" >> $batchname
	echo "artData = <" >> $batchname
	hexdump  -v -e ' "" 4/1 "%02x" " "'  ${tmpdir}/arch.sh >> $batchname
	echo ">;" >> $batchname
	echo "isExecutable = YES;" >> $batchname
	echo "};" >> $batchname
fi
if [ "$minos" != "0" ]; then
	echo "\"osscore\" = {" >> $batchname
	echo "artData = <" >> $batchname
	hexdump  -v -e ' "" 4/1 "%02x" " "'  ${tmpdir}/minos.sh >> $batchname
	echo ">;" >> $batchname
	echo "isExecutable = YES;" >> $batchname
	echo "};" >> $batchname
fi

if [ $force == 1 ]; then
	echo "#!/bin/bash" > ${tmpdir}/force.sh
	echo "comname=\`/usr/sbin/system_profiler SPSoftwareDataType | /usr/bin/grep 'Computer Name' | /usr/bin/awk '{ for (i=3; i<=NF; i++) printf(\"%s \", \$i) }'\`" >> ${tmpdir}/force.sh
	echo 'if [ "$comname" == "'${hints[1]}'" -o "$comname" == "'${hints[1]}' " ]; then' >> ${tmpdir}/force.sh
	echo 'echo 1' >> ${tmpdir}/force.sh
	if [ $nmanhints -gt 1 ]; then	
		for (( i=2; i<=$nmanhints; i++ )); do
			echo 'elif [ "$comname" == "'${hints[$i]}'" -o "$comname" == "'${hints[$i]}' " ]; then' >> ${tmpdir}/force.sh
			echo 'echo 1' >> ${tmpdir}/force.sh
		done
	fi	
	echo 'else' >> ${tmpdir}/force.sh
	echo 'echo 0' >> ${tmpdir}/force.sh
	echo 'fi' >> ${tmpdir}/force.sh
	echo "exit 0" >> ${tmpdir}/force.sh
	
	echo "\"forcescore\" = {" >> $batchname
	echo "artData = <" >> $batchname
	hexdump  -v -e ' "" 4/1 "%02x" " "'  ${tmpdir}/force.sh >> $batchname
	echo ">;" >> $batchname
	echo "isExecutable = YES;" >> $batchname
	echo "};" >> $batchname
fi
echo "};" >> $batchname

echo "name = \"$thejobname\";" >> $batchname
if [ "$email" != "" ];then
	echo "notificationEmail=\"$email\";" >> $batchname 
fi

if [ "$depjobs" != "" ]; then
	echo "\"schedulerParameters\" = {" >> $batchname
	# tasksMustStartSimultaneously and minimumTaskCount also belong in here
	echo "\"dependsOnJobs\" = (" >> $batchname
	echo "$depjobs" >> $batchname
	echo ");" >> $batchname
	echo "};" >> $batchname
fi


echo "inputFiles = {" >> $batchname

if [ $copycmd == 1 ]; then
	
	echo "\"$cmd\" = {" >> $batchname
	echo "fileData = <" >> $batchname
	hexdump  -v -e ' "" 4/1 "%02x" " "'  $cmd >> $batchname
	echo ">;" >> $batchname
	echo "isExecutable = YES;" >> $batchname
	echo "};" >> $batchname
	
	
	cmd="./$cmd"
fi

if [ ! -d "$indir" ]; then
	if [ "$indir" != "" ]; then
		if [ $waitafter == 1 ]; then
			echo "[Specified input/output directoy doesn't exist; creating a new directory]"	
			mkdir $indir
		else
			echo "Error:  specified input directory doesn't exist" >&2
			exit $EX_USAGE		
		fi
	fi
fi

if [ "$indir" != "" ]; then
	cwd=`pwd`
	cd "$indir"
	
	temp2=${tmpdir}/temp2

	find . ! -name '.*' > $temp2

	# could also use find `ls`  ! -name '.*' > tempfile.$$ here to create a list of files/dirs without the ./ prefix
	# then I would need to change every instance of ${line:2} to just $line in the loop below

	ORIGIFS=$IFS
	#save the current value of ifs

	IFS=$(echo -en "\n\b\r")
	#reset ifs to end of line stuff

	exec 3<&0
	#save current value of stdin

	exec 0< $temp2
	#set stdin to read from the temporary file

	while read line
	do
		if [ -d ${line:2} ]; then
			if [ ! "$(ls ${line:2})" ]; then
				##  EMPTY EXCEPT FOR INVISIBLE FILES
				echo "\"${line:2}/holdingfile.$$.hex\" = {" >> $batchname
				echo "fileData = <686f6c64 0a>;" >> $batchname
				echo "isExecutable = NO;" >> $batchname
				echo "};" >> $batchname
			fi
		else
			if [ ${line:2} != $batchname ]; then
				echo "\"${line:2}\" = {" >> $batchname
				# ignore the './' at the start of the filename
				echo "fileData = <" >> $batchname
				hexdump  -v -e ' "" 4/1 "%02x" " "'  $line >> $batchname
				echo ">;" >> $batchname
				if [ -x ${line:2} ]; then
					echo "isExecutable = YES;" >> $batchname
				else
					echo "isExecutable = NO;" >> $batchname
				fi
				echo "};" >> $batchname
			fi
		fi
	done

	exec 0<&3 3<&-
	# restore stdin and free up fd#6 for other processes to use
 
	IFS=$ORIGIFS
	# restore $IFS which was used to determine what the field separators are

	rm $temp2
	
	if [ "$stdin" != "" ]; then
		if [ ! -f $stdin ]; then
			echo "\"$stdin\" = {" >> $batchname
			echo "fileData = <" >> $batchname
			hexdump  -v -e ' "" 4/1 "%02x" " "'  $stdin >> $batchname
			echo ">;" >> $batchname
			if [ -x $stdin ]; then
				echo "isExecutable = YES;" >> $batchname
			else
				echo "isExecutable = NO;" >> $batchname
			fi
			echo "};" >> $batchname
		fi
	fi
	cd "$cwd"
else
	if [ "$stdin" != "" ]; then
		echo "\"$stdin\" = {" >> $batchname
		echo "fileData = <" >> $batchname
		hexdump  -v -e ' "" 4/1 "%02x" " "'  $stdin >> $batchname
		echo ">;" >> $batchname
		if [ -x $stdin ]; then
			echo "isExecutable = YES;" >> $batchname
		else
			echo "isExecutable = NO;" >> $batchname
		fi
		echo "};" >> $batchname
	fi
fi
echo "};" >> $batchname


if [ $force == 0 -a $manhint == 1 ]; then
	echo "schedulerHints = {" >> $batchname

	for (( i = 1; i <= $tasks; i++ ))
	do
		if [ ! "${hints[$i]}" == "" ]; then
			echo $i"=\""${hints[$i]}"\";" >> $batchname
		fi
	done
	echo "};" >> $batchname
fi

echo "taskSpecifications = {" >> $batchname

for (( i = 1; i <= $tasks; i++ ))
do
	echo "$i = {" >> $batchname
	#if [ $i != 1 ];then
	#	echo "dependsOnTasks="$[$i-1]";" >> $batchname
	#fi
	echo "command = \"$cmd\";" >> $batchname
	
#	if [ -f ./$cmd ]; then 
#		cmdisfile=1
#	else
#		cmdisfile=0
#	fi
	
	if [ $nenv_var -gt 0 ]; then
		echo "environment = {" >> $batchname
		for (( j=1; j<=$nenv_var; j++ )); do
			var=`echo ${env_vars[$j]} | awk -F "=" {'print $1'}`
			val=`echo ${env_vars[$j]} | awk -F "=" {'print $2'}`
			echo $var" = "$val";" >> $batchname
		done
		echo "};" >> $batchname
	fi
		
	echo "arguments = (" >> $batchname
	if [ $nargs -ge 1 ]; then
		for (( j=1; j<=$nargs; j++ )); do
			grepf=`echo ${args[$j]} | grep '$task'`
			if [ "$grepf" == "" ]; then
				echo -n "\"${args[$j]}\"" >> $batchname
			else
				task=$i
				echo -n "\""`eval echo ${args[$j]}`"\"" >> $batchname
			fi
			if [ $j -lt $nargs ]; then
				echo "," >> $batchname
			fi
		done
	fi
	echo ");" >> $batchname
	if [ "$stdin" != "" ]; then
		echo "inputStream = \"$stdin\";" >> $batchname
	fi
	
	echo "};" >> $batchname
done

echo "};" >> $batchname
echo "};" >> $batchname
echo "}" >> $batchname

if [ "`plutil -s $batchname`" != "" ]; then
	echo "An error caused the job submission to be aborted, sorry.  plutil is complaining about the batch file I created - please send the 'batch-failed.txt' file to me at matthewdenwood@mac.com"  | fold -s >&2
	cp $batchname 'batch-failed.txt'
	exit $EX_SOFTWARE
fi
	
if [ "$batch" != "" ]; then
	echo "Creating batch file '$batch' for your job in the local directory..."
	mv $batchname $batch
else
	echo "Submitting your job to xgrid..."
	
	plutil -convert binary1 $batchname
	xgrid -job batch $batchname >& $tfile1 && success=1 || success=0
	if [ $success == 1 ]; then
		
		jobstring=`tail +2 $tfile1 | head -n 1`
		cutstr=${jobstring:20}
		len=$(( ${#cutstr}-1 ))
		jobnum=${cutstr:0:len}
		
		# This is just for runjags support:
		if [ -f starter.sh -o -f scriptlauncher.sh ]; then
			echo $jobnum > jobid.txt
			echo "Your job ID is "$jobnum
		else
			echo "Job submission successful, your job ID is "$jobnum
		fi
		
		rm $tfile1
		
	else
		echo "An error occured while submitting the job to xgrid.  The following was returned:"  | fold -s >&2
		cat < $tfile1 >&2
		rm $tfile1
		exit $EX_UNAVAILABLE
	fi
	
fi



if [ "$batch" == "" -a $waitafter == 1 ]; then
	
	# Set up a temporary file I use:
	tfile1=${tmpdir}/temp1
	tfile2=${tmpdir}/temp2
	tempscores=${tmpdir}/tempscores
	prettyprofile=${tmpdir}/pprof
	touch $tfile1
	touch $tfile2
	touch $tempscores
	touch $prettyprofile
	
	echo "Waiting for xgrid job to finish..."
		
	finished=0
		
	while [ $finished == 0 ]; do
		xgrid -job wait -id $jobnum >& $tfile1 && success=1 || success=0	
		if [ $success == 0 ]; then
			echo "An error occured while contacting xgrid.  The following was returned:"  | fold -s >&2
			cat < $tfile1 >&2
			rm $tfile1
			rm $tfile2
			rm $tempscores
			rm $prettyprofile
			exit $EX_UNAVAILABLE
		fi
			
		xgrid -job attributes -id $jobnum >& $tfile1 && success=1 || success=0
		if [ $success == 1 ]; then
				
			errorout=`cat $tfile1 | grep "error = InvalidJobIdentifier"`
			if [ "$errorout" != "" ]; then
				echo "Error:  The job has been deleted" >&2
				exit $EX_USAGE
			fi
				
			grepstring=`cat $tfile1 | grep "jobStatus"`
			status=`echo $grepstring | awk '{print $3}'`
	
			if [ "$status" == "Finished;" ]; then
				finished=1
			else
				if [ "$status" == "Running;" ]; then
					echo "Job is now running"
				elif [ "$status" == "Pending;" ]; then
					grepstring=`cat $tfile1 | grep "suspended"`
					if [ "`echo $grepstring`" == "" ]; then
						echo "Job is now pending"
					else
						echo "Job has been suspended - waiting for restart"
					fi
				elif [ "$status" == "Canceled;" ]; then
					echo "Error:  The job has been cancelled" >&2
					exit $EX_USAGE
				elif [ "$status" == "Failed;" ]; then
					echo "Note:  The job status is showing as failed; results may not be complete (job has not been deleted)" | fold -s 
					deleteafter=0
					finished=1
				else
					echo "Note:  Unknown status returned by job, results may not be returned correctly (job has not been deleted)" | fold -s 
					deleteafter=0
					finished=1
				fi
			fi
		else
			echo "An error occured while attempting to retrieve information about the job - the following was returned:"  | fold -s >&2
			cat $tfile1 >&2
			exit $EX_SOFTWARE
		fi
			
	done	
	
	if [ $profile == 1 ]; then
	
		echo "Retrieving job..."
		xgrid -job log -id $jobnum | grep '\\\" returned score' > $tempscores && success=1 || success=0
	
		if [ $success == 1 ]; then
		
			typesp=${tmpdir}/typesp
			nodesp=${tmpdir}/nodesp		
		
			ORIGIFS=$IFS
			#save the current value of ifs

			IFS=$(echo -en "\n\b\r")
			#reset ifs to end of line stuff

			exec 3<&0
			#save current value of stdin

			exec 0<$tempscores
			#set stdin to read from the temporary file
		
			nentries=0
			while read line
			do
				nentries=$[$nentries+1]
				thescore[$nentries]=`echo $line | awk -F '"' '{print $5}'`
				scoretype[$nentries]=`echo $line | awk -F '"' '{print $3}'`
				thenode[$nentries]=`echo $line | awk -F '"' '{print $7}'`
				echo ${scoretype[$nentries]} >> $typesp
				echo ${thenode[$nentries]} >> $nodesp
				echo "${thenode[$nentries]}|${scoretype[$nentries]}|${thescore[$nentries]}" >> $prettyprofile
			done
		
			# Sort is consistent so I can rely on these to be in the same order:
			cat $prettyprofile | sort > $tfile2		
			cat $typesp | sort | uniq > $prettyprofile
			mv $prettyprofile $typesp
			cat $nodesp | sort | uniq > $prettyprofile
			mv $prettyprofile $nodesp
		
			exec 0<$typesp
			i=0
			while read line
			do
				i=$[$i+1]
				types[$i]=`echo $line`
			done
			ntypes=$i
		
			exec 0<$nodesp
			i=0
			while read line
			do
				i=$[$i+1]
				nodes[$i]=`echo $line`
			done
			nnodes=$i
		
			exec 0<$tfile2
			i=1
			j=0
			while read line
			do
				j=$[$j+1]
				if [ $j == 1 ]; then
					totalscore[$i]=1
				fi
				newscore=`echo $line | awk -F '|' '{print $3}'`
				totalscore[$i]=$[${totalscore[$i]}*$newscore]
				scores[i]=`echo "${scores[$i]}|$newscore|"`
				
				if [ $j == $ntypes ]; then
					echo ${totalscore[$i]}" "$i >> $prettyprofile
					i=$[$i+1]
					j=0
				fi
			done		
				
			# Calculate the original order for feedback at the end:
			cat $prettyprofile | sort -nr -o $tempscores
		
			exec 0<$tempscores
			rank=1
			while read line
			do
				order[$rank]=`echo $line | awk '{print $2}'`
				rank=$[$rank+1]
			done

			exec 0<&3 3<&-
			# restore stdin and free up fd#6 for other processes to use
			IFS=$ORIGIFS
			# restore $IFS which was used to determine what the field separators are
		
		
			basepresent=0
			custompresent=0
			archpresent=0
			serverpresent=0
			forcepresent=0
		
			# Produce prettified output:
			printf "Name" > $prettyprofile
			if [ $ntypes == 1 ]; then
				printf "|Score" >> $prettyprofile
			else
				printf "|TotalScore" >> $prettyprofile
				for (( j=1; j<=$ntypes; j++ )); do
					thetype=${types[$j]}
					if [ "$thetype" == "basescore" ]; then
						thetype="BaseScore"
					fi
					if [ "$thetype" == "archscore" ]; then
						thetype="CPU(arch)"
					fi
					if [ "$thetype" == "serverscore" ]; then
						thetype="Server"
					fi
					if [ "$thetype" == "osscore" ]; then
						thetype="OSverScore"
					fi
					if [ "$thetype" == "customscore" ]; then
						thetype="CustomART"
					fi					
					printf "|$thetype" >> $prettyprofile
				done
			fi
			printf "\n" >> $prettyprofile
			for (( i=1; i<=$nnodes; i++ )); do
				current=${order[$i]}
				printf "${nodes[$current]}|${totalscore[$current]}" >> $prettyprofile
				if [ $ntypes -gt 1 ]; then
					printf "${scores[$current]}" >> $prettyprofile
				fi
				printf "\n" >> $prettyprofile								
			done
		
		else
			echo "An error occured while contacting xgrid.  The following was returned:"  | fold -s >&2
			cat < $tfile1 >&2
		fi
	
		echo "The following ranking scores were obtained:"
		echo ""
		cat $prettyprofile | column -t -s "|"
		echo ""
	fi	

	# outfile is either /dev/null, a file specified as an argument, or a temporary file to dump results to (in which case catoutfile=1):
	if [ "$outfile" != "/dev/null" ]; then
		echo "Retrieving results"
		echo ""
	
		if [ "$indir" == "" ]; then
			xgrid -job results -id $jobnum >& $outfile && success=1 || success=0
		else
			xgrid -job results -id $jobnum -out $indir >& $outfile && success=1 || success=0
		fi
		if [ $success == 0 ]; then
#			if [ "$outfile" == "" ]; then
#				echo ""
#			fi		
			echo "An error occured while retrieving the job results from xgrid.  The job has not been deleted."  | fold -s >&2
			deleteafter=0
		fi
	fi	
	
	if [ $deleteafter == 1 ]; then
#		if [ "$outfile" == "" ]; then
#			echo ""
#		fi
		echo "Deleting job from Xgrid"
		xgrid -job delete -id $jobnum >& $tfile1 && success=1 || success=0
	
		if [ $success == 0 ]; then
			echo "An error occured while deleting the profiling job from xgrid.  The following was returned:"  | fold -s >&2
			cat < $tfile1 >&2
			rm $tfile1
			rm $tfile2
			rm $tempscores
			rm $prettyprofile
			exit $EX_UNAVAILABLE
		fi
	fi
	
	if [ $catoutfile == 1  -a "$outfile" != "/dev/null" ]; then
		# Checks to see if outfile is empty:
		if [ -s $outfile ]; then 
			echo ""			
			echo "Job output is displayed below:"
			echo ""			
			cat $outfile
			echo ""			
		else
			echo ""
		fi
	else
		echo ""
		echo "Finished"
		echo ""
	fi

	rm $tfile1
	rm $tfile2
	rm $tempscores
	rm $prettyprofile
	
else

	echo "Finished"
	echo ""
fi

exit $EX_OK
