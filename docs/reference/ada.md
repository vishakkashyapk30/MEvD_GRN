== Overview ==
Ada cluster consists of ninety two Boston SYS-7048GR-TR nodes equipped with dual Intel Xeon E5-2640 v4 processors, providing 40 virtual cores per node, 128 GB of 2400MT/s DDR4 ECC RAM and four Nvidia GeForce GTX 1080 Ti GPUs, providing  14336 CUDA cores, and 44 GB of GDDR5X VRAM or four Nvidia GeForce RTX 2080 Ti GPUs providing 17408 cores, and 44 GB of GDDR6 VRAM.  The nodes are connected to each other via a Gigabit Ethernet network. All compute nodes have a 1.8 TB local scratch and a 960 GB local SSD scratch. The compute nodes are running Ubuntu 18.04 LTS. [https://slurm.schedmd.com/ SLURM] software is used as a job scheduler and resource manager. The aggregate theoretical peak performance of Ada is 70.66 TFLOPS (CPU) + 4588 TFLOPS (FP32 GPU).

==Applying for an account==
An Ada account is available to IIIT faculty, research staff, and research students. To apply for a new account, please send an email to [mailto:hpc.admin@iiit.ac.in hpc.admin@iiit.ac.in] with the following information. 

*Name
*Roll Number / Employee ID
*Research Center
*Faculty Advisor
*Preferred Login ID

For a CVIT associated account, fill out [https://forms.gle/L4XdmyHvYQidQged9 this form]. The form requires that you sign-in using your Gmail account. If you do not get a response about your account creation within 72 hrs of filling this form, please mail to [mailto:cvit-sudo@googlegroups.com cvit-sudo@googlegroups.com] with the same details. If you are a CVIT student, please watch this [https://www.youtube.com/watch?v=U3_pPJgs2Fg tutorial video]. Strict actions would be taken if requested for a CVIT account without being affiliated to any of the CVIT faculty.

== Accessing Ada ==

===Logging in===
You can log in to Ada by SSH from IIIT LAN. For accessing Ada from off-campus, [https://vpn.iiit.ac.in/ IIIT VPN] must be used.
<source lang="bash">
$ ssh -X user_name@ada.iiit.ac.in
</source>

===Changing initial password===
When you are prompted to enter the password, enter the initial password. In  '''(current) UNIX password:''' , enter the initial password. In '''New password:''' and '''Retype new password:''', enter a new password of  your choice.

[[File:Password2.png]]

===Moving files to and from Ada===

To move a directory from local machine to Ada:
<source lang="bash">
$ scp -r local_directory user_name@ada.iiit.ac.in:
(or)
$ rsync -avz local_directory user_name@ada.iiit.ac.in:
</source>

To move a directory from Ada to local machine:
<source lang="bash">
$ scp -r user_name@ada.iiit.ac.in:remote_directory .
(or)
$ rsync -avz user_name@ada.iiit.ac.in:remote_directory .
</source>

Both '''scp''' and '''rsync''' transfers files locally or over network. If the transfer is interrupted, '''rsync''' has the ability to continue from where it left off when invoked again.

== Partitions, Account, and QoS ==

A partition can be considered as a collection of nodes. There are two partitions in the cluster. The '''short''' partition has two nodes and is for compiling/debugging codes. The '''long''' partition is for serial/parallel jobs that need to run for longer than 6 hours. Node 01-40 contains 4 GeForce GTX 1080 Ti GPUs each while nodes 43-92 contain 4 GeForce RTX 2080 Ti GPUs each. 

{|class="wikitable"
|-
|Partition
|Nodes
|DefMemPerCPU
|MaxMemPerCPU
|Gres
|Maxtime
|Priority
|-
|short
|gnodes[01-02]
|1024 MB
|3000 MB
|gpu:4
|6:00:00
|100
|-
|long
|gnodes[03-92]
|1024 MB
|3000 MB
|gpu:4
|Infinite
|100
|-
|ihub
|gnode[92-112]
|1024
|5000 MB
|gpu:4
|Infinite
|100
|}

A SLURM account is like a bank account, and all users belong to at least one account. The allocated resources to a job are charged to the job's specified account. All users have access to '''research''' account. The accounts '''cvit''', '''nlp''', '''ccnsb''', '''mll''', '''hai''', '''irel''', '''dma''', '''adsac''', '''rrc''', '''plafnet''', and '''biosona''' are accessible only to users in projects/centres that have contributed hardware to the cluster.

{|class="wikitable"
|-
|Account
|Access
|GrpCPUs
|GrpTRES=gres/gpu
|GrpJobs
|GrpSubmitJobs
|Allowed QoS
|-
|research
|ALL
|1640
|164
|820
|1640
|medium
|-
|cvit
|CVIT 
|1080
|108
|540
|1080
|normal
|-
|mll
|MLL
|40
|4
|20
|40
|normal
|-
|nlp
|NLP 
|120
|12
|60
|120
|normal
|-
|ccnsb
|CCNSB
|80
|8
|40
|80
|normal
|-
|cesp
|CESP
|40
|4
|20
|40
|normal
|-
|hai
|HAI
|40
|4
|20
|40
|normal
|-
|irel
|IREL
|120
|12
|60
|120
|normal
|-
|dma
|DMA
|40
|4
|20
|40
|normal
|-
|adsac
|ADSAC
|40
|4
|20
|40
|normal
|-
|rrc
|RRC
|40
|4
|20
|40
|normal
|-
|plafnet
|PLAfNET
|40
|4
|20
|40
|normal
|-
|biosona
|BioSonA
|240
|24
|120
|240
|normal
|-
|}
It is recommended to specify a Quality of Service (QOS) for each job submitted to SLURM. The default QOS for research account is '''medium''' and it is '''normal''' for cvit, mll, nlp, cesp, ccnsb, hai, irel, dma, adsac, rrc, plafnet and biosona accounts.

{|class="wikitable"
|-
|QoS
|MaxCPUsPerUser
|MaxTRESPerUser
|MaxJobsPerUser
|MaxSubmitJobsPerUser
|MaxWall
|MaxTresPerJob
|Priority
|-
|low
|10
|gres/gpu=1
|1
|4
|4-00:00:00
|gres/gpu=1
|0 
|-
|medium
|40
|gres/gpu=4
|4
|8
|4-00:00:00
|gres/gpu=4
|10
|-
|normal
|Account limits
|Account limits
|Account limits
|Account limits
|Infinite
|gres/gpu=4
|0
|-
|}

The following command can be used to list the allowed Accounts and QoSes.
<source lang="bash">
$ sacctmgr show assoc user=$USER format=Account,QOS,DefaultQOS
</source>

=== cvit account ===
To submit jobs to the CVIT account, users should specify SLURM job directive '''-A $USER'''. 
e.g. <source lang="bash">$ sinteractive -c 2 -g 1 -A $USER</source>
The following command will show the associations for an account. 
<source lang="bash">$ sacctmgr show assoc account=$USER</source>

sinteractive is non-standard and restrictive in terms of options and documentation. We recommend '''sbatch''' and '''srun''' for the same [https://slurm.schedmd.com/sbatch.html] . The following is an example of '''srun''' command where we are requesting 1 GPU, 10 CPUs with 2GB memory per CPU (20GB in total) with the '''cvit''' account in the long partition. Please maintain a 1:10 ratio for the number of gpu : number of cpu. The CVIT admins may kill jobs in case this is not followed without warning.  
<source lang="bash">
$ srun --pty --partition=long -A $USER --gres=gpu:1 --mem-per-cpu=2G -c 10 bash -l
</source>
In case, you want a specific node please use the '''--nodelist gnodeXX''' parameter. However, please note that if the requested node is unavailable due to some reason (other users running jobs in the same node, node is under maintenance etc.) then your job will go into the '''pending''' state.  

Sometimes we reserve nodes and give reservation names for dedicated access (in very special cases like conference deadlines and when users need to handle huge amounts of data). Use a reservation argument to gain access. For example we had used a reservation with reservation-name = '''cvit-trial'''
<source lang="bash">
$ srun --reservation <reservation-name> ...
for e.g. 
$ srun --reservation cvit-trial --pty --partition=long -A $USER --gres=gpu:2 --mem-per-cpu=2G -c 10 bash -l

</source>
'''srun''' jobs (also called bash jobs) has a max time limit of 6 hrs. '''srun''' command will allocate a terminal to the user and allow interactive access. This interactive access is meant to be utilized for the general coding and debugging stage. For jobs requiring more than 6 hrs of runtime, the '''sbatch''' command has to be used. An example of this is given in - 

Sample script: batch_script1.sh   
<source lang="bash">
#!/bin/bash
#SBATCH -A research/username (if you want to use the cvit account put username)
#SBATCH -n 10
#SBATCH --gres=gpu:1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=4-00:00:00
#SBATCH --output=op_file.txt

module load u18/cuda/10.1
module load u18/cudnn/7-cuda-10.1

python .....
</source>
We run this script by,
<source>
$ sbatch batch_script1.sh
</source>
The output from the python file (anything that is present in the stdout buffer) will be stored in the op_file.txt. Please use the argument '''flush=True''' in the print statements in your python code to flush the buffer(otherwise it can lead to irregular updates in the file) regularly. If you want to continually check the updates to the file, please use the following command - 
<source>
$ tail -f op_file.txt
</source>
 
CVIT admins control GPUMinutes and the highest number of GPUs allocated to a user. By default on account creation, 1 GPU and 600 GPUMinutes is allowed to a user. Please mail to [mailto:cvit-sudo@googlegroups.com cvit-sudo@googlegroups.com] for changing these limits as and when required.

=== sub account ===

A new SLURM account has been created on Ada to utilize the idle 
nodes/cores/GPUs. The details are as follows: 

<source>
SLURM account: sub  (#SBATCH -A sub)
QOS:     sub        (#SBATCH --qos=sub)
Wall time: 6 hours  (#SBATCH --time=6:00:00)
MaxCPUs: 40
Max GPUs: 4
</source>

The jobs submitted using this account will have low priority and will 
run only when there are no high priority (QOS: normal,medium) pending 
jobs.
t

== File Systems ==
The cluster provides four types of file storage to users. They are referred here as /home, /share1, /scratch and /ssd_scratch. For CVIT users, there is /share3 for long term storage.

The /home located at /home/$USER is an NFS storage and has a disk quota of 25 GB. This space can be used to store your source code and to build your executables. The data store on /home is backed up every day. The /share1 is RAID6 storage available on the master node and has a quota of 100 GB (CVIT users have a group quota of 6 TB). The 6 TB space for CVIT users is used to store public datasets in /share1/dataset directory. If you are planning to use a public dataset, please search it in the said location; you may find it is already present there. Further, you can contact the CVIT admins to move a public dataset to this directory if it is not present already. Each CVIT user also has access to a separate storage space in /share3. This space can be used for long-term storage and for transferring large data files to and from compute node /scratch. By default, each user gets 50 GB of storage space. This can be increased to 200 GB depending on the circumstances. Please mail to [mailto:cvit-sudo@googlegroups.com cvit-sudo@googlegroups.com] for increasing /share3 space. The /scratch is for storing temporary files created during job run time. The /ssd_scratch is for storing temporary files that required high-speed disk I/O access. The files older than 10 days are purged from /scratch and /ssd_scratch. 

{|class="wikitable"
|-
|Space
|Purpose
|Visibility
|Backup
|Quota
|Total Size
|File Deletion Policy
|-
|/home
|Sofware installation space, storing codes and small files
|Master and compute nodes
|Yes
|25 GB
|9.8 TB
| None
|-
|/share1
|Long-term storage 
|Master node only
|No
|100 GB / 6 TB
|20 TB
|None
|-
|/share2
|Long-term storage 
|Master node only
|No
| 
|13 TB
|None
|-
|/share3
|Long-term storage for only CVIT users 
|Master node only
|No
|50 GB/ 100 GB / 150 GB /200 GB
|22 TB
|None
|-
|-
|/scratch
|Temporary storage for large files
|Local disk attached to each compute node
|No
|None
|2.0 TB
| 7 days<sup>#</sup>
|-
|/ssd_scratch
|Temporary storage for jobs that require fast I/O
|Local disk attached to each compute node
|No
|None
|960 GB
|7 days<sup>#</sup>
|}

<sup>#</sup> File deletion based on creation time (ctime).

==Policies==
===Account Policies===
i) Access to HPC account is open to faculty members, post-doctoral researchers
and research students (MS and PhD). Under certain circumstances, a non-research student may be granted access to an account on the HPC system with the endorsement of their faculty advisor.

ii) To apply for an account, an email with the following details should be sent
to hpc.admin@iiit.ac.in

Name,  Roll number / Employee ID, Research Center, Faculty Advisor and Preferred
login name.

iii) Users are required to change their password every six months. They will receive alerts upon login starting one month before the password expiration date.

iv) Accounts will be removed / deleted when a student has been issued no-dues
certificate by IT office and when an employee leaves the institute. Account that
are not accesseed for six months will be locked.

v) Sharing accounts is strictly forbidden and will result in the account being locked for a period of three days to a week.

vi) The account of a student will be locked upon the request of their faculty advisor, while the account of an employee will be locked upon the request of a competent authority.

===Software Policies===
i) Custom Python and R environments must be installed in users' home directories.

ii)  The majority of commonly used software has been installed system-wide and is
available as modules. If there is a need for a specific software to be installed, a request should be sent to hpc.admin@iiit.ac.in .

iii) Installation of any unlicensed proprietary software on the HPC nodes is strictly prohibited.

===Security Policies===
i) Sharing of accounts, even among users working on the same project, is strictly prohibited.

ii) Users are advised to set up passwordless SSH access to their accounts for convenience and ease of use.

iii) Users should not leave their terminal unattended while logged in.

iv) Any suspicious activities and security problems should be reported to hpc.admin@iiit.ac.in immediately.

===Backup Policies===
i) Users' home directories are backed up once every 24 hours, and teh backups
are retained for up to three months.  Users are strongly advised to perform frequent backups of important data to alternate locations.

ii)  The backups are accessible as read-only on the login node.

==Environment Modules ==

The environment module allows users to set shell environmental variables needed for the software. 

To view the list of current loaded modules:
<source lang="bash">
[parithi@ada ~]$ module list

Currently Loaded Modulefiles:
  1) u18/namd/2.14       2) u18/openmpi/4.1.2

</source>

To list the installed modules:
<source lang="bash">
[parithi@ada ~]$ module avail

-------------------------------------------------------- /opt/Modules/versions ---------------------------------------------------------
3.2.10

--------------------------------------------------- /opt/Modules/3.2.10/modulefiles ----------------------------------------------------
u18/amber/18                                       u18/gromacs/2021.4-plumed2
u18/amber/18-plumed                                u18/gromacs/2021.4-plumed2-avx512
u18/ambertools/20                                  u18/gromacs/2021.5-plumed-avx512-ambertools-mmPBSA
u18/ambertools/21                                  u18/gromacs/2022
u18/aria2/1.35.0                                   u18/gromacs/2022-avx512
u18/cmake/3.23.3                                   u18/leptonica/1.82.0
u18/cuda/10.0                                      u18/matlab/R2022a
u18/cuda/10.2                                      u18/namd/2.12
u18/cuda/11.6                                      u18/namd/2.14
u18/cuda/8.0                                       u18/namd/3.0_alpha10
u18/cuda/9.0                                       u18/openblas/0.3.17
u18/cuda/9.2                                       u18/openmpi/2.1.1-gcc5
u18/cudnn/7.6.5-cuda-10.2                          u18/openmpi/3.1.6
u18/cudnn/7.6.5-cuda-9.2                           u18/openmpi/4.0.1-cuda10
u18/cudnn/7.6-cuda-10.0                            u18/openmpi/4.0.7
u18/cudnn/8.3.2-cuda-11.5                          u18/openmpi/4.1.2
u18/cudnn/8.3.3-cuda-10.2                          u18/openmpi/4.1.2-cuda11
u18/cudnn/8.4.0-cuda-11.6                          u18/plumed2/2.6.1
u18/ffmpeg/5.0.1                                   u18/plumed2/2.8.0
u18/freesurfer/7.2.0                               u18/python/3.10.2
u18/fsl/6.0.5.2                                    u18/python/3.11.2
u18/gcc/8.4.0                                      u18/python/3.7.4
u18/gcc/9.4.0                                      u18/python/3.8.3
u18/gflags/2.2.2                                   u18/R/4.1.3
u18/glog/0.5.0                                     u18/singularity-ce/3.9.6
u18/go/1.17.8                                      u18/TensorRT/8.4.0.6
u18/gromacs/2016.3                                 u18/tesseract/5.2.0
u18/gromacs/2016.3-plumed                          u18/use.own
u18/gromacs/2016.3-v2

</source>

To load a module:

<source lang="bash">
[parithi@ada ~]$ module load u18/openmpi/4.1.2
</source>

To remove/unload a module.

<source lang="bash">
[parithi@ada ~]$ module unload u18/openmpi/4.1.2
</source>

To display the changes made by a module to the user shell environment: 

<source lang="bash">
[parithi@ada ~]$ module disp u18/cuda/11.6 
-------------------------------------------------------------------
/opt/Modules/3.2.10/modulefiles/u18/cuda/11.6:

module-whatis	 adds CUDA-11.6 to your environment variable 
prepend-path	 PATH /usr/local/cuda-11.6/bin 
prepend-path	 LD_LIBRARY_PATH /usr/local/cuda-11.6/lib64 
-------------------------------------------------------------------

</source> 

===List of Available Modules ===
{|class="wikitable"
|-
|Software
|Module
|Remarks
|-
|CUDA 9.0
|cuda/9.0
|
|-
|CUDA 9.1
|cuda/9.1
|
|-
|CUDA 10.0
|cuda/10.0
|
|-
|cuDNN 7
|cudnn/7-cuda-9.0
|CUDA 9.0
|-
|cuDNN 7.6.4
|cudnn/7.6.4-cuda-9.0
|CUDA 9.0
|-
|cuDNN 7.1
|cudnn/7.1-cuda-9.1
|CUDA 9.1
|-
|cuDNN 7
|cudnn/7-cuda-10.0
|CUDA 10.0
|-
|cuDNN 7.3
|cudnn/7.3-cuda-10.0
|CUDA 10.0
|-
|cuDNN 7.6
|cudnn/7.6-cuda-10.0
|CUDA 10.0
|-
|OpenBLAS
|openblas/0.3.6
|Haswell
|-
|OpenCV
|opencv/3.3.0
|Python 2.7
|-
|[[OpenMPI]]
|openmpi/2.1.1
|
|-
|OpenMPI
|openmpi/3.1.0
|
|-
|OpenMPI
|openmpi/4.0.0
|
|-
|OpenMPI
|openmpi/4.0.1-cuda10
|CUDA aware
|-
|[[NAMD 2.13]]
|namd/2.13
|
|-
|PLUMED 2.5.2
|plumed/2.5.2
|
|-
|[[GROMACS 2019]]
|gromacs/2019
|
|-
|GROMACS 2019.3
|gromacs/2019.3-plumed
|plumed patched
|-
|[[AMBER 16]]
|amber/16
|Compiled with GCC-4.7
|-
|gflags
|gflags/2.2.1
|
|-
|gflags
|gflags/2.2.2
|
|-
|glog
|glog/0.3.5
|
|-
|glog
|glog/0.4.0
|
|-
|[[MATLAB]]
|matlab/R2019b
|
|-
|MPICH
|mpich-3.2
|Compiled with GCC-5
|-
|Intel TBB
|tbb/2018u1-release
|Compiled with macro TBB_USE_DEBUG=0
|-
|Intel TBB
|tbb/2018u1-debug
|Compiled with macro TBB_USE_DEBUG=1
|-
|FFMPEG
|ffmpeg/3.4
|cuda 9.0, cuvid, nvenc, nonfree, libnpp 
|-
|FFMPEG
|ffmpeg/4.01
|
|}

==SLURM Commands==

{|class="wikitable"
|-
|Description
|Command
|Example
|-
|Submit a batch job
|[https://slurm.schedmd.com/sbatch.html sbatch]
|<pre>sbatch job_script.sh</pre>
|-
|Submit an interactive job
|[[sinteractive]]
|<pre>sinteractive -c 2 -g 1 </pre>
|-
|Cancel a job
|[https://slurm.schedmd.com/scancel.html scancel]
|<pre>scancel job_id</pre>
|-
|List all current jobs for an user
|[https://slurm.schedmd.com/squeue.html squeue]
|<pre>squeue -u username</pre>
|-
|Statistics of a completed job
|[https://slurm.schedmd.com/sacct.html sacct]
|<pre>sacct -j job_id --format=user,jobid,jobname,partition,state,time,start,end,elapsed,alloctres,ncpus,nodelist</pre>
|-
|Pause a job
|[https://slurm.schedmd.com/scontrol.html scontrol]
|<pre>scontrol hold job_id</pre>
|-
|Resume a job
|[https://slurm.schedmd.com/scontrol.html scontrol]
|<pre>scontrol resume job_id</pre>
|-
|Modify attributes of a submitted job
|[https://slurm.schedmd.com/scontrol.html scontrol]
|<pre>scontrol update jobid=job_id TimeLimit=4-00:00:00</pre>
|-
|Display a job's characteristics
|[https://slurm.schedmd.com/scontrol.html scontrol]
|<pre>scontrol show job job_id</pre>
|-
|Display information about nodes and partitions
|[https://slurm.schedmd.com/sinfo.html sinfo]
|<pre>sinfo -a</pre>
|}
==SLURM Job Pending Reasons==
    * None -  The job is in the queue, but the SLURM controller has not determined whether it is pending due to priority, resource availability or another constraints. A pending reason will be assigned shortly. 

    * Priority - The job is waiting for its turn based on priority.

    * ReqNodeNotAvail - A required node is unavailable (e.g., down, drained, or in maintenance).

    * Resources – There are not enough free resources (CPUs, memory, GPUs, etc.).

    * QOSMaxJobsLimit - The job is exceeding the allowed number of jobs per user/account in QoS policy.

    * QOSMaxCpuPerUserLimit - The job is exceeding the allowed number of CPU core per user in QoS policy.

    * AssocGrpGRES - The jobs is execeeding the group resources (CPU, GPU etc)

    * JobHeldAdmin - The job is on hold due to an admin action.

    * JobHeldUser - The job is on hold due to a user action.

    * PartitionDown - The assigned partition is unavailable.

    * PartitionNodeLimit - The job requests more nodes than allowed in the partition.

    * AssocGrpCpuLimit - The job's account or group has hit its CPU limit.

    *AssocGrpMemLimit - The job's account or group has hit its memory limit.

== Job Submission ==
==== Interactive Jobs ====
Cores = 10, partition = long, Account = research, GPU = 1
<source lang="bash"> 
parithi@ada ~]$ sinteractive -c 10 -p long -A research -g 1 
salloc: Granted job allocation 141
parithi@gnode03:/home/parithi$ 
</source>

==== Batch Jobs ====
Sample script: NAMD

<source lang="bash">
#!/bin/bash
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -n 20
#SBATCH --gres=gpu:2
#SBATCH --mem-per-cpu=2048
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=END

module add namd/2.12

scp ada:/share1/$USER/file.tar

charmrun +p$SLURM_NPROCS namd2 +idlepoll +devices $CUDA_VISIBLE_DEVICES apoa1.namd > output-gpu1.out
</source>

Sample script: Python  
<source lang="bash">
#!/bin/bash
#SBATCH -A $USER
#SBATCH -n 40
#SBATCH --gres=gpu:4
#SBATCH --mem-per-cpu=2048
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=END

module add cuda/8.0
module add cudnn/7-cuda-8.0

scp ada:/share1/$USER/file.tar /scratch

tar xf /scratch/file.tar
python .....

scp /scratch/output  ada:/share1/$USER

</source>

The variable CUDA_VISIBLE_DEVICES holds the ids of assigned GPUs.

== Storage and usage tutorial video==
   The below video can help you to get important details about the usage of Ada Machine.
   https://www.youtube.com/watch?v=U3_pPJgs2Fg

== Running JupyterLab ==

=== 1. Start an Interactive Session ===
After logging into the Ada master node, start an interactive session using the following command:

<pre>
sinteractive -c 10 -g 1 -A research
</pre>

This requests '''10 CPU cores''' and '''1 GPU'''.  
Replace <code>research</code> with the Slurm account you have been granted access to.

=== 2. Install JupyterLab and Notebook ===
Install the required packages:

<pre>
pip install jupyterlab
pip install notebook
</pre>

=== 3. Set a Jupyter Password ===
Set a password for Jupyter:

<pre>
jupyter lab password
</pre>

=== 4. Start Jupyter Notebook on the Allocated Compute Node ===
Run Jupyter on a specific port:

<pre>
jupyter notebook --port port1
</pre>

Example:

<pre>
jupyter notebook --port 9977
</pre>

=== 5. Create an SSH Tunnel from Your Local Machine ===
On your local computer (laptop/desktop), start an SSH tunnel to the allocated compute node:

<pre>
ssh -L <port1>:localhost:<port2> -J <username>@ada.iiit.ac.in <username>@gnodexxx
</pre>

* <code>port1</code>: Any available port on your local machine (e.g., 8888)  
* <code>port2</code>: The port where JupyterLab is running on the compute node
*  Replace <code>gnodexxx</code> with the hostname of the allocated compute node.

=== 6. Access JupyterLab in Your Browser ===
Once the tunnel is established, open your web browser on your local machine and go to:

<pre>
http://localhost:<port1>
</pre>

You should now be able to access your JupyterLab interface.

====Note====
* An SSH server is required on the local machine for this method. On Windows, you can use: 
** [https://www.bitvise.com/ssh-serve bitvise]
*  Alternatively, you can configure a reverse SSH tunnel in your <code>.ssh/config</code> file: [https://unix.stackexchange.com/questions/162093/reverse-ssh-tunnel-in-config Link]
* '''A similar approach can be used to run TensorBoard.'''

==Windows Users==
===Connecting to Ada===
You need to connect to [https://vpn.iiit.ac.in/ IIIT VPN] first. After that, open Command Prompt and enter <code>ssh -X username@ada.iiit.ac.in</code>. Then you'll have to enter your Ada credentials. Once done, you're into the head-node of Ada, inside your <code>/home2/username/</code> directory.

===Mounting Ada to Local===
Windows users can mount their Ada drive to their local computer using [https://www.nsoftware.com/sftp/drive/ SFTP Drive V2], which is free for personal use. Doing this allows you to edit your code and upload/modify other files from your local computer.

# The first step in mounting a remote server with SFTP Drive is to create a new drive. Click the <code>New...</code> button in the Drives tab and configure the drive to add it to the drive list. <br> Use these values for configuration:
#:* Remote Host: <code>ada.iiit.ac.in</code>
#:* Username and Password: Your Ada credentials
#:* Root folder on Server: <code>/home2/username</code>
#:* You can leave other fields with their default values. <br>
# In the top left of the application window click the <code>Start</code> button to mount the enabled drives. Alternatively, right click the System Tray icon and select <code>Start</code>.
# Once running, the options to <code>Stop</code> or <code>Restart</code> the application are enabled. Unmount all drives by clicking the <code>Stop</code> button. The buttons in the main application window or the context menu in the System Tray may be used to manage the application.

Detailed steps and alternative methods are listed in the <code>nsoftware.SFTPDrive.htm</code> document in the <code>help</code> folder of the software. Usually, it will be at <code>C:/Program Files/nsoftware/SFTP Drive V2/help</code>.

== Accessing compute nodes via VSCode==

1. Install the '''RemoteSSH''' extension in the VS Code application on your
local machine.

2. Launch an interactive or batch job on Ada and record the compute node that is assigned to your job.

3. On your local machine, edit the ~/.ssh/config file to include the following configuration entry:

(Replace gnodexxx with the name of your assigned compute node and your_Ada_user_id with your Ada username.)

<source>
Host ada
  HostName ada.iiit.ac.in
  User your_Ada_user_id

Host gnode
  HostName gnodexxx
  User your_Ada_user_id

ProxyCommand ssh -W %h:%p ada
</source>

4. In the VSCode, use the '''RemoteSSH''' extension to connect directly to '''gnode'''.