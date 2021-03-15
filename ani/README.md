# NAMD + ANI an illustrative story
----------------------------------

Extended information on how to get [this example](https://github.com/RowleyGroup/NNP-MM) to work:

1. setting up the torch environment

```conda create -n torchani
conda activate torchani
conda install pytorch torchvision -c pytorch
pip install torchani

```

You can test this with starting python:

```
python 
>>> import torchani

```
If you get no errors you are good to go. 


Next we need NAMD installed:
1. Download `NAMD_2.14_Linux-x86_64-multicore.tar.gz`
2. `tar xvzf NAMD_2.14_Linux-x86_64-multicore.tar.gz`
3. Use the output of `pwd` to set the right Alias in your bashrc mine looks like this: 
   `alias namd2='/home/ppxasjsm/Projects/namd_ani/NAMD_2.14_Linux-x86_64-multicore/namd2'` 
4. Make sure everything is ready to go: `source ~/.bashrc`, you may need to reactivate your conda environment. Check that running `namd2` gives this error:
   
   ```
 namd2
Charm++> No provisioning arguments specified. Running with a single PE.
         Use +auto-provision to fully subscribe resources or +p1 to silence this message.
Charm++: standalone mode (not using charmrun)
Charm++> Running in Multicore mode: 1 threads (PEs)
Charm++> Using recursive bisection (scheme 3) for topology aware partitions
Converse/Charm++ Commit ID: v6.10.2-0-g7bf00fa-namd-charm-6.10.2-build-2020-Aug-05-556
Warning> Randomization of virtual memory (ASLR) is turned on in the kernel, thread migration may not work! Run 'echo 0 > /proc/sys/kernel/randomize_va_space' as root to disable it, or try running with '+isomalloc_sync'.
CharmLB> Load balancer assumes all CPUs are same.
Charm++> Running on 1 hosts (1 sockets x 10 cores x 2 PUs = 20-way SMP)
Charm++> cpu topology info is gathered in 0.000 seconds.
Info: NAMD 2.14 for Linux-x86_64-multicore
Info: 
Info: Please visit http://www.ks.uiuc.edu/Research/namd/
Info: for updates, documentation, and support information.
Info: 
Info: Please cite Phillips et al., J. Chem. Phys. 153:044130 (2020) doi:10.1063/5.0014475
Info: in all publications reporting results obtained with NAMD.
Info: 
Info: Based on Charm++/Converse 61002 for multicore-linux-x86_64-iccstatic
Info: Built Mon Aug 24 10:08:36 CDT 2020 by jim on belfast.ks.uiuc.edu
Info: 1 NAMD  2.14  Linux-x86_64-multicore  1    azuma  ppxasjsm
Info: Running on 1 processors, 1 nodes, 1 physical nodes.
Info: CPU topology information available.
Info: Charm++/Converse parallel runtime startup completed at 0.00137437 s
CkLoopLib is used in SMP with simple dynamic scheduling (converse-level notification)
FATAL ERROR: No simulation config file specified on command line.
FATAL ERROR: No simulation config file specified on command line.
[Partition 0][Node 0] End of program
   ```

Next we need to get the example to work. The instructions say install the server and client script. I just added the directory to the python path temporarily by doing this:

```
export PYHTONPATH='/home/ppxasjsm/Projects/namd_ani/NNP-MM'

```
This directory contains the `server.py` and `client.py` script. 

To the get the example to work I followed the instructions:

```
python server.py >& server.01.log&
```

wait a little bit (10s or so), then run:

```
namd2  +p4 md.conf > namd_output
```

This is done in the examples directory. My `md.conf` file looks like this:

```
qmforces on
qmParamPDB qmmm.pdb
qmSoftware custom

# XXX change this to absolute path of socket client script
qmexecpath /home/ppxasjsm/Projects/namd_ani/NNP-MM/client.py

# this should be writeable absolute path on a fast file system (e.g., a RAM disk)
qmBaseDir  /dev/shm
QMColumn occ
qmChargeMode none
qmElecEmbed off

```

The example simulation then runs fine. 

# Getting Torch Ani to work for a protein ligand system with AMBER files
------------------------------------------------------------------------

1. Where is the ANI treated region information stored:

```
qmParamPDB
Description: Name of an optional secondary PDB file where the OCCupancy or BETA column has the indications for QM or MM atoms. QM atoms should have an integer bigger than zero (0) and MM atoms should have zero as the beta or occupancy field. The same file may have indications for bonds between a QM atom and an MM atom. This should be provided in case the PDB file passed to the "coordinates" keyword already has data on its beta or occupancy columns, such as when a SMD simulations is being performed.
Acceptable Values: pdb file
```

Modifying the occurpancy column with pdb-tools:

```
pdb_occ -0.00 qmmm.pdb > mold_qmmm.pdb
```

Make sure that the MOL part is still 1.00 though!
