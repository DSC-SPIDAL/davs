davs
====

Deterministic Annealing Vector Sponge is a scalable and parallel clustering program that operate on both vector space and non-vector spaces. GitHub page for this project is available at http://dsc-spidal.github.io/davs/ 

Publications
-----
Fox, G., Mani, D. R. and Pyne, S. Parallel deterministic annealing
clustering and its application to LC-MS data analysis. IEEE, City, 2013.
Available at http://grids.ucs.indiana.edu/ptliupages/publications/DAVS_IEEE.pdf

Prerequisites
-----
1. Operating System
  * This program is extensively tested and known to work on,
    *  Red Hat Enterprise Linux Server release 5.10 (Tikanga)
    *  Ubuntu 12.04.3 LTS
    *  Ubuntu 12.10
  * This may work in Windows systems depending on the ability to setup OpenMPI properly, however, this has not been tested and we recommend choosing a Linux based operating system instead.
 
2. Java
  * Download Oracle JDK 8 from http://www.oracle.com/technetwork/java/javase/downloads/index.html
  * Extract the archive to a folder named `jdk1.8.0`
  * Set the following environment variables.
  ```
    JAVA_HOME=<path-to-jdk1.8.0-directory>
    PATH=$JAVA_HOME/bin:$PATH
    export JAVA_HOME PATH
  ```
3. Apache Maven
  * Download latest Maven release from http://maven.apache.org/download.cgi
  * Extract it to some folder and set the following environment variables.
  ```
    MVN_HOME=<path-to-Maven-folder>
    $PATH=$MVN_HOME/bin:$PATH
    export MVN_HOME PATH
  ```
4. Habanero Java (HJ) Library
  * Download HJ-lib jar file from http://www.cs.rice.edu/~vs3/hjlib/code/habanero-java-lib-0.1.1.jar
  * Execute the following shell command from the directory containing the HJ-lib jar to install it as a Maven artifact
  ```sh
    mvn install:install-file -DcreateChecksum=true -Dpackaging=jar -Dfile=habanero-java-lib-0.1.1.jar -DgroupId=habanero-java-lib -DartifactId=habanero-java-lib -Dversion=0.1.1;
  ```
  * Find more information about HJ-lib at https://wiki.rice.edu/confluence/display/PARPROG/Download+and+Set+Up
5. OpenMPI
  * We recommend using `OpenMPI 1.8.1` although it works with the previous 1.7 versions. The Java binding is not available in versions prior to 1.7, hence are not recommended. Note, if using a version other than 1.8.1 please remember to set Maven dependency appropriately in the `pom.xml`.
  * Download OpenMPI 1.8.1 from http://www.open-mpi.org/software/ompi/v1.8/downloads/openmpi-1.8.1.tar.gz
  * Extract the archive to a folder named `openmpi-1.8.1`
  * Also create a directory named `build` in some location. We will use this to install OpenMPI
  * Set the following environment variables
  ```
    BUILD=<path-to-build-directory>
    OMPI_181=<path-to-openmpi-1.8.1-directory>
    PATH=$BUILD/bin:$PATH
    LD_LIBRARY_PATH=$BUILD/lib:$LD_LIBRARY_PATH
    export BUILD OMPI_181 PATH LD_LIBRARY_PATH
  ```
  * The instructions to build OpenMPI depend on the platform. Therefore, we highly recommend looking into the `$OMPI_181/INSTALL` file. Platform specific build files are available in `$OMPI_181/contrib/platform` directory.
  * In general, please specify `--prefix=$BUILD` and `--enable-mpi-java` as arguments to `configure` script. If Infiniband is available (highly recommended) specify `--with-verbs=<path-to-verbs-installation>`. In summary, the following commands will build OpenMPI for a Linux system.
  ```
    cd $OMPI_181
    ./configure --prefix=$BUILD --enable-mpi-java
    make;make install
  ```
  * If everything goes well `mpirun --version` will show `mpirun (Open MPI) 1.8.1`. Execute the following command to instal `$OMPI_181/ompi/mpi/java/java/mpi.jar` as a Maven artifact.
  ```
    mvn install:install-file -DcreateChecksum=true -Dpackaging=jar -Dfile=$OMPI_181/ompi/mpi/java/java/mpi.jar -DgroupId=ompi -DartifactId=ompijavabinding -Dversion=1.8.1;
  ```
  * Few examples are available in `$OMPI_181/examples`. Please use `mpijavac` with other parameters similar to `javac` command to compile OpenMPI Java programs. Once compiled `mpirun [options] java -cp <classpath> class-name arguments` command with proper values set as arguments will run the program. 

Building davs
-----
* Check all prerequisites are satisfied before building davs
* Clone this git repository from `git@github.com:DSC-SPIDAL/davs.git` Let's call this directory `davshome`
* Once above two steps are completed, building davs requires only one command, `mvn install`, issued within `davshome`.

Running davs
-----
The following shell script may be used with necessary modifications to run the program.
```sh
#!/bin/bash

# Java classpath. This should include paths to davs dependent jar files and the davs-1.0-ompi1.8.1.jar
# The dependent jar files may be obtained by running mvn dependency:build-classpath command within davshome
cp=<classpath>

# Obtain working directory
wd=`pwd`
# Character x as a variable
x='x'

# A text file listing available nodes
hosts=<path-to-hostfile>
# Number of nodes
nodes=<num-nodes>
# Number of cores per node
corespernode=8

# Options for Java runtime
jopts="-Xms64M -Xmx64M"

# Number of threads to use within one davs process
tpn=<threads-per-process>
# Number of processes per node
ppn=<processes-per-node>
# Total parallelism expressed as a pattern TxPxN
# where T is number of threads per process, P is processes per node, and N is total nodes
pat=$tpn$x$ppn$x$nodes

echo "Running $pat on `date`" >> status.txt
# Invoke MPI to run davs
mpirun --report-bindings --mca btl ^tcp --hostfile $hostfile --map-by node:PE=$(($corespernode / $ppn)) -np $(($nodes*$ppn)) java $jopts -cp $cp edu.indiana.soic.spidal.davs.Program -c config$pat.properties -n $nodes -t $tpn | tee $pat/davs-out.txt
echo "Finished $pat on `date`" >> status.txt
```
The arguments listed in the `mpirun` command fall into three categories.
* OpenMPI Runtime Parameters
  * `--report-bindings` requests OpenMPI runtime to output how processes are mapped to processing elements (cores) in the allocated nodes.
  * `--mca btl ^tcp` instructs to enable transports other than tcp, which is useful when running on Infiniband.
  * `--hostfile` indicates the file listing available nodes. Each node has to be a in a separate line.
  * `--map-by node:PE=$(($corespernode / $ppn))` controls process mapping and binding. This is a topic on its own right, but the specific values in this example requests processes to be mapped by node while binding each to `corespernode/ppn` number of processing elements. A good set of slides on this topic is available at http://www.slideshare.net/jsquyres/open-mpi-explorations-in-process-affinity-eurompi13-presentation
  * `-np $(($nodes*$ppn))` determines the total number of processes to run and in this case it is equal to `nodes*ppn`
* Java Runtime Parameters
  * `$jopts` in this case lists initial and maximum heap sizes for a JVM instance. 
  * `-cp` indicates paths to find required classes where each entry is separated by a `:` (in Linux)
* Program (davs) Parameters
  * `-c` points to the configuration file. This is a Java properties files listing values for each parameter that davs requires. Details on these parameters will follow in a later section.
  * `-n` indicates the total number of nodes
  * `-t` denotes the number of threads to use within one instance of davs

Configuring davs
-----
The following table summarizes the parameters used in davs.



Acknowledgement
-----
We like to express our sincere gratitude to Prof. Vivek Sarkar 
and his team at Rice University for giving us access and 
continuous support for HJ library. We are equally thankful to Prof. Guillermo L�pez Taboada for giving us free unrestricted access to commercially available FastMPJ MPI library, which we evaluated in an earlier internal version. We are also 
thankful to FutureGrid project and its support team for their support with HPC systems. Last but not least OpenMPI community deserves equal recognition for their valuable support.


