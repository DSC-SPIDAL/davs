package edu.indiana.soic.spidal.davs;

import edu.rice.hj.api.SuspendableException;
import mpi.MPI;
import mpi.MPIException;
import edu.indiana.soic.spidal.general.Box;
import edu.indiana.soic.spidal.mpi.MPIPacket;

import static edu.rice.hj.Module1.*;

public class ClusteringSolution {
    public boolean DistributedExecutionMode = false; // If True, run in distributed mode; if false only Global Clusters
    public double[][] M_alpha_kpointer_; // Probability that point in Cluster pointed to by kpointer
    public int[][] Map_alpha_PointertoCreatedIndex; // Map Pointer to Cluster Created Index for each point
    public int[] NumClusters_alpha_; // Number of Clusters associated with point alpha

    public double[] DiffMsummed_k_; // Change in Malpha_k_ summed over alpha
    public double[] C_k_; // Summation of C(k) over parallel threads
    public double[] ClusterScaledSquaredWidth_k_; // Mean Scaled (by Sigma_k_i_) Squared Width

    public double[] FreezingMeasure_k_; // Freezing Measure for Cluster
    public double[] P_k_; // P(k) used in Continuous Clustering
    public int[] Splittable_k_; // If -1 used already; 0 ignored, =1 Try to find eigenvalue; =2 Eigenvalue > 0; = 3 Eigenvalue < 0
    public double[] Eigenvalue_k; // Eigenvalue in Y_k_i_/SQRT(Sigma_k_i_) space
    public int[] OccupationCounts_k_; // Occupation Counts for Clusters
    public int[] SplitPriority_k_; // -1 Infinite priority; 0 zero priority; 1 child just removed otherwise priority increases as value increases

    public double[][] Y_k_i_; // Y(k,i) is center of Cluster k where i runs over vector coordinates (NOT used for Sponge)
    public double[][] YPrevious_k_i_; // YPrevious(k,i) is center of Cluster k where i runs over vector coordinates at previous iteration
    public boolean[] YPreviousActuallySet; // If true this Y has YPrevious set
    public double[][][] Correlation_k_i_j; // Correlation that is second term in second derivation -- Symmetric -- Diagonal Sum is Square of Cluster Width
    public double[][] DC_k_DY_k_i_; // Derivative of C[k] used in splitting
    public double[][] Sigma_k_i_; // Squared Standard Deviation of Cluster k calculated as controlled by Program.SigmaMethod
    public double[][] Eigenvector_k_i; // Splitting vector for cluster k in Y_k_/SQRT(Sigma_k_i_) -- only set for clusters that are candidates

    public int ClustertoSplit = -1; // Current Cluster being split

    public int SpongeCluster = -1; // Label of Sponge
    public double PairwiseHammy = 0.0; // Value of Hamiltonian
    public double OldHammy = 0.0; //  Previous value of Hamiltonian
    public int Ncent_ThisNode = 1; //  the current number of clusters stored in this node
    public int Ncent_Global = 1; //  the current number of clusters summed over all nodes
    public double ActualCoolingFactor = Program.InitialCoolingFactor; // Actual Cooling Factor
    public double Temperature; // The current temperature
    public int IterationSetAt = 0; // Iteration set at
    public boolean Eigenvectorset = false; // If True Eigenvector set
    public int YPreviousSet = -1; // If -1 Y not set, = 0 Y but not YPrevious set, = 1 Y and YPrevious set
    public boolean SolutionSet = false; // If True Solution Set
    public boolean CorrelationsSet = false; // If true Correlation Set
    public int DiffMalpha_k_Set = -1; // -1 initial, 0 set Previous_Malpha_k_, 1 Set
    public double[] AverageWidth; // Average widths over Clusters in each dimension
    public double TotaloverVectorIndicesAverageWidth; // Average Width summed over all dimensions

    public int[] LocalCreatedIndex; // Created Index of Local Cluster
    public int[] LocalStatus; // Status -2 Moved -1 deleted 0 Sponge 1 Global 2 Distributed Stored here 3 Distributed Stored Elsewhere (Not Possible)
    public int[] LocalSplitCreatedIndex; // 1 + CreatedIndex of Parent or -1 -CreatedIndex of child if Split UNTIL SPLIT PROPAGATED

    //  Note quantities below are static and NOT saved in backup. However they change during job
    //  Backup must be followed by a master synchronization setting these variables
    public static int NumberLocalActiveClusters; // Number of Active Clusters stored here -- arrays later on are defined for this range NOT NumberAvailableActiveClusters
    public static int NumberAvailableActiveClusters; // Number of Active Clusters stored here or transported
    public static int NumberGlobalClusters; // Number of status 0 and 1 Clusters
    public static int[] RealClusterIndices; // Maps Active Cluster Indices into Real cluster list
    public static int[] ActiveClusterIndices; // Cluster Numbers of Active Clusters -- maps  index of main arrays into ActiveCluster index
    public static int[] LocalNodeAccPosition; // Current Node Accumulation Positions for this cluster labelled by ActiveCluster Index
    public static int[][] LocalThreadAccPosition; // Current Thread Accumulation Position for this cluster labelled by ActiveCluster Index
    public static int[] LocalHost; // =0 unless Distributed when Packed Host index in form H + MULTIPLIER (H1 + MULTIPLIER * H2 ) where H must be current host labelled by ActiveCluster Index

    //  These quantities are really static. Set at Start and do not change
    public static int NumberofPointsinProcess = -1; // Number of Points in Process
    public static int MaximumCreatedClusters; // Number of allowed created clusters
    public static int MaximumNumberClusterspernode = 0; // Maximum Number of Centers per node summed over points
    public static int MaximumNumberClustersTotal = 0; // Maximum Number of Centers per node summed over points
    public static int TargetClustersperPoint = 0; // Target maximum number of clusters per point
    public static int TargetMinimumClustersperPoint = 0; // Target minimum number of clusters per point
    public static int MaximumClustersperPoint = 0; // Actual maximum number of clusters per point
    public static double ExponentArgumentCut = 20.0; // Include all clusters with (Y(cluster)-X(point))^2 / (2 * Temperature) < ExpArgumentCut in Point Collection
    public static int cachelinesize = 0; // Cacheline increment

    public static int CurrentIteration; // Current Iteration Number to monitor which clusters stale
    public static ClusterIndirection[] UniversalMapping; // Maps CreatedIndex into location -- in distributed model any one node only has partial information
    public static int CenterMaxforCreatedIndex = 0; // Center count to use for CreatedIndex
    public static FullSolution TotalClusterSummary; // Summary of Clusters

    public static int PACKINGMASK; // 2^PACKINGSHIFT - 1
    public static int PACKINGSHIFT;
    public static int PACKINGMULTIPLIER; // 2^PACKINGSHIFT

    public static int ClustersDeleted = 0; // Clusters Deleted this iteration
    public static int ClustersMoved = 0; // Clusters Moved this iteration
    public static int ClustersSplit = 0; // Clusters Split this iteration

    // public static int DEBUGPrint1 = 0;

    public static void SetParameters(int NumberofPointsinProcessINPUT, int MaximumCreatedClustersINPUT, int MaximumNumberClustersTotalINPUT, int MaximumNumberClusterspernodeINPUT, int cachelinesizeINPUT, int TargetClustersperPointINPUT, int MinimumNumberClusterspernodeINPUT, int MaximumClustersperPointINPUT, double ExponentArgumentCutINPUT) {
        if (NumberofPointsinProcess > 0) {
            return;
        }
        NumberofPointsinProcess = NumberofPointsinProcessINPUT;
        MaximumCreatedClusters = MaximumCreatedClustersINPUT;
        MaximumNumberClusterspernode = MaximumNumberClusterspernodeINPUT;
        TargetMinimumClustersperPoint = MinimumNumberClusterspernodeINPUT;
        MaximumNumberClustersTotal = MaximumNumberClustersTotalINPUT;
        cachelinesize = cachelinesizeINPUT;
        TargetClustersperPoint = TargetClustersperPointINPUT;
        MaximumClustersperPoint = MaximumClustersperPointINPUT;
        ExponentArgumentCut = ExponentArgumentCutINPUT;

        LocalNodeAccPosition = new int[MaximumNumberClusterspernode];
        LocalThreadAccPosition = new int[MaximumNumberClusterspernode][];
        LocalHost = new int[MaximumNumberClusterspernode];
        for (int ClusterIndex = 0; ClusterIndex < MaximumNumberClusterspernode; ClusterIndex++) {
            LocalThreadAccPosition[ClusterIndex] = new int[DAVectorUtility.ThreadCount];
        }

        CurrentIteration = 0;

        PACKINGSHIFT = 0;
        PACKINGMASK = 1;
        while (PACKINGMASK <= DAVectorUtility.MPI_Size) {
            PACKINGMASK = 2 * PACKINGMASK;
            PACKINGSHIFT++;
        }
        PACKINGMULTIPLIER = PACKINGMASK;
        PACKINGMASK--;


        // Set Arrays that Map Created Indices to Local Storage
        int UniversalSize = (1 + MaximumCreatedClusters) * PACKINGMULTIPLIER;
        DAVectorUtility.SALSAPrint(0, "Create Universal Mapping " + UniversalSize);
        UniversalMapping = new ClusterIndirection[UniversalSize];
        DAVectorUtility.SALSAPrint(0, "Set up Full Solution " + MaximumNumberClustersTotal);
        TotalClusterSummary = new FullSolution(MaximumNumberClustersTotal);
        DAVectorUtility.SALSAPrint(0, "End ClusteringSolution");
    }

    public static int SetCreatedIndex(int RealLocalClusterIndex) {
        int host = 0;
        if (ParallelClustering.runningSolution.DistributedExecutionMode && ((ParallelClustering.runningSolution.SpongeCluster < 0) || (ParallelClustering.runningSolution.SpongeCluster != RealLocalClusterIndex))) {
            host = DAVectorUtility.MPI_Rank;
        }
        ++CenterMaxforCreatedIndex;
        if (CenterMaxforCreatedIndex > ClusteringSolution.MaximumCreatedClusters) {
            DAVectorUtility.printAndThrowRuntimeException(" Center Index too large " + ClusteringSolution.CenterMaxforCreatedIndex + " " + ParallelClustering.runningSolution.Ncent_ThisNode + " " + ClusteringSolution.MaximumNumberClustersTotal + " " + ClusteringSolution.MaximumCreatedClusters);

        }
        int CreatedIndex = host + 1 + (CenterMaxforCreatedIndex << ClusteringSolution.PACKINGSHIFT);
        UniversalMapping[CreatedIndex] = new ClusterIndirection(CurrentIteration, RealLocalClusterIndex + 1);
        ParallelClustering.runningSolution.LocalCreatedIndex[RealLocalClusterIndex] = CreatedIndex;
        return CreatedIndex;

    } // End SetCreatedIndex()

    public static int MapActivetoCreatedIndex(int ActiveClusterIndex, ClusteringSolution Solution) {
        if (Solution.DistributedExecutionMode && (ActiveClusterIndex >= NumberLocalActiveClusters)) {
            return DistributedClusteringSolution.StorageforTransportedClusters.TotalTransportedCreatedIndex[ActiveClusterIndex - NumberLocalActiveClusters];
        } else {
            return Solution.LocalCreatedIndex[RealClusterIndices[ActiveClusterIndex]];
        }


    } // End MapActivetoCreatedIndex

    public ClusteringSolution(boolean thereisasponge) throws MPIException {
        if (NumberofPointsinProcess < 0) {
            DAVectorUtility.printAndThrowRuntimeException("NumberofPointsinProcess Unset");

        }
        this.Ncent_ThisNode = 1;
        this.SpongeCluster = -1;
        int FirstRealCluster = 0;
        this.DistributedExecutionMode = false;

        if (thereisasponge) {
            this.SpongeCluster = 0;
            this.Ncent_ThisNode = 2;
            FirstRealCluster = 1;
        }
        this.Ncent_Global = this.Ncent_ThisNode;

        //  Point Arrays
        DAVectorUtility.SALSAPrint(0, "Large Arrays Started " + NumberofPointsinProcess + " Clusters per Point " + MaximumClustersperPoint);
        M_alpha_kpointer_ = new double[NumberofPointsinProcess][];
        Map_alpha_PointertoCreatedIndex = new int[NumberofPointsinProcess][];
        NumClusters_alpha_ = new int[NumberofPointsinProcess];

        // Note - parallel for
        final int FirstRealClusterLoopVar = FirstRealCluster;
        try {
            forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) ->
            {
                int indexlen = DAVectorUtility.PointsperThread[threadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[threadIndex] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++) {
                    M_alpha_kpointer_[alpha] = new double[MaximumClustersperPoint];
                    Map_alpha_PointertoCreatedIndex[alpha] = new int[MaximumClustersperPoint];
                    Map_alpha_PointertoCreatedIndex[alpha][0] = 0;
                    Map_alpha_PointertoCreatedIndex[alpha][FirstRealClusterLoopVar] = FirstRealClusterLoopVar;
                    NumClusters_alpha_[alpha] = 1 + FirstRealClusterLoopVar;
                }

            });
        } catch (SuspendableException e) {
            DAVectorUtility.printAndThrowRuntimeException(e.getMessage());
        }

        DAVectorUtility.SALSAPrint(0, "Large Arrays Created " + NumberofPointsinProcess);

        DiffMsummed_k_ = new double[MaximumNumberClusterspernode + cachelinesize];
        C_k_ = new double[MaximumNumberClusterspernode + cachelinesize];
        ClusterScaledSquaredWidth_k_ = new double[MaximumNumberClusterspernode + cachelinesize];
        P_k_ = new double[MaximumNumberClusterspernode + cachelinesize];
        FreezingMeasure_k_ = new double[MaximumNumberClusterspernode + cachelinesize];
        OccupationCounts_k_ = new int[MaximumNumberClusterspernode + cachelinesize];
        Splittable_k_ = new int[MaximumNumberClusterspernode + cachelinesize];
        SplitPriority_k_ = new int[MaximumNumberClusterspernode + cachelinesize];
        Eigenvalue_k = new double[MaximumNumberClusterspernode + cachelinesize];
        ClustertoSplit = -1;

        Y_k_i_ = new double[MaximumNumberClusterspernode][]; // Y(k,i) is center of Cluster k where i runs over vector coordinates
        YPrevious_k_i_ = new double[MaximumNumberClusterspernode][];
        YPreviousActuallySet = new boolean[MaximumNumberClusterspernode];
        Eigenvector_k_i = new double[MaximumNumberClusterspernode][];
        DC_k_DY_k_i_ = new double[MaximumNumberClusterspernode][];
        Sigma_k_i_ = new double[MaximumNumberClusterspernode][]; // Standard Deviation Squared of Cluster k calculated as controlled by Program.SigmaMethod
        Correlation_k_i_j = new double[MaximumNumberClusterspernode][][]; // Correlation that is second term in second derivation -- Symmetric -- Diagonal Sum is Square of Cluster Width

        RealClusterIndices = new int[MaximumNumberClusterspernode];
        ActiveClusterIndices = new int[MaximumNumberClusterspernode];
        LocalCreatedIndex = new int[MaximumNumberClusterspernode];
        LocalStatus = new int[MaximumNumberClusterspernode];
        LocalSplitCreatedIndex = new int[MaximumNumberClusterspernode];

        for (int ClusterIndex = 0; ClusterIndex < MaximumNumberClusterspernode; ClusterIndex++) {
            YPreviousActuallySet[ClusterIndex] = false;
            Y_k_i_[ClusterIndex] = new double[Program.ParameterVectorDimension];
            YPrevious_k_i_[ClusterIndex] = new double[Program.ParameterVectorDimension];
            Eigenvector_k_i[ClusterIndex] = new double[Program.ParameterVectorDimension];
            DC_k_DY_k_i_[ClusterIndex] = new double[Program.ParameterVectorDimension];
            Sigma_k_i_[ClusterIndex] = new double[Program.ParameterVectorDimension];
            int correlsize = 1;
            if (Program.CalculateCorrelationMatrix) {
                correlsize = Program.ParameterVectorDimension;
            }
            Correlation_k_i_j[ClusterIndex] = new double[correlsize][correlsize];
        }

        AverageWidth = new double[Program.ParameterVectorDimension];
        Eigenvectorset = false;
        SolutionSet = false;
        CorrelationsSet = false;
        DiffMalpha_k_Set = -1;


        if (thereisasponge) {
            SplitPriority_k_[this.SpongeCluster] = 0;
            Splittable_k_[this.SpongeCluster] = 0;
            LocalStatus[this.SpongeCluster] = 0;
        }
        SetActiveClusters();
        DAVectorUtility.SALSAPrint(0, "Clusters set up");

    } // End ClusteringSolution

    // Set Active Cluster Array and count
    public final void SetActiveClusters() throws MPIException { // Update Ncent_Global and basic static arrays for locally controlled clusters
        //  Need to set accumulation positions and Host index separately
        //  Need to set remote clusters transported locally separately in
        int ActiveCount = 0;
        int ActiveCount2 = 0;
        for (int RealClusterIndex = 0; RealClusterIndex < this.Ncent_ThisNode; RealClusterIndex++) {
            ActiveClusterIndices[RealClusterIndex] = -1;
            if (this.LocalStatus[RealClusterIndex] < 0 || this.LocalStatus[RealClusterIndex] > 2) {
                continue;
            }
            if (this.LocalStatus[RealClusterIndex] == 2) {
                ++ActiveCount2;
            }
            ActiveClusterIndices[RealClusterIndex] = ActiveCount;
            RealClusterIndices[ActiveCount] = RealClusterIndex;
            ++ActiveCount;
        }
        NumberLocalActiveClusters = ActiveCount;
        NumberAvailableActiveClusters = ActiveCount; // This is incremented on Major Synchronizations
        int ActiveCount01 = ActiveCount - ActiveCount2;
        if (this.DistributedExecutionMode) {
            DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming2);
            // Note - MPI Call - Allreduce - int - sum
            if (DAVectorUtility.MPI_Size > 1){
//                ActiveCount2 = DAVectorUtility.MPI_communicator.<Integer>Allreduce(ActiveCount2, Operation < Integer >.Add);
                ActiveCount2 = DAVectorUtility.mpiOps.allReduce(ActiveCount2, MPI.SUM);
            }
            DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming2);
        }

        //   Check agreement on Global Count
        NumberGlobalClusters = ActiveCount01;
        DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming2);

        // Note - MPI Call - Allreduce - int - max
        if (DAVectorUtility.MPI_Size > 1){
//            NumberGlobalClusters = DAVectorUtility.MPI_communicator.<Integer>Allreduce(NumberGlobalClusters, Operation < Integer >.Max);
            NumberGlobalClusters = DAVectorUtility.mpiOps.allReduce(NumberGlobalClusters,MPI.MAX);
        }

        DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming2);
        if (ActiveCount01 != NumberGlobalClusters) {
            DAVectorUtility.printAndThrowRuntimeException("Error in Global Count Max " + NumberGlobalClusters + " Actual " + ActiveCount01);

        }
        this.Ncent_Global = ActiveCount2 + ActiveCount01;

    } // End SetActiveClusters()

    //  These widths are divided by Sigma_k_i_ and NOT Square Rooted
    public final void SetClusterWidths() throws MPIException {
        DAVectorUtility.StartSubTimer(4);
        final GlobalReductions.FindIndirectVectorDoubleSum FindWidthsGlobal;
        final DistributedReductions.FindIndirectMultiVectorDoubleSum FindWidthsinDistributedMode;
        final int BeginFindWidthsinDistributedMode;

        if (this.DistributedExecutionMode) {
            FindWidthsinDistributedMode = new DistributedReductions.FindIndirectMultiVectorDoubleSum();
            BeginFindWidthsinDistributedMode = FindWidthsinDistributedMode.AddComponents(1);
            FindWidthsinDistributedMode.NodeInitialize();
            FindWidthsGlobal = null;
        } else {
            FindWidthsGlobal = new GlobalReductions.FindIndirectVectorDoubleSum(DAVectorUtility.ThreadCount, ClusteringSolution.NumberLocalActiveClusters);
            FindWidthsinDistributedMode = null;
            BeginFindWidthsinDistributedMode = -1;
        }

        GlobalReductions.FindVectorDoubleSum FindAverageWidth = new GlobalReductions.FindVectorDoubleSum(DAVectorUtility.ThreadCount, Program.ParameterVectorDimension);

        // Note - parallel for
        try {
            forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) -> {
                if (this.DistributedExecutionMode) {
                    FindWidthsinDistributedMode.ThreadInitialize(threadIndex);
                } else {
                    FindWidthsGlobal.startthread(threadIndex);
                }
                int ThreadStorePosition = -1;
                double[] Sigma_Pointer;
                double[] Y_Pointer;
                double[] WorkSpace_Avg = new double[Program.ParameterVectorDimension];
                int ArraySize = Math.min(ClusteringSolution.NumberAvailableActiveClusters,
                        ClusteringSolution.MaximumClustersperPoint);
                int[] ActiveClustersperPoint = new int[ArraySize];
                double[] WorkSpace_k_ = new double[ArraySize];
                int indexlen = DAVectorUtility.PointsperThread[threadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[threadIndex] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++) {
                    int IndirectSize = this.NumClusters_alpha_[alpha];
                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++) {
                        int RealClusterIndex = -1;
                        int RemoteIndex = -1;
                        int ActiveClusterIndex = -1;
                        Box<Integer> tempRef_RealClusterIndex = new Box<>(
                                RealClusterIndex);
                        Box<Integer> tempRef_ActiveClusterIndex = new Box<>(
                                ActiveClusterIndex);
                        Box<Integer> tempRef_RemoteIndex = new Box<>(RemoteIndex);
                        VectorAnnealIterate.ClusterPointersforaPoint(alpha, IndirectClusterIndex, tempRef_RealClusterIndex,
                                tempRef_ActiveClusterIndex, tempRef_RemoteIndex);
                        RealClusterIndex = tempRef_RealClusterIndex.content;
                        ActiveClusterIndex = tempRef_ActiveClusterIndex.content;
                        RemoteIndex = tempRef_RemoteIndex.content;
                        if (RemoteIndex < 0) {
                            Y_Pointer = this.Y_k_i_[RealClusterIndex];
                            Sigma_Pointer = this.Sigma_k_i_[RealClusterIndex];
                            if (this.DistributedExecutionMode) {
                                ThreadStorePosition = ClusteringSolution.LocalThreadAccPosition[ActiveClusterIndex][threadIndex];
                            }
                        } else {
                            Y_Pointer = DistributedClusteringSolution.StorageforTransportedClusters.TotalTransportedY_t_i[RemoteIndex];
                            Sigma_Pointer = DistributedClusteringSolution.StorageforTransportedClusters.TotalTransportedSigma_t_i[RemoteIndex];
                            ThreadStorePosition = DistributedClusteringSolution.StorageforTransportedClusters.TransportedThreadAccPosition[RemoteIndex][threadIndex];
                        }
                        ActiveClustersperPoint[IndirectClusterIndex] = ActiveClusterIndex;
                        double wgt = this.M_alpha_kpointer_[alpha][IndirectClusterIndex];
                        if ((SpongeCluster < 0) || (RealClusterIndex != SpongeCluster)) {
                            WorkSpace_k_[IndirectClusterIndex] = wgt * DAVectorParallelism.getSquaredScaledDistancePointActiveCluster(
                                    alpha, ActiveClusterIndex, this);
                            for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++) {
                                double tmp = Program.PointPosition[alpha][VectorIndex] - Y_Pointer[VectorIndex];
                                WorkSpace_Avg[VectorIndex] = wgt * tmp * tmp / Sigma_Pointer[VectorIndex];
                            }
                            FindAverageWidth.addapoint(threadIndex, WorkSpace_Avg);
                        } else {
                            WorkSpace_k_[IndirectClusterIndex] = wgt * Program.SpongeFactor * Program.SpongeFactor;
                        }
                        if (this.DistributedExecutionMode) {
                            FindWidthsinDistributedMode.addapoint(threadIndex, ThreadStorePosition,
                                    BeginFindWidthsinDistributedMode, WorkSpace_k_[IndirectClusterIndex]);
                        }
                    }
                    if (!this.DistributedExecutionMode) {
                        FindWidthsGlobal.addapoint(threadIndex, IndirectSize, ActiveClustersperPoint, WorkSpace_k_);
                    }
                }

            });
        } catch (SuspendableException e) {
            DAVectorUtility.printAndThrowRuntimeException(e.getMessage());
        }

        if (this.DistributedExecutionMode) {
            FindWidthsinDistributedMode.sumoverthreadsandmpi();
        } else {
            FindWidthsGlobal.sumoverthreadsandmpi();
        }

        // Broadcast zero size information for global clusters
        boolean[] ZeroSizeClusters = new boolean[ClusteringSolution.NumberGlobalClusters];
        int CountGlobalClusters = 0;
        for (int ActiveClusterIndex = 0; ActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; ActiveClusterIndex++) {
            int RealClusterIndex = ClusteringSolution.RealClusterIndices[ActiveClusterIndex];
            if ((ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] < 0) || (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] > 1)) {
                continue;
            }
            boolean zerosizecluster = false;
            if (RealClusterIndex != ParallelClustering.runningSolution.SpongeCluster) {
                zerosizecluster = ParallelClustering.runningSolution.C_k_[RealClusterIndex] <= Program.CountforCluster_C_ktobezero;
            }
            ZeroSizeClusters[CountGlobalClusters] = zerosizecluster;
            ++CountGlobalClusters;
        }
        DAVectorUtility.StartSubTimer(DAVectorUtility.MPIBROADCASTTiming);

        // Note - MPI Call - Broadcast - boolean[]
        if (DAVectorUtility.MPI_Size > 1){
//            DAVectorUtility.MPI_communicator.<Boolean>Broadcast(ZeroSizeClusters, 0);
            DAVectorUtility.mpiOps.broadcast(ZeroSizeClusters,0);
        }
        DAVectorUtility.StopSubTimer(DAVectorUtility.MPIBROADCASTTiming);

        CountGlobalClusters = 0;
        for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++) {
            int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
            double Result;
            if (this.DistributedExecutionMode) {
                Result = FindWidthsinDistributedMode.TotalVectorSum[ClusteringSolution.LocalNodeAccPosition[LocalActiveClusterIndex]][BeginFindWidthsinDistributedMode];
            } else {
                Result = FindWidthsGlobal.TotalVectorSum[LocalActiveClusterIndex];
            }

            boolean zerosizecluster = false;
            if (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] < 2) {
                zerosizecluster = ZeroSizeClusters[CountGlobalClusters];
                ++CountGlobalClusters;
            } else if (ParallelClustering.runningSolution.DistributedExecutionMode) {
                zerosizecluster = ParallelClustering.runningSolution.C_k_[RealClusterIndex] <= Program.CountforCluster_C_ktobezero;
            }

            if (RealClusterIndex == ParallelClustering.runningSolution.SpongeCluster) {
                this.ClusterScaledSquaredWidth_k_[RealClusterIndex] = Program.SpongeFactor * Program.SpongeFactor;
                continue;
            }

            if (zerosizecluster) {
                this.ClusterScaledSquaredWidth_k_[RealClusterIndex] = 0.0;
            } else {
                this.ClusterScaledSquaredWidth_k_[RealClusterIndex] = Result / this.C_k_[RealClusterIndex];
            }
        }

        double C_k_Sum1not0 = 0.0;
        double C_k_Sum2 = 0.0;
        for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++) {
            int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
            if (this.LocalStatus[RealClusterIndex] == 2) {
                C_k_Sum2 += this.C_k_[RealClusterIndex];
            } else if (this.LocalStatus[RealClusterIndex] == 1) {
                C_k_Sum1not0 += this.C_k_[RealClusterIndex];
            }
        }
        if (this.DistributedExecutionMode) {
            DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming2);
            // Note - MPI Call - Allreduce - double - sum
            if (DAVectorUtility.MPI_Size > 1){
//                C_k_Sum2 = DAVectorUtility.MPI_communicator.<Double>Allreduce(C_k_Sum2, Operation < Double >.Add);
                C_k_Sum2 = DAVectorUtility.mpiOps.allReduce(C_k_Sum2, MPI.SUM);
            }
            DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming2);
        }
        double C_Sum = C_k_Sum1not0 + C_k_Sum2;

        //   No contribution from Sponge here
        FindAverageWidth.sumoverthreadsandmpi();
        this.TotaloverVectorIndicesAverageWidth = 0.0;
        for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++) {
            this.AverageWidth[VectorIndex] = FindAverageWidth.TotalVectorSum[VectorIndex] / C_Sum;
            this.TotaloverVectorIndicesAverageWidth += this.AverageWidth[VectorIndex];
        }
        DAVectorUtility.StopSubTimer(4);

    } // End SetClusterWidths

    public final double SetClusterRadius(int RealClusterIndex) { // Note this definition -- Sqrt(squared width) is different from Maximum distance of Points in Cluster from Center as for DA we don't assign rigourously points to clusters

        return Math.sqrt(this.ClusterScaledSquaredWidth_k_[RealClusterIndex]);
    }

    //  Set the C_k_ for Clusters in this Solution
    public final void SetClusterSizes() throws MPIException {
        DAVectorUtility.StartSubTimer(4);
        final GlobalReductions.FindIndirectVectorDoubleSum FindSizes;
        final GlobalReductions.FindIndirectVectorDoubleSum FindFreezing;
        final DistributedReductions.FindIndirectMultiVectorDoubleSum Find2Components;
        final int BeginFindSizes;
        final int BeginFindFreezing;

        if (this.DistributedExecutionMode) {
            Find2Components = new DistributedReductions.FindIndirectMultiVectorDoubleSum();
            BeginFindSizes = Find2Components.AddComponents(1);
            BeginFindFreezing = Find2Components.AddComponents(1);
            Find2Components.NodeInitialize();
            FindSizes = null;
            FindFreezing = null;
        } else {
            FindSizes = new GlobalReductions.FindIndirectVectorDoubleSum(DAVectorUtility.ThreadCount, ClusteringSolution.NumberLocalActiveClusters);
            FindFreezing = new GlobalReductions.FindIndirectVectorDoubleSum(DAVectorUtility.ThreadCount, ClusteringSolution.NumberLocalActiveClusters);
            Find2Components = null;
            BeginFindSizes = -1;
            BeginFindFreezing = -1;
        }

        // Note - parallel for
        try {
            forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) -> {
                if (this.DistributedExecutionMode) {
                    Find2Components.ThreadInitialize(threadIndex);
                } else {
                    FindSizes.startthread(threadIndex);
                    FindFreezing.startthread(threadIndex);
                }
                int ThreadStorePosition = -1;
                int ArraySize = Math.min(ClusteringSolution.NumberAvailableActiveClusters,
                        ClusteringSolution.MaximumClustersperPoint);
                int[] ActiveClustersperPoint = new int[ArraySize];
                double[] WorkSpace_k_ = new double[ArraySize];
                int indexlen = DAVectorUtility.PointsperThread[threadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[threadIndex] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++) {
                    int IndirectSize = this.NumClusters_alpha_[alpha];
                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++) {
                        int RealClusterIndex = -1;
                        int RemoteIndex = -1;
                        int ActiveClusterIndex = -1;
                        Box<Integer> tempRef_RealClusterIndex = new Box<>(
                                RealClusterIndex);
                        Box<Integer> tempRef_ActiveClusterIndex = new Box<>(
                                ActiveClusterIndex);
                        Box<Integer> tempRef_RemoteIndex = new Box<>(RemoteIndex);
                        VectorAnnealIterate.ClusterPointersforaPoint(alpha, IndirectClusterIndex, tempRef_RealClusterIndex,
                                tempRef_ActiveClusterIndex, tempRef_RemoteIndex);
                        RealClusterIndex = tempRef_RealClusterIndex.content;
                        ActiveClusterIndex = tempRef_ActiveClusterIndex.content;
                        RemoteIndex = tempRef_RemoteIndex.content;
                        if (RemoteIndex < 0) {
                            if (this.DistributedExecutionMode) {
                                ThreadStorePosition = ClusteringSolution.LocalThreadAccPosition[ActiveClusterIndex][threadIndex];
                            }
                        } else {
                            ThreadStorePosition = DistributedClusteringSolution.StorageforTransportedClusters.TransportedThreadAccPosition[RemoteIndex][threadIndex];
                        }
                        ActiveClustersperPoint[IndirectClusterIndex] = ActiveClusterIndex;
                        WorkSpace_k_[IndirectClusterIndex] = this.M_alpha_kpointer_[alpha][IndirectClusterIndex] * (1.0 - this.M_alpha_kpointer_[alpha][IndirectClusterIndex]);
                        if (this.DistributedExecutionMode) {
                            if ((ThreadStorePosition < 0) || (ThreadStorePosition >= DistributedClusteringSolution.NodeAccMetaData.NumberofPointsperThread[threadIndex])) {
                                DAVectorUtility.printAndThrowRuntimeException(
                                        " Find Sizes Point " + (alpha + DAVectorUtility.PointStart_Process) + " Cluster Index " + (IndirectClusterIndex) + " Illegal Position " + (ThreadStorePosition) + " Max " + (DistributedClusteringSolution.NodeAccMetaData.NumberofPointsperThread[threadIndex]) + " CreatedIndex " + (ParallelClustering.runningSolution.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex]) + " Node Total " + (DistributedClusteringSolution.NodeAccMetaData.NumberofPointsperNode));

                            }
                            Find2Components.addapoint(threadIndex, ThreadStorePosition, BeginFindSizes,
                                    this.M_alpha_kpointer_[alpha][IndirectClusterIndex]);
                            Find2Components.addapoint(threadIndex, ThreadStorePosition, BeginFindFreezing,
                                    WorkSpace_k_[IndirectClusterIndex]);
                        }
                    }
                    if (!this.DistributedExecutionMode) {
                        FindSizes.addapoint(threadIndex, IndirectSize, ActiveClustersperPoint,
                                this.M_alpha_kpointer_[alpha]);
                        FindFreezing.addapoint(threadIndex, IndirectSize, ActiveClustersperPoint, WorkSpace_k_);
                    }
                }

            });
        } catch (SuspendableException e) {
            DAVectorUtility.printAndThrowRuntimeException(e.getMessage());
        }

        if (this.DistributedExecutionMode) {
            Find2Components.sumoverthreadsandmpi();
        } else {
            FindSizes.sumoverthreadsandmpi();
            FindFreezing.sumoverthreadsandmpi();
        }

        for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++) {
            int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
            if (this.DistributedExecutionMode) {
                this.C_k_[RealClusterIndex] = Find2Components.TotalVectorSum[ClusteringSolution.LocalNodeAccPosition[LocalActiveClusterIndex]][BeginFindSizes];
            } else {
                this.C_k_[RealClusterIndex] = FindSizes.TotalVectorSum[LocalActiveClusterIndex];
            }
        }

        // Broadcast zero size information for global clusters
        boolean[] ZeroSizeClusters = new boolean[ClusteringSolution.NumberGlobalClusters];
        int CountGlobalClusters = 0;
        for (int ActiveClusterIndex = 0; ActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; ActiveClusterIndex++) {
            int RealClusterIndex = ClusteringSolution.RealClusterIndices[ActiveClusterIndex];
            if ((ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] < 0) || (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] > 1)) {
                continue;
            }
            boolean zerosizecluster = false;
            if (RealClusterIndex != ParallelClustering.runningSolution.SpongeCluster) {
                zerosizecluster = ParallelClustering.runningSolution.C_k_[RealClusterIndex] <= Program.CountforCluster_C_ktobezero;
            }
            ZeroSizeClusters[CountGlobalClusters] = zerosizecluster;
            ++CountGlobalClusters;
        }
        DAVectorUtility.StartSubTimer(DAVectorUtility.MPIBROADCASTTiming);

        // Note - MPI Call - Broadcast - boolean[]
        if (DAVectorUtility.MPI_Size > 1){
//            DAVectorUtility.MPI_communicator.<Boolean>Broadcast(ZeroSizeClusters, 0);
            DAVectorUtility.mpiOps.broadcast(ZeroSizeClusters, 0);
        }
        DAVectorUtility.StopSubTimer(DAVectorUtility.MPIBROADCASTTiming);
        CountGlobalClusters = 0;
        for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++) {
            int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
            if (this.DistributedExecutionMode) {
                this.FreezingMeasure_k_[RealClusterIndex] = Find2Components.TotalVectorSum[ClusteringSolution.LocalNodeAccPosition[LocalActiveClusterIndex]][BeginFindFreezing];
            } else {
                this.FreezingMeasure_k_[RealClusterIndex] = FindFreezing.TotalVectorSum[LocalActiveClusterIndex];
            }

            boolean zerosizecluster = false;
            if (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] < 2) {
                zerosizecluster = ZeroSizeClusters[CountGlobalClusters];
                ++CountGlobalClusters;
            } else if (ParallelClustering.runningSolution.DistributedExecutionMode) {
                zerosizecluster = ParallelClustering.runningSolution.C_k_[RealClusterIndex] <= Program.CountforCluster_C_ktobezero;
            }

            if (!zerosizecluster) {
                this.FreezingMeasure_k_[RealClusterIndex] = this.FreezingMeasure_k_[RealClusterIndex] / this.C_k_[RealClusterIndex];
            }
        }

        DAVectorUtility.StopSubTimer(4);

    } // End SetClusterSizes


    //  Set the Correlation_k_i_j for Clusters in this Solution
    // Note NOT Divided by C_k_ for a true correlation as division not wanted in eigencalculation
    public final void SetClusterCorrelations() throws MPIException {
        DAVectorUtility.StartSubTimer(4);
        final int NumberofCorrelationPositions;
        if (Program.CalculateCorrelationMatrix) {
            NumberofCorrelationPositions = (Program.ParameterVectorDimension * (Program.ParameterVectorDimension + 1)) / 2;
        } else {
            NumberofCorrelationPositions = 0;
        }

        final GlobalReductions.FindIndirectVectorDoubleSum[] FindDC_k_DY_k_;
        final GlobalReductions.FindIndirectVectorDoubleSum[] FindCorrelations;
        final DistributedReductions.FindIndirectMultiVectorDoubleSum Find2Components;
        final int BeginDC_k_DY_k_;
        final int BeginCorrelations;

        if (this.DistributedExecutionMode) {
            Find2Components = new DistributedReductions.FindIndirectMultiVectorDoubleSum();
            BeginDC_k_DY_k_ = Find2Components.AddComponents(Program.ParameterVectorDimension);
            if (Program.CalculateCorrelationMatrix) {
                BeginCorrelations = Find2Components.AddComponents(NumberofCorrelationPositions);
            } else {
                BeginCorrelations = -1;
            }
            Find2Components.NodeInitialize();
            FindDC_k_DY_k_ = null;
            FindCorrelations = null;
        } else {
            FindDC_k_DY_k_ = new GlobalReductions.FindIndirectVectorDoubleSum[Program.ParameterVectorDimension];
            for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++) {
                FindDC_k_DY_k_[VectorIndex] = new GlobalReductions.FindIndirectVectorDoubleSum(DAVectorUtility.ThreadCount, ClusteringSolution.NumberLocalActiveClusters);
            }
            if (Program.CalculateCorrelationMatrix) {
                FindCorrelations = new GlobalReductions.FindIndirectVectorDoubleSum[NumberofCorrelationPositions];
                for (int CorrelationPositions = 0; CorrelationPositions < NumberofCorrelationPositions; CorrelationPositions++) {
                    FindCorrelations[CorrelationPositions] = new GlobalReductions.FindIndirectVectorDoubleSum(DAVectorUtility.ThreadCount, ClusteringSolution.NumberLocalActiveClusters);
                }
            } else {
                FindCorrelations = null;
            }
            Find2Components = null;
            BeginDC_k_DY_k_ = -1;
            BeginCorrelations = -1;
        }

        // Note - parallel for
        try {
            forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) -> {
                //	0.5 as cluster halved
                if (this.DistributedExecutionMode) {
                    Find2Components.ThreadInitialize(threadIndex);
                } else {
                    if (Program.CalculateCorrelationMatrix) {
                        for (int CorrelationPositions = 0; CorrelationPositions < NumberofCorrelationPositions; CorrelationPositions++) {
                            FindCorrelations[CorrelationPositions].startthread(threadIndex);
                        }
                    }
                    for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++) {
                        FindDC_k_DY_k_[VectorIndex].startthread(threadIndex);
                    }
                }
                double[] Sigma_Pointer;
                double[] Y_Pointer;
                int ThreadStorePosition = -1;
                int ArraySize = Math.min(ClusteringSolution.NumberAvailableActiveClusters,
                        ClusteringSolution.MaximumClustersperPoint);
                int[] ActiveClustersperPoint = new int[ArraySize];
                double[] DC_k_DY_k_Values = new double[Program.ParameterVectorDimension];
                double[][] DC_k_DY_k_Values_ClusterIndex = new double[Program.ParameterVectorDimension][];
                for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++) {
                    DC_k_DY_k_Values_ClusterIndex[VectorIndex] = new double[ArraySize];
                }
                int Maynotneed = Math.max(1, NumberofCorrelationPositions);
                double[][] WorkSpace_Sym_i_j_ClusterIndex = new double[Maynotneed][];
                double[] WorkSpace_Sym_i_j_ = new double[Maynotneed];
                if (Program.CalculateCorrelationMatrix) {
                    for (int CorrelationPositions = 0; CorrelationPositions < NumberofCorrelationPositions; CorrelationPositions++) {
                        WorkSpace_Sym_i_j_ClusterIndex[CorrelationPositions] = new double[ArraySize];
                    }
                }
                int indexlen = DAVectorUtility.PointsperThread[threadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[threadIndex] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++) {
                    int IndirectSize = this.NumClusters_alpha_[alpha];
                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++) {
                        int RealClusterIndex = -1;
                        int RemoteIndex = -1;
                        int ActiveClusterIndex = -1;
                        Box<Integer> tempRef_RealClusterIndex = new Box<>(
                                RealClusterIndex);
                        Box<Integer> tempRef_ActiveClusterIndex = new Box<>(
                                ActiveClusterIndex);
                        Box<Integer> tempRef_RemoteIndex = new Box<>(RemoteIndex);
                        VectorAnnealIterate.ClusterPointersforaPoint(alpha, IndirectClusterIndex, tempRef_RealClusterIndex,
                                tempRef_ActiveClusterIndex, tempRef_RemoteIndex);
                        RealClusterIndex = tempRef_RealClusterIndex.content;
                        ActiveClusterIndex = tempRef_ActiveClusterIndex.content;
                        RemoteIndex = tempRef_RemoteIndex.content;
                        ActiveClustersperPoint[IndirectClusterIndex] = ActiveClusterIndex;
                        if (RemoteIndex < 0) {
                            Y_Pointer = this.Y_k_i_[RealClusterIndex];
                            Sigma_Pointer = this.Sigma_k_i_[RealClusterIndex];
                            if (this.DistributedExecutionMode) {
                                ThreadStorePosition = ClusteringSolution.LocalThreadAccPosition[ActiveClusterIndex][threadIndex];
                            }
                        } else {
                            Y_Pointer = DistributedClusteringSolution.StorageforTransportedClusters.TotalTransportedY_t_i[RemoteIndex];
                            Sigma_Pointer = DistributedClusteringSolution.StorageforTransportedClusters.TotalTransportedSigma_t_i[RemoteIndex];
                            ThreadStorePosition = DistributedClusteringSolution.StorageforTransportedClusters.TransportedThreadAccPosition[RemoteIndex][threadIndex];
                        }
                        if ((RealClusterIndex == SpongeCluster) && (SpongeCluster >= 0)) {
                            for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++) {
                                DC_k_DY_k_Values_ClusterIndex[VectorIndex][IndirectClusterIndex] = 0.0;
                            }
                            if (Program.CalculateCorrelationMatrix) {
                                for (int CorrelationPositions = 0; CorrelationPositions < NumberofCorrelationPositions; CorrelationPositions++) {
                                    WorkSpace_Sym_i_j_ClusterIndex[CorrelationPositions][IndirectClusterIndex] = 0.0;
                                }
                            }
                            continue;
                        }
                        double CurrentMvalue = this.M_alpha_kpointer_[alpha][IndirectClusterIndex] * 0.5;
                        int countdoubleindex = 0;
                        for (int VectorIndex1 = 0; VectorIndex1 < Program.ParameterVectorDimension; VectorIndex1++) {
                            double fudge = 1.0 / (Math.sqrt(Sigma_Pointer[VectorIndex1]) * this.Temperature);
                            DC_k_DY_k_Values[VectorIndex1] = -(Y_Pointer[VectorIndex1] - Program.PointPosition[alpha][VectorIndex1]) * fudge * CurrentMvalue * (1.0 - CurrentMvalue);
                            DC_k_DY_k_Values_ClusterIndex[VectorIndex1][IndirectClusterIndex] = DC_k_DY_k_Values[VectorIndex1];
                            if (Program.CalculateCorrelationMatrix) {
                                double vecdiff1 = Program.PointPosition[alpha][VectorIndex1] - Y_Pointer[VectorIndex1];
                                for (int VectorIndex2 = VectorIndex1; VectorIndex2 < Program.ParameterVectorDimension; VectorIndex2++) {
                                    double vecdiff2 = Program.PointPosition[alpha][VectorIndex2] - Y_Pointer[VectorIndex2];
                                    WorkSpace_Sym_i_j_[countdoubleindex] = CurrentMvalue * vecdiff1 * vecdiff2;
                                    WorkSpace_Sym_i_j_ClusterIndex[countdoubleindex][IndirectClusterIndex] = WorkSpace_Sym_i_j_[countdoubleindex];
                                    ++countdoubleindex;
                                }
                            }
                        }
                        if (this.DistributedExecutionMode) {
                            if ((ThreadStorePosition < 0) || (ThreadStorePosition >=
                                    DistributedClusteringSolution.NodeAccMetaData.NumberofPointsperThread[threadIndex])) {
                                DAVectorUtility.printAndThrowRuntimeException(
                                        " Find Correlations Point " + (alpha + DAVectorUtility.PointStart_Process) +
                                                " Cluster Index " + IndirectClusterIndex + " Illegal Position " +
                                                ThreadStorePosition + " Max " +
                                                DistributedClusteringSolution.NodeAccMetaData
                                                        .NumberofPointsperThread[threadIndex] +
                                                " CreatedIndex " +
                                                ParallelClustering.runningSolution
                                                        .Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex] +
                                                " Node Total " +
                                                DistributedClusteringSolution.NodeAccMetaData.NumberofPointsperNode
                                );

                            }
                            if (Program.CalculateCorrelationMatrix) {
                                Find2Components.addapoint(threadIndex, ThreadStorePosition, BeginCorrelations,
                                        NumberofCorrelationPositions, WorkSpace_Sym_i_j_);
                            }
                            Find2Components.addapoint(threadIndex, ThreadStorePosition, BeginDC_k_DY_k_,
                                    Program.ParameterVectorDimension, DC_k_DY_k_Values);
                        }
                    }
                    if (!this.DistributedExecutionMode) {
                        if (Program.CalculateCorrelationMatrix) {
                            for (int CorrelationPositions = 0; CorrelationPositions < NumberofCorrelationPositions; CorrelationPositions++) {
                                FindCorrelations[CorrelationPositions].addapoint(threadIndex, IndirectSize,
                                        ActiveClustersperPoint, WorkSpace_Sym_i_j_ClusterIndex[CorrelationPositions]);
                            }
                        }
                        for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++) {
                            FindDC_k_DY_k_[VectorIndex].addapoint(threadIndex, IndirectSize, ActiveClustersperPoint,
                                    DC_k_DY_k_Values_ClusterIndex[VectorIndex]);
                        }
                    }
                }

            });
        } catch (SuspendableException e) {
            DAVectorUtility.printAndThrowRuntimeException(e.getMessage());
        }

        if (this.DistributedExecutionMode) {
            Find2Components.sumoverthreadsandmpi();
        } else {
            if (Program.CalculateCorrelationMatrix) {
                for (int CorrelationPositions = 0; CorrelationPositions < NumberofCorrelationPositions; CorrelationPositions++) {
                    FindCorrelations[CorrelationPositions].sumoverthreadsandmpi();
                }
            }
            for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++) {
                FindDC_k_DY_k_[VectorIndex].sumoverthreadsandmpi();
            }
        }


        int AccumulationPosition = -1;
        double result = 0.0;
        for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++) {
            int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
            if (RealClusterIndex == SpongeCluster) {
                continue;
            }
            if (this.DistributedExecutionMode) {
                AccumulationPosition = ClusteringSolution.LocalNodeAccPosition[LocalActiveClusterIndex];
            }

            int countdoubleindex = 0;
            for (int VectorIndex1 = 0; VectorIndex1 < Program.ParameterVectorDimension; VectorIndex1++) {
                if (this.DistributedExecutionMode) {
                    this.DC_k_DY_k_i_[RealClusterIndex][VectorIndex1] = Find2Components.TotalVectorSum[AccumulationPosition][BeginDC_k_DY_k_ + VectorIndex1];
                } else {
                    this.DC_k_DY_k_i_[RealClusterIndex][VectorIndex1] = FindDC_k_DY_k_[VectorIndex1].TotalVectorSum[LocalActiveClusterIndex];
                }
                if (Program.CalculateCorrelationMatrix) {
                    for (int VectorIndex2 = VectorIndex1; VectorIndex2 < Program.ParameterVectorDimension; VectorIndex2++) {
                        if (this.DistributedExecutionMode) {
                            result = Find2Components.TotalVectorSum[AccumulationPosition][BeginCorrelations + countdoubleindex];
                        } else {
                            result = FindCorrelations[countdoubleindex].TotalVectorSum[LocalActiveClusterIndex];
                        }
                        this.Correlation_k_i_j[RealClusterIndex][VectorIndex1][VectorIndex2] = result / (this.Temperature * Math.sqrt(this.Sigma_k_i_[RealClusterIndex][VectorIndex1] * this.Sigma_k_i_[RealClusterIndex][VectorIndex2]));
                        if (VectorIndex1 != VectorIndex2) {
                            this.Correlation_k_i_j[LocalActiveClusterIndex][VectorIndex2][VectorIndex1] = this.Correlation_k_i_j[LocalActiveClusterIndex][VectorIndex1][VectorIndex2];
                        }
                        ++countdoubleindex;
                    }
                }
            }
        }
        this.CorrelationsSet = true;
        DAVectorUtility.StopSubTimer(4);

    } // End SetClusterCorrelations()

    public final void FindOccupationCounts() throws MPIException {
        DAVectorUtility.StartSubTimer(4);
        final GlobalReductions.FindDoubleArraySum GetOccupationCounts;
        final DistributedReductions.FindIndirectMultiVectorDoubleSum FindOccupationCounts;
        final int BeginOccupationCounts;

        if (this.DistributedExecutionMode) {
            FindOccupationCounts = new DistributedReductions.FindIndirectMultiVectorDoubleSum();
            BeginOccupationCounts = FindOccupationCounts.AddComponents(1);
            FindOccupationCounts.NodeInitialize();
            GetOccupationCounts = null;
        } else {
            GetOccupationCounts = new GlobalReductions.FindDoubleArraySum(DAVectorUtility.ThreadCount, ClusteringSolution.NumberLocalActiveClusters);
            FindOccupationCounts = null;
            BeginOccupationCounts = -1;
        }

        //  Parallel Section setting cluster occupation counts
        // Note - parallel for
        try {
            forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) -> {
                if (this.DistributedExecutionMode) {
                    FindOccupationCounts.ThreadInitialize(threadIndex);
                } else {
                    GetOccupationCounts.startthread(threadIndex);
                }
                int ThreadStorePosition = -1;
                int indexlen = DAVectorUtility.PointsperThread[threadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[threadIndex] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++) {
                    double MaxMvalue = -1.0;
                    int NearestClusterIndirectIndex = -1;
                    int IndirectSize = this.NumClusters_alpha_[alpha];
                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++) {
                        LegalCluster(alpha, IndirectClusterIndex);
                        if (this.M_alpha_kpointer_[alpha][IndirectClusterIndex] > MaxMvalue) {
                            MaxMvalue = this.M_alpha_kpointer_[alpha][IndirectClusterIndex];
                            NearestClusterIndirectIndex = IndirectClusterIndex;
                        }
                    }
                    if (NearestClusterIndirectIndex < 0) {
                        String msg = "";
                        for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++) {
                            msg += " " + String.format("%1$5.4f", this.M_alpha_kpointer_[alpha][IndirectClusterIndex]);
                        }
                        DAVectorUtility.printAndThrowRuntimeException(
                                "Error due to no nearby Cluster for Point " + (alpha + DAVectorUtility.PointStart_Process) + " " + IndirectSize + " Process " + DAVectorUtility.MPI_Rank + msg);

                    }
                    int RealClusterIndex = -1;
                    int RemoteIndex = -1;
                    int ActiveClusterIndex = -1;
                    Box<Integer> tempRef_RealClusterIndex = new Box<>(
                            RealClusterIndex);
                    Box<Integer> tempRef_ActiveClusterIndex = new Box<>(
                            ActiveClusterIndex);
                    Box<Integer> tempRef_RemoteIndex = new Box<>(RemoteIndex);
                    VectorAnnealIterate.ClusterPointersforaPoint(alpha, NearestClusterIndirectIndex,
                            tempRef_RealClusterIndex, tempRef_ActiveClusterIndex, tempRef_RemoteIndex);
                    RealClusterIndex = tempRef_RealClusterIndex.content;
                    ActiveClusterIndex = tempRef_ActiveClusterIndex.content;
                    RemoteIndex = tempRef_RemoteIndex.content;
                    if (RemoteIndex < 0) {
                        if (this.DistributedExecutionMode) {
                            ThreadStorePosition = ClusteringSolution.LocalThreadAccPosition[ActiveClusterIndex][threadIndex];
                        }
                    } else {
                        ThreadStorePosition = DistributedClusteringSolution.StorageforTransportedClusters.TransportedThreadAccPosition[RemoteIndex][threadIndex];
                    }
                    if (this.DistributedExecutionMode) {
                        FindOccupationCounts.addapoint(threadIndex, ThreadStorePosition, 1.0);
                    } else {
                        GetOccupationCounts.addapoint(threadIndex, ActiveClusterIndex);
                    }
                }

            });
        } catch (SuspendableException e) {
            DAVectorUtility.printAndThrowRuntimeException(e.getMessage());
        }

        double Result = 0.0;
        if (this.DistributedExecutionMode) {
            FindOccupationCounts.sumoverthreadsandmpi();
        } else {
            GetOccupationCounts.sumoverthreadsandmpi();
        }

        for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++) {
            if (this.DistributedExecutionMode) {
                Result = FindOccupationCounts.TotalVectorSum[ClusteringSolution.LocalNodeAccPosition[LocalActiveClusterIndex]][0];
            } else {
                Result = GetOccupationCounts.TotalSum[LocalActiveClusterIndex];
            }
            this.OccupationCounts_k_[ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex]] = (int) Math.floor(Result + 0.5);
        }
        DAVectorUtility.StopSubTimer(4);

    } // End FindOccupationCounts()

    public static void SetGlobalClusterNumbers() throws MPIException { // Find Cluster Numbers starting at 0 in Host 0 and incrementing through hosts systematically
        // Set TotalClusterSummary

        if (TotalClusterSummary.IterationSetAt == ParallelClustering.runningSolution.IterationSetAt) {
            return;
        }
        TotalClusterSummary.IterationSetAt = ParallelClustering.runningSolution.IterationSetAt;
        ParallelClustering.runningSolution.FindOccupationCounts();

        int[] LocalClusterRealIndixes = new int[ClusteringSolution.NumberLocalActiveClusters];
        int CountLocal = 0;
        int CountGlobal = 0;
        for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++) {
            int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
            if (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] == 0) {
                TotalClusterSummary.SpongeCluster = CountGlobal;
            }
            if (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] == 2) {
                LocalClusterRealIndixes[CountLocal] = RealClusterIndex;
                ++CountLocal;
            } else {
                UniversalMapping[ParallelClustering.runningSolution.LocalCreatedIndex[RealClusterIndex]].GlobalClusterNumber = CountGlobal;
                TotalClusterSummary.CreatedIndex[CountGlobal] = ParallelClustering.runningSolution.LocalCreatedIndex[RealClusterIndex];
                System.arraycopy(ParallelClustering.runningSolution.Y_k_i_[RealClusterIndex], 0,
                        TotalClusterSummary.CenterPosition[CountGlobal], 0, Program.ParameterVectorDimension);
                TotalClusterSummary.CurrentNode[CountGlobal] = 0;
                TotalClusterSummary.OccupationCount[CountGlobal] = ParallelClustering.runningSolution.OccupationCounts_k_[RealClusterIndex];
                ++CountGlobal;
            }
        }
        int BaseCountGlobal = CountGlobal;
        DAVectorUtility.StartSubTimer(DAVectorUtility.MPIBROADCASTTiming);

        // Note - MPI Call - Broadcast - int
        if (DAVectorUtility.MPI_Size > 1){
//            DAVectorUtility.MPI_communicator.<Integer>Broadcast(CountGlobal, 0);
            CountGlobal = DAVectorUtility.mpiOps.broadcast(CountGlobal, 0);
        }

        DAVectorUtility.StopSubTimer(DAVectorUtility.MPIBROADCASTTiming);
        if (!ParallelClustering.runningSolution.DistributedExecutionMode) {
            TotalClusterSummary.NumberofCenters = CountGlobal;
            return;
        }

        int Maxsize = CountLocal;
        DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming2);
        // Note - MPI Call - Allreduce - int - max
        if (DAVectorUtility.MPI_Size > 1){
//            Maxsize = DAVectorUtility.MPI_communicator.<Integer>Allreduce(Maxsize, Operation < Integer >.Max);
            Maxsize = DAVectorUtility.mpiOps.allReduce(Maxsize, MPI.MAX);
        }
        DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming2);

        MPIPacket FromAfarInt = MPIPacket.newIntegerPacket(Maxsize);
        MPIPacket FromAfarDouble = MPIPacket.newDoublePacket(Maxsize);

        for (int SourceNode = 0; SourceNode < DAVectorUtility.MPI_Size; SourceNode++) {
            int SaveCountGlobal = CountGlobal;

            // Transport Created Index
            if (SourceNode == DAVectorUtility.MPI_Rank) {
                FromAfarInt.setNumberOfPoints(CountLocal);
                FromAfarInt.setFirstPoint(0);
                int numberOfPoints = FromAfarInt.getNumberOfPoints();
                for (int FromAfarIndex = 0; FromAfarIndex < numberOfPoints; FromAfarIndex++) {
                    FromAfarInt.setMArrayIntAt(FromAfarIndex,ParallelClustering.runningSolution.LocalCreatedIndex[LocalClusterRealIndixes[FromAfarIndex]]);
                }
            }
            DAVectorUtility.StartSubTimer(DAVectorUtility.MPIBROADCASTTiming);

            // Note - MPI Call - Broadcast - MPIPacket<Integer>
            if (DAVectorUtility.MPI_Size > 1){
//                DAVectorUtility.MPI_communicator.<MPIPacket<Integer>>Broadcast(FromAfarInt, SourceNode);
                FromAfarInt = DAVectorUtility.mpiOps.broadcast(FromAfarInt, SourceNode);
            }

            DAVectorUtility.StopSubTimer(DAVectorUtility.MPIBROADCASTTiming);
            int numberOfPoints = FromAfarInt.getNumberOfPoints();
            for (int FromAfarIndex = 0; FromAfarIndex < numberOfPoints; FromAfarIndex++) {
                int CreatedIndex = FromAfarInt.getMArrayIntAt(FromAfarIndex);
                if (ClusteringSolution.UniversalMapping[CreatedIndex] == null) {
                    ClusteringSolution.UniversalMapping[CreatedIndex] = new ClusterIndirection(ClusteringSolution.CurrentIteration, 0);
                }

                UniversalMapping[CreatedIndex].GlobalClusterNumber = CountGlobal;
                TotalClusterSummary.CreatedIndex[CountGlobal] = CreatedIndex;
                TotalClusterSummary.CurrentNode[CountGlobal] = SourceNode;
                ++CountGlobal;
            }

            //  Transport Occupation Counts
            if (SourceNode == DAVectorUtility.MPI_Rank) {
                FromAfarInt.setNumberOfPoints(CountLocal);
                FromAfarInt.setFirstPoint(0);
                for (int FromAfarIndex = 0; FromAfarIndex < CountLocal; FromAfarIndex++) {
                    FromAfarInt.setMArrayIntAt(FromAfarIndex,ParallelClustering.runningSolution.OccupationCounts_k_[LocalClusterRealIndixes[FromAfarIndex]]);
                }
            }
            DAVectorUtility.StartSubTimer(DAVectorUtility.MPIBROADCASTTiming);

            // Note - MPI Call - Broadcast - MPIPacket<Integer>
            if (DAVectorUtility.MPI_Size > 1){
//                DAVectorUtility.MPI_communicator.<MPIPacket<Integer>>Broadcast(FromAfarInt, SourceNode);
                FromAfarInt = DAVectorUtility.mpiOps.broadcast(FromAfarInt, SourceNode);
            }


            DAVectorUtility.StopSubTimer(DAVectorUtility.MPIBROADCASTTiming);
            numberOfPoints = FromAfarInt.getNumberOfPoints();
            for (int FromAfarIndex = 0; FromAfarIndex < numberOfPoints; FromAfarIndex++) {
                TotalClusterSummary.OccupationCount[SaveCountGlobal + FromAfarIndex] = FromAfarInt.getMArrayIntAt(FromAfarIndex);
            }

            //  Transport Center positions
            for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++) {
                if (SourceNode == DAVectorUtility.MPI_Rank) {
                    FromAfarDouble.setNumberOfPoints(CountLocal);
                    FromAfarDouble.setFirstPoint(0);
                    for (int FromAfarIndex = 0; FromAfarIndex < CountLocal; FromAfarIndex++) {
                        FromAfarDouble.setMArrayDoubleAt(FromAfarIndex,ParallelClustering.runningSolution.Y_k_i_[LocalClusterRealIndixes[FromAfarIndex]][VectorIndex]);
                    }
                }
                DAVectorUtility.StartSubTimer(DAVectorUtility.MPIBROADCASTTiming);

                // Note - MPI Call - Broadcast - MPIPacket<Double>
                if (DAVectorUtility.MPI_Size > 1){
//                    DAVectorUtility.MPI_communicator.<MPIPacket<Double>>Broadcast(FromAfarDouble, SourceNode);
                    FromAfarDouble = DAVectorUtility.mpiOps.broadcast(FromAfarDouble, SourceNode);
                }

                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIBROADCASTTiming);
                numberOfPoints = FromAfarDouble.getNumberOfPoints();
                for (int FromAfarIndex = 0; FromAfarIndex < numberOfPoints; FromAfarIndex++) {
                    TotalClusterSummary.CenterPosition[SaveCountGlobal + FromAfarIndex][VectorIndex] = FromAfarDouble.getMArrayDoubleAt(FromAfarIndex);
                }
            }
        }
        TotalClusterSummary.NumberofCenters = CountGlobal;

    } // End SetGlobalClusterNumbers()

    //    Return Pointer associate with given Active Cluster Index associated with given point
    //    Return -1 if Point not associated with this cluster
    public final int MapClusterToIndirect(int LocalPointPosition, int ActiveClusterIndex) {
        int CreatedIndex = ClusteringSolution.MapActivetoCreatedIndex(ActiveClusterIndex, this);

        int IndirectSize = this.NumClusters_alpha_[LocalPointPosition];
        for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++) {
            if (CreatedIndex == this.Map_alpha_PointertoCreatedIndex[LocalPointPosition][IndirectClusterIndex]) {
                return IndirectClusterIndex;
            }
        }
        return -1;

    } // End  MapClusterToIndirect(int LocalPointPosition)

    //   Set Clusters for a given point. Tries to preserve existing Malpha -- not always good idea
    public final void SetClustersforaPoint(int LocalPointPosition) {

        int IndirectClusterSize = this.NumClusters_alpha_[LocalPointPosition];

//          Following test removed as cases where clusters removed need Malpha reset even if Clusters perfect
        if ((!DistributedExecutionMode) && (IndirectClusterSize == this.Ncent_ThisNode) && (IndirectClusterSize <= TargetClustersperPoint)) {
            return;
        }

        int TargetNumber1 = Math.min(ClusteringSolution.NumberAvailableActiveClusters, TargetClustersperPoint);
        int TargetNumber = TargetNumber1; // Number of Target non-sponge clusters
        if (this.SpongeCluster >= 0) {
            TargetNumber--;
        }

        // Note multiple threads (in parallel loops) call this method leading to exceptions if StartSubTimer is called
        /*DAVectorUtility.StartSubTimer(13);*/
        int[] ProposedClusters = new int[TargetNumber];
        double[] ClusterDistances = new double[TargetNumber];

        int PositionofWorstCluster = -1;
        for (int listindex = 0; listindex < TargetNumber; listindex++) {
            ProposedClusters[listindex] = -1;
            ClusterDistances[listindex] = 0.0;
        }
        int NearestClusterIndex = -1; // Cluster Index of Nearest Cluster
        double NearestClusterDistance = 0.0;
        int ActiveSpongeIndex = -1;
        if (this.SpongeCluster >= 0) {
            ActiveSpongeIndex = ClusteringSolution.ActiveClusterIndices[this.SpongeCluster];
        }

        //  Note ActiveClusterIndex runs over local and remote clusters
        for (int ActiveClusterIndex = 0; ActiveClusterIndex < ClusteringSolution.NumberAvailableActiveClusters; ActiveClusterIndex++) {
            if (ActiveClusterIndex == ActiveSpongeIndex) {
                continue; // We will add back Sponge Cluster later
            }

            double CurrentClusterScaledDistance = DAVectorParallelism.getSquaredScaledDistancePointActiveCluster(LocalPointPosition, ActiveClusterIndex, this);
            if ((NearestClusterIndex < 0) || (CurrentClusterScaledDistance < NearestClusterDistance)) {
                NearestClusterIndex = ActiveClusterIndex;
                NearestClusterDistance = CurrentClusterScaledDistance;
            }
            Box<Integer> tempRef_PositionofWorstCluster = new Box<>(PositionofWorstCluster);
            GlobalReductions.FindMinimumSet(CurrentClusterScaledDistance, ActiveClusterIndex, tempRef_PositionofWorstCluster, ClusterDistances, ProposedClusters, TargetNumber);
            PositionofWorstCluster = tempRef_PositionofWorstCluster.content;
        }

        //  Make of clusters in order -- starting with sponge if exists
        int[] ClusterIndicesTobeUsed = new int[TargetNumber1];
        int NumberofClustersTobeUsed = 0;
        if (ActiveSpongeIndex >= 0) {
            ClusterIndicesTobeUsed[0] = ActiveSpongeIndex;
            NumberofClustersTobeUsed = 1;
        }

        //  Add in Non sponge clusters
        int Initial_ActualNumberofClusters = NumberofClustersTobeUsed;
        for (int LoopOverMinimumSet = 0; LoopOverMinimumSet < TargetNumber; LoopOverMinimumSet++) {
            int ActiveClusterIndex = ProposedClusters[LoopOverMinimumSet];
            if (ActiveClusterIndex < 0) {
                continue;
            }
            if ((ClusteringSolution.ExponentArgumentCut > 0.0) && ((ClusterDistances[LoopOverMinimumSet] - NearestClusterDistance) > (2.0 * ExponentArgumentCut * this.Temperature))) {
                continue;
            }
            ClusterIndicesTobeUsed[NumberofClustersTobeUsed] = ActiveClusterIndex;
            NumberofClustersTobeUsed++;
        }

        // Add Nearest Cluster if none (other than sponge) selected
        if (NumberofClustersTobeUsed == Initial_ActualNumberofClusters) {
            if ((NumberofClustersTobeUsed < 0) || (NumberofClustersTobeUsed >= TargetNumber1) || (NumberofClustersTobeUsed >= ClusterIndicesTobeUsed.length)) {
                DAVectorUtility.printAndThrowRuntimeException("Error " + NumberofClustersTobeUsed + " " + TargetNumber1 + " Old Number " + IndirectClusterSize + " " + Initial_ActualNumberofClusters + " " + TargetClustersperPoint + " " + ClusterIndicesTobeUsed.length + " " + ClusteringSolution.NumberAvailableActiveClusters + DistributedClusteringSolution.StorageforTransportedClusters.SizeOfTransportedArray);

            }
            ClusterIndicesTobeUsed[NumberofClustersTobeUsed] = NearestClusterIndex;
            ++NumberofClustersTobeUsed;
        }

        //  Now make list for storing back
        double[] Mvalues = new double[NumberofClustersTobeUsed];
        double Msum = 0.0;
        int NearestClusterIndirectIndex = -1;
        for (int IndirectClusterIndex = 0; IndirectClusterIndex < NumberofClustersTobeUsed; IndirectClusterIndex++) {
            Mvalues[IndirectClusterIndex] = 0.0;
            int ActiveClusterIndex = ClusterIndicesTobeUsed[IndirectClusterIndex];
            if (ActiveClusterIndex == NearestClusterIndex) {
                NearestClusterIndirectIndex = IndirectClusterIndex;
            }
            int OldIndirectClusterIndex = this.MapClusterToIndirect(LocalPointPosition, ActiveClusterIndex);
            if (OldIndirectClusterIndex < 0) {
                continue;
            }
            Mvalues[IndirectClusterIndex] = this.M_alpha_kpointer_[LocalPointPosition][OldIndirectClusterIndex];
            Msum += Mvalues[IndirectClusterIndex];
        }
        if (Msum < 0.2) {
            for (int IndirectClusterIndex = 0; IndirectClusterIndex < NumberofClustersTobeUsed; IndirectClusterIndex++) {
                this.M_alpha_kpointer_[LocalPointPosition][IndirectClusterIndex] = 1.0 / (double) NumberofClustersTobeUsed;
            }
            Msum = 1.0;
        }

        this.NumClusters_alpha_[LocalPointPosition] = NumberofClustersTobeUsed;
        if (NumberofClustersTobeUsed <= 0) {
            DAVectorUtility.printAndThrowRuntimeException("Error due to zero New Number of Clusters Point " + (LocalPointPosition + DAVectorUtility.PointStart_Process) + " Old Number " + IndirectClusterSize);

        }

        for (int IndirectClusterIndex = 0; IndirectClusterIndex < NumberofClustersTobeUsed; IndirectClusterIndex++) {
            this.M_alpha_kpointer_[LocalPointPosition][IndirectClusterIndex] = Mvalues[IndirectClusterIndex] / Msum;
            this.Map_alpha_PointertoCreatedIndex[LocalPointPosition][IndirectClusterIndex] = ClusteringSolution.MapActivetoCreatedIndex(ClusterIndicesTobeUsed[IndirectClusterIndex], this);
        }

        this.DiffMalpha_k_Set = -1;

        // Note multiple threads (in parallel loops) call this method leading to exceptions if StartSubTimer is called
        /*DAVectorUtility.StopSubTimer(13);*/
    } // End SetClustersforaPoint(int LocalPointPosition)

    public static void CopySolution(ClusteringSolution From, ClusteringSolution To) {
        To.DistributedExecutionMode = From.DistributedExecutionMode;
        To.SpongeCluster = From.SpongeCluster;
        To.PairwiseHammy = From.PairwiseHammy; // Value of Hamiltonian
        To.OldHammy = From.OldHammy; //  Previous value of Hamiltonian
        To.Ncent_ThisNode = From.Ncent_ThisNode; //the current number of clusters in node
        To.Ncent_Global = From.Ncent_Global; //the current total number of clusters
        To.ActualCoolingFactor = From.ActualCoolingFactor; // Actual Cooling Factor
        To.Temperature = From.Temperature; //the current temperature
        To.Eigenvectorset = From.Eigenvectorset;
        To.SolutionSet = From.SolutionSet;
        To.IterationSetAt = From.IterationSetAt;
        To.CorrelationsSet = From.CorrelationsSet;
        To.DiffMalpha_k_Set = From.DiffMalpha_k_Set;
        To.YPreviousSet = From.YPreviousSet;
        To.ClustertoSplit = From.ClustertoSplit;
        To.TotaloverVectorIndicesAverageWidth = From.TotaloverVectorIndicesAverageWidth;
        System.arraycopy(From.AverageWidth, 0, To.AverageWidth, 0, Program.ParameterVectorDimension);

        int NumberClusters = From.Ncent_ThisNode;
        for (int ClusterIndex = 0; ClusterIndex < NumberClusters; ClusterIndex++) {
            To.C_k_[ClusterIndex] = From.C_k_[ClusterIndex];
            To.ClusterScaledSquaredWidth_k_[ClusterIndex] = From.ClusterScaledSquaredWidth_k_[ClusterIndex];
            To.P_k_[ClusterIndex] = From.P_k_[ClusterIndex];
            To.FreezingMeasure_k_[ClusterIndex] = From.FreezingMeasure_k_[ClusterIndex];
            To.Splittable_k_[ClusterIndex] = From.Splittable_k_[ClusterIndex];
            To.SplitPriority_k_[ClusterIndex] = From.SplitPriority_k_[ClusterIndex];
            To.Eigenvalue_k[ClusterIndex] = From.Eigenvalue_k[ClusterIndex];
            To.DiffMsummed_k_[ClusterIndex] = From.DiffMsummed_k_[ClusterIndex];
            To.OccupationCounts_k_[ClusterIndex] = From.OccupationCounts_k_[ClusterIndex];
            To.LocalCreatedIndex[ClusterIndex] = From.LocalCreatedIndex[ClusterIndex];
            To.LocalStatus[ClusterIndex] = From.LocalStatus[ClusterIndex];
            To.LocalSplitCreatedIndex[ClusterIndex] = From.LocalSplitCreatedIndex[ClusterIndex];
            To.YPreviousActuallySet[ClusterIndex] = From.YPreviousActuallySet[ClusterIndex];

            for (int VectorIndex1 = 0; VectorIndex1 < Program.ParameterVectorDimension; VectorIndex1++) {
                To.Y_k_i_[ClusterIndex][VectorIndex1] = From.Y_k_i_[ClusterIndex][VectorIndex1];
                To.YPrevious_k_i_[ClusterIndex][VectorIndex1] = From.YPrevious_k_i_[ClusterIndex][VectorIndex1];
                To.Eigenvector_k_i[ClusterIndex][VectorIndex1] = From.Eigenvector_k_i[ClusterIndex][VectorIndex1];
                To.DC_k_DY_k_i_[ClusterIndex][VectorIndex1] = From.DC_k_DY_k_i_[ClusterIndex][VectorIndex1];
                To.Sigma_k_i_[ClusterIndex][VectorIndex1] = From.Sigma_k_i_[ClusterIndex][VectorIndex1];

                if (From.CorrelationsSet) {
                    System.arraycopy(To.Correlation_k_i_j[ClusterIndex][VectorIndex1], 0,
                                     From.Correlation_k_i_j[ClusterIndex][VectorIndex1], 0,
                                     Program.ParameterVectorDimension);
                }
            }
        }

        // Note - parallel for
        try {
            forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) -> {
                int indexlen = DAVectorUtility.PointsperThread[threadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[threadIndex] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++) {
                    int IndirectSize = From.NumClusters_alpha_[alpha];
                    To.NumClusters_alpha_[alpha] = IndirectSize;
                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++) {
                        To.M_alpha_kpointer_[alpha][IndirectClusterIndex] = From.M_alpha_kpointer_[alpha][IndirectClusterIndex];
                        To.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex] = From.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex];
                    }
                }

            });
        } catch (SuspendableException e) {
            DAVectorUtility.printAndThrowRuntimeException(e.getMessage());
        }
    } // End CopySolution

    public final void CleanupClusters() { // Must follow by SetActiveClusters()

        int NewClusterIndex = 0;
        this.ClustertoSplit = -1;
        int[] ClusterMapping = new int[this.Ncent_ThisNode];

        for (int OldClusterIndex = 0; OldClusterIndex < this.Ncent_ThisNode; OldClusterIndex++) {
            int status = this.LocalStatus[OldClusterIndex];
            if ((status < 0) || (status > 2)) {
                if (this.SpongeCluster == OldClusterIndex) {
                    DAVectorUtility.printAndThrowRuntimeException(" Attempt to Remove Sponge Cluster " + OldClusterIndex);

                }
                ClusterMapping[OldClusterIndex] = -1;
                continue;
            }

            ClusterMapping[OldClusterIndex] = NewClusterIndex;
            if (this.SpongeCluster == OldClusterIndex) {
                this.SpongeCluster = NewClusterIndex;
            }
            this.DiffMsummed_k_[NewClusterIndex] = this.DiffMsummed_k_[OldClusterIndex];
            this.C_k_[NewClusterIndex] = this.C_k_[OldClusterIndex];
            this.ClusterScaledSquaredWidth_k_[NewClusterIndex] = this.ClusterScaledSquaredWidth_k_[OldClusterIndex];
            this.P_k_[NewClusterIndex] = this.P_k_[OldClusterIndex];
            this.FreezingMeasure_k_[NewClusterIndex] = this.FreezingMeasure_k_[OldClusterIndex];
            this.Splittable_k_[NewClusterIndex] = this.Splittable_k_[OldClusterIndex];
            this.SplitPriority_k_[NewClusterIndex] = this.SplitPriority_k_[OldClusterIndex];
            this.Eigenvalue_k[NewClusterIndex] = this.Eigenvalue_k[OldClusterIndex];
            this.OccupationCounts_k_[NewClusterIndex] = this.OccupationCounts_k_[OldClusterIndex];
            this.LocalCreatedIndex[NewClusterIndex] = this.LocalCreatedIndex[OldClusterIndex];
            this.LocalStatus[NewClusterIndex] = this.LocalStatus[OldClusterIndex];
            this.LocalSplitCreatedIndex[NewClusterIndex] = this.LocalSplitCreatedIndex[OldClusterIndex];
            this.YPreviousActuallySet[NewClusterIndex] = this.YPreviousActuallySet[OldClusterIndex];

            for (int VectorIndex1 = 0; VectorIndex1 < Program.ParameterVectorDimension; VectorIndex1++) {
                this.Y_k_i_[NewClusterIndex][VectorIndex1] = this.Y_k_i_[OldClusterIndex][VectorIndex1];
                this.YPrevious_k_i_[NewClusterIndex][VectorIndex1] = this.YPrevious_k_i_[OldClusterIndex][VectorIndex1];
                this.Sigma_k_i_[NewClusterIndex][VectorIndex1] = this.Sigma_k_i_[OldClusterIndex][VectorIndex1];
                this.Eigenvector_k_i[NewClusterIndex][VectorIndex1] = this.Eigenvector_k_i[OldClusterIndex][VectorIndex1];
                this.DC_k_DY_k_i_[NewClusterIndex][VectorIndex1] = this.DC_k_DY_k_i_[OldClusterIndex][VectorIndex1];

                if (this.CorrelationsSet) {
                    System.arraycopy(this.Correlation_k_i_j[OldClusterIndex][VectorIndex1], 0,
                            this.Correlation_k_i_j[NewClusterIndex][VectorIndex1], 0,
                            Program.ParameterVectorDimension);
                }
            }
            int CreatedIndex = this.LocalCreatedIndex[OldClusterIndex];
            UniversalMapping[CreatedIndex].Availability = 1 + NewClusterIndex;
            ++NewClusterIndex;
        }
        this.Ncent_ThisNode = NewClusterIndex;

    } // End CleanupClusters()

    public final void RemoveCluster(int RemovedIndex) throws MPIException {
        //  Process Child and Parent Clusters before we change numbers
        if (this.DistributedExecutionMode) {
            this.LocalStatus[RemovedIndex] = -1;
            return;
        }
        this.DiffMalpha_k_Set = -1;

        // Reduce Number of Clusters
        //  Global Ncent set later in call to SetActiveClusters
        --this.Ncent_ThisNode;
        if (this.SpongeCluster >= 0) {
            if (this.SpongeCluster == RemovedIndex) {
                DAVectorUtility.printAndThrowRuntimeException(" Attempt to Remove Sponge Cluster " + RemovedIndex);

            }
            if (this.SpongeCluster > RemovedIndex) {
                --this.SpongeCluster;
            }
        }
        final int DeletedCreatedIndex = this.LocalCreatedIndex[RemovedIndex];
        UniversalMapping[DeletedCreatedIndex].Availability = 0;
        this.ClustertoSplit = -1;


        for (int ClusterIndex = RemovedIndex; ClusterIndex < this.Ncent_ThisNode; ClusterIndex++) {
            this.DiffMsummed_k_[ClusterIndex] = this.DiffMsummed_k_[ClusterIndex + 1];
            this.C_k_[ClusterIndex] = this.C_k_[ClusterIndex + 1];
            this.ClusterScaledSquaredWidth_k_[ClusterIndex] = this.ClusterScaledSquaredWidth_k_[ClusterIndex + 1];
            this.P_k_[ClusterIndex] = this.P_k_[ClusterIndex + 1];
            this.FreezingMeasure_k_[ClusterIndex] = this.FreezingMeasure_k_[ClusterIndex + 1];
            this.Splittable_k_[ClusterIndex] = this.Splittable_k_[ClusterIndex + 1];
            this.SplitPriority_k_[ClusterIndex] = this.SplitPriority_k_[ClusterIndex + 1];
            this.Eigenvalue_k[ClusterIndex] = this.Eigenvalue_k[ClusterIndex + 1];
            this.OccupationCounts_k_[ClusterIndex] = this.OccupationCounts_k_[ClusterIndex + 1];
            int CurrentCreatedIndex = this.LocalCreatedIndex[ClusterIndex + 1];
            this.LocalCreatedIndex[ClusterIndex] = CurrentCreatedIndex;
            UniversalMapping[CurrentCreatedIndex].Availability = 1 + ClusterIndex;
            this.LocalStatus[ClusterIndex] = this.LocalStatus[ClusterIndex + 1];
            this.LocalSplitCreatedIndex[ClusterIndex] = this.LocalSplitCreatedIndex[ClusterIndex + 1];
            this.YPreviousActuallySet[ClusterIndex] = this.YPreviousActuallySet[ClusterIndex + 1];

            for (int VectorIndex1 = 0; VectorIndex1 < Program.ParameterVectorDimension; VectorIndex1++) {
                this.Y_k_i_[ClusterIndex][VectorIndex1] = this.Y_k_i_[ClusterIndex + 1][VectorIndex1];
                this.YPrevious_k_i_[ClusterIndex][VectorIndex1] = this.YPrevious_k_i_[ClusterIndex + 1][VectorIndex1];
                this.Sigma_k_i_[ClusterIndex][VectorIndex1] = this.Sigma_k_i_[ClusterIndex + 1][VectorIndex1];
                this.Eigenvector_k_i[ClusterIndex][VectorIndex1] = this.Eigenvector_k_i[ClusterIndex + 1][VectorIndex1];
                this.DC_k_DY_k_i_[ClusterIndex][VectorIndex1] = this.DC_k_DY_k_i_[ClusterIndex + 1][VectorIndex1];

                if (this.CorrelationsSet) {
                    System.arraycopy(this.Correlation_k_i_j[ClusterIndex + 1][VectorIndex1], 0,
                            this.Correlation_k_i_j[ClusterIndex][VectorIndex1], 0,
                            Program.ParameterVectorDimension);
                }
            }
        }

        //  Record changes of count in ActiveCluster and Ncent_Global
        this.YPreviousActuallySet[this.Ncent_ThisNode] = false;
        ParallelClustering.runningSolution.SetActiveClusters();

        if (Program.UseTriangleInequality_DA > 0) {
            DATriangleInequality.DeleteCenter(RemovedIndex, DeletedCreatedIndex);
        }

        final GlobalReductions.FindDoubleArraySum OldClusterNumberHistogram = new GlobalReductions.FindDoubleArraySum(DAVectorUtility.ThreadCount, this.Ncent_ThisNode + 2);
        final GlobalReductions.FindDoubleArraySum NewClusterNumberHistogram = new GlobalReductions.FindDoubleArraySum(DAVectorUtility.ThreadCount, this.Ncent_ThisNode + 1);

        // Note - parallel for
        try {
            forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) -> {
                //  If Program.UseTriangleInequality_DA > 0, zero cluster count will be fixed in DistributedClusteringSolution.ManageMajorSynchronization(true) called later
                if (Program.RemovalDiagnosticPrint) {
                    OldClusterNumberHistogram.startthread(threadIndex);
                    NewClusterNumberHistogram.startthread(threadIndex);
                }
                int indexlen = DAVectorUtility.PointsperThread[threadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[threadIndex] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++) {
                    int IndirectSize = this.NumClusters_alpha_[alpha];
                    int decrement = 0;
                    double Msum = 0.0;
                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++) {
                        int CreatedIndex = this.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex];
                        if (CreatedIndex == DeletedCreatedIndex) {
                            if (decrement > 0) {
                                DAVectorUtility.printAndThrowRuntimeException(
                                        " Double Cluster Entry " + RemovedIndex + " Point " + (alpha + DAVectorUtility.PointStart_Process));

                            }
                            decrement = 1;
                            continue;
                        }
                        if (CreatedIndex <= 0) {
                            DAVectorUtility.printAndThrowRuntimeException(
                                    " Illegal Created Index in Point List " + CreatedIndex + " " + IndirectClusterIndex + " Point " + (alpha + DAVectorUtility.PointStart_Process));

                        }
                        int listmember = UniversalMapping[CreatedIndex].Availability - 1;
                        if (UniversalMapping[CreatedIndex].IterationSet != CurrentIteration) {
                            DAVectorUtility.printAndThrowRuntimeException(
                                    " Out of Date Cluster in Point List " + listmember + " Iterations " + UniversalMapping[CreatedIndex].IterationSet + " " + CurrentIteration + " Point " + (alpha + DAVectorUtility.PointStart_Process));

                        }
                        if (listmember >= this.Ncent_ThisNode) {
                            DAVectorUtility.printAndThrowRuntimeException(
                                    " Bad List Member " + listmember + " Position " + IndirectClusterIndex + " Decrement " + decrement + " Point " + (alpha + DAVectorUtility.PointStart_Process));

                        }
                        if (decrement == 1) {
                            this.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex - decrement] = CreatedIndex;
                            this.M_alpha_kpointer_[alpha][IndirectClusterIndex - decrement] = this.M_alpha_kpointer_[alpha][IndirectClusterIndex];
                        }
                        Msum += this.M_alpha_kpointer_[alpha][IndirectClusterIndex - decrement];
                    }
                    this.NumClusters_alpha_[alpha] = IndirectSize - decrement;
                    int NewIndirectSize = this.NumClusters_alpha_[alpha];
                    if ((NewIndirectSize <= 0) && (Program.UseTriangleInequality_DA == 0)) {
                        DAVectorUtility.printAndThrowRuntimeException(
                                " Zero Cluster Count after removing cluster for Point " + (alpha + DAVectorUtility.PointStart_Process) + " Originally " + IndirectSize);

                    }
                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < NewIndirectSize; IndirectClusterIndex++) {
                        int CreatedIndex = this.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex];
                        if (CreatedIndex <= 0) {
                            DAVectorUtility.printAndThrowRuntimeException(
                                    " Illegal Created Index in Corrected Point List " + CreatedIndex + " " + IndirectClusterIndex + " Point " + (alpha + DAVectorUtility.PointStart_Process));

                        }
                        int listmember = UniversalMapping[CreatedIndex].Availability - 1;
                        if ((listmember >= this.Ncent_ThisNode) || (listmember < 0)) {
                            DAVectorUtility.printAndThrowRuntimeException(
                                    " Bad Corrected List Member " + listmember + " Position " + IndirectClusterIndex + " Point " + (alpha + DAVectorUtility.PointStart_Process));

                        }
                    }
                    if (this.NumClusters_alpha_[alpha] < 1) {
                        if (Program.UseTriangleInequality_DA == 0) {
                            this.SetClustersforaPoint(alpha);
                        }
                    } else {
                        if (Msum < 0.2) {
                            for (int IndirectClusterIndex = 0; IndirectClusterIndex < NewIndirectSize; IndirectClusterIndex++) {
                                this.M_alpha_kpointer_[alpha][IndirectClusterIndex] = 1.0 / (double) NewIndirectSize;
                            }
                        } else {
                            for (int IndirectClusterIndex = 0; IndirectClusterIndex < NewIndirectSize; IndirectClusterIndex++) {
                                this.M_alpha_kpointer_[alpha][IndirectClusterIndex] = this.M_alpha_kpointer_[alpha][IndirectClusterIndex] / Msum;
                            }
                        }
                    }
                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < NewIndirectSize; IndirectClusterIndex++) {
                        ParallelClustering.runningSolution.LegalCluster(alpha, IndirectClusterIndex);
                    }
                    if (Program.RemovalDiagnosticPrint) {
                        OldClusterNumberHistogram.addapoint(threadIndex, IndirectSize);
                        NewClusterNumberHistogram.addapoint(threadIndex, NewIndirectSize);
                    }
                }

            });
        } catch (SuspendableException e) {
            DAVectorUtility.printAndThrowRuntimeException(e.getMessage());
        }

        if (Program.RemovalDiagnosticPrint) {
            OldClusterNumberHistogram.sumoverthreadsandmpi();
            NewClusterNumberHistogram.sumoverthreadsandmpi();

            String message = "\nOld Counts ";
            int ArraySize = 2 + this.Ncent_ThisNode;
            for (int ClusterIndex = 0; ClusterIndex < ArraySize; ClusterIndex++) {
                message += (new Double(OldClusterNumberHistogram.TotalSum[ClusterIndex])) + " ";
            }
            ArraySize--;
            message += "\nNew Counts ";
            for (int ClusterIndex = 0; ClusterIndex < ArraySize; ClusterIndex++) {
                message += (new Double(NewClusterNumberHistogram.TotalSum[ClusterIndex])) + " ";
            }
            DAVectorUtility.SALSAPrint(1, "Cluster " + RemovedIndex + " Removed" + message);
        }

    } // End RemoveCluster

    public final boolean LegalCluster(int alpha, int IndirectClusterIndex) {
        int IndirectSize = ParallelClustering.runningSolution.NumClusters_alpha_[alpha];
        int CreatedIndex = this.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex];
        int MappedClusterIndex = ClusteringSolution.UniversalMapping[CreatedIndex].Availability;
        int ClusterIterationNo = ClusteringSolution.UniversalMapping[CreatedIndex].IterationSet;
        if ((ClusterIterationNo < ClusteringSolution.CurrentIteration) || (MappedClusterIndex == 0)) {
            String errormessage = " Point " + (alpha + DAVectorUtility.PointStart_Process) + " Cluster Index " + IndirectClusterIndex + " Created Index " + CreatedIndex + " Actual Iteration " + ClusteringSolution.CurrentIteration + " Cluster Iteration " + ClusterIterationNo + " Mapped Index " + MappedClusterIndex + " Full Set Created Indices ";
            for (int errorloop = 0; errorloop < IndirectSize; errorloop++) {
                errormessage += ParallelClustering.runningSolution.Map_alpha_PointertoCreatedIndex[alpha][errorloop] + " ";
            }
            DAVectorUtility.printAndThrowRuntimeException(errormessage);

        }
        double Mvalue = this.M_alpha_kpointer_[alpha][IndirectClusterIndex];
        if (Double.isNaN(Mvalue) || Double.isInfinite(Mvalue)) {
            String message = "";
            for (int loop = 0; loop < IndirectSize; loop++) {
                message += this.Map_alpha_PointertoCreatedIndex[alpha][loop] + " Malpha " + String.format("%1$6.5f", ParallelClustering.runningSolution.M_alpha_kpointer_[alpha][loop]) + " * ";
            }
            DAVectorUtility.printAndThrowRuntimeException("Arithmetic Error Point " + (alpha + DAVectorUtility.PointStart_Process) + " Number of Clusters " + IndirectSize + " Indirect Index " + IndirectClusterIndex + " Created Index " + this.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex] + "\n" + message);

        }
        return true;

    } // End LegalCluster

} // End ClusteringSolution