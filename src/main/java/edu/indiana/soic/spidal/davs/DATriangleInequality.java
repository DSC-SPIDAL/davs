package edu.indiana.soic.spidal.davs;

import edu.indiana.soic.spidal.general.Box;
import mpi.MPIException;

import java.util.Arrays;
import java.util.function.IntToDoubleFunction;

import static edu.rice.hj.Module0.launchHabaneroApp;
import static edu.rice.hj.Module1.forallChunked;


public class DATriangleInequality
{
    //  Deterministic Annealing version of Kmeans Triangle Inequality Speed up package
    //  Note uses square of distance not distance and finds all clusters to associate with a point. Not just the best
    // i.e. all clusters within exponential cut of best
    //  Only used to set associated clusters for each point and NOT to estimate distances
    //  Associated Clusters only reset each Major Synchronization and so many (potentially) sub iterations in between
    //  Uses Active Indices for Clusters and has special actions for deleting and splitting clusters

    //  Result of this is list of centers in order of distance from point stored in LB_Current_alpha_ClusterPointer
    //  Need to extract those you want
    //  Utilities cope with deleted and split clusters

    public static IntToDoubleFunction GetRadiusperCenter;

    //  The following are pointers to structures outside DATriangleInequality

    public static double[][] PointPosition; // Distributed Point Positions

    //  The following point dependent structures are internal to DATriangleInequality
    public static int[] NearestCentertoPoint_alpha; // Nearest Cluster Center to each Point
    public static double[] Distance_NearestCentertoPoint_alpha; // Distance from Point to Nearest Cluster Center

    public static PointData[] LB_Old_alpha_ClusterPointer; //  "Old" Lower Bound data structure
    public static PointData[] LB_Last_alpha_ClusterPointer; //  "Last" Lower Bound data structure
    public static PointData[] LB_Current_alpha_ClusterPointer; //  "Current" Lower Bound data structure

    public static PointData[] LB_Buffer1_alpha_ClusterPointer; //  Buffer1 for "Current" and "Last" Lower Bound data structure
    public static PointData[] LB_Buffer2_alpha_ClusterPointer; //  Buffer2 for "Current" and "Last" Lower Bound data structure

    //  The following are pointers to structures outside DATriangleInequality and need not be changed for delete or split
    public static int[] CenterStatus; // < 0 deleted 0 Sponge 1 Global >1 Distributed
    public static int[] RealClusterIndices = null; // Real index values associated with Active Index
    public static ClusteringSolution Solution = null; // Clustering Solution needed to calculate Sigma
    public static double[][] CenterY_k_i_Current; // Defined and stored outside TriangleInequality. This is a Pointer
    public static double[][] CenterSigma_k_i_Current; // Standard Deviation of Cluster

    //  The following are local to DATriangleInequality and must be changed for delete or split
    public static double[][] CenterOld; // This is only used in TriangleInequality and is "old" list of centers
    public static double[][] CenterLast; // This is only used in TriangleInequality and is "last" list of centers
    public static int[] IncumbentCreatedIndices; // Created Indices of current centers

    public static CenterData[] CurrentInterCenter; // Center and Center-Center Data stored in all nodes
    public static double[][] DistributedCenterDifference; // Center-Center data stored for "old centers

    public static int ActiveSpongexIndex = -1; // Active Cluster Index for Sponge
    public static int RealSpongeIndex = -1; // Real Cluster Index for Sponge
    public static double SpongeFactor; // Sponge is SpongeFactor/2T

    public static int Ncent_TotalAvailable = 0; // Total Number of Clusters Available
    public static int Ncent_RealTotal = 0; // Number of Clusters stored here
    public static double ExponentialxTemperature = 0.0; // Upper Limit on Exponential Times Temperature

    public static boolean StartingOut = true; // True when starting
    public static int ParameterVectorDimension = 2; // Vector Dimension of Points
    public static int MaxNcent_Global = 50; // Maximum number of centers
    public static int Ncent_Global_Parallel = 50; // Use Parallelism if total number of centers greater than this
    public static int OldCenterOption = 0; // -1 Don't use, 0 Use with incremental update, >0 Refresh every OldCenterOption iterations

    public static int FakeCenterCount_Process = 0; // Total number of Centers summed over all threads calculated in this process
    public static int FakeCenterCount_Largest = 0; // Largest number of points in all processes
    public static int FakeCenterStart_Process = 0; //    First data point in this process

    //	Within a process, data points will be divided into ThreadCount segments( each thread a segment), the array keep the size of each segment
    public static int[] FullParallel_CentersperThread = null; //how many data points each thread will take care for thread+process parallelism
    public static int[] FullParallel_StartCenterperThread = null; //the starting point that a thread will take care for thread+process parallelism

    public static int[] LocalParallel_CentersperThread = null; //how many data points each thread will take care for thread+process parallelism
    public static int[] LocalParallel_StartCenterperThread = null; //the starting point that a thread will take care for thread+process parallelism

    public static boolean UseParallelismoverCenters = false; // If true major center loops run in parallel
    //  Even if centers calculated in parallel, they are stored globally per process

    //  Specify parameters for controlling use of triangle inequality
    public static int UseTriangleInequality = 0; // if 0 do NOT use trianglew inequality; > 0 use it in a way specified by integer
    public static double TriangleInequality_Delta1_old = 0.1; // Test for center change and old lower bounds (normalized by radius)
    public static double TriangleInequality_Delta1_current = 0.1; // Test for Center change and current lower bounds (normalized by radius)
    public static int MaxClusterLBsperPoint = 50; // Maximum number of Lower Bound values
    public static int MaxCentersperCenter = 49; // Maximum number of Centers in Center Difference Array

    public static GlobalReductions.FindVectorDoubleSum FindDiagnosticSums_Center; // Accumulate Diagnostic Sums in a parallel Environment
    public static int NumberDiagnosticSums_Center = 0; //   Number of Center Diagnostic Sums
    public static int SaveMPI_Size = 0; //  Save DAVectorUtility.MPI_Size
    public static GlobalReductions.FindVectorDoubleSum FindDiagnosticSums_Points; // Accumulate Diagnostic Sums in a parallel Environment
    public static int NumberDiagnosticSums_Points = 0; //   Number of Diagnostic Sums

    //  Specify parameters for monitoring use of triangle inequality
    public static double NumberFullDistancesCalculatedCC = 0.0; // Count Number of Actual Full Distance Computations for center-center cases
    public static double NumberFullDistancesCalculatedCC_off = 0.0; // Count Number of Actual Full Distance Computations for center-center cases with mixed update status
    public static double NumberFullDistancesCalculatedCC_CleanReject = 0.0; // Count Number of Actual Full Distance Computations for center-center cases cleanly rejected
    public static double NumberFullDistancesCalculatedCC_NotStored = 0.0; // Count Number of Actual Full Distance Computations for center-center cases Not Stored
    public static double NumberFullDistancesCalculatedPC = 0.0; // Count Number of Actual Full Distance Computations for point-center cases
    public static double CalculationCenterIterations = 0.0; // Count Center times Iterations
    public static double CountToobigAShift = 0.0; // Count Centers with large shift between iteration
    public static double AccumulateRadius = 0.0; // Accumulate Radius
    public static double AccumulateDmax = 0.0; // Accumulate Dmax
    public static double CountCCEntries = 0.0; // Accumulate Number of Center-Center distances
    public static double AccumulateCCvalues = 0.0; // Average value for Center-Center distances
    public static double CalculationPointIterations = 0.0; // Accumulate Product of Points and Iterations
    public static double AccumulateOLD_LB = 0.0; // OLD Lower Bound used successfully
    public static double AccumulateCURRENT_LB = 0.0; // CURRENT Lower Bound used successfully
    public static double AccumulateOLD_LBFail = 0.0; // OLD Lower Bound used unsuccessfully
    public static double AccumulateCURRENT_LBFail = 0.0; // CURRENT Lower Bound used unsuccessfully
    public static double CountLBEntries = 0.0; //  Number in Lower Bound Array
    public static double AccumulateBestScore = 0.0; // Best Score
    public static double AccumulateMAX_LB = 0.0; // Maximum Lower Bound
    public static double AccumulateBreakReason7 = 0.0; // Break reason 7 -- while ends
    public static double AccumulateBreakReason8 = 0.0; // Break reason 8 -- Center Scan Dmax
    public static double AccumulateBreakReason9 = 0.0; // Break reason 9 -- Center Scan Middle
    public static double AccumulateBreakReason10 = 0.0; // Break reason 10 -- Center Scan Start
    public static double NumberFullDistancesCalculatedPC_Old = 0.0; // Count Number of Actual Full Distance Computations for point-center cases to refresh Old Centers
    public static double NumberFullDistancesCalculatedPC_Last = 0.0; // Count Number of Actual Full Distance Computations for point-center cases to refresh Last Centers
    public static double NumberFullDistancesCalculatedLoop = 0.0; // Number of Full Distances calculated in main lower bound loop
    public static double PossibleNumberFullDistancesCalculated = 0.0; // Possible Number of Full Distances calculated in main lower bound loop
    public static double NumberCentersinBestSet = 0.0; // Number of centers examined in best loop Method = -2
    public static double NumberCentersinSimpleSet = 0.0; // Number of centers examined in Method = -1 loop
    public static double NumberCentersinCenterbasedset = 0.0; // Number of centers examined in method >=0 loop
    public static double AverageCentersLookedAt = 0.0; // Average number of centers NOT looked at
    public static double NumberCentersRemovedBackwards = 0.0; // Centers Removed in Backwards test
    public static double NumberBackwardTests = 0.0; // Number of Backwards test
    public static double NumberForwardTests = 0.0; // Number of Forwards test

    public static double RadiusTest1 = 0.0; // Number of Radius Tests
    public static double RadiusTest2 = 0.0; // Zero Radius
    public static double RadiusTest3 = 0.0; // Change < 0.05 Radius
    public static double RadiusTest4 = 0.0; // Change < 0.1 Radius
    public static double RadiusTest5 = 0.0; // Change < 0.2 Radius

    public static int CalculationIterations = 0; // Divide NumberFullDistancesCalculated to get Number per position and per iteration


    public static class PointData
    {
        public double[] DistceLowerBound; // Lower Bound on Point Center Distance
        public int[] CenterIndices; // Centers whose Point-Center distances stored in DistceLowerBound
        public int NumCenters; // Number of Centers in list
        public int AssociatedClusterIndex; // Center associated with this point
        public double NearestDistance; // Distance to center

        public PointData(int ArraySize)
        {
            DistceLowerBound = new double[ArraySize];
            CenterIndices = new int[ArraySize];
            NumCenters = 0;
            AssociatedClusterIndex = -1;
            NearestDistance = 0.0;
        }

    } // End PointData for TriangleInequality

    public static class CenterData
    { // Only gathered for current Centers i.e. Clusters
        // Center positions available for current and last iteration

        public static int NumberPackingInts = 2;
        public static int NumberPackingDoubles = 4;
        public int RecalculateStatus; // If 1 Recalculate all distances for this center at start of each x loop; if 2 force 1; if 0 do not recalculate
        public double ClusterRadius; // Radius of Cluster
        public int NearbyClusters; // Number of Nearby Clusters Stored
        public double Dmax; // Minimum of distances NOT in CennterDifference
        public double LastCenterChange; // Value of distance between current and last center position
        public double OldCenterChange; // Value of distance between current and old center position
        public double[] CenterDifference; // Distances between centers
        public int[] CenterDifferenceIndices; // Indices of Centers in CenterDifference

        public CenterData()
        {
            CenterDifference = new double[MaxCentersperCenter];
            CenterDifferenceIndices = new int[MaxCentersperCenter];
        } // End constructor CenterData()

        public final void serialize(int[] intarray, double[] doublearray, int intShift, int doubleShift)
        {
            intarray[intShift] = this.NearbyClusters;
            intarray[1 + intShift] = this.RecalculateStatus;
            doublearray[doubleShift] = this.ClusterRadius;
            doublearray[1 + doubleShift] = this.Dmax;
            doublearray[2 + doubleShift] = this.LastCenterChange;
            doublearray[3 + doubleShift] = this.OldCenterChange;

            for (int loopindex = 0; loopindex < this.NearbyClusters; loopindex++)
            {
                intarray[NumberPackingInts + loopindex + intShift] = this.CenterDifferenceIndices[loopindex];
                doublearray[NumberPackingDoubles + loopindex + doubleShift] = this.CenterDifference[loopindex];
            }

        } // End serialize

        public final void deserialize(int[] intarray, double[] doublearray, int intShift, int doubleShift)
        {
            this.NearbyClusters = intarray[intShift];
            this.RecalculateStatus = intarray[1 + intShift];
            this.ClusterRadius = doublearray[doubleShift];
            this.Dmax = doublearray[1 + doubleShift];
            this.LastCenterChange = doublearray[2 + doubleShift];
            this.OldCenterChange = doublearray[3 + doubleShift];

            for (int loopindex = 0; loopindex < this.NearbyClusters; loopindex++)
            {
                this.CenterDifferenceIndices[loopindex] = intarray[NumberPackingInts + loopindex + intShift];
                this.CenterDifference[loopindex] = doublearray[NumberPackingDoubles + loopindex + doubleShift];
            }

        } // End deserialize

    } // End CenterData for TriangleInequality

    public static void SetExternalFunctions(IntToDoubleFunction GetRadiusperCenter)
    {
        DATriangleInequality.GetRadiusperCenter = GetRadiusperCenter;

    } // End SetExternalFunctions

    public static void SetTriangleInequalityParameters(int UseTriangleInequalityINPUT, int MaxClusterLBsperPointINPUT, int MaxCentersperCenterINPUT, double TriangleInequality_Delta1_oldINPUT, double TriangleInequality_Delta1_currentINPUT, int OldCenterOptionINPUT)
    {
        DATriangleInequality.UseTriangleInequality = UseTriangleInequalityINPUT; // if 0 do NOT use triangle inequality; > 0 use it in a way specified by integer
        DATriangleInequality.MaxClusterLBsperPoint = MaxClusterLBsperPointINPUT; // Maximum number of Lower Bound values
        DATriangleInequality.MaxCentersperCenter = MaxCentersperCenterINPUT; // Maximum number of Centers in Center Difference Array
        DATriangleInequality.TriangleInequality_Delta1_old = TriangleInequality_Delta1_oldINPUT; // Test for center change and old lower bounds (normalized by radius)
        DATriangleInequality.TriangleInequality_Delta1_current = TriangleInequality_Delta1_currentINPUT; // Test for Center change and current lower bounds (normalized by radius)
        DATriangleInequality.OldCenterOption = OldCenterOptionINPUT;

    } // End SetTriangleInequalityParameters


    public static void InitializeTriangleInequality(double[][] PointPositionINPUT, double[][] CenterCurrentINPUT, double[][] CenterSigma_k_i_CurrentINPUT, int[] CenterStatusINPUT, int[] RealClusterIndicesINPUT, ClusteringSolution SolutionINPUT, int MaxNcent_GlobalINPUT, int Ncent_Global_ParallelINPUT, int ParameterVectorDimensionINPUT)
    {
        DATriangleInequality.PointPosition = PointPositionINPUT;
        DATriangleInequality.SaveMPI_Size = DAVectorUtility.MPI_Size;

        DATriangleInequality.CenterY_k_i_Current = CenterCurrentINPUT; // Array for current centers
        DATriangleInequality.CenterSigma_k_i_Current = CenterSigma_k_i_CurrentINPUT;
        DATriangleInequality.CenterStatus = CenterStatusINPUT; // Status of Arrays
        DATriangleInequality.RealClusterIndices = RealClusterIndicesINPUT; // Map Active to Real Indices
        DATriangleInequality.Solution = SolutionINPUT;

        DATriangleInequality.MaxNcent_Global = MaxNcent_GlobalINPUT;
        DATriangleInequality.Ncent_Global_Parallel = Ncent_Global_ParallelINPUT;
        DATriangleInequality.ParameterVectorDimension = ParameterVectorDimensionINPUT;

        //  Set Initial Parallelism based on DATriangleInequality.MaxNcent_Global NOT DATriangleInequality.Ncent_Global
        DATriangleInequality.SetParallelCenterDecomposition();

        //  Center Related Data
        DATriangleInequality.CenterOld = new double[DATriangleInequality.MaxNcent_Global][];
        DATriangleInequality.CenterLast = new double[DATriangleInequality.MaxNcent_Global][];
        DATriangleInequality.IncumbentCreatedIndices = new int[DATriangleInequality.MaxNcent_Global];
        DATriangleInequality.CurrentInterCenter = new CenterData[DATriangleInequality.MaxNcent_Global];
        DATriangleInequality.DistributedCenterDifference = new double[DATriangleInequality.FakeCenterCount_Process][];
        for (int LocalCenterIndex = 0; LocalCenterIndex < DATriangleInequality.FakeCenterCount_Process; LocalCenterIndex++)
        {
            DATriangleInequality.DistributedCenterDifference[LocalCenterIndex] = new double[DATriangleInequality.MaxNcent_Global];
        }

        if (DATriangleInequality.UseParallelismoverCenters)
        { // Centers Parallel over Threads NOT nodes
            // Note - parallel for
            launchHabaneroApp(() -> {
                forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) -> // End Sum over Threads
                {
                    int NumberCentersperThread = DATriangleInequality.LocalParallel_CentersperThread[threadIndex];
                    int BeginCenter = DATriangleInequality.LocalParallel_StartCenterperThread[threadIndex];
                    for (int CenterIndex = BeginCenter; CenterIndex < NumberCentersperThread + BeginCenter; CenterIndex++) {
                        DATriangleInequality.CenterOld[CenterIndex] = new double[DATriangleInequality.ParameterVectorDimension];
                        DATriangleInequality.CenterLast[CenterIndex] = new double[DATriangleInequality.ParameterVectorDimension];
                        DATriangleInequality.CurrentInterCenter[CenterIndex] = new CenterData();
                    }

                });
            });
        }
        else
        { // Centers Sequential
            for (int CenterIndex = 0; CenterIndex < MaxNcent_Global; CenterIndex++)
            {
                DATriangleInequality.CenterOld[CenterIndex] = new double[DATriangleInequality.ParameterVectorDimension];
                DATriangleInequality.CenterLast[CenterIndex] = new double[DATriangleInequality.ParameterVectorDimension];
                DATriangleInequality.CurrentInterCenter[CenterIndex] = new CenterData();
            }
        }

        // Point Related Data
        DATriangleInequality.NearestCentertoPoint_alpha = new int[DAVectorUtility.PointCount_Process];
        DATriangleInequality.Distance_NearestCentertoPoint_alpha = new double[DAVectorUtility.PointCount_Process];
        if (DATriangleInequality.OldCenterOption >= 0)
        {
            DATriangleInequality.LB_Old_alpha_ClusterPointer = new PointData[DAVectorUtility.PointCount_Process];
        }
        DATriangleInequality.LB_Buffer1_alpha_ClusterPointer = new PointData[DAVectorUtility.PointCount_Process];
        DATriangleInequality.LB_Buffer2_alpha_ClusterPointer = new PointData[DAVectorUtility.PointCount_Process];

        // Note - parallel for
        launchHabaneroApp(() -> {
            forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) -> // End Loop over Threads for points
            {
                int indexlen = DAVectorUtility.PointsperThread[threadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[threadIndex] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++) {
                    if (DATriangleInequality.OldCenterOption >= 0) {
                        DATriangleInequality.LB_Old_alpha_ClusterPointer[alpha] = new PointData(
                                DATriangleInequality.MaxClusterLBsperPoint);
                    }
                    DATriangleInequality.LB_Buffer1_alpha_ClusterPointer[alpha] = new PointData(
                            DATriangleInequality.MaxClusterLBsperPoint);
                    DATriangleInequality.LB_Buffer2_alpha_ClusterPointer[alpha] = new PointData(
                            DATriangleInequality.MaxClusterLBsperPoint);
                }

            });
        });

        DATriangleInequality.LB_Current_alpha_ClusterPointer = DATriangleInequality.LB_Buffer1_alpha_ClusterPointer;
        DATriangleInequality.LB_Last_alpha_ClusterPointer = null;

    } // End Initialize()

    public static void SetParallelCenterDecomposition()
    {
        DATriangleInequality.UseParallelismoverCenters = false;
        if (DATriangleInequality.MaxNcent_Global <= DATriangleInequality.Ncent_Global_Parallel)
        {
            DATriangleInequality.FakeCenterStart_Process = 0;
            DATriangleInequality.FakeCenterCount_Process = DATriangleInequality.MaxNcent_Global;
            DATriangleInequality.FakeCenterCount_Largest = DATriangleInequality.MaxNcent_Global;
            return;
        }
        DATriangleInequality.UseParallelismoverCenters = true;
        DAVectorUtility.SALSAPrint(0, "Center Parallelism used as well as Point Parallelism in Triangle Inequality");

        //  Disable Node Parallelism for Centers as this will be inherited from Program
        DATriangleInequality.FakeCenterStart_Process = 0;
        DATriangleInequality.FakeCenterCount_Process = DATriangleInequality.MaxNcent_Global;
        DATriangleInequality.FakeCenterCount_Largest = DATriangleInequality.MaxNcent_Global;

        //	Now divide centers among threads for this process assuming Maximum Number of Centers
        DATriangleInequality.SetParallelCenterThreadsDecomposition(DATriangleInequality.FakeCenterCount_Process);

    } // End SetParallelCenterDecomposition()

    public static void SetParallelCenterThreadsDecomposition(int NumbertoDivide)
    {
        //  Process only Thread Center Parallelism
        Range[] Local_ThreadRanges = RangePartitioner.Partition(NumbertoDivide, DAVectorUtility.ThreadCount);
        DATriangleInequality.LocalParallel_CentersperThread = new int[DAVectorUtility.ThreadCount];
        DATriangleInequality.LocalParallel_StartCenterperThread = new int[DAVectorUtility.ThreadCount];

        for (int j = 0; j < DAVectorUtility.ThreadCount; j++)
        {
            DATriangleInequality.LocalParallel_CentersperThread[j] = Local_ThreadRanges[j].getLength();
            DATriangleInequality.LocalParallel_StartCenterperThread[j] = Local_ThreadRanges[j].getStartIndex();
        }

        DATriangleInequality.FullParallel_CentersperThread = DATriangleInequality.LocalParallel_CentersperThread;
        DATriangleInequality.FullParallel_StartCenterperThread = DATriangleInequality.LocalParallel_StartCenterperThread;
    }

    public static void NextIteration() throws MPIException { //  Move to next iteration

        DATriangleInequality.SetUpLinktoDAClustering();

        //  Set Current Number of Clusters and Parallelism based on DATriangleInequality.Ncent_Global NOT DATriangleInequality.MaxNcent_Global
        DATriangleInequality.UseParallelismoverCenters = false;
        if (DATriangleInequality.Ncent_TotalAvailable > DATriangleInequality.Ncent_Global_Parallel)
        {
            DATriangleInequality.UseParallelismoverCenters = true;
            DATriangleInequality.SetParallelCenterThreadsDecomposition(DATriangleInequality.Ncent_TotalAvailable);
        } // Thread Parallelism over centers set

        if (DATriangleInequality.UseTriangleInequality != 0)
        {
            DATriangleInequality.LB_Last_alpha_ClusterPointer = DATriangleInequality.LB_Current_alpha_ClusterPointer;
            if (DATriangleInequality.LB_Last_alpha_ClusterPointer.equals(DATriangleInequality.LB_Buffer1_alpha_ClusterPointer))
            {
                DATriangleInequality.LB_Current_alpha_ClusterPointer = DATriangleInequality.LB_Buffer2_alpha_ClusterPointer;
            }
            else
            {
                DATriangleInequality.LB_Current_alpha_ClusterPointer = DATriangleInequality.LB_Buffer1_alpha_ClusterPointer;
            }
            if ((!DATriangleInequality.StartingOut) && (DATriangleInequality.LB_Last_alpha_ClusterPointer == null))
            { // Last not set properly
                DAVectorUtility.printAndThrowRuntimeException("Incorrect Last Lower bound " + " Iter " + DATriangleInequality.CalculationIterations);

            }

        }

        //  Set up Diagnostics
        DATriangleInequality.NumberDiagnosticSums_Center = 15;
        DATriangleInequality.FindDiagnosticSums_Center = new GlobalReductions.FindVectorDoubleSum(DAVectorUtility.ThreadCount, DATriangleInequality.NumberDiagnosticSums_Center);
        DATriangleInequality.NumberDiagnosticSums_Points = 24;
        DATriangleInequality.FindDiagnosticSums_Points = new GlobalReductions.FindVectorDoubleSum(DAVectorUtility.ThreadCount, DATriangleInequality.NumberDiagnosticSums_Points);

        // Analyze over centers
        DAVectorUtility.StartSubTimer(11);
        DATriangleInequality.CenterFacingAnalysis();
        DAVectorUtility.StopSubTimer(11);

        // Analyze over points
        DAVectorUtility.StartSubTimer(12);
        if (DATriangleInequality.StartingOut)
        {
            DATriangleInequality.InitialPointAnalysis();
        }
        else
        {
            DATriangleInequality.PointDependentAnalysis();
        }
        DAVectorUtility.StopSubTimer(12);

        //  Process Diagnostics

        DAVectorUtility.MPI_Size = 1;
        DATriangleInequality.FindDiagnosticSums_Center.sumoverthreadsandmpi();
        DAVectorUtility.MPI_Size = DATriangleInequality.SaveMPI_Size;
        DATriangleInequality.FindDiagnosticSums_Points.sumoverthreadsandmpi();

        DATriangleInequality.NumberFullDistancesCalculatedCC += DATriangleInequality.FindDiagnosticSums_Center.TotalVectorSum[0];
        DATriangleInequality.CalculationCenterIterations += DATriangleInequality.FindDiagnosticSums_Center.TotalVectorSum[1];
        DATriangleInequality.CountToobigAShift += DATriangleInequality.FindDiagnosticSums_Center.TotalVectorSum[2];
        DATriangleInequality.AccumulateRadius += DATriangleInequality.FindDiagnosticSums_Center.TotalVectorSum[3];
        DATriangleInequality.AccumulateDmax += DATriangleInequality.FindDiagnosticSums_Center.TotalVectorSum[4];
        DATriangleInequality.CountCCEntries += DATriangleInequality.FindDiagnosticSums_Center.TotalVectorSum[5];
        DATriangleInequality.AccumulateCCvalues += DATriangleInequality.FindDiagnosticSums_Center.TotalVectorSum[6];
        DATriangleInequality.NumberFullDistancesCalculatedCC_off += DATriangleInequality.FindDiagnosticSums_Center.TotalVectorSum[7];
        DATriangleInequality.NumberFullDistancesCalculatedCC_CleanReject += DATriangleInequality.FindDiagnosticSums_Center.TotalVectorSum[8];
        DATriangleInequality.NumberFullDistancesCalculatedCC_NotStored += DATriangleInequality.FindDiagnosticSums_Center.TotalVectorSum[9];

        DATriangleInequality.RadiusTest1 += DATriangleInequality.FindDiagnosticSums_Center.TotalVectorSum[10];
        DATriangleInequality.RadiusTest2 += DATriangleInequality.FindDiagnosticSums_Center.TotalVectorSum[11];
        DATriangleInequality.RadiusTest3 += DATriangleInequality.FindDiagnosticSums_Center.TotalVectorSum[12];
        DATriangleInequality.RadiusTest4 += DATriangleInequality.FindDiagnosticSums_Center.TotalVectorSum[13];
        DATriangleInequality.RadiusTest5 += DATriangleInequality.FindDiagnosticSums_Center.TotalVectorSum[14];

        DATriangleInequality.NumberFullDistancesCalculatedPC += DATriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[0]; // Total
        DATriangleInequality.CalculationPointIterations += DATriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[1]; // Accumulate Product of Points and Iterations
        DATriangleInequality.AccumulateOLD_LB += DATriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[2]; // OLD Lower Bound used successfully
        DATriangleInequality.AccumulateCURRENT_LB += DATriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[3]; // CURRENT Lower Bound used successfully
        DATriangleInequality.CountLBEntries += DATriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[4]; //  Number in Lower Bound Array
        DATriangleInequality.AccumulateBestScore += DATriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[5]; // Best Score
        DATriangleInequality.AccumulateMAX_LB += DATriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[6]; // Maximum Lower Bound
        DATriangleInequality.AccumulateBreakReason7 += DATriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[7]; // Break in LB while as all centers done
        DATriangleInequality.AccumulateBreakReason8 += DATriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[8]; // Break in LB while as Center Scan Dmax
        DATriangleInequality.AccumulateBreakReason9 += DATriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[9]; // Break in LB while as Center Scan middle
        DATriangleInequality.AccumulateBreakReason10 += DATriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[10]; // Break in LB while Center Scan start
        DATriangleInequality.NumberFullDistancesCalculatedPC_Old += DATriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[11]; // Total PC calculation in Old Lower Bound Update
        DATriangleInequality.NumberFullDistancesCalculatedPC_Last += DATriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[12]; // Total PC calculation in Last Lower Bound Update
        DATriangleInequality.AccumulateOLD_LBFail += DATriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[13]; // OLD Lower Bound used UNsuccessfully
        DATriangleInequality.AccumulateCURRENT_LBFail += DATriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[14]; // CURRENT Lower Bound used UNsuccessfully
        DATriangleInequality.NumberFullDistancesCalculatedLoop += DATriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[15]; // Number of Full Distances calculated in main lowwer bound loop
        DATriangleInequality.NumberCentersinBestSet += DATriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[16]; // Number of centers examined in best loop Method = -2
        DATriangleInequality.NumberCentersinSimpleSet += DATriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[17]; // Number of centers examined in Method = -1 loop
        DATriangleInequality.NumberCentersinCenterbasedset += DATriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[18]; // Number of centers examined in method >=0 loop
        DATriangleInequality.AverageCentersLookedAt += DATriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[19]; // Average number of centers NOT looked at
        DATriangleInequality.NumberCentersRemovedBackwards += DATriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[20]; // Centers Removed in Backwards test
        DATriangleInequality.NumberBackwardTests += DATriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[21]; // Number of backward tests
        DATriangleInequality.NumberForwardTests += DATriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[22]; // Number of forward tests
        DATriangleInequality.PossibleNumberFullDistancesCalculated += DATriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[23]; // Possible Number of Full Distances calculated in main lowwer bound loop

        //  End Iteration
        DATriangleInequality.EndIteration();

    } // End NextIteration()

    public static void SetUpLinktoDAClustering()
    { // Set up DA Clustering wrt to Triangle Equality

        DATriangleInequality.Ncent_TotalAvailable = ClusteringSolution.NumberAvailableActiveClusters;
        DATriangleInequality.ExponentialxTemperature = Program.ExpArgumentCut4 * DATriangleInequality.Solution.Temperature;
        DATriangleInequality.ActiveSpongexIndex = -1;
        DATriangleInequality.SpongeFactor = 0.0;
        if (DATriangleInequality.Ncent_RealTotal != DATriangleInequality.Solution.Ncent_ThisNode)
        {
            DAVectorUtility.printAndThrowRuntimeException("Illegal Number of Centers called " + DATriangleInequality.Solution.Ncent_ThisNode + " in DATriangleInequality " + DATriangleInequality.Ncent_RealTotal);

        }
        for (int ActiveIndices = 0; ActiveIndices < DATriangleInequality.Ncent_TotalAvailable; ActiveIndices++)
        {
            int RealClusterIndex = DATriangleInequality.RealClusterIndices[ActiveIndices];
            if (DATriangleInequality.CenterStatus[RealClusterIndex] == 0)
            {
                DATriangleInequality.ActiveSpongexIndex = ActiveIndices;
                DATriangleInequality.RealSpongeIndex = RealClusterIndex;
                DATriangleInequality.SpongeFactor = Program.SpongeFactor * Program.SpongeFactor;
            }
            if (DATriangleInequality.CenterStatus[RealClusterIndex] < 0)
            { // Deleted Clusters are not allowed
                DAVectorUtility.printAndThrowRuntimeException("Illegal deleted Cluster " + RealClusterIndex + " in DATriangleInequality");

            }
            if (DATriangleInequality.IncumbentCreatedIndices[RealClusterIndex] != DATriangleInequality.Solution.LocalCreatedIndex[RealClusterIndex])
            { // Created Indices MUST match
                DAVectorUtility.printAndThrowRuntimeException("Inconsistent Created Indices " + RealClusterIndex + " Actual " + DATriangleInequality.Solution.LocalCreatedIndex[RealClusterIndex] + " Local " + DATriangleInequality.IncumbentCreatedIndices[RealClusterIndex] + " in DATriangleInequality");

            }
        }

    } // End SetUpLinktoDAClustering()

    public static void AddCenter(int CreatedIndex)
    { // Needed at start when centers are created without splitting

        DATriangleInequality.IncumbentCreatedIndices[DATriangleInequality.Ncent_RealTotal] = CreatedIndex;
        ++DATriangleInequality.Ncent_RealTotal;
        ++DATriangleInequality.Ncent_TotalAvailable;

    } // End AddCenter(int CreatedIndex)

    public static void DeleteCenter(int RealIndextoDelete, int CreatedIndextoDelete)
    { //  Cope with Centers being deleted; delete just one

        //  Some structures are dealt with outside DATriangleInequality and so need no action
        //  DATriangleInequality.RealClusterIndices
        //  DATriangleInequality.CenterStatus
        //  DATriangleInequality.CenterY_k_i_Current
        //  DATriangleInequality.CenterSigma_k_i_Current
        //
        //  DATriangleInequality.Solution

        //  Check out any errors in input data
        if (RealIndextoDelete < 0)
        {
            return;
        }
        if (RealIndextoDelete == DATriangleInequality.RealSpongeIndex)
        { // Cannot delete Sponge
            DAVectorUtility.printAndThrowRuntimeException("Illegal deleted Cluster " + RealIndextoDelete + " in DATriangleInequality is Sponge");

        }
        if (RealIndextoDelete >= DATriangleInequality.Ncent_RealTotal)
        { // Index too large
            DAVectorUtility.printAndThrowRuntimeException("Illegal deleted Cluster " + RealIndextoDelete + " in DATriangleInequality is too big " + DATriangleInequality.Ncent_RealTotal);

        }
        if (CreatedIndextoDelete != DATriangleInequality.IncumbentCreatedIndices[RealIndextoDelete])
        { // Inconsistent Created Indices
            DAVectorUtility.printAndThrowRuntimeException("Illegal deleted Cluster " + RealIndextoDelete + " in DATriangleInequality has inconsistent Created Indices " + CreatedIndextoDelete + " " + DATriangleInequality.IncumbentCreatedIndices[RealIndextoDelete]);

        }

        //  ClusterMapping will tell us how to map kept centers
        final int[] ClusterMapping = new int[DATriangleInequality.Ncent_RealTotal];
        int NewRealIndex = 0;
        for (int RealClusterIndex = 0; RealClusterIndex < DATriangleInequality.Ncent_RealTotal; RealClusterIndex++)
        {
            if (RealClusterIndex == RealIndextoDelete)
            {
                ClusterMapping[RealClusterIndex] = -1;
                continue;
            }
            ClusterMapping[RealClusterIndex] = NewRealIndex;
            ++NewRealIndex;
        }

        --DATriangleInequality.Ncent_RealTotal;
        --DATriangleInequality.Ncent_TotalAvailable;
        if (DATriangleInequality.RealSpongeIndex >= 0)
        {
            if (DATriangleInequality.RealSpongeIndex > RealIndextoDelete)
            {
                --DATriangleInequality.RealSpongeIndex;
                --DATriangleInequality.ActiveSpongexIndex;
            }
        }

        //  Loop over Center indexed arrays
        for (int RealClusterIndex = RealIndextoDelete; RealClusterIndex < DATriangleInequality.Ncent_RealTotal; RealClusterIndex++)
        { // Just start at deleted center
            DATriangleInequality.IncumbentCreatedIndices[RealClusterIndex] = DATriangleInequality.IncumbentCreatedIndices[RealClusterIndex + 1];

            //  These are arrays of size ParameterVectorDimension
            //DATriangleInequality.CenterOld[RealClusterIndex] = DATriangleInequality.CenterOld[RealClusterIndex + 1];
            //DATriangleInequality.CenterLast[RealClusterIndex] = DATriangleInequality.CenterLast[RealClusterIndex + 1];
            for (int VectorIndex = 0; VectorIndex < DATriangleInequality.ParameterVectorDimension; VectorIndex++)
            {
                DATriangleInequality.CenterOld[RealClusterIndex][VectorIndex] = DATriangleInequality.CenterOld[RealClusterIndex + 1][VectorIndex];
                DATriangleInequality.CenterLast[RealClusterIndex][VectorIndex] = DATriangleInequality.CenterLast[RealClusterIndex + 1][VectorIndex];
            }
        }

        for (int RowRealClusterIndex = 0; RowRealClusterIndex <= DATriangleInequality.Ncent_RealTotal; RowRealClusterIndex++)
        { // Consider all centers
            int RowMapped = ClusterMapping[RowRealClusterIndex];
            if (RowMapped < 0)
            {
                continue;
            }

            for (int ColumnRealClusterIndex = 0; ColumnRealClusterIndex <= DATriangleInequality.Ncent_RealTotal; ColumnRealClusterIndex++)
            {
                int ColumnMapped = ClusterMapping[ColumnRealClusterIndex];
                if (ColumnMapped < 0)
                {
                    continue;
                }
                DATriangleInequality.DistributedCenterDifference[RowMapped][ColumnMapped] = DATriangleInequality.DistributedCenterDifference[RowRealClusterIndex][ColumnRealClusterIndex];
            }

            if (RowMapped != RowRealClusterIndex)
            {
                DATriangleInequality.CurrentInterCenter[RowMapped].RecalculateStatus = DATriangleInequality.CurrentInterCenter[RowRealClusterIndex].RecalculateStatus;
                DATriangleInequality.CurrentInterCenter[RowMapped].ClusterRadius = DATriangleInequality.CurrentInterCenter[RowRealClusterIndex].ClusterRadius;
                DATriangleInequality.CurrentInterCenter[RowMapped].Dmax = DATriangleInequality.CurrentInterCenter[RowRealClusterIndex].Dmax;
                DATriangleInequality.CurrentInterCenter[RowMapped].LastCenterChange = DATriangleInequality.CurrentInterCenter[RowRealClusterIndex].LastCenterChange;
                DATriangleInequality.CurrentInterCenter[RowMapped].OldCenterChange = DATriangleInequality.CurrentInterCenter[RowRealClusterIndex].OldCenterChange;
            }
            int DifferenceNumber = DATriangleInequality.CurrentInterCenter[RowRealClusterIndex].NearbyClusters;
            int NewLocalIndex = 0;
            for (int NearbyClusterLocalIndex = 0; NearbyClusterLocalIndex < DifferenceNumber; NearbyClusterLocalIndex++)
            {
                int OldLocalIndex = DATriangleInequality.CurrentInterCenter[RowRealClusterIndex].CenterDifferenceIndices[NearbyClusterLocalIndex];
                int MappedOldLocalIndex = ClusterMapping[OldLocalIndex];
                if (MappedOldLocalIndex < 0)
                {
                    continue;
                }
                DATriangleInequality.CurrentInterCenter[RowMapped].CenterDifferenceIndices[NewLocalIndex] = MappedOldLocalIndex;
                DATriangleInequality.CurrentInterCenter[RowMapped].CenterDifference[NewLocalIndex] = DATriangleInequality.CurrentInterCenter[RowRealClusterIndex].CenterDifference[NearbyClusterLocalIndex];
                ++NewLocalIndex;
            }
            DATriangleInequality.CurrentInterCenter[RowMapped].NearbyClusters = NewLocalIndex;

        } // End loop over Center Indexed arrays moving down one

        //  Point Dependent Lower Bounds
        // Note - parallel for
        launchHabaneroApp(() -> {
            forallChunked(0, DAVectorUtility.ThreadCount - 1,
                    (threadIndex) -> // End Threads over Points -   End Loop over Points
                    {
                        int indexlen = DAVectorUtility.PointsperThread[threadIndex];
                        int beginpoint = DAVectorUtility.StartPointperThread[threadIndex] - DAVectorUtility.PointStart_Process;
                        for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++) {
                            if (DATriangleInequality.OldCenterOption >= 0) {
                                DATriangleInequality.CleanupLowerBoundList(
                                        DATriangleInequality.LB_Old_alpha_ClusterPointer[alpha], ClusterMapping);
                            }
                            DATriangleInequality.CleanupLowerBoundList(
                                    DATriangleInequality.LB_Current_alpha_ClusterPointer[alpha], ClusterMapping);
                        }

                    }
            );
        });

    } // End DeleteCenter

    //  Clean up a point dependent lower bound list
    public static void CleanupLowerBoundList(DATriangleInequality.PointData LowerBoundStructure, int[] ClusterMapping)
    {
        int OldAssociatedCenter = LowerBoundStructure.AssociatedClusterIndex;
        if (OldAssociatedCenter >= 0)
        {
            LowerBoundStructure.AssociatedClusterIndex = ClusterMapping[OldAssociatedCenter];
            if (ClusterMapping[OldAssociatedCenter] < 0)
            {
                LowerBoundStructure.NearestDistance = -1.0;
            }
        }
        else
        {
            LowerBoundStructure.AssociatedClusterIndex = -1;
            LowerBoundStructure.NearestDistance = -1.0;
        }

        int OldNumLowerBounds = LowerBoundStructure.NumCenters;
        if (OldNumLowerBounds <= 0)
        {
            return;
        }
        int NewNumLowerBounds = 0;
        for (int LowerBoundLoop = 0; LowerBoundLoop < OldNumLowerBounds; LowerBoundLoop++)
        {
            int CenterIndex = LowerBoundStructure.CenterIndices[LowerBoundLoop];
            if (CenterIndex < 0)
            {
                continue;
            }
            int Mapped = ClusterMapping[CenterIndex];
            if (Mapped < 0)
            {
                continue;
            }
            LowerBoundStructure.CenterIndices[NewNumLowerBounds] = Mapped;
            LowerBoundStructure.DistceLowerBound[NewNumLowerBounds] = LowerBoundStructure.DistceLowerBound[LowerBoundLoop];
            ++NewNumLowerBounds;
        }
        LowerBoundStructure.NumCenters = NewNumLowerBounds;

    } // End CleanupLowerBoundList

    public static void SplitCenter(int RealIndextoSplit, int CreatedIndextoSplit, int NewRealIndex, int NewCreatedIndex)
    { // Add a new center due to a split

        //  Some structures are dealt with outside DATriangleInequality as pointers to outside structures
        //  DATriangleInequality.RealClusterIndices
        //  DATriangleInequality.CenterStatus
        //  DATriangleInequality.Solution
        //  DATriangleInequality.CenterY_k_i_Current
        //  DATriangleInequality.CenterSigma_k_i_Current

        if (RealIndextoSplit == DATriangleInequality.RealSpongeIndex)
        { // Cannot split Sponge
            DAVectorUtility.printAndThrowRuntimeException("Illegal Split Cluster " + RealIndextoSplit + " in DATriangleInequality is Sponge");

        }
        if (RealIndextoSplit >= DATriangleInequality.Ncent_RealTotal)
        { // Index too large
            DAVectorUtility.printAndThrowRuntimeException("Illegal Split Cluster " + RealIndextoSplit + " in DATriangleInequality is too big " + DATriangleInequality.Ncent_RealTotal);

        }
        if (CreatedIndextoSplit != DATriangleInequality.IncumbentCreatedIndices[RealIndextoSplit])
        { // Inconsistent Created Indices
            DAVectorUtility.printAndThrowRuntimeException("Illegal Split Cluster " + RealIndextoSplit + " in DATriangleInequality has inconsistent Created Indices " + CreatedIndextoSplit + " " + DATriangleInequality.IncumbentCreatedIndices[RealIndextoSplit]);

        }

        if (NewRealIndex != DATriangleInequality.Ncent_RealTotal)
        { // Inconsistent Total and new Split Index
            DAVectorUtility.printAndThrowRuntimeException("Illegal Split Cluster " + NewRealIndex + " in DATriangleInequality has inconsistent Totals " + DATriangleInequality.Ncent_RealTotal);

        }
        if (DATriangleInequality.Ncent_RealTotal >= DATriangleInequality.MaxNcent_Global)
        { // Inconsistent Total and new Split Index
            DAVectorUtility.printAndThrowRuntimeException("Illegal Split Cluster " + NewRealIndex + " in DATriangleInequality has inconsistent Totals " + DATriangleInequality.Ncent_RealTotal + " Max is " + DATriangleInequality.MaxNcent_Global);

        }

        ++DATriangleInequality.Ncent_RealTotal;
        ++DATriangleInequality.Ncent_TotalAvailable;

        //  Center indexed arrays
        //  No change needed for existing centers; just add a new center as a value copy
        for (int VectorIndex = 0; VectorIndex < DATriangleInequality.ParameterVectorDimension; VectorIndex++)
        {
            DATriangleInequality.CenterOld[NewRealIndex][VectorIndex] = DATriangleInequality.CenterOld[RealIndextoSplit][VectorIndex];
            DATriangleInequality.CenterLast[NewRealIndex][VectorIndex] = DATriangleInequality.CenterLast[RealIndextoSplit][VectorIndex];
        }
        //  Set incumbent created index
        DATriangleInequality.IncumbentCreatedIndices[NewRealIndex] = NewCreatedIndex;

        //  Add Final Column to DATriangleInequality.DistributedCenterDifference
        for (int RowRealClusterIndex = 0; RowRealClusterIndex < DATriangleInequality.Ncent_RealTotal - 1; RowRealClusterIndex++)
        {
            DATriangleInequality.DistributedCenterDifference[RowRealClusterIndex][NewRealIndex] = DATriangleInequality.DistributedCenterDifference[RowRealClusterIndex][RealIndextoSplit];
        }

        //  Add Final Row to DATriangleInequality.DistributedCenterDifference
        System.arraycopy(DATriangleInequality.DistributedCenterDifference[RealIndextoSplit], 0,
                         DATriangleInequality.DistributedCenterDifference[NewRealIndex], 0,
                         DATriangleInequality.Ncent_RealTotal);

        // Add a new row to DATriangleInequality.CurrentInterCenter
        DATriangleInequality.CurrentInterCenter[NewRealIndex].RecalculateStatus = 2; // Special Flag to force recalculation
        DATriangleInequality.CurrentInterCenter[NewRealIndex].ClusterRadius = DATriangleInequality.CurrentInterCenter[RealIndextoSplit].ClusterRadius;
        DATriangleInequality.CurrentInterCenter[NewRealIndex].Dmax = DATriangleInequality.CurrentInterCenter[RealIndextoSplit].Dmax;
        DATriangleInequality.CurrentInterCenter[NewRealIndex].LastCenterChange = DATriangleInequality.CurrentInterCenter[RealIndextoSplit].LastCenterChange;
        DATriangleInequality.CurrentInterCenter[NewRealIndex].OldCenterChange = DATriangleInequality.CurrentInterCenter[RealIndextoSplit].OldCenterChange;

        int DifferenceNumber = DATriangleInequality.CurrentInterCenter[RealIndextoSplit].NearbyClusters;
        System.arraycopy(DATriangleInequality.CurrentInterCenter[RealIndextoSplit].CenterDifference, 0,
                         DATriangleInequality.CurrentInterCenter[NewRealIndex].CenterDifference, 0, DifferenceNumber);
        DATriangleInequality.CurrentInterCenter[NewRealIndex].NearbyClusters = DifferenceNumber;

        DATriangleInequality.CurrentInterCenter[RealIndextoSplit].RecalculateStatus = 2; // Special Flag to force recalculation

        //  Point Dependent Lower Bounds
        //  These are not touched as new cluster will be added automatically
        //  OLD Point Dependent Lower Bounds however need to have cluster added
        if (DATriangleInequality.OldCenterOption >= 0)
        {
            // Note - parallel for
            launchHabaneroApp(() -> {
                forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) -> {
                    int indexlen = DAVectorUtility.PointsperThread[threadIndex];
                    int beginpoint = DAVectorUtility.StartPointperThread[threadIndex] - DAVectorUtility.PointStart_Process;
                    for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++) {
                        int NumOldbounds = DATriangleInequality.LB_Old_alpha_ClusterPointer[alpha].NumCenters;
                        double MaxVal = -1.0;
                        int MaxValIndex = -1;
                        int SplitIndex = -1;
                        for (int loopindex = 0; loopindex < NumOldbounds; loopindex++) {
                            int SaveIndex = DATriangleInequality.LB_Old_alpha_ClusterPointer[alpha].CenterIndices[loopindex];
                            if (SaveIndex == RealIndextoSplit) {
                                SplitIndex = loopindex;
                            }
                            if ((MaxVal < 0.0) || (DATriangleInequality.LB_Old_alpha_ClusterPointer[alpha].DistceLowerBound[loopindex] > MaxVal)) {
                                MaxValIndex = loopindex;
                                MaxVal = DATriangleInequality.LB_Old_alpha_ClusterPointer[alpha].DistceLowerBound[loopindex];
                            }
                        }
                        if (SplitIndex < 0) {
                            continue;
                        }
                        if (NumOldbounds < DATriangleInequality.MaxClusterLBsperPoint) {
                            DATriangleInequality.LB_Old_alpha_ClusterPointer[alpha].DistceLowerBound[NumOldbounds] = DATriangleInequality.LB_Old_alpha_ClusterPointer[alpha].DistceLowerBound[SplitIndex];
                            DATriangleInequality.LB_Old_alpha_ClusterPointer[alpha].CenterIndices[NumOldbounds] = NewRealIndex;
                            ++DATriangleInequality.LB_Old_alpha_ClusterPointer[alpha].NumCenters;
                        } else {
                            DATriangleInequality.LB_Old_alpha_ClusterPointer[alpha].DistceLowerBound[MaxValIndex] = DATriangleInequality.LB_Old_alpha_ClusterPointer[alpha].DistceLowerBound[SplitIndex];
                            DATriangleInequality.LB_Old_alpha_ClusterPointer[alpha].CenterIndices[MaxValIndex] = NewRealIndex;
                        }
                    }

                });
            });
        } // End skipping of Old centers
    } // End SplitCenter

    public static void EndIteration() throws MPIException {
        if (DATriangleInequality.UseParallelismoverCenters)
        { // Centers Parallel over Threads NOT nodes
            // Note - parallel for
            launchHabaneroApp(() -> {
                forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) -> // End Sum over Threads
                {
                    int NumberCentersperThread = DATriangleInequality.LocalParallel_CentersperThread[threadIndex];
                    int BeginCenter = DATriangleInequality.LocalParallel_StartCenterperThread[threadIndex];
                    for (int ActiveCenterIndex = BeginCenter; ActiveCenterIndex < NumberCentersperThread + BeginCenter; ActiveCenterIndex++) {
                        int RealCenterIndex = DATriangleInequality.RealClusterIndices[ActiveCenterIndex];
                        System.arraycopy(DATriangleInequality.CenterY_k_i_Current[RealCenterIndex], 0,
                                DATriangleInequality.CenterLast[RealCenterIndex], 0,
                                DATriangleInequality.ParameterVectorDimension);
                    }

                });
            });
        }
        else
        { // Centers Sequential
            for (int ActiveCenterIndex = 0; ActiveCenterIndex < DATriangleInequality.Ncent_TotalAvailable; ActiveCenterIndex++)
            {
                int RealCenterIndex = DATriangleInequality.RealClusterIndices[ActiveCenterIndex];
                System.arraycopy(DATriangleInequality.CenterY_k_i_Current[RealCenterIndex], 0,
                                 DATriangleInequality.CenterLast[RealCenterIndex], 0,
                                 DATriangleInequality.ParameterVectorDimension);
            }
        }
        DATriangleInequality.StartingOut = false;

        //  Iteration Count
        DATriangleInequality.CalculationIterations++;

        // Note - MPI Call - Barrier
        if (DAVectorUtility.MPI_Size > 1){
            DAVectorUtility.mpiOps.barrier(); // Make certain all processes have finished iteration
        }

        // Reset Cluster Centers done externally

    } // End EndIteration()

    public static void CenterFacingAnalysis() throws MPIException { // This does NOT get factor of 2 speed up from symmetry of distance computation

        //  First set O(Center) changes using thread not distributed parallelism
        // Distributed version could be used
        if (DATriangleInequality.UseParallelismoverCenters)
        { // Centers Parallel over Threads NOT nodes

            double[] AccumulateCC = new double[DAVectorUtility.ThreadCount];
            double[] AccumulateTooBig = new double[DAVectorUtility.ThreadCount];
            double[] AccumulateRadius = new double[DAVectorUtility.ThreadCount];
            for (int ThreadIndex = 0; ThreadIndex < DAVectorUtility.ThreadCount; ThreadIndex++)
            {
                AccumulateCC[ThreadIndex] = 0.0;
                AccumulateTooBig[ThreadIndex] = 0.0;
                AccumulateRadius[ThreadIndex] = 0.0;
            }

            // Note - parallel for
            launchHabaneroApp(() -> {
                forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) -> // End Sum over Threads
                {
                    int NumberCentersperThread = DATriangleInequality.LocalParallel_CentersperThread[threadIndex];
                    int BeginCenter = DATriangleInequality.LocalParallel_StartCenterperThread[threadIndex];
                    for (int ActiveCenterIndex = BeginCenter; ActiveCenterIndex < NumberCentersperThread + BeginCenter; ActiveCenterIndex++) {
                        int RealCenterIndex = DATriangleInequality.RealClusterIndices[ActiveCenterIndex];
                        Box<Double> tempRef_Object = new Box<>(
                                AccumulateCC[threadIndex]);
                        Box<Double> tempRef_Object2 = new Box<>(
                                AccumulateTooBig[threadIndex]);
                        Box<Double> tempRef_Object3 = new Box<>(
                                AccumulateRadius[threadIndex]);
                        DATriangleInequality.CFA_IndividualCenterCalculation(RealCenterIndex, tempRef_Object,
                                tempRef_Object2, tempRef_Object3);
                        AccumulateCC[threadIndex] = tempRef_Object.content;
                        AccumulateTooBig[threadIndex] = tempRef_Object2.content;
                        AccumulateRadius[threadIndex] = tempRef_Object3.content;
                    }

                });
            });

            for (int ThreadIndex = 0; ThreadIndex < DAVectorUtility.ThreadCount; ThreadIndex++)
            {
                DATriangleInequality.NumberFullDistancesCalculatedCC += AccumulateCC[ThreadIndex];
                DATriangleInequality.CountToobigAShift += AccumulateTooBig[ThreadIndex];
                DATriangleInequality.AccumulateRadius += AccumulateRadius[ThreadIndex];
            }

        } // End Center Parallel over threads

        else
        { // Centers Sequential
            for (int ActiveCenterIndex = 0; ActiveCenterIndex < Ncent_TotalAvailable; ActiveCenterIndex++)
            {
                int RealCenterIndex = DATriangleInequality.RealClusterIndices[ActiveCenterIndex];
                Box<Double> tempRef_NumberFullDistancesCalculatedCC = new Box<>(DATriangleInequality.NumberFullDistancesCalculatedCC);
                Box<Double> tempRef_CountToobigAShift = new Box<>(DATriangleInequality.CountToobigAShift);
                Box<Double> tempRef_AccumulateRadius = new Box<>(DATriangleInequality.AccumulateRadius);
                DATriangleInequality.CFA_IndividualCenterCalculation(RealCenterIndex, tempRef_NumberFullDistancesCalculatedCC, tempRef_CountToobigAShift, tempRef_AccumulateRadius);
                DATriangleInequality.NumberFullDistancesCalculatedCC = tempRef_NumberFullDistancesCalculatedCC.content;
                DATriangleInequality.CountToobigAShift = tempRef_CountToobigAShift.content;
                DATriangleInequality.AccumulateRadius = tempRef_AccumulateRadius.content;
            }

        } // End Sequential center analysis doing things O(Number of Centers)


        //  Reconcile TriangleInequality.CountToobigAShift and TriangleInequality.CurrentInterCenter[CenterIndex].RecalculateStatus if MPI Parallel
		if (DAVectorUtility.MPI_Size > 1)
		{
			DAVectorUtility.StartSubTimer(DAVectorUtility.MPIBROADCASTTiming);
            // Note - MPI call - Broadcast - double
			DATriangleInequality.CountToobigAShift = DAVectorUtility.mpiOps.broadcast(DATriangleInequality.CountToobigAShift,0);

			int[] TransmitRecalculateStatus = new int[DATriangleInequality.Ncent_RealTotal];
			for (int ActiveCenterIndex = 0; ActiveCenterIndex < DATriangleInequality.Ncent_TotalAvailable; ActiveCenterIndex++)
			{
				int RealCenterIndex = DATriangleInequality.RealClusterIndices[ActiveCenterIndex];
				TransmitRecalculateStatus[RealCenterIndex] = DATriangleInequality.CurrentInterCenter[RealCenterIndex].RecalculateStatus;
			}
            // Note - MPI call - Broadcast - int[]
			DAVectorUtility.mpiOps.broadcast(TransmitRecalculateStatus,0);
			for (int ActiveCenterIndex = 0; ActiveCenterIndex < DATriangleInequality.Ncent_TotalAvailable; ActiveCenterIndex++)
			{
				int RealCenterIndex = DATriangleInequality.RealClusterIndices[ActiveCenterIndex];
				DATriangleInequality.CurrentInterCenter[RealCenterIndex].RecalculateStatus = TransmitRecalculateStatus[RealCenterIndex];
			}
			DAVectorUtility.StopSubTimer(DAVectorUtility.MPIBROADCASTTiming);
		}

        //  Now do full analysis with O(Number of Centers ^2)
        if (DATriangleInequality.UseParallelismoverCenters)
        { // Centers Parallel over Threads NOT Nodes

            //  Parallel Center Processing
            // Note - parallel for
            launchHabaneroApp(() -> {
                forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) -> // End Loop over Threads
                {
                    int NumberCentersToProcess = DATriangleInequality.FullParallel_CentersperThread[threadIndex];
                    int BeginCenterThread = DATriangleInequality.FullParallel_StartCenterperThread[threadIndex];
                    for (int ActiveCenterIndex = BeginCenterThread; ActiveCenterIndex < NumberCentersToProcess + BeginCenterThread; ActiveCenterIndex++) {
                        int RealCenterIndex = DATriangleInequality.RealClusterIndices[ActiveCenterIndex];
                        DATriangleInequality.CFA_CenterCenterCalculation(threadIndex,
                                DATriangleInequality.FakeCenterStart_Process, RealCenterIndex);
                    }

                });
            });
        } // End Parallel Centers analysis of O(Number of Centers ^2) calculation

        else
        { // Centers Sequential

            for (int ActiveCenterIndex = 0; ActiveCenterIndex < DATriangleInequality.Ncent_TotalAvailable; ActiveCenterIndex++)
            {
                DATriangleInequality.CFA_CenterCenterCalculation(0, 0, DATriangleInequality.RealClusterIndices[ActiveCenterIndex]);
            }

        } // End Sequential Centers analysis of O(Number of Centers ^2) calculation

        //  Loop over all centers addressing TriangleInequality.CurrentInterCenter[CenterIndex].RecalculateStatus == 1
        if (DATriangleInequality.StartingOut)
        {
            return;
        }
        for (int ActiveCenterIndex = 0; ActiveCenterIndex < Ncent_TotalAvailable; ActiveCenterIndex++)
        {
            int RealCenterIndex = DATriangleInequality.RealClusterIndices[ActiveCenterIndex];
            if (DATriangleInequality.CurrentInterCenter[RealCenterIndex].RecalculateStatus == 0)
            {
                continue;
            }
            DATriangleInequality.CurrentInterCenter[RealCenterIndex].LastCenterChange = 0.0;
            DATriangleInequality.CurrentInterCenter[RealCenterIndex].OldCenterChange = 0.0;
            for (int VectorIndex = 0; VectorIndex < DATriangleInequality.ParameterVectorDimension; VectorIndex++)
            {
                DATriangleInequality.CenterOld[RealCenterIndex][VectorIndex] = DATriangleInequality.CenterY_k_i_Current[RealCenterIndex][VectorIndex];
                DATriangleInequality.CenterLast[RealCenterIndex][VectorIndex] = DATriangleInequality.CenterY_k_i_Current[RealCenterIndex][VectorIndex];
            }
        }

    } // End CenterFacingAnalysis()

    //  Initial O(Number of Centers) Computation
    public static void CFA_IndividualCenterCalculation(int RealCenterIndex, Box<Double> ThreadCCAccumulate, Box<Double> ThreadTooBigAccumulate, Box<Double> ThreadRadiusAccumulate)
    {
        double CurrentLastDifference = 0.0;
        double CurrentOldDifference = 0.0;
        double Radius = 0.0;
        boolean TooBigaShift = false;

        if (RealCenterIndex != DATriangleInequality.RealSpongeIndex)
        { // Skip Center Computation if Sponge
            if (DATriangleInequality.StartingOut)
            { // Set Old Centers
                for (int VectorIndex = 0; VectorIndex < DATriangleInequality.ParameterVectorDimension; VectorIndex++)
                {
                    DATriangleInequality.CenterOld[RealCenterIndex][VectorIndex] = DATriangleInequality.CenterY_k_i_Current[RealCenterIndex][VectorIndex];
                    DATriangleInequality.CenterLast[RealCenterIndex][VectorIndex] = DATriangleInequality.CenterY_k_i_Current[RealCenterIndex][VectorIndex];
                }
            }

            else
            {
                Radius = DATriangleInequality.GetRadiusperCenter.applyAsDouble(RealCenterIndex);
                CurrentLastDifference = DAVectorParallelism.getNOTSquaredScaledDistancebetweenVectors(DATriangleInequality.CenterY_k_i_Current[RealCenterIndex], DATriangleInequality.CenterLast[RealCenterIndex], DATriangleInequality.CenterSigma_k_i_Current[RealCenterIndex]);
                CurrentOldDifference = DAVectorParallelism.getNOTSquaredScaledDistancebetweenVectors(DATriangleInequality.CenterY_k_i_Current[RealCenterIndex], DATriangleInequality.CenterOld[RealCenterIndex], DATriangleInequality.CenterSigma_k_i_Current[RealCenterIndex]);
                ThreadCCAccumulate.content += 2.0;

                if ((CurrentLastDifference > DATriangleInequality.TriangleInequality_Delta1_current * Radius) || (CurrentOldDifference > DATriangleInequality.TriangleInequality_Delta1_old * Radius))
                {
                    TooBigaShift = true;
                    ThreadTooBigAccumulate.content += 1.0;
                }
            }
        }

        //  =2 Means affected by Split Update so force recalculate
        if (DATriangleInequality.CurrentInterCenter[RealCenterIndex].RecalculateStatus == 2)
        {
            TooBigaShift = true;
        }
        else
        {
            DATriangleInequality.CurrentInterCenter[RealCenterIndex].RecalculateStatus = 0;
        }

        if (TooBigaShift)
        {
            DATriangleInequality.CurrentInterCenter[RealCenterIndex].RecalculateStatus = 1;
        }

        DATriangleInequality.CurrentInterCenter[RealCenterIndex].ClusterRadius = Radius;
        ThreadRadiusAccumulate.content += Radius;
        DATriangleInequality.CurrentInterCenter[RealCenterIndex].LastCenterChange = CurrentLastDifference;
        DATriangleInequality.CurrentInterCenter[RealCenterIndex].OldCenterChange = CurrentOldDifference;

    } // End CFA_IndividualCenterCalculation

    public static void CFA_CenterCenterCalculation(int ThreadIndex, int CenterStart, int RealCenterIndex)
    { // Do calculation for one center in O(Number of Centers ^2)  part of CenterFacingAnalysis()

        if (RealCenterIndex == DATriangleInequality.RealSpongeIndex)
        {
            DATriangleInequality.CurrentInterCenter[RealCenterIndex].NearbyClusters = 0;
            DATriangleInequality.CurrentInterCenter[RealCenterIndex].Dmax = 0.0;
            return;
        }
        DATriangleInequality.FindDiagnosticSums_Center.addapoint(ThreadIndex, 1.0, 1);
        boolean MainChanged = false;
        if (DATriangleInequality.CurrentInterCenter[RealCenterIndex].RecalculateStatus > 0)
        {
            MainChanged = true;
        }
        double CurrentOldDifference = DATriangleInequality.CurrentInterCenter[RealCenterIndex].OldCenterChange; // Change Old-Current for current center
        double CurrentLastDifference = DATriangleInequality.CurrentInterCenter[RealCenterIndex].LastCenterChange; // Change Last-Current for current cente

        if (!StartingOut)
        {
            DATriangleInequality.FindDiagnosticSums_Center.addapoint(ThreadIndex, 1.0, 10);
            double Radius = DATriangleInequality.CurrentInterCenter[RealCenterIndex].ClusterRadius;
            if (Radius < 0.000001)
            {
                DATriangleInequality.FindDiagnosticSums_Center.addapoint(ThreadIndex, 1.0, 11);
            }
            else
            {
                if (CurrentLastDifference < 0.05 * Radius)
                {
                    DATriangleInequality.FindDiagnosticSums_Center.addapoint(ThreadIndex, 1.0, 12);
                }
                if (CurrentLastDifference < 0.1 * Radius)
                {
                    DATriangleInequality.FindDiagnosticSums_Center.addapoint(ThreadIndex, 1.0, 13);
                }
                if (CurrentLastDifference < 0.2 * Radius)
                {
                    DATriangleInequality.FindDiagnosticSums_Center.addapoint(ThreadIndex, 1.0, 14);
                }
            }
        }

        //  Set up ordered matrix of inter-center distances
        int currentcut = -1;
        double Minrejected = -1.0;
        int ActualNumber = 0;
        int NumCentersinOrderedList = Math.min(DATriangleInequality.Ncent_RealTotal - 1, DATriangleInequality.MaxCentersperCenter);
        Smaller [] smallers = new Smaller[NumCentersinOrderedList];

        //  Prepare Dynamic Looping
        boolean[] LookedAt = new boolean[DATriangleInequality.Ncent_RealTotal];
        for (int CenterIndex1 = 0; CenterIndex1 < DATriangleInequality.Ncent_RealTotal; CenterIndex1++)
        {
            LookedAt[CenterIndex1] = false;
        }
        int NumCentersToGo = DATriangleInequality.Ncent_RealTotal; // Number of Centers still to look at -- end when <= 0

        int CurrentCandidateCenter = -1; // Current Center being looked at in Method = -1 approach
        int Method = -2; // Method = -2 Start with Previous list, Method = -1 scan List in order of centers
        int Methodminus2Size = 0;
        if (DATriangleInequality.StartingOut)
        {
            Method = -1;
        }
        else
        {
            Methodminus2Size = DATriangleInequality.CurrentInterCenter[RealCenterIndex].NearbyClusters;
        }
        int MethodPosition = -1; // Position in Center List of Center Method = -2
        int RealCenterIndexPrime = -1;
        while (NumCentersToGo > 0)
        {
            if (Method == -1)
            {
                ++CurrentCandidateCenter;
                RealCenterIndexPrime = CurrentCandidateCenter;
            }
            else // Method = -2
            {
                ++MethodPosition;
                if (MethodPosition >= Methodminus2Size)
                {
                    Method = -1;
                    continue;
                }
                RealCenterIndexPrime = DATriangleInequality.CurrentInterCenter[RealCenterIndex].CenterDifferenceIndices[MethodPosition];
            }
            if ((RealCenterIndexPrime < 0) || (RealCenterIndexPrime >= DATriangleInequality.Ncent_RealTotal))
            {
                DAVectorUtility.printAndThrowRuntimeException("Center " + RealCenterIndex + " Method " + Method + " Incorrect Center Index in InterCenter " + RealCenterIndexPrime + " Total " + DATriangleInequality.Ncent_RealTotal + " Thread " + ThreadIndex + " Rank " + DAVectorUtility.MPI_Rank);

            }

            if (LookedAt[RealCenterIndexPrime])
            {
                continue;
            }
            LookedAt[RealCenterIndexPrime] = true;
            NumCentersToGo--;

            if ((RealCenterIndexPrime == RealCenterIndex) || (RealCenterIndexPrime == DATriangleInequality.RealSpongeIndex))
            {
                DATriangleInequality.DistributedCenterDifference[RealCenterIndex - CenterStart][RealCenterIndexPrime] = 0.0;
                continue;
            }
            boolean PrimeChanged = DATriangleInequality.CurrentInterCenter[RealCenterIndexPrime].RecalculateStatus > 0;
            if ((!DATriangleInequality.StartingOut) && (!MainChanged) && (!PrimeChanged) && (currentcut >= 0) && (ActualNumber == NumCentersinOrderedList))
            {
                double AtLeast = DATriangleInequality.DistributedCenterDifference[RealCenterIndex - CenterStart][RealCenterIndexPrime] - CurrentOldDifference - DATriangleInequality.CurrentInterCenter[RealCenterIndexPrime].OldCenterChange;
                if (AtLeast >= smallers[currentcut].value)
                {
                    if (Minrejected < 0.0)
                    {
                        Minrejected = AtLeast;
                    }
                    else
                    {
                        if (Minrejected <= AtLeast)
                        {
                            DATriangleInequality.FindDiagnosticSums_Center.addapoint(ThreadIndex, 1.0, 8);
                            continue;
                        }
                        Minrejected = AtLeast;
                    }
                    continue;
                }
            }
            double Centerdistce = DAVectorParallelism.getNOTSquaredScaledDistancebetweenVectors(DATriangleInequality.CenterY_k_i_Current[RealCenterIndex], DATriangleInequality.CenterY_k_i_Current[RealCenterIndexPrime], DATriangleInequality.CenterSigma_k_i_Current[RealCenterIndex]);
            DATriangleInequality.FindDiagnosticSums_Center.addapoint(ThreadIndex, 1.0, 0);
            Box<Integer> tempRef_currentcut = new Box<>(currentcut);
            Box<Integer> tempRef_ActualNumber = new Box<>(ActualNumber);
            Box<Double> tempRef_Minrejected = new Box<>(Minrejected);
            DATriangleInequality.FindMinimumSetwithRemainder(Centerdistce, RealCenterIndexPrime, tempRef_currentcut, smallers, NumCentersinOrderedList, tempRef_ActualNumber, tempRef_Minrejected);
            currentcut = tempRef_currentcut.content;
            ActualNumber = tempRef_ActualNumber.content;
            Minrejected = tempRef_Minrejected.content;

            if (DATriangleInequality.StartingOut || (MainChanged && PrimeChanged))
            {
                DATriangleInequality.DistributedCenterDifference[RealCenterIndex - CenterStart][RealCenterIndexPrime] = Centerdistce;
            }
            else
            {
                if ((!MainChanged) && PrimeChanged)
                {
                    DATriangleInequality.DistributedCenterDifference[RealCenterIndex - CenterStart][RealCenterIndexPrime] = DAVectorParallelism.getNOTSquaredScaledDistancebetweenVectors(DATriangleInequality.CenterY_k_i_Current[RealCenterIndexPrime], DATriangleInequality.CenterOld[RealCenterIndex], DATriangleInequality.CenterSigma_k_i_Current[RealCenterIndex]);
                    DATriangleInequality.FindDiagnosticSums_Center.addapoint(ThreadIndex, 1.0, 7);
                }
                if (MainChanged && (!PrimeChanged))
                {
                    DATriangleInequality.DistributedCenterDifference[RealCenterIndex - CenterStart][RealCenterIndexPrime] = DAVectorParallelism.getNOTSquaredScaledDistancebetweenVectors(DATriangleInequality.CenterY_k_i_Current[RealCenterIndex], DATriangleInequality.CenterOld[RealCenterIndexPrime], DATriangleInequality.CenterSigma_k_i_Current[RealCenterIndex]);
                    DATriangleInequality.FindDiagnosticSums_Center.addapoint(ThreadIndex, 1.0, 7);
                }
            }

        } // End While loop over centers

        Arrays.sort(smallers, 0, ActualNumber); // last parameter is actually the toIndex(exclusive) 0+ActualNumber

        //  Set Data Structure
        DATriangleInequality.CurrentInterCenter[RealCenterIndex].NearbyClusters = ActualNumber;
        DATriangleInequality.CurrentInterCenter[RealCenterIndex].Dmax = Minrejected;
        DATriangleInequality.FindDiagnosticSums_Center.addapoint(ThreadIndex, Minrejected, 4);
        DATriangleInequality.FindDiagnosticSums_Center.addapoint(ThreadIndex, (double) ActualNumber, 5);

        for (int loopindex = 0; loopindex < ActualNumber; loopindex++)
        {
            DATriangleInequality.CurrentInterCenter[RealCenterIndex].CenterDifference[loopindex] = smallers[loopindex].value;
            DATriangleInequality.FindDiagnosticSums_Center.addapoint(ThreadIndex, smallers[loopindex].value, 6);
            DATriangleInequality.CurrentInterCenter[RealCenterIndex].CenterDifferenceIndices[loopindex] = smallers[loopindex].index;
        }
        double NotStored = Math.max(0.0, (double)(DATriangleInequality.Ncent_RealTotal - 1 - DATriangleInequality.MaxCentersperCenter));

    } // End CFA_CenterCenterCalculation

    public static void InitialPointAnalysis()
    {
        //  Set Associated Cluster and initial lower bounds before any previous point data available
        //  It loops through centers from beginning to end
        //  When a new best center is found, it changes to loop over associated centers
        //  If best center loop exhausted it reverts to ploughing through center list
        // Note - parallel for
        launchHabaneroApp(() -> {
            forallChunked(0, DAVectorUtility.ThreadCount - 1,
                    (threadIndex) -> // End Sum over Threads -  End loop over Points
                    {
                        //  Loop over candidates
                        //  Start with Sponge Cluster if it exists for Method = -1
                        //  Need to look at explicit distance
                        //  New best Center discovered
                        //  BestCenter and BestScore are not updated and not used
                        int indexlen = DAVectorUtility.PointsperThread[threadIndex];
                        int beginpoint = DAVectorUtility.StartPointperThread[threadIndex] - DAVectorUtility.PointStart_Process;
                        for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++) {
                            DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex, 1.0, 1);
                            int NumCentersToGo = DATriangleInequality.Ncent_RealTotal;
                            boolean[] LookedAt = new boolean[DATriangleInequality.Ncent_RealTotal];
                            for (int CenterIndex = 0; CenterIndex < DATriangleInequality.Ncent_RealTotal; CenterIndex++) {
                                LookedAt[CenterIndex] = false;
                            }
                            int CurrentCandidateCenter = -2;
                            if (DATriangleInequality.RealSpongeIndex < 0) {
                                CurrentCandidateCenter = -1;
                            }
                            int Method = -1;
                            int MethodPosition = -1;
                            int BestCenter = -2;
                            double BestScore = -1.0;
                            int ExaminedCenter = -1;
                            int BreakReason = 7;
                            int Num_LBValuesforPoint = DATriangleInequality.MaxClusterLBsperPoint;
                            int ActualNumber = 0;
                            int currentcut = -1;
                            double Minrejected = -1.0;
                            int ActuallyUse = DATriangleInequality.Ncent_RealTotal;
                            Smaller[] smallers = new Smaller[ActuallyUse];
                            while (NumCentersToGo > 0) {
                                if (Method == -1) {
                                    CurrentCandidateCenter++;
                                    if (CurrentCandidateCenter == -1) {
                                        ExaminedCenter = DATriangleInequality.RealSpongeIndex;
                                        if (LookedAt[ExaminedCenter]) {
                                            continue;
                                        }
                                    } else {
                                        if (LookedAt[CurrentCandidateCenter]) {
                                            continue;
                                        }
                                        ExaminedCenter = CurrentCandidateCenter;
                                    }
                                } else {
                                    ++MethodPosition;
                                    double TestScore = BestScore + Math.sqrt(
                                            BestScore * BestScore + DATriangleInequality.ExponentialxTemperature);
                                    if (MethodPosition >= DATriangleInequality.CurrentInterCenter[Method].NearbyClusters) {
                                        if (DATriangleInequality.CurrentInterCenter[Method].Dmax > TestScore) {
                                            BreakReason = 8;
                                            break;
                                        }
                                        Method = -1;
                                        DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex, 1.0, 17);
                                        continue;
                                    }
                                    if (DATriangleInequality.CurrentInterCenter[Method].CenterDifference[MethodPosition] >= TestScore) {
                                        BreakReason = 9;
                                        break;
                                    }
                                    ExaminedCenter = DATriangleInequality.CurrentInterCenter[Method].CenterDifferenceIndices[MethodPosition];
                                    if (LookedAt[ExaminedCenter]) {
                                        continue;
                                    }
                                    DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex, 1.0, 18);
                                }
                                double PointCenterDistce;
                                if (ExaminedCenter == DATriangleInequality.RealSpongeIndex) {
                                    PointCenterDistce = DATriangleInequality.SpongeFactor;
                                } else {
                                    PointCenterDistce = DAVectorParallelism.getNOTSquaredScaledDistancebetweenVectors(
                                            DATriangleInequality.CenterY_k_i_Current[ExaminedCenter],
                                            DATriangleInequality.PointPosition[alpha],
                                            DATriangleInequality.CenterSigma_k_i_Current[ExaminedCenter]);
                                    DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex, 1.0, 0);
                                    DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex, 1.0, 15);
                                }
                                Box<Integer> tempRef_currentcut = new Box<>(
                                        currentcut);
                                Box<Integer> tempRef_ActualNumber = new Box<>(
                                        ActualNumber);
                                Box<Double> tempRef_Minrejected = new Box<>(
                                        Minrejected);
                                DATriangleInequality.FindMinimumSetwithRemainder(PointCenterDistce, ExaminedCenter,
                                        tempRef_currentcut, smallers, Num_LBValuesforPoint, tempRef_ActualNumber,
                                        tempRef_Minrejected);
                                currentcut = tempRef_currentcut.content;
                                ActualNumber = tempRef_ActualNumber.content;
                                Minrejected = tempRef_Minrejected.content;
                                if ((BestCenter >= 0) && (PointCenterDistce > BestScore)) {
                                    --NumCentersToGo;
                                    LookedAt[ExaminedCenter] = true;
                                    continue;
                                }
                                BestScore = PointCenterDistce;
                                BestCenter = ExaminedCenter;
                                if (ExaminedCenter != DATriangleInequality.RealSpongeIndex) {
                                    if (DATriangleInequality.CurrentInterCenter[ExaminedCenter].CenterDifference[0] >= (BestScore + Math.sqrt(
                                            BestScore * BestScore + DATriangleInequality.ExponentialxTemperature))) {
                                        BreakReason = 10;
                                        break;
                                    }
                                    Method = ExaminedCenter;
                                    MethodPosition = -1;
                                }
                                --NumCentersToGo;
                                LookedAt[ExaminedCenter] = true;
                            }
                            DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex, 1.0, BreakReason);
                            DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex,
                                    (double) NumCentersToGo,
                                    19);
                            Arrays.sort(smallers, 0, ActualNumber);
                            if (Num_LBValuesforPoint < ActualNumber) {
                                Minrejected = smallers[Num_LBValuesforPoint].value;
                                ActualNumber = Num_LBValuesforPoint;
                            }
                            if (BestCenter != smallers[0].index) {
                                DAVectorUtility.printAndThrowRuntimeException(
                                        "Inconsistent Best Centers 1 " + smallers[0].index + " " + BestCenter);

                            }
                            if (Math.abs(BestScore - smallers[0].value) > 0.001 * BestScore) {
                                DAVectorUtility.printAndThrowRuntimeException(
                                        "Incorrect Score " + String.format("%1$5.4f", BestScore) + " " + String.format(
                                                "%1$5.4f", smallers[0].value)
                                );

                            }
                            if (BestCenter < 0) {
                                DAVectorUtility.printAndThrowRuntimeException(
                                        "Illegal Best Cluster in Initial Analysis " + BestCenter + " Cluster Count " + DATriangleInequality.Ncent_RealTotal + " Point " + (alpha + DAVectorUtility.PointStart_Process) + " Rank " + DAVectorUtility.MPI_Rank);

                            }
                            if (ActualNumber > DATriangleInequality.Ncent_RealTotal) {
                                DAVectorUtility.printAndThrowRuntimeException(
                                        "Too many centers in Initial Analysis " + ActualNumber + " " + DATriangleInequality.Ncent_RealTotal + " Point " + (alpha + DAVectorUtility.PointStart_Process) + " Rank " + DAVectorUtility.MPI_Rank);

                            }
                            DATriangleInequality.LB_Current_alpha_ClusterPointer[alpha].AssociatedClusterIndex = BestCenter;
                            DATriangleInequality.LB_Current_alpha_ClusterPointer[alpha].NearestDistance = BestScore;
                            DATriangleInequality.LB_Current_alpha_ClusterPointer[alpha].NumCenters = ActualNumber;
                            for (int loopindex = 0; loopindex < ActualNumber; loopindex++) {
                                DATriangleInequality.LB_Current_alpha_ClusterPointer[alpha].DistceLowerBound[loopindex] = smallers[loopindex].value;
                                DATriangleInequality.LB_Current_alpha_ClusterPointer[alpha].CenterIndices[loopindex] = smallers[loopindex].index;
                            }
                            if (DATriangleInequality.OldCenterOption >= 0) {
                                DATriangleInequality.LB_Old_alpha_ClusterPointer[alpha].AssociatedClusterIndex = BestCenter;
                                DATriangleInequality.LB_Old_alpha_ClusterPointer[alpha].NearestDistance = BestScore;
                                DATriangleInequality.LB_Old_alpha_ClusterPointer[alpha].NumCenters = ActualNumber;
                                for (int loopindex = 0; loopindex < ActualNumber; loopindex++) {
                                    DATriangleInequality.LB_Old_alpha_ClusterPointer[alpha].DistceLowerBound[loopindex] = smallers[loopindex].value;
                                    DATriangleInequality.LB_Old_alpha_ClusterPointer[alpha].CenterIndices[loopindex] = smallers[loopindex].index;
                                }
                            }
                            DATriangleInequality.NearestCentertoPoint_alpha[alpha] = BestCenter;
                            DATriangleInequality.Distance_NearestCentertoPoint_alpha[alpha] = BestScore;
                            DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex, (double) ActualNumber,
                                    4);
                            DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex, BestScore, 5);
                            DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex,
                                    DATriangleInequality.LB_Current_alpha_ClusterPointer[alpha].DistceLowerBound[ActualNumber - 1],
                                    6);
                            DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex,
                                    DATriangleInequality.MaxNcent_Global, 23);
                        }

                    }
            );
        });

    } // End InitialPointAnalysis()

    public static void PointDependentAnalysis()
    { // Use on all but first iteration
        //  Use lower bounds and find interesting centers for points

        //  Set Associated Cluster and reset lower bounds following previous iterations
        //  It loops through centers in order of centers near its previous center
        //  When a new best center is found, it changes to loop over its associated centers
        //  If best center loop exhausted it reverts to ploughing through center list
        // Note - parallel for
        launchHabaneroApp(() -> {
            forallChunked(0, DAVectorUtility.ThreadCount - 1,
                    (threadIndex) -> { // End loop over Threads -  End loop over Points
                        // Set up Look up arrays
                        //  Update bounds for recalculated centers in Old arrays
                        //  Update bounds for recalculated centers in Last arrays
                        //  Arrays for finding minimum set
                        // Note temporary arrays in this function are as big as possible
                        //  stored values are restricted by DATriangleInequality.MaxClusterLBsperPoint
                        //  Set up parameters of Lower Bound Finding
                        //  Set up Scan for best center
                        //  Loop over candidates
                        //  Start with Sponge Cluster if it exists for Method = -1
                        //  See if bound exists
                        //  Calculate Explicit Distance
                        //  Note here BestScore is score wrt current best center BestCenter
                        //  PointCenterDistceis new  center distance associated with ExaminedCenter which is NOT new best center
                        //  Next logic purges ExaminedCenter
                        //  New best Center discovered
                        //  Note here BestScore is score wrt current best center BestCenter
                        //  PointCenterDistceis new best center distance associated with ExaminedCenter
                        //  Next logic purges old best center
                        //  Combine Last and Current Lower Bounds (stored above in Small values/indices and store back into Current
                        //  For DA always need sorted list as it picks off initial centers for use by point
                        //  Store LB_Current_alpha_ClusterPointer
                        //  Diagnostics

                        int indexlen = DAVectorUtility.PointsperThread[threadIndex];
                        int beginpoint = DAVectorUtility.StartPointperThread[threadIndex] - DAVectorUtility.PointStart_Process;
                        for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++) {
                            DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex, 1.0, 1);
                            boolean[] LookedAt = new boolean[DATriangleInequality.Ncent_RealTotal];
                            int[] LastBoundLookup = new int[DATriangleInequality.Ncent_RealTotal];
                            int Oldsize = 1;
                            if (DATriangleInequality.OldCenterOption >= 0) {
                                Oldsize = DATriangleInequality.Ncent_RealTotal;
                            }
                            int[] OldBoundLookup = new int[Oldsize];
                            for (int CenterIndex = 0; CenterIndex < DATriangleInequality.Ncent_RealTotal; CenterIndex++) {
                                LookedAt[CenterIndex] = false;
                                if (DATriangleInequality.OldCenterOption >= 0) {
                                    OldBoundLookup[CenterIndex] = -1;
                                }
                                LastBoundLookup[CenterIndex] = -1;
                            }
                            if (DATriangleInequality.OldCenterOption >= 0) {
                                if (DATriangleInequality.LB_Old_alpha_ClusterPointer[alpha].NumCenters > DATriangleInequality.Ncent_RealTotal) {
                                    DAVectorUtility.printAndThrowRuntimeException(
                                            "Incorrect Center Count in Old Bound " + DATriangleInequality.LB_Old_alpha_ClusterPointer[alpha].NumCenters + " " + DATriangleInequality.Ncent_RealTotal + " Point " + (alpha + DAVectorUtility.PointStart_Process) + " Rank " + DAVectorUtility.MPI_Rank);

                                }
                            }
                            if (DATriangleInequality.LB_Last_alpha_ClusterPointer[alpha].NumCenters > DATriangleInequality.Ncent_RealTotal) {
                                DAVectorUtility.printAndThrowRuntimeException(
                                        "Incorrect Center Count in Last Bound " + DATriangleInequality.LB_Last_alpha_ClusterPointer[alpha].NumCenters + " " + DATriangleInequality.Ncent_RealTotal + " Point " + (alpha + DAVectorUtility.PointStart_Process) + " Rank " + DAVectorUtility.MPI_Rank);

                            }
                            if (DATriangleInequality.OldCenterOption >= 0) {
                                for (int loopindex = 0; loopindex < DATriangleInequality.LB_Old_alpha_ClusterPointer[alpha].NumCenters; loopindex++) {
                                    OldBoundLookup[DATriangleInequality.LB_Old_alpha_ClusterPointer[alpha].CenterIndices[loopindex]] = loopindex;
                                }
                            }
                            for (int loopindex = 0; loopindex < DATriangleInequality.LB_Last_alpha_ClusterPointer[alpha].NumCenters; loopindex++) {
                                LastBoundLookup[DATriangleInequality.LB_Last_alpha_ClusterPointer[alpha].CenterIndices[loopindex]] = loopindex;
                            }
                            if (DATriangleInequality.OldCenterOption >= 0) {
                                for (int loopindex = 0; loopindex < DATriangleInequality.LB_Old_alpha_ClusterPointer[alpha].NumCenters; loopindex++) {
                                    int CenterIndex = DATriangleInequality.LB_Old_alpha_ClusterPointer[alpha].CenterIndices[loopindex];
                                    if (DATriangleInequality.CurrentInterCenter[CenterIndex].RecalculateStatus > 0) {
                                        double lowerbound = DAVectorParallelism.getNOTSquaredScaledDistancebetweenVectors(
                                                DATriangleInequality.CenterY_k_i_Current[CenterIndex],
                                                DATriangleInequality.PointPosition[alpha],
                                                DATriangleInequality.CenterSigma_k_i_Current[CenterIndex]);
                                        DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex, 1.0,
                                                11);
                                        DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex, 1.0, 0);
                                        DATriangleInequality.LB_Old_alpha_ClusterPointer[alpha].DistceLowerBound[loopindex] = lowerbound;
                                        if (LastBoundLookup[CenterIndex] < 0) {
                                            continue;
                                        }
                                        DATriangleInequality.LB_Last_alpha_ClusterPointer[alpha].DistceLowerBound[LastBoundLookup[CenterIndex]] = lowerbound;
                                    }
                                }
                            }
                            for (int loopindex = 0; loopindex < DATriangleInequality.LB_Last_alpha_ClusterPointer[alpha].NumCenters; loopindex++) {
                                double lowerbound;
                                int CenterIndex = DATriangleInequality.LB_Last_alpha_ClusterPointer[alpha].CenterIndices[loopindex];
                                if (DATriangleInequality.OldCenterOption >= 0) {
                                    if (OldBoundLookup[CenterIndex] >= 0) {
                                        continue;
                                    }
                                }
                                if (DATriangleInequality.CurrentInterCenter[CenterIndex].RecalculateStatus > 0) {
                                    lowerbound = DAVectorParallelism.getNOTSquaredScaledDistancebetweenVectors(
                                            DATriangleInequality.CenterY_k_i_Current[CenterIndex],
                                            DATriangleInequality.PointPosition[alpha],
                                            DATriangleInequality.CenterSigma_k_i_Current[CenterIndex]);
                                    DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex, 1.0, 12);
                                    DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex, 1.0, 0);
                                    DATriangleInequality.LB_Last_alpha_ClusterPointer[alpha].DistceLowerBound[loopindex] = lowerbound;
                                }
                            }
                            int Num_LBValuesforPoint = DATriangleInequality.MaxClusterLBsperPoint;
                            int ActuallyUse = DATriangleInequality.Ncent_RealTotal;
                            Smaller[] smallers = new Smaller[ActuallyUse];
                            int ActualNumber = 0;
                            int currentcut = -1;
                            double Minrejected = -1.0;
                            int NumCentersToGo = DATriangleInequality.Ncent_RealTotal;
                            int CurrentCandidateCenter = -2;
                            if (DATriangleInequality.RealSpongeIndex < 0) {
                                CurrentCandidateCenter = -1;
                            }
                            int Method = -2;
                            int Bestposition = -2;
                            int MethodPosition = -1;
                            int BestCenter = DATriangleInequality.LB_Last_alpha_ClusterPointer[alpha].AssociatedClusterIndex;
                            if (BestCenter < 0) {
                                Method = -1;
                            }
                            double BestScore = -1.0;
                            int ExaminedCenter = BestCenter;
                            int BreakReason = 7;
                            double BreakValue = -1.0;
                            while (NumCentersToGo > 0) {
                                if (Method == -2) {
                                    Bestposition++;
                                    if (Bestposition == -1) {
                                        ExaminedCenter = BestCenter;
                                        DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex, 1.0, 16);
                                    } else {
                                        if (Bestposition < DATriangleInequality.LB_Last_alpha_ClusterPointer[alpha].NumCenters) {
                                            ExaminedCenter = DATriangleInequality.LB_Last_alpha_ClusterPointer[alpha].CenterIndices[Bestposition];
                                            if (LookedAt[ExaminedCenter]) {
                                                continue;
                                            }
                                            DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex, 1.0,
                                                    16);
                                        } else {
                                            Method = -1;
                                        }
                                    }
                                }
                                if (Method == -1) {
                                    CurrentCandidateCenter++;
                                    if (CurrentCandidateCenter == -1) {
                                        ExaminedCenter = DATriangleInequality.RealSpongeIndex;
                                        if (LookedAt[ExaminedCenter]) {
                                            continue;
                                        }
                                        DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex, 1.0, 17);
                                    } else {
                                        if (LookedAt[CurrentCandidateCenter]) {
                                            continue;
                                        }
                                        ExaminedCenter = CurrentCandidateCenter;
                                        DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex, 1.0, 17);
                                    }
                                } else if (Method >= 0) {
                                    ++MethodPosition;
                                    double TestScore = BestScore + Math.sqrt(
                                            BestScore * BestScore + DATriangleInequality.ExponentialxTemperature);
                                    if (MethodPosition >= DATriangleInequality.CurrentInterCenter[Method].NearbyClusters) {
                                        if (DATriangleInequality.CurrentInterCenter[Method].Dmax > TestScore) {
                                            BreakReason = 8;
                                            BreakValue = DATriangleInequality.CurrentInterCenter[Method].Dmax - BestScore;
                                            break;
                                        }
                                        Method = -2;
                                        continue;
                                    }
                                    if (DATriangleInequality.CurrentInterCenter[Method].CenterDifference[MethodPosition] >= TestScore) {
                                        BreakReason = 9;
                                        BreakValue = DATriangleInequality.CurrentInterCenter[Method].CenterDifference[MethodPosition] - BestScore;
                                        break;
                                    }
                                    ExaminedCenter = DATriangleInequality.CurrentInterCenter[Method].CenterDifferenceIndices[MethodPosition];
                                    if (LookedAt[ExaminedCenter]) {
                                        continue;
                                    }
                                    DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex, 1.0, 18);
                                }
                                if (BestScore > -0.5) {
                                    if (DATriangleInequality.OldCenterOption >= 0) {
                                        int status_old = OldBoundLookup[ExaminedCenter];
                                        if (status_old >= 0) {
                                            if ((DATriangleInequality.LB_Old_alpha_ClusterPointer[alpha].DistceLowerBound[status_old] - CurrentInterCenter[ExaminedCenter].OldCenterChange) > Math.sqrt(
                                                    BestScore * BestScore + DATriangleInequality.ExponentialxTemperature)) {
                                                --NumCentersToGo;
                                                LookedAt[ExaminedCenter] = true;
                                                DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex,
                                                        1.0,
                                                        2);
                                                continue;
                                            } else {
                                                DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex,
                                                        1.0,
                                                        13);
                                            }
                                        }
                                    }
                                    int status_last = LastBoundLookup[ExaminedCenter];
                                    if (status_last >= 0) {
                                        if ((DATriangleInequality.LB_Last_alpha_ClusterPointer[alpha].DistceLowerBound[status_last] - CurrentInterCenter[ExaminedCenter].LastCenterChange) > Math.sqrt(
                                                BestScore * BestScore + DATriangleInequality.ExponentialxTemperature)) {
                                            --NumCentersToGo;
                                            LookedAt[ExaminedCenter] = true;
                                            DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex, 1.0,
                                                    3);
                                            continue;
                                        } else {
                                            DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex, 1.0,
                                                    14);
                                        }
                                    }
                                }
                                if (LastBoundLookup[ExaminedCenter] >= 0) {
                                    LastBoundLookup[ExaminedCenter] = -1;
                                }
                                double PointCenterDistce;
                                if (ExaminedCenter == DATriangleInequality.RealSpongeIndex) {
                                    PointCenterDistce = DATriangleInequality.SpongeFactor;
                                } else {
                                    PointCenterDistce = DAVectorParallelism.getNOTSquaredScaledDistancebetweenVectors(
                                            CenterY_k_i_Current[ExaminedCenter],
                                            DATriangleInequality.PointPosition[alpha],
                                            DATriangleInequality.CenterSigma_k_i_Current[ExaminedCenter]);
                                    DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex, 1.0, 0);
                                    DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex, 1.0, 15);
                                }
                                Box<Integer> tempRef_currentcut = new Box<>(currentcut);
                                Box<Integer> tempRef_ActualNumber = new Box<>(ActualNumber);
                                Box<Double> tempRef_Minrejected = new Box<>(Minrejected);
                                DATriangleInequality.FindMinimumSetwithRemainder(PointCenterDistce, ExaminedCenter,
                                        tempRef_currentcut, smallers, Num_LBValuesforPoint, tempRef_ActualNumber,
                                        tempRef_Minrejected);
                                currentcut = tempRef_currentcut.content;
                                ActualNumber = tempRef_ActualNumber.content;
                                Minrejected = tempRef_Minrejected.content;
                                if ((BestScore > -0.5) && (PointCenterDistce > BestScore)) {
                                    if (Method < 0) {
                                        Method = ExaminedCenter;
                                        MethodPosition = -1;
                                    }
                                    --NumCentersToGo;
                                    LookedAt[ExaminedCenter] = true;
                                    if (ExaminedCenter != DATriangleInequality.RealSpongeIndex) {
                                        Box<Integer> tempRef_NumCentersToGo = new Box<>(
                                                NumCentersToGo);
                                        DATriangleInequality.BackwardTest(threadIndex,
                                                DATriangleInequality.CurrentInterCenter[ExaminedCenter],
                                                tempRef_NumCentersToGo, LookedAt, PointCenterDistce + Math.sqrt(
                                                        BestScore * BestScore + DATriangleInequality.ExponentialxTemperature)
                                        );
                                        NumCentersToGo = tempRef_NumCentersToGo.content;
                                    }
                                    continue;
                                }
                                if ((BestScore > -0.5) && (BestCenter != DATriangleInequality.RealSpongeIndex)) {
                                    Box<Integer> tempRef_NumCentersToGo2 = new Box<>(
                                            NumCentersToGo);
                                    DATriangleInequality.BackwardTest(threadIndex,
                                            DATriangleInequality.CurrentInterCenter[BestCenter],
                                            tempRef_NumCentersToGo2,
                                            LookedAt, BestScore + Math.sqrt(
                                                    PointCenterDistce * PointCenterDistce + DATriangleInequality.ExponentialxTemperature)
                                    );
                                    NumCentersToGo = tempRef_NumCentersToGo2.content;
                                }
                                BestScore = PointCenterDistce;
                                BestCenter = ExaminedCenter;
                                if (ExaminedCenter != DATriangleInequality.RealSpongeIndex) {
                                    if (DATriangleInequality.CurrentInterCenter[ExaminedCenter].CenterDifference[0] >= (BestScore + Math.sqrt(
                                            BestScore * BestScore + DATriangleInequality.ExponentialxTemperature))) {
                                        BreakReason = 10;
                                        BreakValue = DATriangleInequality.CurrentInterCenter[ExaminedCenter].CenterDifference[0] - BestScore;
                                        break;
                                    }
                                    Method = ExaminedCenter;
                                    MethodPosition = -1;
                                }
                                --NumCentersToGo;
                                DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex, 1.0, 22);
                                LookedAt[ExaminedCenter] = true;
                            }
                            DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex, 1.0, BreakReason);
                            DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex,
                                    (double) NumCentersToGo,
                                    19);
                            for (int loopindex = 0; loopindex < DATriangleInequality.LB_Last_alpha_ClusterPointer[alpha].NumCenters; loopindex++) {
                                int CenterIndex = DATriangleInequality.LB_Last_alpha_ClusterPointer[alpha].CenterIndices[loopindex];
                                if (LastBoundLookup[CenterIndex] < 0) {
                                    continue;
                                }
                                double lowerbound = DATriangleInequality.LB_Last_alpha_ClusterPointer[alpha].DistceLowerBound[loopindex] - CurrentInterCenter[CenterIndex].LastCenterChange;
                                if (BreakValue > -0.5) {
                                    lowerbound = Math.max(lowerbound, BreakValue);
                                }
                                Box<Integer> tempRef_currentcut2 = new Box<>(
                                        currentcut);
                                Box<Integer> tempRef_ActualNumber2 = new Box<>(
                                        ActualNumber);
                                Box<Double> tempRef_Minrejected2 = new Box<>(
                                        Minrejected);
                                DATriangleInequality.FindMinimumSetwithRemainder(lowerbound, CenterIndex,
                                        tempRef_currentcut2, smallers, Num_LBValuesforPoint, tempRef_ActualNumber2,
                                        tempRef_Minrejected2);
                                currentcut = tempRef_currentcut2.content;
                                ActualNumber = tempRef_ActualNumber2.content;
                                Minrejected = tempRef_Minrejected2.content;
                            }
                            if (BestCenter < 0) {
                                DAVectorUtility.printAndThrowRuntimeException(
                                        "Illegal Best Cluster in Ongoing Analysis " + BestCenter + " Cluster Count " + DATriangleInequality.Ncent_RealTotal + " Point " + (alpha + DAVectorUtility.PointStart_Process) + " Rank " + DAVectorUtility.MPI_Rank);

                            }
                            if (ActualNumber > DATriangleInequality.Ncent_RealTotal) {
                                DAVectorUtility.printAndThrowRuntimeException(
                                        "Too many centers in Ongoing Analysis " + ActualNumber + " " + DATriangleInequality.Ncent_RealTotal + " Point " + (alpha + DAVectorUtility.PointStart_Process) + " Rank " + DAVectorUtility.MPI_Rank);

                            }
                            Arrays.sort(smallers, 0, ActualNumber);
                            if (Num_LBValuesforPoint < ActualNumber) {
                                Minrejected = smallers[Num_LBValuesforPoint].value;
                                ActualNumber = Num_LBValuesforPoint;
                            }
                            DATriangleInequality.LB_Current_alpha_ClusterPointer[alpha].AssociatedClusterIndex = BestCenter;
                            DATriangleInequality.LB_Current_alpha_ClusterPointer[alpha].NearestDistance = BestScore;
                            DATriangleInequality.LB_Current_alpha_ClusterPointer[alpha].NumCenters = ActualNumber;
                            for (int loopindex = 0; loopindex < ActualNumber; loopindex++) {
                                DATriangleInequality.LB_Current_alpha_ClusterPointer[alpha].DistceLowerBound[loopindex] = smallers[loopindex].value;
                                DATriangleInequality.LB_Current_alpha_ClusterPointer[alpha].CenterIndices[loopindex] = smallers[loopindex].index;
                            }
                            DATriangleInequality.NearestCentertoPoint_alpha[alpha] = BestCenter;
                            DATriangleInequality.Distance_NearestCentertoPoint_alpha[alpha] = BestScore;
                            DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex,
                                    (double) ActualNumber, 4);
                            DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex, BestScore, 5);
                            DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex,
                                    DATriangleInequality.LB_Current_alpha_ClusterPointer[alpha].DistceLowerBound[ActualNumber - 1],
                                    6);
                            DATriangleInequality.FindDiagnosticSums_Points.addapoint(threadIndex,
                                    DATriangleInequality.MaxNcent_Global, 23);
                        }

                    }
            );
        });

    } // End PointDependentAnalysis()


    public static void BackwardTest(int ThreadIndex, CenterData CurrentCenterData, Box<Integer> NumberCentersLeft, boolean[] LookedAtArray, double SumofScores)
    { // Test backwards in Center-Center Array
        DATriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 21);

        for (int CenterIndexPointer = CurrentCenterData.NearbyClusters - 1; CenterIndexPointer >= 0; CenterIndexPointer--)
        {
            if (CurrentCenterData.CenterDifference[CenterIndexPointer] < SumofScores)
            {
                return;
            }
            int CenterIndex = CurrentCenterData.CenterDifferenceIndices[CenterIndexPointer];
            if (LookedAtArray[CenterIndex])
            {
                continue;
            }
            LookedAtArray[CenterIndex] = true;
            --NumberCentersLeft.content;
            DATriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 20);
        }
    }

    //  Return Active Cluster Indices and Number of Clusters
    //  Error code -1 Fatal no clusters
    //  1 Less than Minimum
    //  2 Less than Target
    //  3 Hit Storage Limit
    //  0 OK

    //  Assumes DATriangleInequality.LB_Current_alpha_ClusterPointer[alpha].DistceLowerBound sorted in increasing value
    public static int SetAssociatedCenters(int alpha, Box<Integer> NumberCenters, int[] ListCenters)
    {
        //            public static int SetAssociatedCenters(int alpha, out int NumberCenters, int[] ListCenters, out double BestDistce, double[] ListDistces)
        //BestDistce = 0.0;
        double BestDistce = 0.0;
        NumberCenters.content = DATriangleInequality.LB_Current_alpha_ClusterPointer[alpha].NumCenters;

        if (NumberCenters.content <= 0)
        {
            NumberCenters.content = 0;
            return -1;
        }

        BestDistce = DATriangleInequality.LB_Current_alpha_ClusterPointer[alpha].NearestDistance;
        double TestDistce = BestDistce * BestDistce + Program.ExpArgumentCut2 * DATriangleInequality.Solution.Temperature;

        boolean toomany = false;
        for (int CentersinLBList = 0; CentersinLBList < NumberCenters.content; CentersinLBList++)
        {
            if (CentersinLBList >= Program.maxNcentperPoint)
            {
                toomany = true;
                NumberCenters.content = CentersinLBList;
                break;
            }
            double distance = DATriangleInequality.LB_Current_alpha_ClusterPointer[alpha].DistceLowerBound[CentersinLBList];
            distance = distance * distance;
            if (distance > TestDistce)
            {
                NumberCenters.content = CentersinLBList;
                break;
            }
            ListCenters[CentersinLBList] = ClusteringSolution.ActiveClusterIndices[DATriangleInequality.LB_Current_alpha_ClusterPointer[alpha].CenterIndices[CentersinLBList]];
            //               ListDistces[CentersinLBList] = DATriangleInequality.LB_Current_alpha_ClusterPointer[alpha].DistceLowerBound[CentersinLBList];
        }
        if (NumberCenters.content == 0)
        {
            DAVectorUtility.printAndThrowRuntimeException("No valid Centers out of " + (DATriangleInequality.Ncent_RealTotal) + " Num point list centers " + (DATriangleInequality.LB_Current_alpha_ClusterPointer[alpha].NumCenters) + " Point " + (alpha + DAVectorUtility.PointStart_Process) + " Rank " + (DAVectorUtility.MPI_Rank));

        }
        if (NumberCenters.content < Program.targetMinimumNcentperPoint)
        {
            return 1;
        }
        if (NumberCenters.content < Program.targetNcentperPoint)
        {
            return 2;
        }
        if (toomany)
        {
            return 3;
        }
        return 0;

    } // End SetAssociatedCenters(int alpha, out int NumberCenters, int[] ListCenters)


    public static void PrintDiagnostics()
    {
        DAVectorUtility.SALSAPrint(0, "\nDeterministic Annealing Triangle Inequality Diagnostics ***** Number of Iterations " + DATriangleInequality.CalculationIterations);

        DAVectorUtility.SALSAPrint(0, "Triangle Inequality Option " + DATriangleInequality.UseTriangleInequality + " Old Center " + DATriangleInequality.OldCenterOption);
        DAVectorUtility.SALSAPrint(0, "Test for center change and old lower bounds (normalized by radius) " + String.format("%1$5.4f", DATriangleInequality.TriangleInequality_Delta1_old));
        DAVectorUtility.SALSAPrint(0, "Test for center change and current lower bounds (normalized by radius) " + String.format("%1$5.4f", DATriangleInequality.TriangleInequality_Delta1_current));
        DAVectorUtility.SALSAPrint(0, "Maximum Allowed # Centers " + DATriangleInequality.MaxNcent_Global + " Cut for Center Parallelism " + DATriangleInequality.Ncent_Global_Parallel + " Parallelism used? " + DATriangleInequality.UseParallelismoverCenters);
        DAVectorUtility.SALSAPrint(0, "Vector Dimension " + DATriangleInequality.ParameterVectorDimension);
        DAVectorUtility.SALSAPrint(0, "Maximum number of Lower Bound values for each Point " + DATriangleInequality.MaxClusterLBsperPoint);
        DAVectorUtility.SALSAPrint(0, "Maximum number of Centers in Center Difference Array " + DATriangleInequality.MaxCentersperCenter);

        double tmp1 = 1.0;
        if (DATriangleInequality.CalculationIterations > 0)
        {
            tmp1 = 1.0 / (double)DATriangleInequality.CalculationIterations;
        }
        double tmp2 = 1.0;
        if (DATriangleInequality.CalculationPointIterations > 0)
        {
            tmp2 = 1.0 / DATriangleInequality.CalculationPointIterations;
        }
        double tmp3 = 1.0;
        if (DATriangleInequality.CalculationCenterIterations > 0)
        {
            tmp3 = 1.0 / DATriangleInequality.CalculationCenterIterations;
        }
        double tmp4 = 1.0;
        if (DATriangleInequality.CountCCEntries > 0)
        {
            tmp4 = 1.0 / DATriangleInequality.CountCCEntries;
        }

        double tmp5 = 1.0;
        if (DATriangleInequality.RadiusTest1 > 0)
        {
            tmp5 = 1.0 / DATriangleInequality.RadiusTest1;
        }

        DAVectorUtility.SALSAPrint(0, "\nPoints per Iteration " + String.format("%1$3.2f", DATriangleInequality.CalculationPointIterations * tmp1));
        DAVectorUtility.SALSAPrint(0, "Centers per Iteration " + String.format("%1$3.2f", DATriangleInequality.CalculationCenterIterations * tmp1));
        DAVectorUtility.SALSAPrint(0, " Total Center-Center Calc " + String.format("%1$5.4E", DATriangleInequality.NumberFullDistancesCalculatedCC) + " Point Center " + String.format("%1$5.4E", DATriangleInequality.NumberFullDistancesCalculatedPC));
        DAVectorUtility.SALSAPrint(0, "Center-Center Distances per Center-Iteration " + String.format("%1$7.6f", DATriangleInequality.NumberFullDistancesCalculatedCC * tmp3));
        DAVectorUtility.SALSAPrint(0, "Center-Center Distances per Center-Iteration " + String.format("%1$7.6f", DATriangleInequality.NumberFullDistancesCalculatedCC_off * tmp3) + " Mixed Update status");
        DAVectorUtility.SALSAPrint(0, "Center-Center Distances per Center-Iteration " + String.format("%1$7.6f", DATriangleInequality.NumberFullDistancesCalculatedCC_CleanReject * tmp3) + " Cleanly Rejected");
        DAVectorUtility.SALSAPrint(0, "Center-Center Distances per Center-Iteration " + String.format("%1$7.6f", DATriangleInequality.NumberFullDistancesCalculatedCC_NotStored * tmp3) + " Not Stored");
        DAVectorUtility.SALSAPrint(0, "Center Moving Too much per Iteration " + String.format("%1$7.6f", DATriangleInequality.CountToobigAShift * tmp1));
        DAVectorUtility.SALSAPrint(0, "Center Zero Radius " + String.format("%1$7.6f", DATriangleInequality.RadiusTest2 * tmp5));
        DAVectorUtility.SALSAPrint(0, "Center Change < 0.05 Radius " + String.format("%1$7.6f", DATriangleInequality.RadiusTest3 * tmp5));
        DAVectorUtility.SALSAPrint(0, "Center Change < 0.1 Radius " + String.format("%1$7.6f", DATriangleInequality.RadiusTest4 * tmp5));
        DAVectorUtility.SALSAPrint(0, "Center Change < 0.2 Radius " + String.format("%1$7.6f", DATriangleInequality.RadiusTest5 * tmp5));

        DAVectorUtility.SALSAPrint(0, "\nAverage Center Radius " + String.format("%1$7.6f", DATriangleInequality.AccumulateRadius * tmp3));
        DAVectorUtility.SALSAPrint(0, "Largest Center Center Distance NOT Recorded (-1 if all done) " + String.format("%1$7.6f", DATriangleInequality.AccumulateDmax * tmp3));
        DAVectorUtility.SALSAPrint(0, "Average number of recorded Center-Center distances " + String.format("%1$7.6f", DATriangleInequality.CountCCEntries * tmp3));
        DAVectorUtility.SALSAPrint(0, "Average value  of recorded Center-Center distances " + String.format("%1$7.6f", DATriangleInequality.AccumulateCCvalues * tmp4));

        DAVectorUtility.SALSAPrint(0, "\nAverage number of OLD Lower Bounds tested Good " + String.format("%1$7.6f", DATriangleInequality.AccumulateOLD_LB * tmp2));
        DAVectorUtility.SALSAPrint(0, "Average number of LAST Lower Bounds tested Good " + String.format("%1$7.6f", DATriangleInequality.AccumulateCURRENT_LB * tmp2));
        DAVectorUtility.SALSAPrint(0, "Average number of OLD Lower Bounds tested Bad " + String.format("%1$7.6f", DATriangleInequality.AccumulateOLD_LBFail * tmp2));
        DAVectorUtility.SALSAPrint(0, "Average number of LAST Lower Bounds tested Bad " + String.format("%1$7.6f", DATriangleInequality.AccumulateCURRENT_LBFail * tmp2));

        DAVectorUtility.SALSAPrint(0, "\nAverage number of CURRENT Lower Bound  distance entries " + String.format("%1$7.6f", DATriangleInequality.CountLBEntries * tmp2));
        DAVectorUtility.SALSAPrint(0, "Average number of CURRENT Best Score " + String.format("%1$7.6f", DATriangleInequality.AccumulateBestScore * tmp2));
        DAVectorUtility.SALSAPrint(0, "Average number of CURRENT Maximum value of Lower Bound " + String.format("%1$7.6f", DATriangleInequality.AccumulateMAX_LB * tmp2));

        DAVectorUtility.SALSAPrint(0, "\nPoint-Centers Distances per Point-Iteration " + String.format("%1$7.6f", DATriangleInequality.PossibleNumberFullDistancesCalculated * tmp2) + " Possible");
        DAVectorUtility.SALSAPrint(0, "Point-Centers Distances per Point-Iteration " + String.format("%1$7.6f", DATriangleInequality.NumberFullDistancesCalculatedPC * tmp2));
        DAVectorUtility.SALSAPrint(0, "Point-Centers Distances per Point-Iteration " + String.format("%1$7.6f", DATriangleInequality.NumberFullDistancesCalculatedPC_Old * tmp2) + " OLD refresh");
        DAVectorUtility.SALSAPrint(0, "Point-Centers Distances per Point-Iteration " + String.format("%1$7.6f", DATriangleInequality.NumberFullDistancesCalculatedPC_Last * tmp2) + " LAST refresh");
        DAVectorUtility.SALSAPrint(0, "Point-Centers Distances per Point-Iteration " + String.format("%1$7.6f", DATriangleInequality.NumberFullDistancesCalculatedLoop * tmp2) + " In while loop");

        DAVectorUtility.SALSAPrint(0, "\nAverage number of Breaks in LB/Center distance test as while limit hit " + String.format("%1$7.6f", DATriangleInequality.AccumulateBreakReason7 * tmp2));
        DAVectorUtility.SALSAPrint(0, "Average number of Breaks in LB/Center distance test as Center Scan Dmax   " + String.format("%1$7.6f", DATriangleInequality.AccumulateBreakReason8 * tmp2));
        DAVectorUtility.SALSAPrint(0, "Average number of Breaks in LB/Center distance test as Center Scan Middle " + String.format("%1$7.6f", DATriangleInequality.AccumulateBreakReason9 * tmp2));
        DAVectorUtility.SALSAPrint(0, "Average number of Breaks in LB/Center distance test as Center Scan Start  " + String.format("%1$7.6f", DATriangleInequality.AccumulateBreakReason10 * tmp2));
        DAVectorUtility.SALSAPrint(0, "Number of Centers NOT looked at " + String.format("%1$7.6f", DATriangleInequality.AverageCentersLookedAt * tmp2) + " After while loop");
        DAVectorUtility.SALSAPrint(0, "Number of Centers removed with Method  -2 (Best) " + String.format("%1$7.6f", DATriangleInequality.NumberCentersinBestSet * tmp2) + " In while loop");
        DAVectorUtility.SALSAPrint(0, "Number of Centers removed with Method  -1 (Linear) " + String.format("%1$7.6f", DATriangleInequality.NumberCentersinSimpleSet * tmp2) + " In while loop");
        DAVectorUtility.SALSAPrint(0, "Number of Centers removed with Method  >=0 (Center facing) " + String.format("%1$7.6f", DATriangleInequality.NumberCentersinCenterbasedset * tmp2) + " In while loop");
        DAVectorUtility.SALSAPrint(0, "Number of such Center based forward Tests " + String.format("%1$7.6f", DATriangleInequality.NumberForwardTests * tmp2) + " In while loop");
        DAVectorUtility.SALSAPrint(0, "Number of Centers removed with Backward Test " + String.format("%1$7.6f", DATriangleInequality.NumberCentersRemovedBackwards * tmp2) + " In while loop");
        DAVectorUtility.SALSAPrint(0, "Number of Center based Backward Tests " + String.format("%1$7.6f", DATriangleInequality.NumberBackwardTests * tmp2) + " In while loop");

    } // End PrintDiagnostics()


    //  Support finding list of minimum values by inserting new point with value newvalue and index newindex into lists SmallValues SmallIndices
    //  In SmallIndices negative values correspond to unset values
    //  NumberSmallOnes is total number wanted
    //  Currentcut is position in SmallValues, SmallIndices of largest min value
    //  Minrejected is Minimum of values not included
    //  Minrejected must be set to -1.0 and currentcut to -1
    public static void FindMinimumSetwithRemainder(double newvalue, int newindex, Box<Integer> currentcut, Smaller[] smallers, int NumberSmallones, Box<Integer> ActualNumber, Box<Double> Minrejected)
    {
        if (currentcut.content < 0)
        {
            currentcut.content = 0;
            smallers[0].value = newvalue;
            smallers[0].index = newindex;
            ActualNumber.content = 1;
            return;
        }
        if (ActualNumber.content < NumberSmallones)
        { // Not all positions are filled so add at next available
            // Reset currentcut if worst
            smallers[ActualNumber.content].value = newvalue;
            smallers[ActualNumber.content].index = newindex;
            if (smallers[ActualNumber.content].value > smallers[currentcut.content].value)
            {
                currentcut.content = ActualNumber.content;
            }
            ++ActualNumber.content;
            return;
        }
        if (newvalue >= smallers[currentcut.content].value)
        {
            if (Minrejected.content < 0.0)
            {
                Minrejected.content = newvalue;
                return;
            }
            Minrejected.content = Math.min(Minrejected.content, newvalue);
            return;
        }

        // Replace currentcut position with new values and Reset new worst position
        Minrejected.content = smallers[currentcut.content].value;
        smallers[currentcut.content].value = newvalue;
        smallers[currentcut.content].index = newindex;
        double maxvalue = -1.0;
        for (int ivalue = 0; ivalue < NumberSmallones; ivalue++)
        {
            if (smallers[ivalue].index < 0)
            {
                continue;
            }
            if (smallers[ivalue].value > maxvalue)
            {
                currentcut.content = ivalue;
                maxvalue = smallers[ivalue].value;
            }
        }

    } // End FindMinimumSetwithRemainder

} // end class TriangleInequality

// End Namespace edu.indiana.soic.spidal.davs