package edu.indiana.soic.spidal.davs;

import com.google.common.base.Strings;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;
import edu.indiana.soic.spidal.general.Box;
import mpi.MPIException;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.regex.Pattern;

import static edu.rice.hj.Module0.launchHabaneroApp;
import static edu.rice.hj.Module1.forallChunked;


public class Kmeans {

    public static double[][] PointPosition; // Point Positions
    public static int[] InitialPointAssignment; // Initial Assignment of Points to Clusters
    public static int MaxNcent_Global = 50; // Maximum number of centers
    public static int ParameterVectorDimension = 2; // Vector Dimension of Points
    public static int Ncent_Global = 1; // Total Number of Clusters
    public static int Ncent_Global_Parallel = 50; // Use Parallelism if total number of centers greater than this
    public static double[][] ClusterCenter; // Cluster Centers stored in every node
    public static int[] ClusterSize; // Size of Cluster (number of points)
    public static double[] ClusterRadius; // Radius of Cluster
    public static double[] ClusterWidth; // Width of Cluster
    //  Center Parallelism over threads and processes
    public static int CenterCount_Process = 0; // Total number of Centers summed over all threads calculated in this
    // process
    public static int CenterCount_Largest = 0; // Largest number of points in all processes
    public static int CenterStart_Process = 0; //    First data point in this process
    //	Within a process, data points will be divided into ThreadCount segments( each thread a segment),
    // the array keep the size of each segment
    public static int[] FullParallel_CentersperThread = null; //how many data points each thread will take care for
    // thread+process parallelism
    public static int[] FullParallel_StartCenterperThread = null; //the starting point that a thread will take care
    // for thread+process parallelism
    //  Center Parallelism over threads withon each process
    public static int[] LocalParallel_CentersperThread = null; //how many data points each thread will take care for
    // thread+process parallelism
    public static int[] LocalParallel_StartCenterperThread = null; //the starting point that a thread will take care
    // for thread+process parallelism
    public static boolean UseParallelismoverCenters = false; // If true major center loops run in parallel over nodes
    public static int CountKmeansIterations = 0; // Count Iterations in Kmeans
    public static int KmeansIterationCut = 1000; // Cut on Iterations to stop
    public static double KmeansCenterChangeCut = 0.001; // Fractional change to stop
    public static double AverageCenterChange = 0.0; //  Average Center Change each iteration
    public static double AverageRadius = 0.0; // Average Radius each iteration
    public static double AverageWidth = 0.0; // Average Width each iteration

    public static void InitializeKmeans(double[][] PointPositionINPUT, String FileName, int ClusterPosition,
                                        int FirstClusterValue, int StartingPosition, int Ncent_GlobalINPUT,
                                        int MaxNcent_GlobalINPUT, int Ncent_Global_ParallelINPUT,
                                        int ParameterVectorDimensionINPUT, double CenterChangeINPUT,
                                        int IterationCutINPUT) {

        Kmeans.PointPosition = PointPositionINPUT;

        Kmeans.Ncent_Global = Ncent_GlobalINPUT;
        Kmeans.MaxNcent_Global = MaxNcent_GlobalINPUT;
        Kmeans.Ncent_Global_Parallel = Ncent_Global_ParallelINPUT;
        Kmeans.KmeansCenterChangeCut = CenterChangeINPUT;
        Kmeans.KmeansIterationCut = IterationCutINPUT;

        Kmeans.ParameterVectorDimension = ParameterVectorDimensionINPUT;

        Kmeans.SetParallelCenterDecomposition();

        Kmeans.InitialPointAssignment = new int[DAVectorUtility.PointCount_Process];
        if (FileName.length() == 0) {
            return;
        }

        //  Read Initial Assignments
        DAVectorUtility.SALSAPrint(0, "Kmeans Read File " + FileName + " Points " + DAVectorUtility.PointCount_Global +
                " Starting at position " + StartingPosition + " Dimension " + Kmeans.ParameterVectorDimension + " " +
                "Cluster Position " + ClusterPosition + " Initial value " + FirstClusterValue);
        Kmeans.ReadDataFromFile(FileName, ClusterPosition, FirstClusterValue, StartingPosition);

    } // End InitializeKmeans

    public static void SetupKmeans(int[] FullAssignment) { // Set up initial assignment from existing array

        // Note - parallel for
        launchHabaneroApp(() -> {
            forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) ->

            {
                int indexlen = DAVectorUtility.PointsperThread[threadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[threadIndex] - DAVectorUtility.PointStart_Process;
                System.arraycopy(FullAssignment, beginpoint + DAVectorUtility.PointStart_Process,
                        Kmeans.InitialPointAssignment, beginpoint,
                        indexlen + beginpoint - beginpoint); //  End loop over Points

            }); // End Sum over Threads
        });

    } // End SetupKmeans

    public static void RunKmeans(double[][] ClusterCenterINPUT, int[] ClusterSizeINPUT, double[] ClusteRadiusINPUT,
                                 Box<Integer> Ncent_GlobalFINAL,
                                 Box<Double> AverageWidthFINAL) throws MPIException {

        java.util.ArrayList KeepPCfractions = new java.util.ArrayList(200);
        java.util.ArrayList KeepCCfractions = new java.util.ArrayList(200);

        //  Inherit Solution arrays
        Kmeans.ClusterCenter = ClusterCenterINPUT;
        Kmeans.ClusterSize = ClusterSizeINPUT;
        Kmeans.ClusterRadius = ClusteRadiusINPUT;
        Kmeans.ClusterWidth = new double[Kmeans.MaxNcent_Global];

        //  Set up TriangleInequality
        KmeansTriangleInequality.SetExternalFunctions(Kmeans::GetClusterRadius, Kmeans::GetClusterCenters, Kmeans::FindClusterCenters);
        KmeansTriangleInequality
                .InitializeTriangleInequality(Kmeans.PointPosition, Kmeans.ClusterCenter, Kmeans.Ncent_Global,
                                              Kmeans.MaxNcent_Global, Kmeans.Ncent_Global_Parallel,
                                              Kmeans.ParameterVectorDimension);

        DAVectorUtility.SALSAPrint(0, "Start Kmeans ****** Number of Centers " + Kmeans.Ncent_Global + " Max Number " +
                "of Centers " + Kmeans.MaxNcent_Global + " Center Limit for Parallelism " +
                Kmeans.Ncent_Global_Parallel + " Vector Dimension " + Kmeans.ParameterVectorDimension);

        Kmeans.FindClusterCenters(true, Kmeans.InitialPointAssignment, null, null);
        Kmeans.CountKmeansIterations = 0;
        boolean StartStop = false;
        int CountStops = 0;
        while (Kmeans.CountKmeansIterations < Kmeans.KmeansIterationCut) {
            double save1 = KmeansTriangleInequality.NumberFullDistancesCalculatedCC;
            double save2 = KmeansTriangleInequality.NumberFullDistancesCalculatedPC;

            KmeansTriangleInequality.NextIteration();
            ++Kmeans.CountKmeansIterations;
            boolean WillStop = false;
            if (!StartStop) {
                if (Kmeans.AverageCenterChange < Kmeans.AverageRadius * Kmeans.KmeansCenterChangeCut) {
                    StartStop = true;
                }
            } else {
                ++CountStops;
                if (CountStops > 10) {
                    WillStop = true;
                }
            }
            double tmp1 = (KmeansTriangleInequality.NumberFullDistancesCalculatedCC - save1) /
                    (double) Kmeans.MaxNcent_Global;
            double tmp2 = (KmeansTriangleInequality.NumberFullDistancesCalculatedPC - save2) /
                    ((double) Kmeans.MaxNcent_Global * (double) DAVectorUtility.PointCount_Global);
            double tmp3 = KmeansTriangleInequality.NumberFullDistancesCalculatedPC / ((double) Kmeans.MaxNcent_Global *
                    (double) (DAVectorUtility.PointCount_Global * Kmeans.CountKmeansIterations));
            double tmp4 = (KmeansTriangleInequality.NumberFullDistancesCalculatedPC +
                    KmeansTriangleInequality.NumberFullDistancesCalculatedCC) / ((double) Kmeans.MaxNcent_Global *
                    (double) (DAVectorUtility.PointCount_Global * Kmeans.CountKmeansIterations));
            DAVectorUtility.SALSAPrint(0, "Iteration " + Kmeans.CountKmeansIterations + " Average Center Change " +
                    String.format("%1$5.4E", Kmeans.AverageCenterChange) + " Average Radius " +
                    String.format("%1$5.4E", Kmeans.AverageRadius) + " Average Width " +
                    String.format("%1$5.4E", Kmeans.AverageWidth) + " CC calcs per C " + String.format("%1$5.4f", tmp1) +
                    " PC calcs per P&C " + String.format("%1$7.6f", tmp2) + " Cumul PC / Max " +
                    String.format("%1$7.6f", tmp3) + " Cumul PC+CC / PC Max " + String.format("%1$7.6f", tmp4));
            KeepPCfractions.add(tmp2);
            KeepCCfractions.add(tmp1 / DAVectorUtility.PointCount_Global);
            if (((Kmeans.CountKmeansIterations % 10) == 1) || WillStop) {
                String message = " Sizes";
                for (int CenterIndex = 0; CenterIndex < Kmeans.Ncent_Global; CenterIndex++) {
                    message += " " + Kmeans.ClusterSize[CenterIndex];
                }
                DAVectorUtility.SALSAPrint(0, message);
            }
            if (WillStop) {
                break;
            }
        }
        DAVectorUtility.SALSAPrint(0, "End Kmeans Iterations " + Kmeans.CountKmeansIterations + " Iteration Cut " +
                Kmeans.KmeansIterationCut + " Average Center Change " +
                String.format("%1$5.4E", Kmeans.AverageCenterChange) + " Average Radius " +
                String.format("%1$5.4E", Kmeans.AverageRadius) + " Average Width " +
                String.format("%1$5.4E", Kmeans.AverageWidth) + " Fractional Cut " +
                String.format("%1$5.4f", Kmeans.KmeansCenterChangeCut));
        KmeansTriangleInequality.PrintDiagnostics();
        String messagePC = "\nPC Calcs per Point iteration";
        String messageCC = "\nCC Calcs per Point iteration";
        int numPC = KeepPCfractions.size();
        for (int linecount = 0; linecount < numPC; linecount++) {
            messagePC += " " + String.format("%1$5.4f", (double) KeepPCfractions.get(linecount)) + ",";
            messageCC += " " + String.format("%1$5.4f", (double) KeepCCfractions.get(linecount)) + ",";
        }
        DAVectorUtility.SALSAPrint(0, messagePC);
        DAVectorUtility.SALSAPrint(0, messageCC);
        Ncent_GlobalFINAL.content = Kmeans.Ncent_Global;
        AverageWidthFINAL.content = Kmeans.AverageWidth;

        //  Print Histograms
        if (KmeansTriangleInequality.UseTriangleInequality != 0) {
            KmeansTriangleInequality.PlotPointHistograms(Math.sqrt(AverageWidthFINAL.content));
            KmeansTriangleInequality.PlotCenterHistograms(Math.sqrt(AverageWidthFINAL.content));
        }

    } // End RunKmeans()

    public static void SetParallelCenterDecomposition() {
        Kmeans.UseParallelismoverCenters = false;
        if (Kmeans.Ncent_Global <= Kmeans.Ncent_Global_Parallel) {
            Kmeans.CenterStart_Process = 0;
            Kmeans.CenterCount_Process = Kmeans.MaxNcent_Global;
            Kmeans.CenterCount_Largest = Kmeans.MaxNcent_Global;
            return;
        }
        Kmeans.UseParallelismoverCenters = true;

        //	First divide centers among processes
        Range[] processRanges = RangePartitioner.Partition(Kmeans.MaxNcent_Global, DAVectorUtility.MPI_Size);
        Range processRange = processRanges[DAVectorUtility.MPI_Rank]; // The answer for this process

        Kmeans.CenterStart_Process = processRange.getStartIndex();
        Kmeans.CenterCount_Process = processRange.getLength();

        Kmeans.CenterCount_Largest = Integer.MIN_VALUE;
        for (Range r : processRanges) {
            Kmeans.CenterCount_Largest = Math.max(r.getLength(), Kmeans.CenterCount_Largest);
        }


        //	Now divide centers among threads for this process
        Range[] Full_ThreadRanges = RangePartitioner.Partition(processRange, DAVectorUtility.ThreadCount);
        Kmeans.FullParallel_CentersperThread = new int[DAVectorUtility.ThreadCount];
        Kmeans.FullParallel_StartCenterperThread = new int[DAVectorUtility.ThreadCount];

        for (int j = 0; j < DAVectorUtility.ThreadCount; j++) {
            Kmeans.FullParallel_CentersperThread[j] = Full_ThreadRanges[j].getLength();
            Kmeans.FullParallel_StartCenterperThread[j] = Full_ThreadRanges[j].getStartIndex();
        }

        //  Process only Center Parallelism
        Range[] Local_ThreadRanges = RangePartitioner.Partition(Kmeans.MaxNcent_Global, DAVectorUtility.ThreadCount);
        Kmeans.LocalParallel_CentersperThread = new int[DAVectorUtility.ThreadCount];
        Kmeans.LocalParallel_StartCenterperThread = new int[DAVectorUtility.ThreadCount];

        for (int j = 0; j < DAVectorUtility.ThreadCount; j++) {
            Kmeans.LocalParallel_CentersperThread[j] = Local_ThreadRanges[j].getLength();
            Kmeans.LocalParallel_StartCenterperThread[j] = Local_ThreadRanges[j].getStartIndex();
        }

    } // End SetParallelCenterDecomposition()

    public static void GetClusterCenters(int CenterIndex,
                                         double[] IndividualClusterCenter) { //  Return Cluster Center vector

        IndividualClusterCenter = Kmeans.ClusterCenter[CenterIndex];

    } // End GetClusterCenters

    public static double GetClusterRadius(int CenterIndex) { // Return Cluster Radius

        return Kmeans.ClusterRadius[CenterIndex];

    } // End GetClusterRadius

    // Use LastClusterCenter if Size 0
    public static void FindClusterCenters(boolean begin, int[] NearestCentertoPoint,
                                          double[] Distance_NearestCentertoPoint,
                                          double[][] LastClusterCenter) throws MPIException { // Calculate Cluster Parameters

        GlobalReductions.FindVectorDoubleSum3 FindCenterVectorSums =
                new GlobalReductions.FindVectorDoubleSum3(DAVectorUtility.ThreadCount, Kmeans.ParameterVectorDimension,
                                                          Kmeans.Ncent_Global);
        GlobalReductions.FindVectorIntSum FindCenterSizeSums =
                new GlobalReductions.FindVectorIntSum(DAVectorUtility.ThreadCount, Kmeans.Ncent_Global);
        GlobalReductions.FindVectorDoubleMax FindClusterRadius =
                new GlobalReductions.FindVectorDoubleMax(DAVectorUtility.ThreadCount, Kmeans.Ncent_Global);
        GlobalReductions.FindVectorDoubleSum FindClusterWidth =
                new GlobalReductions.FindVectorDoubleSum(DAVectorUtility.ThreadCount, Kmeans.Ncent_Global);

        // Note - parallel for
        launchHabaneroApp(() -> {
            forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) ->
            {
                FindCenterVectorSums.startthread(threadIndex);
                FindCenterSizeSums.startthread(threadIndex);
                int indexlen = DAVectorUtility.PointsperThread[threadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[threadIndex] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++) {
                    int ClusterIndex = NearestCentertoPoint[alpha];
                    if ((ClusterIndex >= Kmeans.Ncent_Global) || (ClusterIndex < 0)) {
                        DAVectorUtility.printAndThrowRuntimeException(
                                "Illegal Cluster Index " + ClusterIndex + " Number " +
                                        "" + Kmeans.Ncent_Global + " Point " +
                                        (alpha + DAVectorUtility.PointStart_Process) +
                                        " " +
                                        "Rank " + DAVectorUtility.MPI_Rank);

                    }
                    FindCenterVectorSums.addapoint(threadIndex, PointPosition[alpha], ClusterIndex);
                    FindCenterSizeSums.addapoint(threadIndex, 1, ClusterIndex);
                    if (!begin) {
                        FindClusterRadius.addapoint(threadIndex, Distance_NearestCentertoPoint[alpha], ClusterIndex);
                        FindClusterWidth.addapoint(threadIndex, Distance_NearestCentertoPoint[alpha] *
                                Distance_NearestCentertoPoint[alpha], ClusterIndex);
                    }
                } // End loop over Points

            }); // End Sum over Threads
        });

        FindCenterVectorSums.sumoverthreadsandmpi();
        FindCenterSizeSums.sumoverthreadsandmpi();
        if (!begin) {
            FindClusterRadius.sumoverthreadsandmpi();
            FindClusterWidth.sumoverthreadsandmpi();
        }
        Kmeans.AverageRadius = 0.0;
        Kmeans.AverageWidth = 0.0;
        Kmeans.AverageCenterChange = 0.0;

        if (Kmeans.UseParallelismoverCenters) { // Centers Parallel over Threads NOT nodes
            double[] AccumulateRadius = new double[DAVectorUtility.ThreadCount];
            double[] AccumulateWidth = new double[DAVectorUtility.ThreadCount];
            double[] AccumulateCenterChange = new double[DAVectorUtility.ThreadCount];
            for (int ThreadIndex = 0; ThreadIndex < DAVectorUtility.ThreadCount; ThreadIndex++) {
                AccumulateRadius[ThreadIndex] = 0.0;
                AccumulateWidth[ThreadIndex] = 0.0;
                AccumulateCenterChange[ThreadIndex] = 0.0;
            }

            // Note - parallel for
            launchHabaneroApp(() -> {
                forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) ->
                {
                    int indexlen = KmeansTriangleInequality.LocalParallel_CentersperThread[threadIndex];
                    int beginpoint = KmeansTriangleInequality.LocalParallel_StartCenterperThread[threadIndex];
                    for (int CenterIndex = beginpoint; CenterIndex < indexlen + beginpoint; CenterIndex++) {
                        Kmeans.ClusterSize[CenterIndex] = FindCenterSizeSums.TotalVectorSum[CenterIndex];
                        if (Kmeans.ClusterSize[CenterIndex] > 0) {
                            double divisor = 1.0 / (double) Kmeans.ClusterSize[CenterIndex];
                            for (int VectorIndex = 0; VectorIndex < Kmeans.ParameterVectorDimension; VectorIndex++) {
                                int totalindex = VectorIndex + CenterIndex * Kmeans.ParameterVectorDimension;
                                Kmeans.ClusterCenter[CenterIndex][VectorIndex] =
                                        FindCenterVectorSums.TotalVectorSum[totalindex] * divisor;
                            }
                            if (!begin) {
                                Kmeans.ClusterRadius[CenterIndex] = FindClusterRadius.TotalVectorMax[CenterIndex];
                                Kmeans.ClusterWidth[CenterIndex] = FindClusterWidth.TotalVectorSum[CenterIndex];
                                AccumulateRadius[threadIndex] += Kmeans.ClusterRadius[CenterIndex];
                                AccumulateWidth[threadIndex] += Kmeans.ClusterWidth[CenterIndex];
                                AccumulateCenterChange[threadIndex] += DAVectorParallelism
                                        .getNOTSquaredUNScaledDistancebetweenVectors(Kmeans.ClusterCenter[CenterIndex],
                                                LastClusterCenter[CenterIndex]);
                            }
                        } else {
                            if (begin) {
                                DAVectorUtility.printAndThrowRuntimeException(
                                        "Empty Input Cluster " + CenterIndex + " " +
                                                "Number " + Kmeans.Ncent_Global +
                                                " Rank " + DAVectorUtility.MPI_Rank);

                            }
                            System.arraycopy(LastClusterCenter[CenterIndex], 0, Kmeans.ClusterCenter[CenterIndex], 0,
                                    Kmeans.ParameterVectorDimension);
                            Kmeans.ClusterRadius[CenterIndex] = 0.0;
                            Kmeans.ClusterWidth[CenterIndex] = 0.0;
                        }
                    }

                }); // End Sum over Threads
            });

            if (!begin) {
                for (int ThreadIndex = 0; ThreadIndex < DAVectorUtility.ThreadCount; ThreadIndex++) {
                    Kmeans.AverageRadius += AccumulateRadius[ThreadIndex];
                    Kmeans.AverageWidth += AccumulateWidth[ThreadIndex];
                    Kmeans.AverageCenterChange += AccumulateCenterChange[ThreadIndex];
                }
            }
        } else { // Centers Sequential

            for (int CenterIndex = 0; CenterIndex < Kmeans.Ncent_Global; CenterIndex++) {
                Kmeans.ClusterSize[CenterIndex] = FindCenterSizeSums.TotalVectorSum[CenterIndex];

                if (Kmeans.ClusterSize[CenterIndex] > 0) {
                    double divisor = 1.0 / (double) Kmeans.ClusterSize[CenterIndex];
                    for (int VectorIndex = 0; VectorIndex < Kmeans.ParameterVectorDimension; VectorIndex++) {
                        int totalindex = VectorIndex + CenterIndex * Kmeans.ParameterVectorDimension;
                        Kmeans.ClusterCenter[CenterIndex][VectorIndex] =
                                FindCenterVectorSums.TotalVectorSum[totalindex] * divisor;
                    }

                    if (!begin) {
                        Kmeans.ClusterRadius[CenterIndex] = FindClusterRadius.TotalVectorMax[CenterIndex];
                        Kmeans.ClusterWidth[CenterIndex] = FindClusterWidth.TotalVectorSum[CenterIndex];
                        Kmeans.AverageRadius += Kmeans.ClusterRadius[CenterIndex];
                        Kmeans.AverageWidth += Kmeans.ClusterWidth[CenterIndex];
                        Kmeans.AverageCenterChange += DAVectorParallelism
                                .getNOTSquaredUNScaledDistancebetweenVectors(Kmeans.ClusterCenter[CenterIndex],
                                                                             LastClusterCenter[CenterIndex]);
                    }
                } else {
                    if (begin) {
                        DAVectorUtility.printAndThrowRuntimeException("Empty Input Cluster " + CenterIndex + " Number" +
                                                                              " " + Kmeans.Ncent_Global + " Rank " +
                                                                              DAVectorUtility.MPI_Rank);

                    }
                    Kmeans.ClusterRadius[CenterIndex] = 0.0;
                    Kmeans.ClusterWidth[CenterIndex] = 0.0;
                    System.arraycopy(LastClusterCenter[CenterIndex], 0, Kmeans.ClusterCenter[CenterIndex], 0,
                            Kmeans.ParameterVectorDimension);
                }

            } // End Sequential Center Loop
        }

        if (begin) {
            return;
        }
        Kmeans.AverageCenterChange = Kmeans.AverageCenterChange / Kmeans.Ncent_Global;
        Kmeans.AverageRadius = Kmeans.AverageRadius / Kmeans.Ncent_Global;
        Kmeans.AverageWidth = Kmeans.AverageWidth / DAVectorUtility.PointCount_Global;

    } // End FindClusterCenters(int[] NearestCentertoPoint, double[][] LastClusterCenter)

    public static void ReadDataFromFile(String fname, int ClusterPosition, int FirstClustervalue,
                                        int StartPointPosition) {
        int FirstPointPosition = 0;
        int TotalNumberPointstoRead = 0;
        FirstPointPosition = DAVectorUtility.PointStart_Process;
        TotalNumberPointstoRead = DAVectorUtility.PointCount_Process;
        java.util.Random RandomObject = new java.util.Random(10101010 + DAVectorUtility.MPI_Rank);
        if (ClusterPosition < 0) {
            DAVectorUtility
                    .SALSAPrint(0, "Random Start 10101010 plus rank ******************* Option " + ClusterPosition);
        }
        int MinSplitSize = ClusterPosition + 1;
        if (StartPointPosition >= 0) {
            MinSplitSize = Math.max(MinSplitSize, StartPointPosition + Kmeans.ParameterVectorDimension);
        } else {
            DAVectorUtility.printAndThrowRuntimeException(
                    "Illegal Start Position on Points file " + fname + " Rank " + DAVectorUtility.MPI_Rank +
                            " POsition " + StartPointPosition + " Number to Read " +
                            TotalNumberPointstoRead);

        }

        boolean success = false;
        String line = " Unset";
        int countLinesinFile = 0;

        if (Strings.isNullOrEmpty(fname)) {
            DAVectorUtility.printAndThrowRuntimeException(new IllegalArgumentException(Constants.ERR_EMPTY_FILE_NAME));
        }

        try (BufferedReader br = Files.newBufferedReader(Paths.get(fname), Charset.defaultCharset())) {
            Pattern pattern = Pattern.compile("[\t ,]");
            while ((line = br.readLine()) != null) {
                if (Strings.isNullOrEmpty(line)) {
                    continue; // continue on empty lines - "while" will break on null anyway;
                }

                String[] splits = pattern.split(line.trim());
                if (splits.length < MinSplitSize) {
                    DAVectorUtility.SALSAPrint(0, "Count " + countLinesinFile + " Illegal data length on Point file " +
                            splits.length + " " + MinSplitSize + " " + line);
                    continue;
                } // Skip header lines

                Double junk = Doubles.tryParse(splits[StartPointPosition]);
                if (junk == null) continue; // Skip header lines

                if (countLinesinFile < FirstPointPosition) {
                    countLinesinFile += 1;
                    continue;
                }

                int ActualPointPosition = countLinesinFile - FirstPointPosition;
                Integer label = 0;

                Kmeans.PointPosition[ActualPointPosition][0] = Double.parseDouble(splits[StartPointPosition]);
                Kmeans.PointPosition[ActualPointPosition][1] = Double.parseDouble(splits[StartPointPosition + 1]);
                if (Kmeans.ParameterVectorDimension > 2) {
                    for (int VectorIndex = 2; VectorIndex < Kmeans.ParameterVectorDimension; VectorIndex++) {
                        Kmeans.PointPosition[ActualPointPosition][VectorIndex] =
                                Double.parseDouble(splits[VectorIndex + StartPointPosition]);
                    }
                }

                if (ClusterPosition >= 0) {
                    label = Ints.tryParse(splits[ClusterPosition]);
                    if (label == null) {
                        label = FirstClustervalue;
                    }
                    Kmeans.InitialPointAssignment[ActualPointPosition] = label - FirstClustervalue;
                } else {
                    Kmeans.InitialPointAssignment[ActualPointPosition] =
                            RandomObject.nextInt(Program.InitialNcent);
                    if (ClusterPosition == -2) { // Force each cluster to have one point
                        if (countLinesinFile < Program.InitialNcent) {
                            Kmeans.InitialPointAssignment[ActualPointPosition] = countLinesinFile;
                        }
                    }
                    if (ClusterPosition == -3) {
                        int divisor = Program.NumberDataPoints / Program.InitialNcent;
                        if (countLinesinFile % divisor == 0) {
                            Kmeans.InitialPointAssignment[ActualPointPosition] = countLinesinFile / divisor;
                        }
                    }
                    if (ClusterPosition == -4) {
                        int divisor = Program.NumberDataPoints / Program.InitialNcent;
                        Kmeans.InitialPointAssignment[ActualPointPosition] = countLinesinFile / divisor;
                    }
                }
                ++ActualPointPosition;
                ++countLinesinFile;
                if (countLinesinFile >= (FirstPointPosition + TotalNumberPointstoRead)) {
                    break;
                }
            }
            if (countLinesinFile != (FirstPointPosition + TotalNumberPointstoRead)) {
                DAVectorUtility.printAndThrowRuntimeException(
                        "Illegal count on Points file " + fname + " Rank " + DAVectorUtility.MPI_Rank +
                                " Lines in File " + countLinesinFile + " Number to Read " + TotalNumberPointstoRead);

            }
            success = true;
            br.close();
        } catch (IOException e) {
            System.err.format("Failed reading Points data " + DAVectorUtility.MPI_Rank + " " + countLinesinFile +
                                      " Start " + FirstPointPosition + " Number " + TotalNumberPointstoRead + " " +
                                      line + ": %s%n", e);

        }
        if (!success) {
            DAVectorUtility.printAndThrowRuntimeException("DA Vector File read error " + fname);

        }
    } // End ReadDataFromFile

} //  End Kmeans
// End Namespace edu.indiana.soic.spidal.davs
