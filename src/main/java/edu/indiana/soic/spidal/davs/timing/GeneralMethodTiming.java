package edu.indiana.soic.spidal.davs.timing;

import com.google.common.base.Stopwatch;
import edu.indiana.soic.spidal.davs.DAVectorUtility;
import mpi.MPIException;

import java.nio.LongBuffer;
import java.util.concurrent.TimeUnit;

public class GeneralMethodTiming {
    public static enum TimingTask{
        SET_CLUSTER_STATISTICS, SET_POINT_STATISTICS, SET_NEARBY_CLUSTERS, OUTPUT_STATUS,SETUP, CLUSTER_COMPARISON, EXPERIMENT_ANALYSIS,
        SET_WIDTH
    }

    private static Stopwatch timerSetClusterStatistics = Stopwatch.createUnstarted();
    private static Stopwatch timerOutputStatus = Stopwatch.createUnstarted();
    private static Stopwatch timerSetup = Stopwatch.createUnstarted();
    private static Stopwatch timerClusterComparison = Stopwatch.createUnstarted();
    private static Stopwatch timerExperimentAnalysis = Stopwatch.createUnstarted();
    private static Stopwatch timerSetPointStatistics = Stopwatch.createUnstarted();
    private static Stopwatch timerSetNearbyClusters = Stopwatch.createUnstarted();
    private static Stopwatch timerSetWidth = Stopwatch.createUnstarted();

    private static long tSetClusterStatistics;
    private static long tOutputStatus;
    private static long tSetup;
    private static long tClusterComparison;
    private static long tExperimentAnalysis;
    private static long tSetPointStatistics;
    private static long tSetNearbyClusters;
    private static long tSetWidth;

    private static long countSetPointStatistics;
    private static long countSetClusterStatistics;
    private static long countSetNearbyClusters;
    private static long countOutputStatus;
    private static long countSetup;
    private static long countClusterComparison;
    private static long countExperimentAnalysis;
    private static long countSetWidth;

    public static void startTiming(TimingTask task){
        switch (task){
            case SET_CLUSTER_STATISTICS:
                timerSetClusterStatistics.start();
                ++countSetClusterStatistics;
                break;
            case SET_POINT_STATISTICS:
                timerSetPointStatistics.start();
                ++countSetPointStatistics;
                break;
            case SET_NEARBY_CLUSTERS:
                timerSetNearbyClusters.start();
                ++countSetNearbyClusters;
                break;
            case OUTPUT_STATUS:
                timerOutputStatus.start();
                ++countOutputStatus;
                break;
            case SETUP:
                timerSetup.start();
                ++countSetup;
                break;
            case CLUSTER_COMPARISON:
                timerClusterComparison.start();
                ++countClusterComparison;
                break;
            case EXPERIMENT_ANALYSIS:
                timerExperimentAnalysis.start();
                ++countExperimentAnalysis;
                break;
            case SET_WIDTH:
                timerSetWidth.start();
                ++countSetWidth;
                break;

        }
    }

    public static void endTiming(TimingTask task){
        switch (task){
            case SET_CLUSTER_STATISTICS:
                timerSetClusterStatistics.stop();
                tSetClusterStatistics += timerSetClusterStatistics.elapsed(TimeUnit.MILLISECONDS);
                timerSetClusterStatistics.reset();
                break;
            case SET_POINT_STATISTICS:
                timerSetPointStatistics.stop();
                tSetPointStatistics += timerSetPointStatistics.elapsed(TimeUnit.MILLISECONDS);
                timerSetPointStatistics.reset();
                break;
            case SET_NEARBY_CLUSTERS:
                timerSetNearbyClusters.stop();
                tSetNearbyClusters += timerSetNearbyClusters.elapsed(TimeUnit.MILLISECONDS);
                timerSetNearbyClusters.reset();
                break;
            case OUTPUT_STATUS:
                timerOutputStatus.stop();
                tOutputStatus += timerOutputStatus.elapsed(TimeUnit.MILLISECONDS);
                timerOutputStatus.reset();
                break;
            case SETUP:
                timerSetup.stop();
                tSetup += timerSetup.elapsed(TimeUnit.MILLISECONDS);
                timerSetup.reset();
                break;
            case CLUSTER_COMPARISON:
                timerClusterComparison.stop();
                tClusterComparison += timerClusterComparison.elapsed(TimeUnit.MILLISECONDS);
                timerClusterComparison.reset();
                break;
            case EXPERIMENT_ANALYSIS:
                timerExperimentAnalysis.stop();
                tExperimentAnalysis += timerExperimentAnalysis.elapsed(TimeUnit.MILLISECONDS);
                timerExperimentAnalysis.reset();
                break;
            case SET_WIDTH:
                timerSetWidth.stop();
                tSetWidth += timerSetWidth.elapsed(TimeUnit.MILLISECONDS);
                timerSetWidth.reset();
                break;
        }
    }

    public static long getTotalTime(TimingTask task){
        switch (task){
            case SET_CLUSTER_STATISTICS:
                return tSetClusterStatistics;
            case SET_POINT_STATISTICS:
                return tSetPointStatistics;
            case SET_NEARBY_CLUSTERS:
                return tSetNearbyClusters;
            case OUTPUT_STATUS:
                return tOutputStatus;
            case SETUP:
                return tSetup;
            case CLUSTER_COMPARISON:
                return tClusterComparison;
            case EXPERIMENT_ANALYSIS:
                return tExperimentAnalysis;
            case SET_WIDTH:
                return tSetWidth;
        }
        return  0;
    }

    public static double getAverageTime(TimingTask task){
        switch (task){
            case SET_CLUSTER_STATISTICS:
                return tSetClusterStatistics *1.0/ countSetClusterStatistics;
            case SET_POINT_STATISTICS:
                return tSetPointStatistics *1.0/ countSetPointStatistics;
            case SET_NEARBY_CLUSTERS:
                return tSetNearbyClusters *1.0/ countSetNearbyClusters;
            case OUTPUT_STATUS:
                return tOutputStatus *1.0/ countOutputStatus;
            case SETUP:
                return tSetup *1.0/ countSetup;
            case CLUSTER_COMPARISON:
                return tClusterComparison *1.0/ countClusterComparison;
            case EXPERIMENT_ANALYSIS:
                return tExperimentAnalysis *1.0/ countExperimentAnalysis;
            case SET_WIDTH:
                return tSetWidth *1.0/ countSetWidth;
        }
        return  0.0;
    }

    public static long[] getCountDistribution(TimingTask task) throws MPIException{
        LongBuffer mpiOnlyTimingBuffer =  DAVectorUtility.mpiOnlyBuffer;
        mpiOnlyTimingBuffer.position(0);
        switch (task){
            case SET_CLUSTER_STATISTICS:
                mpiOnlyTimingBuffer.put(countSetClusterStatistics);
                break;
            case SET_POINT_STATISTICS:
                mpiOnlyTimingBuffer.put(countSetPointStatistics);
                break;
            case SET_NEARBY_CLUSTERS:
                mpiOnlyTimingBuffer.put(countSetNearbyClusters);
                break;
            case OUTPUT_STATUS:
                mpiOnlyTimingBuffer.put(countOutputStatus);
                break;
            case SETUP:
                mpiOnlyTimingBuffer.put(countSetup);
                break;
            case CLUSTER_COMPARISON:
                mpiOnlyTimingBuffer.put(countClusterComparison);
                break;
            case EXPERIMENT_ANALYSIS:
                mpiOnlyTimingBuffer.put(countExperimentAnalysis);
                break;
            case SET_WIDTH:
                mpiOnlyTimingBuffer.put(countSetWidth);
                break;
        }
        long [] mpiOnlyTimingArray = new long[DAVectorUtility.MPI_Size];
        DAVectorUtility.gather(mpiOnlyTimingBuffer, 1, 0);
        mpiOnlyTimingBuffer.position(0);
        mpiOnlyTimingBuffer.get(mpiOnlyTimingArray);
        return mpiOnlyTimingArray;
    }

}
