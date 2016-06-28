package edu.indiana.soic.spidal.davs.timing;

import com.google.common.base.Stopwatch;
import edu.indiana.soic.spidal.davs.DAVectorUtility;
import mpi.MPIException;

import java.nio.LongBuffer;
import java.util.concurrent.TimeUnit;

public class GeneralTiming {
    public static enum TimingTask{
        KMEANS,LCMS,DA
    }

    private static Stopwatch timerKmeans = Stopwatch.createUnstarted();
    private static Stopwatch timerLCMS = Stopwatch.createUnstarted();
    private static Stopwatch timerDA = Stopwatch.createUnstarted();

    private static long tKmeans;
    private static long tLCMS;
    private static long tDA;

    private static long countLCMS;
    private static long countKmeans;
    private static long countDA;

    public static void startTiming(TimingTask task){
        switch (task){
            case KMEANS:
                timerKmeans.start();
                ++countKmeans;
                break;
            case LCMS:
                timerLCMS.start();
                ++countLCMS;
                break;
            case DA:
                timerDA.start();
                ++countDA;
                break;
        }
    }

    public static void endTiming(TimingTask task){
        switch (task){
            case KMEANS:
                timerKmeans.stop();
                tKmeans += timerKmeans.elapsed(TimeUnit.MILLISECONDS);
                timerKmeans.reset();
                break;
            case LCMS:
                timerLCMS.stop();
                tLCMS += timerLCMS.elapsed(TimeUnit.MILLISECONDS);
                timerLCMS.reset();
                break;
            case DA:
                timerDA.stop();
                tDA += timerDA.elapsed(TimeUnit.MILLISECONDS);
                timerDA.reset();
                break;
        }
    }

    public static double getTotalTime(TimingTask task){
        switch (task){
            case KMEANS:
                return tKmeans;
            case LCMS:
                return tLCMS;
            case DA:
                return tDA;
        }
        return  0.0;
    }

    public static double getAverageTime(TimingTask task){
        switch (task){
            case KMEANS:
                return tKmeans *1.0/ countKmeans;
            case LCMS:
                return tLCMS *1.0/ countLCMS;
            case DA:
                return tDA *1.0/ countDA;
        }
        return  0.0;
    }

    public static long[] getCountDistribution(TimingTask task) throws MPIException{
        LongBuffer mpiOnlyTimingBuffer =  DAVectorUtility.mpiOnlyBuffer;
        mpiOnlyTimingBuffer.position(0);
        switch (task){
            case KMEANS:
                mpiOnlyTimingBuffer.put(countKmeans);
                break;
            case LCMS:
                mpiOnlyTimingBuffer.put(countLCMS);
                break;
            case DA:
                mpiOnlyTimingBuffer.put(countDA);
                break;
        }
        long [] mpiOnlyTimingArray = new long[DAVectorUtility.MPI_Size];
        DAVectorUtility.gather(mpiOnlyTimingBuffer, 1, 0);
        mpiOnlyTimingBuffer.position(0);
        mpiOnlyTimingBuffer.get(mpiOnlyTimingArray);
        return mpiOnlyTimingArray;
    }

}
