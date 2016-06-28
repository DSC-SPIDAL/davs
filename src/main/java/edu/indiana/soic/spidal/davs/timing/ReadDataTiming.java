package edu.indiana.soic.spidal.davs.timing;

import com.google.common.base.Stopwatch;
import edu.indiana.soic.spidal.davs.DAVectorUtility;
import mpi.MPIException;

import java.nio.LongBuffer;
import java.util.concurrent.TimeUnit;

public class ReadDataTiming {
    public static enum TimingTask{
        TOTAL_READ
    }

    private static Stopwatch timerTotalRead = Stopwatch.createUnstarted();

    private static long tTotalRead;

    private static long countTotalRead;

    public static void startTiming(TimingTask task){
        switch (task){
            case TOTAL_READ:
                timerTotalRead.start();
                ++countTotalRead;
                break;
        }
    }

    public static void endTiming(TimingTask task){
        switch (task){
            case TOTAL_READ:
                timerTotalRead.stop();
                tTotalRead += timerTotalRead.elapsed(TimeUnit.MILLISECONDS);
                timerTotalRead.reset();
                break;
        }
    }

    public static double getTotalTime(TimingTask task){
        switch (task){
            case TOTAL_READ:
                return tTotalRead;
        }
        return  0.0;
    }

    public static double getAverageTime(TimingTask task){
        switch (task){
            case TOTAL_READ:
                return tTotalRead *1.0/ countTotalRead;
        }
        return  0.0;
    }

    public static long[] getCountDistribution(TimingTask task) throws MPIException{
        LongBuffer mpiOnlyTimingBuffer =  DAVectorUtility.mpiOnlyBuffer;
        mpiOnlyTimingBuffer.position(0);
        switch (task){
            case TOTAL_READ:
                mpiOnlyTimingBuffer.put(countTotalRead);
                break;
        }
        long [] mpiOnlyTimingArray = new long[DAVectorUtility.MPI_Size];
        DAVectorUtility.gather(mpiOnlyTimingBuffer, 1, 0);
        mpiOnlyTimingBuffer.position(0);
        mpiOnlyTimingBuffer.get(mpiOnlyTimingArray);
        return mpiOnlyTimingArray;
    }

}
