package edu.indiana.soic.spidal.davs.timing;

import com.google.common.base.Stopwatch;
import edu.indiana.soic.spidal.davs.DAVectorUtility;
import mpi.MPIException;

import java.nio.LongBuffer;
import java.util.concurrent.TimeUnit;

public class SetupTiming {
    public static enum TimingTask{
        SETUP
    }

    private static Stopwatch timerSetup = Stopwatch.createUnstarted();

    private static long tSetup;

    private static long countSetup;

    public static void startTiming(TimingTask task){
        switch (task){
            case SETUP:
                timerSetup.start();
                ++countSetup;
                break;
        }
    }

    public static void endTiming(TimingTask task){
        switch (task){
            case SETUP:
                timerSetup.stop();
                tSetup += timerSetup.elapsed(TimeUnit.MILLISECONDS);
                timerSetup.reset();
                break;
        }
    }

    public static double getTotalTime(TimingTask task){
        switch (task){
            case SETUP:
                return tSetup;
        }
        return  0.0;
    }

    public static double getAverageTime(TimingTask task){
        switch (task){
            case SETUP:
                return tSetup *1.0/ countSetup;
        }
        return  0.0;
    }

    public static long[] getCountDistribution(TimingTask task) throws MPIException{
        LongBuffer mpiOnlyTimingBuffer =  DAVectorUtility.mpiOnlyBuffer;
        mpiOnlyTimingBuffer.position(0);
        switch (task){
            case SETUP:
                mpiOnlyTimingBuffer.put(countSetup);
                break;
        }
        long [] mpiOnlyTimingArray = new long[DAVectorUtility.MPI_Size];
        DAVectorUtility.gather(mpiOnlyTimingBuffer, 1, 0);
        mpiOnlyTimingBuffer.position(0);
        mpiOnlyTimingBuffer.get(mpiOnlyTimingArray);
        return mpiOnlyTimingArray;
    }

}
