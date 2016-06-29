package edu.indiana.soic.spidal.davs.timing;

import com.google.common.base.Stopwatch;
import edu.indiana.soic.spidal.davs.DAVectorUtility;
import mpi.MPIException;

import java.nio.LongBuffer;
import java.util.concurrent.TimeUnit;

public class SectionTiming {
    public static enum TimingTask{
        SC1,SC2,SC3
    }

    private static Stopwatch timerSection1 = Stopwatch.createUnstarted();
    private static Stopwatch timerSection2 = Stopwatch.createUnstarted();
    private static Stopwatch timerSection3 = Stopwatch.createUnstarted();

    private static long tSection1;
    private static long tSection2;
    private static long tSection3;

    private static long countSection1;
    private static long countSection2;
    private static long countSection3;

    public static void startTiming(TimingTask task){
        switch (task){
            case SC1:
                timerSection1.start();
                ++countSection2;
                break;
            case SC2:
                timerSection2.start();
                ++countSection1;
                break;
            case SC3:
                timerSection3.start();
                ++countSection3;
                break;
        }
    }

    public static void endTiming(TimingTask task){
        switch (task){
            case SC1:
                timerSection1.stop();
                tSection1 += timerSection1.elapsed(TimeUnit.MILLISECONDS);
                timerSection1.reset();
                break;
            case SC2:
                timerSection2.stop();
                tSection2 += timerSection2.elapsed(TimeUnit.MILLISECONDS);
                timerSection2.reset();
                break;
            case SC3:
                timerSection3.stop();
                tSection3 += timerSection3.elapsed(TimeUnit.MILLISECONDS);
                timerSection3.reset();
                break;
        }
    }

    public static long getTotalTime(TimingTask task){
        switch (task){
            case SC1:
                return tSection1;
            case SC2:
                return tSection2;
            case SC3:
                return tSection3;
        }
        return  0;
    }

    public static double getAverageTime(TimingTask task){
        switch (task){
            case SC1:
                return tSection1 *1.0/ countSection2;
            case SC2:
                return tSection2 *1.0/ countSection1;
            case SC3:
                return tSection3 *1.0/ countSection3;
        }
        return  0.0;
    }

    public static long[] getCountDistribution(TimingTask task) throws MPIException{
        LongBuffer mpiOnlyTimingBuffer =  DAVectorUtility.mpiOnlyBuffer;
        mpiOnlyTimingBuffer.position(0);
        switch (task){
            case SC1:
                mpiOnlyTimingBuffer.put(countSection2);
                break;
            case SC2:
                mpiOnlyTimingBuffer.put(countSection1);
                break;
            case SC3:
                mpiOnlyTimingBuffer.put(countSection3);
                break;
        }
        long [] mpiOnlyTimingArray = new long[DAVectorUtility.MPI_Size];
        DAVectorUtility.gather(mpiOnlyTimingBuffer, 1, 0);
        mpiOnlyTimingBuffer.position(0);
        mpiOnlyTimingBuffer.get(mpiOnlyTimingArray);
        return mpiOnlyTimingArray;
    }

}
