package edu.indiana.soic.spidal.davs;

import com.google.common.base.Stopwatch;
import com.google.common.base.Strings;
import mpi.Intracomm;
import mpi.MPIException;
import edu.indiana.soic.spidal.general.Box;
import edu.indiana.soic.spidal.mpi.MpiOps;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.Charset;
import java.nio.file.*;
import java.util.ArrayList;
import java.util.concurrent.TimeUnit;

public class DAVectorUtility
{
	public static int PointCount_Global = 0; // Total number of points summed over all threads and processes
	public static int PointCount_Process = 0; // Total number of points summed over all threads in this process
	public static int PointCount_Largest = 0; // Largest number of points in all processes
	public static int PointStart_Process = 0; //    First data point in this process

	// Parallel Parameters
	public static int MPI_Rank = 0; // Rank of process
	public static int MPI_Size = 1; // Number of MPI Processes

	public static Intracomm MPI_communicator = null; //MPI communicator

    public static MpiOps mpiOps = null; // MPI operations

	public static String ParallelPattern = " "; // Pattern of parallel execution (e.g. 8x1x8 indicates [threads/process][processes/node][nodes])
	public static String PatternLabel = " "; // Title line for print

	//  Within a job data points will be divided into MPI_Size parts -- each part is assigned to a separate MPI Process
	public static int[] PointsperProcess = null; //how many data points each process will take care
	public static int[][] PointsperThreadperProcess = null; // Number of data points in each process-thread

	//	Within a process, data points will be divided into ThreadCount segments( each thread a segment), the array keep the size of each segment
	public static int[] PointsperThread = null; //how many data points each thread will take care
	public static int[] StartPointperThread = null; //the starting point that a thread will take care

	public static int ThreadCount = 1; // maximum number of parallel threads in a process
	public static int NodeCount = 1; // maximum number of separate nodes in run
	public static int MPIperNodeCount = 1; // Number of MPI processes per node

    // Timing Parameters
    public static Stopwatch mainTimer;
	public static Stopwatch PreciseTimer; //    Hold Precise Timing
	public static int NumberofSubTimings = 0; // Number of subtimings
	public static Stopwatch[] SubTimers; // Timing Objects
	public static long mainDuration = 0; // duration measured by mainTimer in milliseconds
	public static double HPDuration = 0.0; // Time with Precision measured in microseconds
	public static double[] SubDurations; // Hold partial timing
	public static String[] SubTimingNames; //  Labels of partial timing
	public static int[] SubTimingCalls; // Number of calls
	public static boolean[] SubTimingEnable;
	public static int[] TimingOutputOrder; // Order of Output
	public static int MPIREDUCETiming1 = -1;
	public static int MPIREDUCETiming2 = -1;
	public static int MPIREDUCETiming3 = -1;
	public static int MPIREDUCETiming4 = -1;
	public static int MPIREDUCETiming5 = -1;
	public static int MPIREDUCETiming6 = -1;
	public static int MPISENDRECEIVETiming = -1;
	public static int MPIGATHERTiming = -1;
	public static int MPIDistributedREDUCETiming = -1;
	public static int MPIBROADCASTTiming = -1;
	public static int MPISynchTiming = -1;
	public static int ThreadTiming = -1; // viewed as MPI timing

	//  These are general parameters for C# codes
	public static ArrayList<String> CosmicOutput = new ArrayList<>(1000); // Monitoring Output
	public static boolean ConsoleDebugOutput = true; // If true send Monitoring output to console
	public static int DebugPrintOption = 2; // Control Printing (= 0 None, ==1 Summary, = 2 Full)

	public static java.util.ArrayList TemperatureValues = new java.util.ArrayList(10000); // Record Temperatures per step
	public static java.util.ArrayList ClusterCountValues = new java.util.ArrayList(10000); // Record Cluster Counts per Step

    public static void printException(Exception e){
        System.out.println("SALSA Error " + e.getMessage());
    }
    public static void printAndThrowRuntimeException(RuntimeException e){
        System.out.println("SALSA Error " + e.getMessage());
        throw e;
    }
	public static void printAndThrowRuntimeException(String message) {
		System.out.println("SALSA Error " + message);
		throw new RuntimeException(message);
	} // end printAndThrowRuntimeException

	// PrintOption = 0 Essential Printout
	// PrintOption = 1 Summary Printout
	// PrintOption = 2 Only if full print out requested
	public static void SALSAPrint(int PrintOption, String StufftoPrint)
	{
		if (MPI_Rank != 0)
		{
			return;
		}
		if (DebugPrintOption < PrintOption)
		{
			return;
		}
		CosmicOutput.add(StufftoPrint);

		if (ConsoleDebugOutput)
		{
			System.out.println(StufftoPrint);
		}

    } // End SALSAPrint

	public static void SALSAFullPrint(int PrintOption, String StufftoPrint)
	{
		if (DebugPrintOption < PrintOption)
		{
			return;
		}
		if (MPI_Rank == 0)
		{
			CosmicOutput.add(StufftoPrint);
		}

		if (ConsoleDebugOutput)
		{
			System.out.println(" Node:" + MPI_Rank + " " + StufftoPrint);
		}

    } // End SALSAFullPrint

	public static void SALSASyncPrint(int PrintOption, String GlobalStufftoPrint, String StufftoPrint) throws MPIException {
		if (DebugPrintOption < PrintOption)
		{
			return;
		}
		String TotalStufftoPrint = "";
		if (StufftoPrint.length() > 0)
		{
			TotalStufftoPrint = " Node:" + MPI_Rank + " " + StufftoPrint;
		}
		DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming3);
        // Note - MPI Call - Allreduce - String - Done with Allgather followed by local concatenation
        if (DAVectorUtility.MPI_Size > 1){
//            TotalStufftoPrint = DAVectorUtility.MPI_communicator.<String>Allreduce(TotalStufftoPrint, Operation<String>.Add);
            TotalStufftoPrint = DAVectorUtility.mpiOps.allReduce(TotalStufftoPrint);
        }
		DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming3);
		TotalStufftoPrint = GlobalStufftoPrint + TotalStufftoPrint;
		if (MPI_Rank != 0)
		{
			return;
		}
		CosmicOutput.add(TotalStufftoPrint);

		if (ConsoleDebugOutput)
		{
			System.out.println(TotalStufftoPrint);
		}

    } // End SALSASyncPrint

	public static String PrintFixedInteger(int data, int digits)
	{
        String returned = String.valueOf(data);
		for (int charloop = returned.length(); charloop < digits; charloop++)
		{
			returned = " " + returned;
		}
		return returned;

	} // End PrintFixedInteger

	public static String PadString(String start, int digits)
	{
		String returned = start;
		for (int charloop = start.length(); charloop < digits; charloop++)
		{
			returned += " ";
		}
		return returned;
	} // End padstring

	public static String LeftPrintFixedInteger(int data, int digits)
	{
        String returned = String.valueOf(data);
		for (int charloop = returned.length(); charloop < digits; charloop++)
		{
			returned = returned + " ";
		}
		return returned;

	} // End LeftPrintFixedInteger

    public static String formatElapsedMillis(long elapsed){
        String format = "%dd:%02dH:%02dM:%02dS:%03dmS";
        short millis = (short)(elapsed % (1000.0));
        elapsed = (elapsed - millis) / 1000; // remaining elapsed in seconds
        byte seconds = (byte)(elapsed % 60.0);
        elapsed = (elapsed - seconds) / 60; // remaining elapsed in minutes
        byte minutes =  (byte)(elapsed % 60.0);
        elapsed = (elapsed - minutes) / 60; // remaining elapsed in hours
        byte hours = (byte)(elapsed % 24.0);
        long days = (elapsed - hours) / 24; // remaining elapsed in days
        return String.format(format, days, hours, minutes,  seconds, millis);
    }


    public static void InitializeTiming(int InputNumberofTimers)
	{
		NumberofSubTimings = InputNumberofTimers;
		SubTimers = new Stopwatch[NumberofSubTimings]; // Timing Objects
		SubDurations = new double[NumberofSubTimings]; // Hold partial timing
		SubTimingEnable = new boolean[NumberofSubTimings];
		SubTimingNames = new String[NumberofSubTimings];
		SubTimingCalls = new int[NumberofSubTimings];
		TimingOutputOrder = new int[NumberofSubTimings];

		for (int itimer = 0; itimer < NumberofSubTimings; itimer++)
		{
			SubTimers[itimer] = Stopwatch.createUnstarted();
			SubDurations[itimer] = 0.0;
			SubTimingEnable[itimer] = true;
			SubTimingCalls[itimer] = 0;
			TimingOutputOrder[itimer] = itimer;
		}
		PreciseTimer = Stopwatch.createStarted();

        mainTimer = Stopwatch.createStarted(); // Start the main timer

	} // End InitializeTiming

	public static void SetUpSubTimer(int TimingIndex, String TimingLabel)
	{
		if (TimingIndex >= NumberofSubTimings)
		{
			printAndThrowRuntimeException("Error in Timing Index " + TimingIndex + " Max " + NumberofSubTimings);
			return;
		}
		SubTimingNames[TimingIndex] = TimingLabel;
		SubTimingEnable[TimingIndex] = true;
		SubDurations[TimingIndex] = 0.0;
	} // End SetUpSubTimer

	public static void SetUpMPISubTimers(int StartTimingIndex, String MPISetLabel)
	{
		if ((StartTimingIndex + 12) >= NumberofSubTimings)
		{
			printAndThrowRuntimeException("Error in  MPI Timing Index " + StartTimingIndex + 12 + " Max " + NumberofSubTimings);
			return;
		}
		SetUpSubTimer(StartTimingIndex, MPISetLabel + "MPI Reduce1(GlobalReductions)");
		SetUpSubTimer(StartTimingIndex + 1, MPISetLabel + "MPI Reduce1(GlobalReductions while Dist)");
		SetUpSubTimer(StartTimingIndex + 2, MPISetLabel + "MPI Reduce2(DAVectorClustering)");
		SetUpSubTimer(StartTimingIndex + 3, MPISetLabel + "MPI Reduce3(DAVectorDist)");
		SetUpSubTimer(StartTimingIndex + 4, MPISetLabel + "MPI Reduce4(Dist Reduce Initial)");
		SetUpSubTimer(StartTimingIndex + 5, MPISetLabel + "MPI Reduce5(Dist Exec Control)");
		SetUpSubTimer(StartTimingIndex + 6, MPISetLabel + "MPI Reduce6(Dist Reduce Transport)");
		SetUpSubTimer(StartTimingIndex + 7, MPISetLabel + "MPI DistributedReduce");
		SetUpSubTimer(StartTimingIndex + 8, MPISetLabel + "MPI Gather");
		SetUpSubTimer(StartTimingIndex + 9, MPISetLabel + "MPI Bcast");
		SetUpSubTimer(StartTimingIndex + 10, MPISetLabel + "MPI SendRecv");
		SetUpSubTimer(StartTimingIndex + 11, MPISetLabel + "MPI Sync");
		SetUpSubTimer(StartTimingIndex + 12, MPISetLabel + "Thread Timing");

		if ((!MPISetLabel.equals("")) && (!MPISetLabel.equals("Lib ")))
		{
			return;
		}
		MPIREDUCETiming1 = StartTimingIndex; // StartTimingIndex+1 also used for different calls in Dist
		MPIREDUCETiming2 = StartTimingIndex + 2;
		MPIREDUCETiming3 = StartTimingIndex + 3;
		MPIREDUCETiming4 = StartTimingIndex + 4;
		MPIREDUCETiming5 = StartTimingIndex + 5;
		MPIREDUCETiming6 = StartTimingIndex + 6;
		MPIBROADCASTTiming = StartTimingIndex + 9;
		MPIGATHERTiming = StartTimingIndex + 8;
		MPIDistributedREDUCETiming = StartTimingIndex + 7;
		MPISENDRECEIVETiming = StartTimingIndex + 10;
		MPISynchTiming = StartTimingIndex + 11;
		ThreadTiming = StartTimingIndex + 12;

	} // End SetUpMPISubTimers

	public static int MPIaddviaGather(int thisnodevalue) throws MPIException {
        if (DAVectorUtility.MPI_Size > 1){
            // Note - MPI Call - Allgather - int
//            Nodevalues = DAVectorUtility.MPI_communicator.<Integer>Allgather(thisnodevalue);
            int [] Nodevalues = DAVectorUtility.mpiOps.allGather(thisnodevalue);
            int maxval = 0;
            for (int nodes = 0; nodes < DAVectorUtility.MPI_Size; nodes++)
            {
                maxval += Nodevalues[nodes];
            }
            return maxval;
        } else {
            return thisnodevalue;
        }
	}


	public static void InterimTiming()
	{
		PreciseTimer.stop();
		HPDuration += PreciseTimer.elapsed(TimeUnit.MICROSECONDS);
		PreciseTimer.reset();PreciseTimer.start();

	} // end InterimTiming

	public static void EndTiming()
	{
		mainTimer.stop();
		PreciseTimer.stop();
		HPDuration += PreciseTimer.elapsed(TimeUnit.MICROSECONDS);
        mainDuration = mainTimer.elapsed(TimeUnit.MILLISECONDS);
        mainTimer.reset();
        PreciseTimer.reset();

	} // end EndTiming

	public static void StartSubTimer(int TimingIndex)
	{
		if (TimingIndex < 0)
		{
			return;
		}

		if (SubTimingEnable[TimingIndex])
		{
			SubTimers[TimingIndex].start();
			++SubTimingCalls[TimingIndex];
		}

	} // End StartSubTimer

	public static void StopSubTimer(int TimingIndex)
	{
		if (TimingIndex < 0)
		{
			return;
		}

		if (SubTimingEnable[TimingIndex])
		{
            SubTimers[TimingIndex].stop();
			SubDurations[TimingIndex] += SubTimers[TimingIndex].elapsed(TimeUnit.MICROSECONDS);
            SubTimers[TimingIndex].reset();
		}

	} // End StopSubTimer

	public static void CopyVector(int[] VectorC, int[] VectorA, int StartIndex, int TotalSize)
	{ // Can be parallelized
        System.arraycopy(VectorA, 0, VectorC, StartIndex, TotalSize);
	}

    public static double SynchronizeMPIDouble(double sync) throws MPIException {
        if (DAVectorUtility.MPI_Size > 1)
        {
            DAVectorUtility.StartSubTimer(DAVectorUtility.MPISynchTiming);
            // Note - MPI Call - Broadcast - double
            sync = DAVectorUtility.mpiOps.broadcast(sync,0);
            DAVectorUtility.StopSubTimer(DAVectorUtility.MPISynchTiming);
        }
        return sync;
    }

    public static int SynchronizeMPIInteger(int sync) throws MPIException {
        if (DAVectorUtility.MPI_Size > 1)
        {
            DAVectorUtility.StartSubTimer(DAVectorUtility.MPISynchTiming);
            // Note - MPI Call - Broadcast - int
            sync = DAVectorUtility.mpiOps.broadcast(sync,0);
            DAVectorUtility.StopSubTimer(DAVectorUtility.MPISynchTiming);
        }
        return sync;
    }

    public static boolean SynchronizeMPIBoolean(boolean sync) throws MPIException {
        if (DAVectorUtility.MPI_Size > 1)
        {
            DAVectorUtility.StartSubTimer(DAVectorUtility.MPISynchTiming);
            // Note - MPI Call - Broadcast - boolean
            sync = DAVectorUtility.mpiOps.broadcast(sync,0);
            DAVectorUtility.StopSubTimer(DAVectorUtility.MPISynchTiming);
        }
        return sync;
    }


	public static void SetUpHistogramRange(int PointsinHistogram, Box<Double> Histmin, Box<Double> Histmax)
	{ // Choose good range to get rounded labels

		if (Histmax.content <= Histmin.content)
		{
			return;
		}
		if (Math.abs(Histmin.content) < 0.000000001)
		{
			Histmin.content = 0.0;
		}
		else
		{
			double newminvalue = NewBinSize(Math.abs(Histmin.content));
			if (Histmin.content < 0.0)
			{
				Histmin.content = -newminvalue;
			}
			else
			{
				Histmin.content = newminvalue;
			}
		}
		double binsize = (Histmax.content - Histmin.content) / PointsinHistogram;
		Histmax.content = Histmin.content + PointsinHistogram * NewBinSize(binsize);

	} // End SetUpHistogramRange

	public static double NewBinSize(double oldBinSize)
	{ // Round Bin size up to a pretty value
		if (oldBinSize <= 0.000000001 || oldBinSize > 10000000000.0)
		{
			return oldBinSize;
		}

		double logvalue = Math.log10(oldBinSize);
		int intlogvalue = (int)Math.floor(logvalue);
		double fudgepower = 1.0 - (double)intlogvalue;
		double fudge = Math.pow(10.0, fudgepower);
		double scaled = fudge * oldBinSize;
		scaled = Math.min(scaled, 100.0);
		double Intversionofnewbinsize = Math.ceil(scaled) / fudge;
		//            SALSAUtility.SALSAPrint(1, "Hist " + oldBinSize.ToString("F4") + " " + intlogvalue + " " + Intversionofnewbinsize.ToString("F4") + " scaled " + scaled.ToString("F2"));
		return Intversionofnewbinsize;
	}

	public static void writeClusterResults(String file, ArrayList<String> lines)
	{
        Path filePath = Paths.get(file);
        OpenOption mode = Files.exists(filePath) ? StandardOpenOption.APPEND : StandardOpenOption.CREATE;
        try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(filePath, Charset.defaultCharset(),
                mode), true)) {
            if (Strings.isNullOrEmpty(file)) {
                DAVectorUtility.printAndThrowRuntimeException(new IllegalArgumentException(Constants
                        .ERR_EMPTY_FILE_NAME));
            }

            writer.println();
            for (String line : lines) {
                writer.println(line);
            }
            writer.close();
        } catch (IOException e) {
            System.err.format("Failed writing cluster results due to I/O exception: %s%n", e);
        }
    } // End writeClusterResults

	public static void writeExperimentalShifts(String file,int experimentID, double MZshift, double MZSD, double RTshift, double RTSD){
		Path filePath = Paths.get(file);
		OpenOption mode = Files.exists(filePath) ? StandardOpenOption.APPEND : StandardOpenOption.CREATE;


		try(PrintWriter writer = new PrintWriter(Files.newBufferedWriter(filePath,Charset.defaultCharset(),mode),true)) {

			if(mode == StandardOpenOption.CREATE){
				writer.println("Experiment\tMZshift\tMZSD\tRTshift\tRTSD");
			}

			writer.println(experimentID + "\t" + MZshift + "\t" + MZSD + "\t" + RTshift + "\t" + RTSD);
			writer.flush();
			writer.close();
		} catch (IOException e) {
			System.err.format("Failed writing Experiment Shifts results due to I/O exception: %s%n", e);
		}
	}
} // End DAVectorUtility
 // End Namespace SALSALibrary
