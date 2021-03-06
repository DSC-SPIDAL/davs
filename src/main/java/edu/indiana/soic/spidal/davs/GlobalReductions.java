package edu.indiana.soic.spidal.davs;

import edu.indiana.soic.spidal.general.Box;
import edu.indiana.soic.spidal.mpi.MPIReducePlusIndex;
import mpi.MPI;
import mpi.MPIException;

import static edu.rice.hj.Module0.launchHabaneroApp;
import static edu.rice.hj.Module1.forallChunked;


public class GlobalReductions
{

	public static class FindBoolOr
	{
		public double[] NumberofPoints;
		public int NumberofThreads;
		public boolean[] Orvalue;
		public double TotalNumberofPoints;
		public boolean TotalOr;

		public FindBoolOr(int NumThreads)
		{
			NumberofThreads = NumThreads;

			NumberofPoints = new double[NumThreads];
			Orvalue = new boolean[NumThreads];

			TotalNumberofPoints = 0.0;
			TotalOr = false;
			for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
			{
				NumberofPoints[ThreadNo] = 0.0;
				Orvalue[ThreadNo] = false;
			}
		}

		public final void addapoint(int ThreadNo, boolean value)
		{
			NumberofPoints[ThreadNo] += 1.0;
			Orvalue[ThreadNo] = Orvalue[ThreadNo] || value;
		}
		public final void sumoverthreadsandmpi() throws MPIException {
			for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
			{
				TotalNumberofPoints += NumberofPoints[ThreadNo];
				TotalOr = Orvalue[ThreadNo] || TotalOr;
			}


			if (DAVectorUtility.MPI_Size > 1)
			{
				DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming1);
                // Note - MPI Call - Allreduce - double - sum
//				TotalNumberofPoints = DAVectorUtility.MPI_communicator.<Double>Allreduce(TotalNumberofPoints, Operation<Double>.Add);
				TotalNumberofPoints = DAVectorUtility.mpiOps.allReduce(TotalNumberofPoints, MPI.SUM);
                // Note - MPI Call - Allreduce - boolean - logical OR
//				TotalOr = DAVectorUtility.MPI_communicator.<Boolean>Allreduce(TotalOr, Operation<Boolean>.LogicalOr);
				TotalOr = DAVectorUtility.mpiOps.allReduce(TotalOr, MPI.LOR);
				DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming1);
			}
		}
	} // End FindBoolOr

	public static class FindIntSum
	{
		private int[] NumberofPoints;
		private int NumberofThreads;
		private int[] Intvalue;
		public int TotalNumberofPoints;
		public int TotalInt;

		public FindIntSum(int NumThreads)
		{
			NumberofThreads = NumThreads;

			NumberofPoints = new int[NumThreads];
			Intvalue = new int[NumThreads];

			TotalNumberofPoints = 0;
			TotalInt = 0;
			for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
			{
				NumberofPoints[ThreadNo] = 0;
				Intvalue[ThreadNo] = 0;
			}
		}

		public final void addapoint(int ThreadNo, int value)
		{
			NumberofPoints[ThreadNo] += 1;
			Intvalue[ThreadNo] += value;
		}
		public final void sumoverthreadsandmpi() throws MPIException {
			for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
			{
				TotalNumberofPoints += NumberofPoints[ThreadNo];
				TotalInt += Intvalue[ThreadNo];
			}
			if (DAVectorUtility.MPI_Size > 1)
			{
				DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming1);
                // Note - MPI Call - Allreduce - int - sum
//				TotalNumberofPoints = DAVectorUtility.MPI_communicator.<Integer>Allreduce(TotalNumberofPoints, Operation<Integer>.Add);
				TotalNumberofPoints = DAVectorUtility.mpiOps.allReduce(TotalNumberofPoints, MPI.SUM);
                // Note - MPI Call - Allreduce - int - sum
//				TotalInt = DAVectorUtility.MPI_communicator.<Integer>Allreduce(TotalInt, Operation<Integer>.Add);
				TotalInt = DAVectorUtility.mpiOps.allReduce(TotalInt, MPI.SUM);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming1);
			}
        }
	} // End FindIntSum

	public static class FindDoubleArraySum
	{ // Used to do histograms
		// Must call startthread method at start of threads
		public double[] NumberofPoints;
		public int NumberofThreads;
		public int NumberinSum;
		public double[][] Sum;
		public double TotalNumberofPoints;
		public double[] TotalSum;

		public FindDoubleArraySum(int NumThreads, int ArraySize)
		{
			NumberofThreads = NumThreads;
			NumberinSum = ArraySize;

			NumberofPoints = new double[NumThreads];
			TotalSum = new double[ArraySize];
			Sum = new double[NumThreads][];

			TotalNumberofPoints = 0.0;
			for (int loop = 0; loop < ArraySize; loop++)
			{
				TotalSum[loop] = 0.0;
			}
		}
		public final void startthread(int ThreadNo)
		{
			NumberofPoints[ThreadNo] = 0.0;
			Sum[ThreadNo] = new double[NumberinSum];
			for (int loop = 0; loop < NumberinSum; loop++)
			{
				Sum[ThreadNo][loop] = 0.0;
			}
		}

		public final void addapoint(int ThreadNo, int loopvalue)
		{
			if ((loopvalue < 0) || (loopvalue >= NumberinSum))
			{
				return;
			}
			NumberofPoints[ThreadNo] += 1.0;
			Sum[ThreadNo][loopvalue] += 1.0;
		}

		public final void sumoverthreadsandmpi() throws MPIException {
			for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
			{
				TotalNumberofPoints += NumberofPoints[ThreadNo];
				for (int loop = 0; loop < NumberinSum; loop++)
				{
					TotalSum[loop] += Sum[ThreadNo][loop];
				}
			}
			if (DAVectorUtility.MPI_Size > 1)
			{
				DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming1);
                // Note - MPI Call - Allreduce - double - sum
//				TotalNumberofPoints = DAVectorUtility.MPI_communicator.<Double>Allreduce(TotalNumberofPoints, Operation<Double>.Add);
				TotalNumberofPoints = DAVectorUtility.mpiOps.allReduce(TotalNumberofPoints, MPI.SUM);
                // Note - MPI Call - Allreduce - double[] - sum
//				TotalSum = DAVectorUtility.MPI_communicator.<Double>Allreduce(TotalSum, Operation<Double>.Add);
				DAVectorUtility.mpiOps.allReduce(TotalSum, MPI.SUM);

				DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming1);
			}
        }
	} // End FindDoubleArraySum

	public static class FindDoubleMax
	{
		public double[] NumberofPoints;
		public int NumberofThreads;
		public double[] Maxvalue;
		public double TotalNumberofPoints;
		public double TotalMax;

		public FindDoubleMax(int NumThreads)
		{
			NumberofThreads = NumThreads;

			NumberofPoints = new double[NumThreads];
			Maxvalue = new double[NumThreads];

			TotalNumberofPoints = 0.0;
			TotalMax = 0.0;
			for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
			{
				NumberofPoints[ThreadNo] = 0.0;
				Maxvalue[ThreadNo] = 0.0;
			}
		}

		public final void addapoint(int ThreadNo, double value)
		{
			NumberofPoints[ThreadNo] += 1.0;
			Maxvalue[ThreadNo] = Math.max(Maxvalue[ThreadNo], value);
		}
		public final void sumoverthreadsandmpi() throws MPIException {
			for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
			{
				TotalNumberofPoints += NumberofPoints[ThreadNo];
				TotalMax = Math.max(TotalMax, Maxvalue[ThreadNo]);
			}
			if (DAVectorUtility.MPI_Size > 1)
			{
				DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming1);
                // Note - MPI Call - Allreduce - double - sum
//				TotalNumberofPoints = DAVectorUtility.MPI_communicator.<Double>Allreduce(TotalNumberofPoints, Operation<Double>.Add);
				TotalNumberofPoints = DAVectorUtility.mpiOps.allReduce(TotalNumberofPoints, MPI.SUM);
                // Note - MPI Call - Allreduce - double - sum
//				TotalMax = DAVectorUtility.MPI_communicator.<Double>Allreduce(TotalMax, Operation<Double>.Max);
				TotalMax = DAVectorUtility.mpiOps.allReduce(TotalMax, MPI.MAX);
				DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming1);
			}
        }
	} // End FindDoubleMax

	public static class FindMeanSigma
	{
		public double[] NumberofPoints;
		public int NumberofThreads;
		public double[] mean;
		public double[] square;
		public double TotalNumberofPoints;
		public double Totalmean;
		public double Totalsquare;
		public double Totalsigma;

		public FindMeanSigma(int NumThreads)
		{
			NumberofThreads = NumThreads;

			NumberofPoints = new double[NumThreads];
			mean = new double[NumThreads];
			square = new double[NumThreads];

			TotalNumberofPoints = 0.0;
			Totalmean = 0.0;
			Totalsquare = 0.0;
			Totalsigma = 0.0;
			for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
			{
				NumberofPoints[ThreadNo] = 0.0;
				mean[ThreadNo] = 0.0;
				square[ThreadNo] = 0.0;
			}
		}

		public final void addapoint(int ThreadNo, double value1)
		{
			NumberofPoints[ThreadNo] += 1.0;
			mean[ThreadNo] += value1;
			square[ThreadNo] += value1 * value1;
		}
		public final void sumoverthreadsandmpi() throws MPIException {
			for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
			{
				TotalNumberofPoints += NumberofPoints[ThreadNo];
				Totalmean += mean[ThreadNo];
				Totalsquare += square[ThreadNo];
			}
			if (DAVectorUtility.MPI_Size > 1)
			{
				DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming1);
                // Note - MPI Call - Allreduce - double - sum
//				TotalNumberofPoints = DAVectorUtility.MPI_communicator.<Double>Allreduce(TotalNumberofPoints, Operation<Double>.Add);
				TotalNumberofPoints = DAVectorUtility.mpiOps.allReduce(TotalNumberofPoints, MPI.SUM);
                // Note - MPI Call - Allreduce - double - sum
//				Totalmean = DAVectorUtility.MPI_communicator.<Double>Allreduce(Totalmean, Operation<Double>.Add);
				Totalmean = DAVectorUtility.mpiOps.allReduce(Totalmean, MPI.SUM);
                // Note - MPI Call - Allreduce - double - sum
//				Totalsquare = DAVectorUtility.MPI_communicator.<Double>Allreduce(Totalsquare, Operation<Double>.Add);
				Totalsquare = DAVectorUtility.mpiOps.allReduce(Totalsquare, MPI.SUM);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming1);
			}

			if (TotalNumberofPoints < 0.5)
			{
				return;
			}

			Totalmean = Totalmean / TotalNumberofPoints;
			Totalsquare = (Totalsquare / TotalNumberofPoints) - Totalmean * Totalmean;
			Totalsigma = Math.sqrt(Math.max(0.0, Totalsquare));
		}
	} // End FindMeanSigma

	public static class FindDoubleSum
	{
		private double[] NumberofPoints;
		private int NumberofThreads;
		private double[] TotalinThread;
		public double TotalNumberofPoints;
		public double Total;

		public FindDoubleSum(int NumThreads)
		{
			NumberofThreads = NumThreads;

			NumberofPoints = new double[NumThreads];
			TotalinThread = new double[NumThreads];

			TotalNumberofPoints = 0.0;
			Total = 0.0;
			for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
			{
				NumberofPoints[ThreadNo] = 0.0;
				TotalinThread[ThreadNo] = 0.0;
			}
		}

		public final void addapoint(int ThreadNo, double value1)
		{
			NumberofPoints[ThreadNo] += 1.0;
			TotalinThread[ThreadNo] += value1;
		}
		public final void sumoverthreadsandmpi() throws MPIException {
			for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
			{
				TotalNumberofPoints += NumberofPoints[ThreadNo];
				Total += TotalinThread[ThreadNo];
			}
			if (DAVectorUtility.MPI_Size > 1)
			{
				DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming1);
                // Note - MPI Call - Allreduce - double - sum
//				TotalNumberofPoints = DAVectorUtility.MPI_communicator.<Double>Allreduce(TotalNumberofPoints, Operation<Double>.Add);
				TotalNumberofPoints = DAVectorUtility.mpiOps.allReduce(TotalNumberofPoints, MPI.SUM);
                // Note - MPI Call - Allreduce - double - sum
//				Total = DAVectorUtility.MPI_communicator.<Double>Allreduce(Total, Operation<Double>.Add);
				Total = DAVectorUtility.mpiOps.allReduce(Total, MPI.SUM);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming1);
			}

		}
	} // End FindDoubleSum

	public static class FindDoubleMean
	{
		public double[] NumberofPoints;
		public int NumberofThreads;
		public double[] mean;
		public double TotalNumberofPoints;
		public double Totalmean;

		public FindDoubleMean(int NumThreads)
		{
			NumberofThreads = NumThreads;

			NumberofPoints = new double[NumThreads];
			mean = new double[NumThreads];

			TotalNumberofPoints = 0.0;
			Totalmean = 0.0;
			for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
			{
				NumberofPoints[ThreadNo] = 0.0;
				mean[ThreadNo] = 0.0;
			}
		}

		public final void addapoint(int ThreadNo, double value1)
		{
			NumberofPoints[ThreadNo] += 1.0;
			mean[ThreadNo] += value1;
		}
		public final void sumoverthreadsandmpi() throws MPIException {
			for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
			{
				TotalNumberofPoints += NumberofPoints[ThreadNo];
				Totalmean += mean[ThreadNo];
			}
			if (DAVectorUtility.MPI_Size > 1)
			{
				DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming1);
                // Note - MPI Call - Allreduce - double -sum
//				TotalNumberofPoints = DAVectorUtility.MPI_communicator.<Double>Allreduce(TotalNumberofPoints, Operation<Double>.Add);
				TotalNumberofPoints = DAVectorUtility.mpiOps.allReduce(TotalNumberofPoints, MPI.SUM);
                // Note - MPI Call - Allreduce - double -sum
//				Totalmean = DAVectorUtility.MPI_communicator.<Double>Allreduce(Totalmean, Operation<Double>.Add);
				Totalmean = DAVectorUtility.mpiOps.allReduce(Totalmean, MPI.SUM);
				DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming1);
			}

			if (TotalNumberofPoints < 0.5)
			{
				return;
			}

			Totalmean = Totalmean / TotalNumberofPoints;
		}
	} // End FindDoubleMean


	public static class FindVectorIntSum
	{
		private int[] NumberofPoints;
		private int NumberofThreads;
		private int[][] VectorSum;
		public int TotalNumberofPoints;
		public int[] TotalVectorSum;
		private int ArraySize;

		public FindVectorIntSum(int NumThreads, int NumberinArray)
		{
			NumberofThreads = NumThreads;
			ArraySize = NumberinArray;

			NumberofPoints = new int[NumThreads];
			TotalVectorSum = new int[ArraySize];
			VectorSum = new int[NumThreads][];

			TotalNumberofPoints = 0;
			for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
			{
				TotalVectorSum[ArrayLoop] = 0;
			}
		}

		public final void startthread(int ThreadNo)
		{
			NumberofPoints[ThreadNo] = 0;
			VectorSum[ThreadNo] = new int[ArraySize];
			for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
			{
				VectorSum[ThreadNo][ArrayLoop] = 0;
			}
		}

		public final void addapoint(int ThreadNo, int[] value1)
		{
			NumberofPoints[ThreadNo] += 1;
			for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
			{
				VectorSum[ThreadNo][ArrayLoop] += value1[ArrayLoop];
			}
		}

		public final void addapoint(int ThreadNo, int value1, int position)
		{
			NumberofPoints[ThreadNo] += 1;
			VectorSum[ThreadNo][position] += value1;
		}

		public final void sumoverthreadsandmpi() throws MPIException {
			for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
			{
				TotalNumberofPoints += NumberofPoints[ThreadNo];
				for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
				{
					TotalVectorSum[ArrayLoop] += VectorSum[ThreadNo][ArrayLoop];
				}
			}
			if (DAVectorUtility.MPI_Size > 1)
			{
				DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming1);
                // Note - MPI Call - Allreduce - int - sum
//				TotalNumberofPoints = DAVectorUtility.MPI_communicator.<Integer>Allreduce(TotalNumberofPoints, Operation<Integer>.Add);
				TotalNumberofPoints = DAVectorUtility.mpiOps.allReduce(TotalNumberofPoints, MPI.SUM);
                // Note - MPI Call - Allreduce - int[] - sum
//				TotalVectorSum = DAVectorUtility.MPI_communicator.<Integer>Allreduce(TotalVectorSum, Operation<Integer>.Add);
				DAVectorUtility.mpiOps.allReduce(TotalVectorSum, MPI.SUM);
				DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming1);
			}
		}
	} // End FindVectorIntSum

	public static class FindVectorDoubleMax
	{
		private double[] NumberofPoints;
		private int NumberofThreads;
		private double[][] VectorMax;
		public double TotalNumberofPoints;
		public double[] TotalVectorMax;
		private int ArraySize;

		public FindVectorDoubleMax(int NumThreads, int NumberinArray)
		{
			NumberofThreads = NumThreads;
			ArraySize = NumberinArray;

			NumberofPoints = new double[NumThreads];
			TotalVectorMax = new double[ArraySize];
			VectorMax = new double[NumThreads][];

			TotalNumberofPoints = 0.0;
			for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
			{
				NumberofPoints[ThreadNo] = 0.0;
				VectorMax[ThreadNo] = new double[ArraySize];
				for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
				{
					VectorMax[ThreadNo][ArrayLoop] = -1.0E10;
				}
			}
			for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
			{
				TotalVectorMax[ArrayLoop] = -1.0E10;
			}
		}

		public final void addapoint(int ThreadNo, double[] value)
		{
			NumberofPoints[ThreadNo] += 1.0;
			for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
			{
				VectorMax[ThreadNo][ArrayLoop] = Math.max(VectorMax[ThreadNo][ArrayLoop], value[ArrayLoop]);
			}
		}

		public final void addapoint(int ThreadNo, double value1, int position)
		{
			NumberofPoints[ThreadNo] += 1.0;
				VectorMax[ThreadNo][position] = Math.max(VectorMax[ThreadNo][position], value1);
		}

		public final void sumoverthreadsandmpi() throws MPIException {
			for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
			{
				TotalNumberofPoints += NumberofPoints[ThreadNo];
				for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
				{
					TotalVectorMax[ArrayLoop] = Math.max(TotalVectorMax[ArrayLoop], VectorMax[ThreadNo][ArrayLoop]);
				}
			}

			if (DAVectorUtility.MPI_Size > 1)
			{
				DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming1);
                // Note - MPI Call - Allreduce - double - sum
//				TotalNumberofPoints = DAVectorUtility.MPI_communicator.<Double>Allreduce(TotalNumberofPoints, Operation<Double>.Add);
				TotalNumberofPoints = DAVectorUtility.mpiOps.allReduce(TotalNumberofPoints, MPI.SUM);
                // Note - MPI Call - Allreduce - double[] - max
//				TotalVectorMax = DAVectorUtility.MPI_communicator.<Double>Allreduce(TotalVectorMax, Operation<Double>.Max);
				DAVectorUtility.mpiOps.allReduce(TotalVectorMax, MPI.MAX);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming1);
			}

		}

	} // End FindVectorDoubleMax

	public static class FindVectorDoubleSum
	{
		private double[] NumberofPoints;
		private int NumberofThreads;
		private double[][] VectorSum;
		public double TotalNumberofPoints;
		public double[] TotalVectorSum;
		private int ArraySize;

		public FindVectorDoubleSum(int NumThreads, int NumberinArray)
		{
			NumberofThreads = NumThreads;
			ArraySize = NumberinArray;

			NumberofPoints = new double[NumThreads];
			TotalVectorSum = new double[ArraySize];
			VectorSum = new double[NumThreads][];

			TotalNumberofPoints = 0.0;
			for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
			{
				NumberofPoints[ThreadNo] = 0.0;
				VectorSum[ThreadNo] = new double[ArraySize];
				for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
				{
					VectorSum[ThreadNo][ArrayLoop] = 0.0;
				}
			}
			for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
			{
				TotalVectorSum[ArrayLoop] = 0.0;
			}
		}

		public final void addapoint(int ThreadNo, double[] value1)
		{
			NumberofPoints[ThreadNo] += 1.0;
			for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
			{
				VectorSum[ThreadNo][ArrayLoop] += value1[ArrayLoop];
			}
		}

		public final void addapoint(int ThreadNo, double value1, int position)
		{
			NumberofPoints[ThreadNo] += 1.0;
			VectorSum[ThreadNo][position] += value1;
		}

		public final void sumoverthreadsandmpi() throws MPIException {
			for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
			{
				TotalNumberofPoints += NumberofPoints[ThreadNo];
				for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
				{
					TotalVectorSum[ArrayLoop] += VectorSum[ThreadNo][ArrayLoop];
				}
			}

			if (DAVectorUtility.MPI_Size > 1)
			{
				DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming1);
                // Note - MPI Call - Allreduce - double - sum
//				TotalNumberofPoints = DAVectorUtility.MPI_communicator.<Double>Allreduce(TotalNumberofPoints, Operation<Double>.Add);
				TotalNumberofPoints = DAVectorUtility.mpiOps.allReduce(TotalNumberofPoints, MPI.SUM);
                // Note - MPI Call - Allreduce - double[] - sum
//				TotalVectorSum = DAVectorUtility.MPI_communicator.<Double>Allreduce(TotalVectorSum, Operation<Double>.Add);
				DAVectorUtility.mpiOps.allReduce(TotalVectorSum, MPI.SUM);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming1);
			}

		}

	} // End FindVectorDoubleSum

	public static class FindVectorDoubleSum2
	{
		private double[] NumberofPoints;
		private int NumberofThreads;
		private double[][] VectorSum;
		public double TotalNumberofPoints;
		public double[] TotalVectorSum;
		private int ArraySize;
		private Range[] ParallelArrayRanges;

		public FindVectorDoubleSum2(int NumThreads, int NumberinArray)
		{
			NumberofThreads = NumThreads;
			ArraySize = NumberinArray;

			NumberofPoints = new double[NumThreads];
			TotalVectorSum = new double[ArraySize];
			VectorSum = new double[NumThreads][];

			TotalNumberofPoints = 0.0;
			for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
			{
				TotalVectorSum[ArrayLoop] = 0.0;
			}

			ParallelArrayRanges = RangePartitioner.Partition(NumberinArray, DAVectorUtility.ThreadCount);
		}

		public final void startthread(int ThreadNo)
		{
			NumberofPoints[ThreadNo] = 0.0;
			VectorSum[ThreadNo] = new double[ArraySize];
			for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
			{
				VectorSum[ThreadNo][ArrayLoop] = 0.0;
			}
		}

		public final void addapoint(int ThreadNo, double[] value1)
		{
			NumberofPoints[ThreadNo] += 1.0;
			for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
			{
				VectorSum[ThreadNo][ArrayLoop] += value1[ArrayLoop];
			}
		}

		public final void sumoverthreadsandmpi() throws MPIException {
            // Note - parallel for
            launchHabaneroApp(() -> {
                forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) ->
                {
                    int beginindex = ParallelArrayRanges[threadIndex].getStartIndex();
                    int indexlength = ParallelArrayRanges[threadIndex].getLength();
                    for (int ArrayLoop = beginindex; ArrayLoop < beginindex + indexlength; ArrayLoop++) {
                        TotalVectorSum[ArrayLoop] = 0.0;
                        for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                            TotalVectorSum[ArrayLoop] += VectorSum[ThreadNo][ArrayLoop];
                        }
                    }

                });
            });

            if (DAVectorUtility.MPI_Size > 1)
			{
				DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming1);
                // Note - MPI Call - Allreduce - double - sum
//				TotalNumberofPoints = DAVectorUtility.MPI_communicator.<Double>Allreduce(TotalNumberofPoints, Operation<Double>.Add);
				TotalNumberofPoints = DAVectorUtility.mpiOps.allReduce(TotalNumberofPoints, MPI.SUM);
                // Note - MPI Call - Allreduce - double[] - sum
//				TotalVectorSum = DAVectorUtility.MPI_communicator.<Double>Allreduce(TotalVectorSum, Operation<Double>.Add);
				DAVectorUtility.mpiOps.allReduce(TotalVectorSum, MPI.SUM);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming1);
			}

		}

	} // End FindVectorDoubleSum2

	public static class FindVectorDoubleSum3
	{
		private double[] NumberofPoints;
		private int NumberofThreads;
		private double[][] VectorSum;
		public double TotalNumberofPoints;
		public double[] TotalVectorSum;
		private int ArraySize;
		private int ArraySize1;
		private int ArraySize2;
		private Range[] ParallelArrayRanges;

		public FindVectorDoubleSum3(int NumThreads, int NumberinArray1, int NumberinArray2)
		{
			NumberofThreads = NumThreads;
			ArraySize = NumberinArray1 * NumberinArray2;
			ArraySize1 = NumberinArray1;
			ArraySize2 = NumberinArray2;

			NumberofPoints = new double[NumThreads];
			TotalVectorSum = new double[ArraySize];

			VectorSum = new double[NumThreads][];

			TotalNumberofPoints = 0.0;

			ParallelArrayRanges = RangePartitioner.Partition(ArraySize, DAVectorUtility.ThreadCount);
		}

		public final void startthread(int ThreadNo)
		{
			NumberofPoints[ThreadNo] = 0.0;
			VectorSum[ThreadNo] = new double[ArraySize];
			for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
			{
				VectorSum[ThreadNo][ArrayLoop] = 0.0;
			}

		}

		public final void addapoint(int ThreadNo, double[] value1, int index2)
		{
			NumberofPoints[ThreadNo] += 1.0;
			int additive = index2 * ArraySize1;
			for (int ArrayLoop1 = 0; ArrayLoop1 < ArraySize1; ArrayLoop1++)
			{
				VectorSum[ThreadNo][ArrayLoop1 + additive] += value1[ArrayLoop1];
			}
		}

		public final void sumoverthreadsandmpi() throws MPIException {
			DAVectorUtility.StartSubTimer(DAVectorUtility.ThreadTiming);
            // Note - parallel for
            launchHabaneroApp(() -> {
                forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) -> {
                    int beginindex = ParallelArrayRanges[threadIndex].getStartIndex();
                    int indexlength = ParallelArrayRanges[threadIndex].getLength();
                    for (int ArrayLoop = beginindex; ArrayLoop < beginindex + indexlength; ArrayLoop++) {
                        double tmp = 0.0;
                        for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                            tmp += VectorSum[ThreadNo][ArrayLoop];
                        }
                        TotalVectorSum[ArrayLoop] = tmp;
                    }

                });
            });

            DAVectorUtility.StopSubTimer(DAVectorUtility.ThreadTiming);

			if (DAVectorUtility.MPI_Size > 1)
			{
				DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming1);
                // Note - MPI Call - Allreduce - double - sum
//				TotalNumberofPoints = DAVectorUtility.MPI_communicator.<Double>Allreduce(TotalNumberofPoints, Operation<Double>.Add);
				TotalNumberofPoints = DAVectorUtility.mpiOps.allReduce(TotalNumberofPoints, MPI.SUM);
                int bigsize = TotalVectorSum.length;
				if (bigsize <= 4096)
				{
                    // Note - MPI Call - Allreduce - double[] - sum
//					TotalVectorSum = DAVectorUtility.MPI_communicator.<Double>Allreduce(TotalVectorSum, Operation<Double>.Add);
					DAVectorUtility.mpiOps.allReduce(TotalVectorSum, MPI.SUM);
                }
				else
				{
					double[] buffer = new double[4096];
					int start = 0;
                    // TODO - This block can be improved without having a separate buffer by using offsets with Allreduce call
					while (start < bigsize)
					{
						int whatsleft = Math.min(bigsize - start, 4096);
                        System.arraycopy(TotalVectorSum, start, buffer, 0, whatsleft);
                        // Note - MPI Call - Allreduce - double[] - sum
//						buffer = DAVectorUtility.MPI_communicator.<Double>Allreduce(buffer, Operation<Double>.Add);
						DAVectorUtility.mpiOps.allReduce(buffer, MPI.SUM);
                        System.arraycopy(buffer, 0, TotalVectorSum, start, whatsleft);
						start += whatsleft;
					}

				}
				DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming1);
			}
		}

	} // End FindVectorDoubleSum3

	public static class FindIndirectVectorDoubleSum
	{
		private double[] NumberofPoints;
		private int NumberofThreads;
		private double[][] VectorSum;
		public double TotalNumberofPoints;
		public double[] TotalVectorSum;
		private int ArraySize;
		private Range[] ParallelArrayRanges;

		public FindIndirectVectorDoubleSum(int NumThreads, int NumberinArray)
		{
			NumberofThreads = NumThreads;
			ArraySize = NumberinArray;

			NumberofPoints = new double[NumThreads];
			TotalVectorSum = new double[ArraySize];
			VectorSum = new double[NumThreads][];

			TotalNumberofPoints = 0.0;

			ParallelArrayRanges = RangePartitioner.Partition(NumberinArray, DAVectorUtility.ThreadCount);
		}

		public final void startthread(int ThreadNo)
		{
			NumberofPoints[ThreadNo] = 0.0;
			VectorSum[ThreadNo] = new double[ArraySize];
			for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
			{
				VectorSum[ThreadNo][ArrayLoop] = 0.0;
			}
		}

		public final void addapoint(int ThreadNo, int NumberLocations, int[] location, double[] value1)
		{
			if (NumberLocations <= 0)
			{
				return;
			}
			NumberofPoints[ThreadNo] += 1.0;
			for (int ArrayLoop = 0; ArrayLoop < NumberLocations; ArrayLoop++)
			{
				VectorSum[ThreadNo][location[ArrayLoop]] += value1[ArrayLoop];
			}
		}
		public final void sumoverthreadsandmpi() throws MPIException {

			for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
			{
				TotalNumberofPoints += NumberofPoints[ThreadNo];
			}

            // Note - parallel for
            launchHabaneroApp(() -> {
                forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) ->
                {
                    int beginindex = ParallelArrayRanges[threadIndex].getStartIndex();
                    int indexlength = ParallelArrayRanges[threadIndex].getLength();
                    for (int ArrayLoop = beginindex; ArrayLoop < beginindex + indexlength; ArrayLoop++) {
                        TotalVectorSum[ArrayLoop] = 0.0;
                        for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                            TotalVectorSum[ArrayLoop] += VectorSum[ThreadNo][ArrayLoop];
                        }
                    }

                });
            });

            if (DAVectorUtility.MPI_Size > 1)
			{
				DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming1);
                // Note - MPI Call - Allreduce - double - sum
//				TotalNumberofPoints = DAVectorUtility.MPI_communicator.<Double>Allreduce(TotalNumberofPoints, Operation<Double>.Add);
				TotalNumberofPoints = DAVectorUtility.mpiOps.allReduce(TotalNumberofPoints, MPI.SUM);
                // Note - MPI Call - Allreduce - double[] - sum
//				TotalVectorSum = DAVectorUtility.MPI_communicator.<Double>Allreduce(TotalVectorSum, Operation<Double>.Add);
				DAVectorUtility.mpiOps.allReduce(TotalVectorSum, MPI.SUM);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming1);
			}

		}

	} // End FindIndirectVectorDoubleSum

	public static class FindArrayMean
	{
		public double[] NumberofPoints;
		public int NumberofThreads;
		public double[][] mean;
		public double TotalNumberofPoints;
		public double[] Totalmean;
		public int ArraySize;

		public FindArrayMean(int NumThreads, int NumberinArray)
		{
			NumberofThreads = NumThreads;
			ArraySize = NumberinArray;

			NumberofPoints = new double[NumThreads];
			Totalmean = new double[ArraySize];
			mean = new double[NumThreads][];

			TotalNumberofPoints = 0.0;
			for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
			{
				NumberofPoints[ThreadNo] = 0.0;
				mean[ThreadNo] = new double[ArraySize];
				for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
				{
					mean[ThreadNo][ArrayLoop] = 0.0;
				}
			}
			for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
			{
				Totalmean[ArrayLoop] = 0.0;
			}
		}

		public final void addapoint(int ThreadNo, double[] value1)
		{
			NumberofPoints[ThreadNo] += 1.0;
			for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
			{
				mean[ThreadNo][ArrayLoop] += value1[ArrayLoop];
			}
		}
		public final void sumoverthreadsandmpi() throws MPIException {
			for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
			{
				TotalNumberofPoints += NumberofPoints[ThreadNo];
				for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
				{
					Totalmean[ArrayLoop] += mean[ThreadNo][ArrayLoop];
				}
			}
			if (DAVectorUtility.MPI_Size > 1)
			{
				DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming1);
                // Note - MPI Call - Allreduce - double - sum
//				TotalNumberofPoints = DAVectorUtility.MPI_communicator.<Double>Allreduce(TotalNumberofPoints, Operation<Double>.Add);
                DAVectorUtility.mpiOps.allReduce(TotalNumberofPoints,MPI.SUM);
                // Note - MPI Call - Allreduce - double[] - sum
//				Totalmean = DAVectorUtility.MPI_communicator.<Double>Allreduce(Totalmean, Operation<Double>.Add);
                DAVectorUtility.mpiOps.allReduce(Totalmean,MPI.SUM);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming1);
			}

			if (TotalNumberofPoints < 0.5)
			{
				return;
			}
			for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
			{
				Totalmean[ArrayLoop] = Totalmean[ArrayLoop] / TotalNumberofPoints;
			}
		}
	} // End FindArrayMean

	public static class FindMinorMaxValuewithIndex
	{
		public double[] NumberofPoints;
		public int NumberofThreads;
		public double[] MaxOrMinvalue;
		public double TotalNumberofPoints;
		public double TotalMaxOrMin;
		public int TotalIndexValue;
		public int[] IndexValue;
		public int MinMaxPointer; // =0 Min = 1 Max

		public FindMinorMaxValuewithIndex(int NumThreads, int UseMax)
		{
			NumberofThreads = NumThreads;
			MinMaxPointer = UseMax;

			NumberofPoints = new double[NumThreads];
			MaxOrMinvalue = new double[NumThreads];
			IndexValue = new int[NumThreads];

			TotalNumberofPoints = 0.0;
			TotalMaxOrMin = 0.0;
			TotalIndexValue = -1;
			for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
			{
				NumberofPoints[ThreadNo] = 0.0;
				MaxOrMinvalue[ThreadNo] = 0.0;
				IndexValue[ThreadNo] = -1;
			}
		}

		public final void addapoint(int ThreadNo, int indexposition, double value)
		{
			NumberofPoints[ThreadNo] += 1.0;
			if (MinMaxPointer != 0)
			{ // Max
				if ((IndexValue[ThreadNo] >= 0) && (MaxOrMinvalue[ThreadNo] > value))
				{
					return;
				}
			}
			else
			{ // Min
				if ((IndexValue[ThreadNo] >= 0) && (MaxOrMinvalue[ThreadNo] <= value))
				{
					return;
				}
			}
			MaxOrMinvalue[ThreadNo] = value;
			IndexValue[ThreadNo] = indexposition;
		}

		public final void sumoverthreadsandmpi() throws MPIException {
			for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
			{
				if (IndexValue[ThreadNo] < 0)
				{
					continue;
				}

				TotalNumberofPoints += NumberofPoints[ThreadNo];
				if (MinMaxPointer != 0)
				{
					if ((TotalIndexValue >= 0) && (TotalMaxOrMin > MaxOrMinvalue[ThreadNo]))
					{
						continue;
					}
				}
				else
				{
					if ((TotalIndexValue >= 0) && (TotalMaxOrMin <= MaxOrMinvalue[ThreadNo]))
					{
						continue;
					}
				}

				TotalMaxOrMin = MaxOrMinvalue[ThreadNo];
				TotalIndexValue = IndexValue[ThreadNo];
			}
			if (DAVectorUtility.MPI_Size > 1)
			{
				DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming1);
				if (MinMaxPointer != 0)
				{
                    // Note - MPI Call - Allreduce - MPIReducePlusIndex - MAX_WITH_INDEX
                    MPIReducePlusIndex result = DAVectorUtility.mpiOps.allReduce(new MPIReducePlusIndex(TotalIndexValue,TotalMaxOrMin), MPIReducePlusIndex.Op.MAX_WITH_INDEX);
                    TotalIndexValue = result.getIndex();
                    TotalMaxOrMin = result.getValue();
				}
				else
				{
                    // Note - MPI Call - Allreduce - MPIReducePlusIndex - MIN_WITH_INDEX
                    MPIReducePlusIndex result = DAVectorUtility.mpiOps.allReduce(new MPIReducePlusIndex(TotalIndexValue,TotalMaxOrMin), MPIReducePlusIndex.Op.MIN_WITH_INDEX);
                    TotalIndexValue = result.getIndex();
                    TotalMaxOrMin = result.getValue();
				}
                // Note - MPI Call - Allreduce - double - sum
//				TotalNumberofPoints = DAVectorUtility.MPI_communicator.<Double>Allreduce(TotalNumberofPoints, Operation<Double>.Add);
                TotalNumberofPoints = DAVectorUtility.mpiOps.allReduce(TotalNumberofPoints,MPI.SUM);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming1);
			}
		}
	} // End FindMinorMaxValuewithIndex

	public static class FindManyMinValuewithIndex
	{ // Finds top LimitNumberStored points by minimizing value given in addapoint
		// Rsults returned in order of figure of merit
		// Store values and indices
		// Uses FindMinimumSet
		// Results TotalNumberofPoints OrderedMinValue OrderedIndexValue

		public double[] NumberofPoints;
		public int NumberofThreads;
		public int Numbertofind;
		public double[][] MinValuebythread;
		public int[][] IndexValuebythread;
		public int[] CurrentWorstbythread;

		public double TotalNumberofPoints;
		public double[] TotalMinValue;
		public int[] TotalIndexValue;
		public int TotalWorst;
		public double[] OrderedMinValue;
		public int[] OrderedIndexValue;

		public FindManyMinValuewithIndex(int NumThreads, int LimitNumberStored)
		{
			NumberofThreads = NumThreads;
			Numbertofind = LimitNumberStored;

			NumberofPoints = new double[NumThreads];
			MinValuebythread = new double[NumThreads][];
			IndexValuebythread = new int[NumThreads][];
			CurrentWorstbythread = new int[NumThreads];

			TotalNumberofPoints = 0.0;
			TotalMinValue = new double[LimitNumberStored];
			TotalIndexValue = new int[LimitNumberStored];
			OrderedMinValue = new double[LimitNumberStored];
			OrderedIndexValue = new int[LimitNumberStored];

			for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
			{
				NumberofPoints[ThreadNo] = 0.0;
				CurrentWorstbythread[ThreadNo] = -1;
				MinValuebythread[ThreadNo] = new double[LimitNumberStored];
				IndexValuebythread[ThreadNo] = new int[LimitNumberStored];
				for (int storeloop = 0; storeloop < LimitNumberStored; storeloop++)
				{
					MinValuebythread[ThreadNo][storeloop] = -1.0;
					IndexValuebythread[ThreadNo][storeloop] = -1;
				}
			}
		}

		public final void addapoint(int ThreadNo, int indexposition, double value)
		{
			NumberofPoints[ThreadNo] += 1.0;
			Box<Integer> tempRef_Object = new Box<>(CurrentWorstbythread[ThreadNo]);
			FindMinimumSet(value, indexposition, tempRef_Object, MinValuebythread[ThreadNo], IndexValuebythread[ThreadNo], Numbertofind);
			CurrentWorstbythread[ThreadNo] = tempRef_Object.content;
		}

		public final void sumoverthreadsandmpi() throws MPIException {

			for (int storeloop = 0; storeloop < Numbertofind; storeloop++)
			{
				TotalMinValue[storeloop] = -1.0;
				TotalIndexValue[storeloop] = -1;
			}
			TotalWorst = -1;
			for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
			{
				for (int storeloop = 0; storeloop < Numbertofind; storeloop++)
				{
					if (IndexValuebythread[ThreadNo][storeloop] < 0)
					{
						continue; // End this thread
					}
					Box<Integer> tempRef_TotalWorst = new Box<>(TotalWorst);
					FindMinimumSet(MinValuebythread[ThreadNo][storeloop], IndexValuebythread[ThreadNo][storeloop], tempRef_TotalWorst, TotalMinValue, TotalIndexValue, Numbertofind);
					TotalWorst = tempRef_TotalWorst.content;
				}
			}
			if (DAVectorUtility.MPI_Size > 1)
			{
				DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming1);
                // Note - MPI Call - Allreduce - double - sum
//				TotalNumberofPoints = DAVectorUtility.MPI_communicator.<Double>Allreduce(TotalNumberofPoints, Operation<Double>.Add);
                TotalNumberofPoints = DAVectorUtility.mpiOps.allReduce(TotalNumberofPoints,MPI.SUM);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming1);
			}
			// Sort in absolute order and accumulate over processes. This takes Numbertofindsteps
			for (int OrderLoop = 0; OrderLoop < Numbertofind; OrderLoop++)
			{
				int localindex = -1; // unset
				double localvalue = -1.0;
				int loopused = -1;
				for (int internalloop = 0; internalloop < Numbertofind; internalloop++)
				{ // Find minimum
					if (TotalIndexValue[internalloop] < 0)
					{
						continue;
					}
					if ((localindex < 0) || (TotalMinValue[internalloop] < localvalue))
					{
						localindex = TotalIndexValue[internalloop];
						localvalue = TotalMinValue[internalloop];
						loopused = internalloop;
					}
				}
				int oldlocalindex = localindex;
				if (DAVectorUtility.MPI_Size > 1)
				{
					DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming1);
                    // Note - MPI Call - Allreduce - MPIReducePlusIndex - MIN_WITH_INDEX
                    MPIReducePlusIndex result = DAVectorUtility.mpiOps.allReduce(new MPIReducePlusIndex(localindex, localvalue), MPIReducePlusIndex.Op.MIN_WITH_INDEX);
					localvalue = result.getValue();
                    localindex = result.getIndex();
                    DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming1);
				}

				OrderedMinValue[OrderLoop] = localvalue;
				OrderedIndexValue[OrderLoop] = localindex;
				if ((oldlocalindex >= 0) && (OrderedIndexValue[OrderLoop] == oldlocalindex))
				{
					TotalIndexValue[loopused] = -1;
					TotalMinValue[loopused] = -1.0;
				}
			} // Loop over Order Loop

        }
	} // End FindManyMinValuewithIndex

	//  Support finding list of minimum values by inserting new point with value newvalue and index newindex into lists SmallValues SmallIndices
	//  In SmallIndices negative values correspond to unset values
	//  NumberSmallOnes is total number wanted
	//  Currentcut is position in SmallValues, SmallIndices of largest min value
	public static void FindMinimumSet(double newvalue, int newindex, Box<Integer> currentcut, double[] SmallValues, int[] SmallIndices, int NumberSmallones)
	{
		if (currentcut.content < 0)
		{
			currentcut.content = 0;
			SmallValues[0] = newvalue;
			SmallIndices[0] = newindex;
			return;
		}
		if (SmallIndices[NumberSmallones - 1] < 0)
		{ // Not all positions are filled so add at next available
			// Reset currentcut if worst
			for (int ivalue = 0; ivalue < NumberSmallones; ivalue++)
			{
				if (SmallIndices[ivalue] < 0)
				{
					SmallValues[ivalue] = newvalue;
					SmallIndices[ivalue] = newindex;
					if (SmallValues[ivalue] > SmallValues[currentcut.content])
					{
						currentcut.content = ivalue;
					}
					return;
				}
			}
		}
		if (newvalue >= SmallValues[currentcut.content])
		{
			return;
		}

		// Replace currentcut position with new values and Reset new worst position
		SmallValues[currentcut.content] = newvalue;
		SmallIndices[currentcut.content] = newindex;
		double maxvalue = -1.0;
		for (int ivalue = 0; ivalue < NumberSmallones; ivalue++)
		{
			if (SmallIndices[ivalue] < 0)
			{
				continue;
			}
			if (SmallValues[ivalue] > maxvalue)
			{
				currentcut.content = ivalue;
				maxvalue = SmallValues[ivalue];
			}
		}

    } // End FindMinimumSet


} // End class GlobalReductions