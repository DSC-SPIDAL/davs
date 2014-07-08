package salsa.davectorsponge;

import edu.rice.hj.api.SuspendableException;
import mpi.MPI;
import mpi.MPIException;
import salsa.general.Box;

import static edu.rice.hj.Module1.forallChunked;

//  These algorithms use distributed cluster formalism with both clusters and points in parallel
//  Instead of AllReduce we use local shifts up and down the decomposition
//  FindIndirectMultiVectorDoubleSum replaces routine by this name in GlobalReductions.cs
//  Current version only does multiple components and is more efficient than simple one in GlobalReductions.cs
public class DistributedReductions
{

	public static class FindIndirectMultiVectorDoubleSum
	{
		private int NumberofThreads;
		private double[][][] VectorSum;
		public double[][] TotalVectorSum;
		private int ArraySize;
		private int[] ThreadArraySize;
		private int NumberDoubleComponents;

		public FindIndirectMultiVectorDoubleSum()
		{
			NumberofThreads = DAVectorUtility.ThreadCount;
			ArraySize = DistributedClusteringSolution.NodeAccMetaData.NumberofPointsperNode;
			ThreadArraySize = new int[NumberofThreads];
			NumberDoubleComponents = 0;
		}

		public final int AddComponents(int Number)
		{
			int old = NumberDoubleComponents;
			NumberDoubleComponents += Number;
			return old;
		}

		public final void NodeInitialize()
		{
			TotalVectorSum = new double[ArraySize][];
			VectorSum = new double[NumberofThreads][][];
			for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
			{
				TotalVectorSum[ArrayLoop] = new double[NumberDoubleComponents];
				for (int ComponentIndex = 0; ComponentIndex < NumberDoubleComponents; ComponentIndex++)
				{
					TotalVectorSum[ArrayLoop][ComponentIndex] = 0.0;
				}
			}
		}

		public final void ThreadInitialize(int ThreadNo)
		{
			ThreadArraySize[ThreadNo] = DistributedClusteringSolution.NodeAccMetaData.NumberofPointsperThread[ThreadNo];
			VectorSum[ThreadNo] = new double[ThreadArraySize[ThreadNo]][];
			for (int ArrayLoop = 0; ArrayLoop < ThreadArraySize[ThreadNo]; ArrayLoop++)
			{
				VectorSum[ThreadNo][ArrayLoop] = new double[NumberDoubleComponents];
				for (int ComponentIndex = 0; ComponentIndex < NumberDoubleComponents; ComponentIndex++)
				{
					VectorSum[ThreadNo][ArrayLoop][ComponentIndex] = 0.0;
				}
			}
		}

		// Location points to position in local Thread Array and value is double to be accumulated for each component
		public final void addapoint(int ThreadNo, int location, double[] value)
		{
			for (int ComponentIndex = 0; ComponentIndex < NumberDoubleComponents; ComponentIndex++)
			{
				VectorSum[ThreadNo][location][ComponentIndex] += value[ComponentIndex];
			}
		}

		public final void addapoint(int ThreadNo, int location, int StartComponent, int LengthComponent, double[] value)
		{
			for (int RestrictedComponentIndex = 0; RestrictedComponentIndex < LengthComponent; RestrictedComponentIndex++)
			{
				VectorSum[ThreadNo][location][RestrictedComponentIndex + StartComponent] += value[RestrictedComponentIndex];
			}
		}

		public final void addapoint(int ThreadNo, int location, int ComponentIndex, double value)
		{
			VectorSum[ThreadNo][location][ComponentIndex] += value;
		}

		public final void addapoint(int ThreadNo, int location, double value)
		{
			VectorSum[ThreadNo][location][0] += value;
		}

		public final void sumoverthreadsandmpi() throws MPIException {
            // Note - parallel for
            try {
                forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) ->
                {
                    int beginindex = DistributedClusteringSolution.ParallelNodeAccumulationRanges[threadIndex].getStartIndex();
                    int indexlength = DistributedClusteringSolution.ParallelNodeAccumulationRanges[threadIndex].getLength();
                    for (int NodeAccumulationIndex = beginindex; NodeAccumulationIndex < beginindex + indexlength; NodeAccumulationIndex++) {
                        for (int ThreadNo = 0; ThreadNo < DAVectorUtility.ThreadCount; ThreadNo++) {
                            int IndexforThread = DistributedClusteringSolution.NodeAccMetaData.AccumulationNodetoThreadClusterAssociations[NodeAccumulationIndex][ThreadNo];
                            if (IndexforThread >= 0) {
                                for (int ComponentIndex = 0; ComponentIndex < NumberDoubleComponents; ComponentIndex++) {
                                    TotalVectorSum[NodeAccumulationIndex][ComponentIndex] += VectorSum[ThreadNo][IndexforThread][ComponentIndex];
                                }
                            }
                        }
                    }

                }); // End Sum over Threads
            } catch (SuspendableException e) {
                DAVectorUtility.printAndThrowRuntimeException(e.getMessage());
            }

            //  Sum over nodes using pipelined transport
			if (DAVectorUtility.MPI_Size > 1)
			{
				// Divide Clusters into 3 types
				// NodeAccumulationClusterStatus = 2 Locally controlled Distributed Cluster Will be Updated
				// NodeAccumulationClusterStatus = 0, 1 Global Cluster
				// NodeAccumulationClusterStatus = 3 Remotely controlled Distributed Cluster 
				// Types 1 and 3 may send data (both) Up and Down

				// Process Global Clusters -- there are same total in all nodes and must be stored in same order
				//  They do not need to be in same absolute position
				int LocalTotal = 0; // Number of Global Clusters
				int NumberNodeAccumulationPoints = DistributedClusteringSolution.NodeAccMetaData.NumberofPointsperNode;
				for (int NodeAccumulationIndex = 0; NodeAccumulationIndex < NumberNodeAccumulationPoints; NodeAccumulationIndex++)
				{
					if (DistributedClusteringSolution.NodeAccMetaData.NodeAccumulationClusterStatus[NodeAccumulationIndex] > 1)
					{
						continue;
					}
					++LocalTotal;
				}
				if (LocalTotal > 0)
				{
					int LocalTotaltimesComponents = LocalTotal * NumberDoubleComponents;
					double[] GlobalClusterComponent = new double[LocalTotaltimesComponents];
					int[] GlobalClusterIndex = new int[LocalTotal];
					int NumberGlobal1 = 0;
					int NumberGlobal2 = 0;
					for (int NodeAccumulationIndex = 0; NodeAccumulationIndex < NumberNodeAccumulationPoints; NodeAccumulationIndex++)
					{
						int LocalStatus = DistributedClusteringSolution.NodeAccMetaData.NodeAccumulationClusterStatus[NodeAccumulationIndex];
						if ((LocalStatus < 0) || (LocalStatus > 1))
						{
							continue;
						}
						for (int ComponentIndex = 0; ComponentIndex < NumberDoubleComponents; ComponentIndex++)
						{
							GlobalClusterComponent[NumberGlobal2] = TotalVectorSum[NodeAccumulationIndex][ComponentIndex];
							++NumberGlobal2;
						}
						GlobalClusterIndex[NumberGlobal1] = NodeAccumulationIndex;
						++NumberGlobal1;
						if (NumberGlobal1 >= LocalTotal)
						{
							break;
						}
					}
					DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming6);
                    // Note - MPI Call - Allreduce - double [] - sum
//					GlobalClusterComponent = DAVectorUtility.MPI_communicator.<Double>Allreduce(GlobalClusterComponent, Operation<Double>.Add);
					DAVectorUtility.mpiOps.allReduce(GlobalClusterComponent, MPI.SUM);

                    DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming6);
					NumberGlobal2 = 0;
					for (int LocalIndex = 0; LocalIndex < LocalTotal; LocalIndex++)
					{
						for (int ComponentIndex = 0; ComponentIndex < NumberDoubleComponents; ComponentIndex++)
						{
							TotalVectorSum[GlobalClusterIndex[LocalIndex]][ComponentIndex] = GlobalClusterComponent[ComponentIndex];
							++NumberGlobal2;
						}
					}
				} // End Case where there are Global Clusters

				//  Distributed Clusters
				DAVectorUtility.StartSubTimer(DAVectorUtility.MPIDistributedREDUCETiming);
				int FinalClusterCount = 0; // Dummy
				DistributedSynchronization.TransportviaPipeline DoDistributedTransfer
                        = new DistributedSynchronization.TransportviaPipeline(
                        1, false, NumberDoubleComponents, 0, NumberNodeAccumulationPoints, 3);
                Box<Integer> tempRef_FinalClusterCount = new Box<>(FinalClusterCount);
				DoDistributedTransfer.PipelineDistributedBroadcast(TotalVectorSum, TotalVectorSum, null, null,
                        DistributedClusteringSolution.NodeAccMetaData.NodeAccumulationCreatedIndices,
                        DistributedClusteringSolution.NodeAccMetaData.NodeAccumulationCreatedIndices,
                        DistributedClusteringSolution.NodeAccMetaData.NodeAccumulationClusterHosts,
                        DistributedClusteringSolution.NodeAccMetaData.NodeAccumulationClusterHosts, tempRef_FinalClusterCount);
				FinalClusterCount = tempRef_FinalClusterCount.content;
				DAVectorUtility.StopSubTimer(DAVectorUtility.MPIDistributedREDUCETiming);

			} // End Case where MPI needed

		} // End sumoverthreadsandmpi()

	} // End FindIndirectMultiVectorDoubleSum

} // End class DistributedReductions