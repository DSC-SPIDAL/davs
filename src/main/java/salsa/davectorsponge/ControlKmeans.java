package salsa.davectorsponge;

import edu.rice.hj.api.SuspendableException;
import mpi.MPIException;
import salsa.general.Box;
import salsa.mpi.MPIPacket;

import static edu.rice.hj.Module1.forallChunked;

// Control Kmeans
public class ControlKmeans
{
	public static boolean initialized = false;

	public ControlKmeans() throws MPIException {

		//allocate memory on first and indeed only call
		if (!initialized)
		{
			if (!Program.DoKmeans)
			{
				DAVectorUtility.printAndThrowRuntimeException(" Invalid Kmeans Request");
			}

			initialized = true;

			Program.UseSponge = false;
			Program.ContinuousClustering = false;
			Program.MaxNumberSplitClusters = Math.max(Program.MaxNumberSplitClusters, 1);

			ParallelClustering.runningSolution = new ClusteringSolution(Program.UseSponge);
			ParallelClustering.savedSolution = ParallelClustering.runningSolution;
			ParallelClustering.bestSolution = ParallelClustering.runningSolution;
			DAVectorUtility.SALSAPrint(0, "Clustering Solutions Created");

		} //end Initialization

		InitializeSolution(ParallelClustering.runningSolution);

		ClusteringSolution.ClustersDeleted = 0;
		ClusteringSolution.ClustersMoved = 0;
		ClusteringSolution.ClustersSplit = 0;

		//  Run Kmeans
		Box<Integer> tempRef_Ncent_Global = new Box<>(ParallelClustering.runningSolution.Ncent_Global);
		Box<Double> tempRef_TotaloverVectorIndicesAverageWidth = new Box<>(ParallelClustering.runningSolution.TotaloverVectorIndicesAverageWidth);
		Kmeans.RunKmeans(ParallelClustering.runningSolution.Y_k_i_, ParallelClustering.runningSolution.OccupationCounts_k_, ParallelClustering.runningSolution.ClusterScaledSquaredWidth_k_, tempRef_Ncent_Global, tempRef_TotaloverVectorIndicesAverageWidth);
		ParallelClustering.runningSolution.Ncent_Global = tempRef_Ncent_Global.content;
		ParallelClustering.runningSolution.TotaloverVectorIndicesAverageWidth = tempRef_TotaloverVectorIndicesAverageWidth.content;

    } // End ControlKmeans

	public static void SetSolution(ClusteringSolution StartSolution, int Ncent_GlobalINPUT) throws MPIException {

		for (int RealClusterIndex = 0; RealClusterIndex < Ncent_GlobalINPUT; RealClusterIndex++)
		{
			StartSolution.P_k_[0] = 1.0;
			StartSolution.SplitPriority_k_[RealClusterIndex] = -1;
			StartSolution.Splittable_k_[RealClusterIndex] = 0;
			StartSolution.LocalSplitCreatedIndex[RealClusterIndex] = 0;
			StartSolution.LocalStatus[RealClusterIndex] = 1;
			int CreatedIndex = ClusteringSolution.SetCreatedIndex(RealClusterIndex);
		}

		StartSolution.Ncent_Global = Ncent_GlobalINPUT;
		StartSolution.Ncent_ThisNode = Ncent_GlobalINPUT;
		StartSolution.SetActiveClusters();
	}

	public static void InitializeSolution(ClusteringSolution StartSolution) throws MPIException {
		// Temperature always 1 for K means
		StartSolution.SpongeCluster = -1;

		StartSolution.Temperature = 1.0;
		Program.ActualStartTemperature = StartSolution.Temperature; // For Output
		Program.TargetEndTemperature = StartSolution.Temperature; // For Output

		DAVectorUtility.SALSAPrint(1, "Points " + DAVectorUtility.PointCount_Global + " Kmeans");
		StartSolution.ActualCoolingFactor = Program.InitialCoolingFactor;

		StartSolution.PairwiseHammy = 0.0;
		ControlKmeans.SetSolution(StartSolution, Program.InitialNcent);

    } // End InitializeSolution

	public static double[][] CaptureKmeans(int[] FullAssignment) throws MPIException { // Save Kmeans assignment and centers

        // Note - parallel for
        try {
            forallChunked(0, DAVectorUtility.ThreadCount - 1,
                    (threadIndex) -> // End Sum over Threads -  End loop over Points
                    {
                        //                    FullAssignment[alpha + DAVectorUtility.PointStart_Process] = Kmeans.InitialPointAssignment[alpha];
                        int indexlen = DAVectorUtility.PointsperThread[threadIndex];
                        int beginpoint = DAVectorUtility.StartPointperThread[threadIndex] - DAVectorUtility.PointStart_Process;
                        System.arraycopy(KmeansTriangleInequality.NearestCentertoPoint, beginpoint, FullAssignment,
                                beginpoint + DAVectorUtility.PointStart_Process, indexlen + beginpoint - beginpoint);

                    });
        } catch (SuspendableException e) {
            DAVectorUtility.printAndThrowRuntimeException(e.getMessage());
        }

        if (DAVectorUtility.MPI_Size > 1)
		{
			int MPItag = 100;
			for (int mpiloop = 1; mpiloop < DAVectorUtility.MPI_Size; mpiloop++)
			{
				if (DAVectorUtility.MPI_Rank == 0)
				{
                    // Note - MPI Call - Receive - MPIPacket<Integer>
//					fromsource = DAVectorUtility.MPI_communicator.<MPIPacket<Integer>>Receive(mpiloop, MPItag);
					MPIPacket fromsource = DAVectorUtility.mpiOps.receive(mpiloop, MPItag, MPIPacket.Type.Integer);
                    int numberOfPoints = fromsource.getNumberOfPoints();
                    for (int index = 0; index < numberOfPoints; index++)
					{
						FullAssignment[index + fromsource.getFirstPoint()] = fromsource.getMArrayIntAt(index);
					}
				}
				else
				{
					if (DAVectorUtility.MPI_Rank == mpiloop)
					{
						MPIPacket tosend = MPIPacket.newIntegerPacket(DAVectorUtility.PointCount_Process);
						tosend.setFirstPoint(DAVectorUtility.PointStart_Process);
						tosend.setNumberOfPoints(DAVectorUtility.PointCount_Process);
						for (int index = 0; index < DAVectorUtility.PointCount_Process; index++)
						{
							tosend.setMArrayIntAt(index,FullAssignment[index + DAVectorUtility.PointStart_Process]);
						}
                        // Note - MPI Call - Send - MPIPacket<Integer>
//						DAVectorUtility.MPI_communicator.<MPIPacket<Integer>>Send(tosend, 0, MPItag);
						DAVectorUtility.mpiOps.send(tosend, 0, MPItag);
                    }
				}
				DAVectorUtility.MPI_communicator.barrier();
			}
		}

		return Kmeans.ClusterCenter;

	} // End CaptureKmeans


} // End ControlKmeans