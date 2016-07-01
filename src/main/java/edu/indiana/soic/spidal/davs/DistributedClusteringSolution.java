package edu.indiana.soic.spidal.davs;

import edu.indiana.soic.spidal.general.Box;
import mpi.MPI;
import mpi.MPIException;

import static edu.rice.hj.Module0.launchHabaneroApp;
import static edu.rice.hj.Module1.forallChunked;


public class DistributedClusteringSolution
{

	public static int MaxMPITransportBuffer; // Maximum size of buffer for input and output transport
	public static int MaxTransportedClusterStorage; // Maximum size of Storage for input and output transport
	public static int MaxNumberAccumulationsperNode; // Maximum size of buffer for Node Accumulations
	public static int MaxDoubleComponents; // Maximum Number of Double Components Allowed
	public static int MaxIntegerComponents; // Maximum Number of Integer Components Allowed

	public static double[] NodeCutUpperLimit; // Upper limit of one dimemsional cut for each node
	public static double[] ThreadCutUpperLimit; // Upper limit of one dimensional cut for each thread in this node; thread cuts in other nodes are NOT known

	public static TransportedClusterStoredInformation StorageforTransportedClusters; // Structure Storing Distributed Clusters hosted elsewhere
	public static TemporaryLocalClusterInfoForTransport TemporaryLocalClustersforSynchronization; // Structure for storing Local Cluster Data for Synchronization

	public static Range[] ParallelNodeAccumulationRanges; // Ranges for parallel loops over Node Accumulations
	public static NodeAccumulationMetaData NodeAccMetaData; // Hold all the metadata for Node Accumulatio

    // TODO - MPI Status
//	public static MPI.CompletedStatus MPISecStatus; // MPI Send Receive Status

	public static String HostLinkageMessage; // Holds Information on Linkage

	public DistributedClusteringSolution(int MaxBuf, int MaxStorage, int MaxAcc, int _MaxDoubleComponents, int _MaxIntegerComponents) throws MPIException { // Initialize Distributed Clusters

		MaxTransportedClusterStorage = MaxStorage;
		MaxNumberAccumulationsperNode = MaxAcc;
		MaxDoubleComponents = _MaxDoubleComponents;
		MaxIntegerComponents = _MaxIntegerComponents;
		MaxMPITransportBuffer = MaxBuf;

		//  This is place on node where information on clusters external to this node are stored
		StorageforTransportedClusters = new TransportedClusterStoredInformation(MaxStorage);

		//  This is temporary storage for local cluster data so it can be pipelined to other nodes
		TemporaryLocalClustersforSynchronization = new TemporaryLocalClusterInfoForTransport(MaxStorage);

		//  Accumulation Arrays to support distributed reductions
		NodeAccMetaData = new NodeAccumulationMetaData(MaxNumberAccumulationsperNode);

		//  Define Responsibity of Nodes and Threads
		SetOneDimensionalUpperLimits();

		//  Array for RECEIVING MPI Information
		DistributedSynchronization.TransportComponent = new MPITransportComponentPacket(MaxMPITransportBuffer, MaxDoubleComponents, MaxIntegerComponents);

	} // End DistributedClusteringSolution

	//  Note that in ambiguous cases (Cluster value = Upper limit, the current host node makes the decision and broadcasts that
	//  It won't matter which decision it makes 
	public static void SetOneDimensionalUpperLimits() throws MPIException {
		double[] NodeCutLowerLimit = new double[DAVectorUtility.MPI_Size];
		NodeCutUpperLimit = new double[DAVectorUtility.MPI_Size];
		ThreadCutUpperLimit = new double[DAVectorUtility.ThreadCount];

		double Lower = Program.PointPosition[0][0];
		double Upper = Program.PointPosition[DAVectorUtility.PointCount_Process - 1][0];

		DAVectorUtility.StartSubTimer(DAVectorUtility.MPIGATHERTiming);

        if (DAVectorUtility.MPI_Size > 1){
            // Note - MPI Call - Allgather double
//            NodeCutLowerLimit = DAVectorUtility.MPI_communicator.<Double>Allgather(Lower);
            DAVectorUtility.mpiOps.allGather(Lower,NodeCutLowerLimit);
            // Note - MPI Call - Allgather double
//            NodeCutUpperLimit = DAVectorUtility.MPI_communicator.<Double>Allgather(Upper);
            DAVectorUtility.mpiOps.allGather(Upper,NodeCutUpperLimit);
        } else {
            NodeCutLowerLimit[0] = Lower;
            NodeCutUpperLimit[0] = Upper;
        }

		DAVectorUtility.StopSubTimer(DAVectorUtility.MPIGATHERTiming);

		String message = "\nPosition Cuts";
		for (int rank = 0; rank < DAVectorUtility.MPI_Size - 1; rank++)
		{
			NodeCutUpperLimit[rank] = 0.5 * (NodeCutUpperLimit[rank] + NodeCutLowerLimit[rank + 1]);
			message += " " + String.format("%1$5.4E", NodeCutUpperLimit[rank]);
		}
		DAVectorUtility.SALSAPrint(1, message);

		for (int ThreadNo = 0; ThreadNo < DAVectorUtility.ThreadCount; ThreadNo++)
		{
			int indextop = DAVectorUtility.PointsperThread[ThreadNo] + DAVectorUtility.StartPointperThread[ThreadNo] - DAVectorUtility.PointStart_Process - 1;
			ThreadCutUpperLimit[ThreadNo] = Program.PointPosition[indextop][0];
			if (ThreadNo != (DAVectorUtility.ThreadCount - 1))
			{
				ThreadCutUpperLimit[ThreadNo] = 0.5 * (ThreadCutUpperLimit[ThreadNo] + Program.PointPosition[indextop + 1][0]);
			}
		}

	} // End SetOneDimensionalUpperLimits()

	//  Set Cluster Host associations needed before any Major Synchronization
	public static void ClusterHostlinkage(boolean useClustertype1) throws MPIException {
		int TestType = 2;
		if (useClustertype1)
		{
			TestType = 1;
		}

		String MoveMessage = "";
		String UpMessage = "";
		String DownMessage = "";
		for (int ActiveClusterIndex = 0; ActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; ActiveClusterIndex++)
		{
			int RealClusterIndex = ClusteringSolution.RealClusterIndices[ActiveClusterIndex];
			int H = DAVectorUtility.MPI_Rank;
			if (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] == 0)
			{
				ClusteringSolution.LocalHost[ActiveClusterIndex] = H | (H | (H << ClusteringSolution.PACKINGSHIFT)) << ClusteringSolution.PACKINGSHIFT;
				ClusteringSolution.UniversalMapping[ParallelClustering.runningSolution.LocalCreatedIndex[RealClusterIndex]].PackedHost = ClusteringSolution.LocalHost[ActiveClusterIndex];
				continue;
			}
			if (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] != TestType)
			{
				continue;
			}

			// New Children MUST be identical to Parent
			int Parent = ParallelClustering.runningSolution.LocalSplitCreatedIndex[RealClusterIndex];
			if (Parent > 0)
			{
				ClusteringSolution.LocalHost[ActiveClusterIndex] = ClusteringSolution.UniversalMapping[Parent - 1].PackedHost;
				ClusteringSolution.UniversalMapping[ParallelClustering.runningSolution.LocalCreatedIndex[RealClusterIndex]].PackedHost = ClusteringSolution.UniversalMapping[Parent - 1].PackedHost;
				continue;
			}
			double TestPosition = ParallelClustering.runningSolution.Y_k_i_[RealClusterIndex][0];

			H = 0;
			while (H < DAVectorUtility.MPI_Size - 1)
			{
				if (TestPosition < NodeCutUpperLimit[H])
				{
					break;
				}
				H++;
			}

			int H1 = H;
			while (H1 >= 1)
			{
				double tmp = TestPosition - NodeCutUpperLimit[H1 - 1];
				if (tmp > 0.0)
				{
					double scaled1Ddistance = tmp * tmp / ParallelClustering.runningSolution.Sigma_k_i_[RealClusterIndex][0];
					if (scaled1Ddistance >= (2.0 * Program.ExpArgumentCut3 * ParallelClustering.runningSolution.Temperature))
					{
						break;
					}
				}
				H1--;
			}

			int H2 = H;
			while (H2 < DAVectorUtility.MPI_Size - 1)
			{
				double tmp = TestPosition - NodeCutUpperLimit[H2];
				if (tmp < 0.0)
				{
					double scaled1Ddistance = tmp * tmp / ParallelClustering.runningSolution.Sigma_k_i_[RealClusterIndex][0];
					if (scaled1Ddistance >= (2.0 * Program.ExpArgumentCut3 * ParallelClustering.runningSolution.Temperature))
					{
						break;
					}
				}
				H2++;
			}
			int CreatedIndex = ParallelClustering.runningSolution.LocalCreatedIndex[RealClusterIndex];
			if (useClustertype1)
			{
				if (H == DAVectorUtility.MPI_Rank)
				{
					if (H1 != DAVectorUtility.MPI_Rank)
					{
						DownMessage += CreatedIndex + " ";
					}
					if (H2 != DAVectorUtility.MPI_Rank)
					{
						UpMessage += CreatedIndex + " ";
					}
				}
			}
			else
			{
				if (H != DAVectorUtility.MPI_Rank)
				{
					MoveMessage += CreatedIndex + " ";
				}
				if (H1 != DAVectorUtility.MPI_Rank)
				{
					DownMessage += CreatedIndex + " ";
				}
				if (H2 != DAVectorUtility.MPI_Rank)
				{
					UpMessage += CreatedIndex + " ";
				}
			}

			H1 = H1 | (H2 << ClusteringSolution.PACKINGSHIFT);
			ClusteringSolution.LocalHost[ActiveClusterIndex] = H | (H1 << ClusteringSolution.PACKINGSHIFT);
			ClusteringSolution.UniversalMapping[CreatedIndex].PackedHost = ClusteringSolution.LocalHost[ActiveClusterIndex];
		}

		HostLinkageMessage = "";
		if (MoveMessage.length() != 0)
		{
			HostLinkageMessage += "Move " + MoveMessage;
		}
		if (DownMessage.length() != 0)
		{
			HostLinkageMessage += "Down " + DownMessage;
		}
		if (UpMessage.length() != 0)
		{
			HostLinkageMessage += "Up " + UpMessage;
		}
		if (!useClustertype1)
		{
			return;
		}

		//  If Processing Status 1 clusters, then take results from Processor 0
		//  If Processing Status 2 clusters, then each Processor decides on utility and placement of its clusters
		DAVectorUtility.StartSubTimer(DAVectorUtility.MPIBROADCASTTiming);

        // Note - MPI Call - Broadcast - int[]
        if (DAVectorUtility.MPI_Size > 1){
//            DAVectorUtility.MPI_communicator.<Integer>Broadcast(ClusteringSolution.LocalHost, 0);
            DAVectorUtility.mpiOps.broadcast(ClusteringSolution.LocalHost, 0);
        }
		DAVectorUtility.StopSubTimer(DAVectorUtility.MPIBROADCASTTiming);

	} // End ClusterHostlinkage()

	//  Return Array pointers associated with any cluster based on CreatedIndex hashing
	public static int IndicesperCluster(int CreatedIndex, int ThreadNo, Box<Integer> LocalClusterIndex, Box<Integer> TransportedClusterIndex_instorage, Box<Integer> NodeAccumulationPosition, Box<Integer> ThreadAccumulationPosition)
	{
		if (ClusteringSolution.UniversalMapping[CreatedIndex].IterationSet < ClusteringSolution.CurrentIteration)
		{
			return -1;
		}
		int position = ClusteringSolution.UniversalMapping[CreatedIndex].Availability;
		if (position == 0)
		{
			return -2;
		}
		if (position > 0)
		{
			TransportedClusterIndex_instorage.content = -1;
			LocalClusterIndex.content = position - 1;
			int ActiveClusterIndex = ClusteringSolution.ActiveClusterIndices[LocalClusterIndex.content];
			NodeAccumulationPosition.content = ClusteringSolution.LocalNodeAccPosition[ActiveClusterIndex];
			if (NodeAccumulationPosition.content == -1)
			{
				return -3;
			}
			ThreadAccumulationPosition.content = -1;
			if (ThreadNo >= 0)
			{
				ThreadAccumulationPosition.content = ClusteringSolution.LocalThreadAccPosition[ActiveClusterIndex][ThreadNo];
			}
			return 0;
		}
		LocalClusterIndex.content = -1;
		TransportedClusterIndex_instorage.content = -position - 1;
		NodeAccumulationPosition.content = StorageforTransportedClusters.TransportedNodeAccPosition[TransportedClusterIndex_instorage.content];
		if (NodeAccumulationPosition.content == -1)
		{
			return -4;
		}
		ThreadAccumulationPosition.content = -1;
		if (ThreadNo >= 0)
		{
			ThreadAccumulationPosition.content = StorageforTransportedClusters.TransportedThreadAccPosition[TransportedClusterIndex_instorage.content][ThreadNo];
		}
		return 1;

	} // End IndicesperCluster

	//  Only called in distributed mode
	//  Perform Minor Synchronization when cluster positions Updated
	public static void MinorSynchronizationTransportDistributedClusterCenters() throws MPIException {
		DAVectorUtility.StartSubTimer(7);
		int NumberofLocalDistributedClusters = 0;
		++Program.NumberMinorSynchs;
		if (VectorAnnealIterate.EMIterationCount % Program.PrintInterval == 0)
		{
			if (ParallelClustering.runningSolution.DistributedExecutionMode)
			{
				DAVectorUtility.SALSAPrint(1, "Minor Synchronization " + Program.NumberMinorSynchs + " Iteration " + ClusteringSolution.CurrentIteration + " Temperature " + String.format("%1$5.4f", ParallelClustering.runningSolution.Temperature) + " Cluster Count " + ParallelClustering.runningSolution.Ncent_Global);
			}
		}
		int H = DAVectorUtility.MPI_Rank;

		for (int ActiveClusterIndex = 0; ActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; ActiveClusterIndex++)
		{
			int RealClusterIndex = ClusteringSolution.RealClusterIndices[ActiveClusterIndex];
			if (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] < 2)
			{
				continue;
			}
			int PackedHost = ClusteringSolution.LocalHost[ActiveClusterIndex];
//C# TO JAVA CONVERTER WARNING: The right shift operator was not replaced by Java's logical right shift operator since the left operand was not confirmed to be of an unsigned type, but you should review whether the logical right shift operator (>>>) is more appropriate:
			int H1 = (PackedHost >> ClusteringSolution.PACKINGSHIFT) & ClusteringSolution.PACKINGMASK;
//C# TO JAVA CONVERTER WARNING: The right shift operator was not replaced by Java's logical right shift operator since the left operand was not confirmed to be of an unsigned type, but you should review whether the logical right shift operator (>>>) is more appropriate:
			int H2 = PackedHost >> (2 * ClusteringSolution.PACKINGSHIFT);
			if ((H1 == H) && (H2 == H))
			{
				if ((PackedHost & ClusteringSolution.PACKINGMASK) != H)
				{
					DAVectorUtility.printAndThrowRuntimeException(" Inconsistent Host Range " + (PackedHost & ClusteringSolution.PACKINGMASK) + " in Node " + H);

				}
				continue;
			}
			TemporaryLocalClustersforSynchronization.totalTransportedOriginalHost[NumberofLocalDistributedClusters] = PackedHost;
			TemporaryLocalClustersforSynchronization.totalTransportedCreatedIndex[NumberofLocalDistributedClusters] = ParallelClustering.runningSolution.LocalCreatedIndex[RealClusterIndex];
            System.arraycopy(ParallelClustering.runningSolution.Y_k_i_[RealClusterIndex], 0,
                    TemporaryLocalClustersforSynchronization.totalTransportedY_t_i[NumberofLocalDistributedClusters], 0,
                    Program.ParameterVectorDimension);
			TemporaryLocalClustersforSynchronization.totalTransportedY_t_i[NumberofLocalDistributedClusters][Program.ParameterVectorDimension] = ParallelClustering.runningSolution.P_k_[RealClusterIndex];
			++NumberofLocalDistributedClusters;
		} // End Loop over LocalClusterIndex

		TemporaryLocalClustersforSynchronization.sizeOfTransportedArray = NumberofLocalDistributedClusters;
		Program.ActualMaxTransportedClusterStorage = Math.max(Program.ActualMaxTransportedClusterStorage, NumberofLocalDistributedClusters);
		DistributedSynchronization.TransportviaPipeline DoSynchronization = new DistributedSynchronization.TransportviaPipeline(0, true, 3, 0, NumberofLocalDistributedClusters, 2);

		int FinalClusterCount = 0; // Dummy
		Box<Integer> tempRef_FinalClusterCount = new Box<>(FinalClusterCount);
		DoSynchronization.PipelineDistributedBroadcast(TemporaryLocalClustersforSynchronization.totalTransportedY_t_i, StorageforTransportedClusters.TotalTransportedY_t_i, null, null, TemporaryLocalClustersforSynchronization.totalTransportedCreatedIndex, StorageforTransportedClusters.TotalTransportedCreatedIndex, TemporaryLocalClustersforSynchronization.totalTransportedOriginalHost, StorageforTransportedClusters.TotalTransportedOriginalHost, tempRef_FinalClusterCount);
		FinalClusterCount = tempRef_FinalClusterCount.content;

		// Load P_t and Reset Sigmas if necessary
		for (int Clustersfromafar = 0; Clustersfromafar < StorageforTransportedClusters.SizeOfTransportedArray; Clustersfromafar++)
		{
			StorageforTransportedClusters.TotalTransportedP_t[Clustersfromafar] = StorageforTransportedClusters.TotalTransportedY_t_i[Clustersfromafar][Program.ParameterVectorDimension];
			if ((Program.SigmaMethod > 1) && (StorageforTransportedClusters.TotalTransportedStatus[Clustersfromafar][0] == 3))
			{
				Box<double[]> tempRef_Object = new Box<>(StorageforTransportedClusters.TotalTransportedSigma_t_i[Clustersfromafar]);
				Program.CalculateSigma(StorageforTransportedClusters.TotalTransportedY_t_i[Clustersfromafar],
                        tempRef_Object);
				StorageforTransportedClusters.TotalTransportedSigma_t_i[Clustersfromafar] = tempRef_Object.content;
			}
		}
		DAVectorUtility.StopSubTimer(7);

	} // End MinorSynchronizationTransportDistributedClusterCenters()

	//  Perform Major Synchronization Data Transport needed when Cluster structure changes
	//  The CurrentIteration value is incremented at such events and information about splits and deletes are propagated
	public static void MajorSynchronizationTransportDistributedClusterCenters() throws MPIException {

		int NumberofLocalDistributedClusters = 0;
		int H = DAVectorUtility.MPI_Rank;
		// string message = "";

		for (int ActiveClusterIndex = 0; ActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; ActiveClusterIndex++)
		{
			int RealClusterIndex = ClusteringSolution.RealClusterIndices[ActiveClusterIndex];
			if (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] < 2)
			{
				continue;
			}
			int PackedHost = ClusteringSolution.LocalHost[ActiveClusterIndex];
//C# TO JAVA CONVERTER WARNING: The right shift operator was not replaced by Java's logical right shift operator since the left operand was not confirmed to be of an unsigned type, but you should review whether the logical right shift operator (>>>) is more appropriate:
			int H1 = (PackedHost >> ClusteringSolution.PACKINGSHIFT) & ClusteringSolution.PACKINGMASK;
//C# TO JAVA CONVERTER WARNING: The right shift operator was not replaced by Java's logical right shift operator since the left operand was not confirmed to be of an unsigned type, but you should review whether the logical right shift operator (>>>) is more appropriate:
			int H2 = PackedHost >> (2 * ClusteringSolution.PACKINGSHIFT);
			if ((H1 == H) && (H2 == H))
			{
				if ((PackedHost & ClusteringSolution.PACKINGMASK) != H)
				{
					DAVectorUtility.printAndThrowRuntimeException(" Inconsistent Host Range " + (PackedHost & ClusteringSolution.PACKINGMASK) + " in Node " + H);

				}
				continue;
			}
			// message += RealClusterIndex.ToString() + "(" + ParallelClustering.RunningDAVectorSolt.LocalCreatedIndex[RealClusterIndex].ToString() + ") ";
			TemporaryLocalClustersforSynchronization.totalTransportedOriginalHost[NumberofLocalDistributedClusters] = PackedHost;
			TemporaryLocalClustersforSynchronization.totalTransportedCreatedIndex[NumberofLocalDistributedClusters] = ParallelClustering.runningSolution.LocalCreatedIndex[RealClusterIndex];
			TemporaryLocalClustersforSynchronization.totalTransportedStatus[NumberofLocalDistributedClusters][0] = 3;
			TemporaryLocalClustersforSynchronization.totalTransportedStatus[NumberofLocalDistributedClusters][1] = ParallelClustering.runningSolution.LocalSplitCreatedIndex[RealClusterIndex];
            System.arraycopy(ParallelClustering.runningSolution.Y_k_i_[RealClusterIndex], 0,
                    TemporaryLocalClustersforSynchronization.totalTransportedY_t_i[NumberofLocalDistributedClusters], 0,
                    Program.ParameterVectorDimension);
			TemporaryLocalClustersforSynchronization.totalTransportedY_t_i[NumberofLocalDistributedClusters][Program.ParameterVectorDimension] = ParallelClustering.runningSolution.P_k_[RealClusterIndex];
			++NumberofLocalDistributedClusters;
		} // End Loop over LocalClusterIndex

		TemporaryLocalClustersforSynchronization.sizeOfTransportedArray = NumberofLocalDistributedClusters;
		Program.ActualMaxTransportedClusterStorage = Math.max(Program.ActualMaxTransportedClusterStorage, NumberofLocalDistributedClusters);
		DistributedSynchronization.TransportviaPipeline DoSynchronization = new DistributedSynchronization.TransportviaPipeline(2, true, 3, 1, NumberofLocalDistributedClusters, 0);

		int FinalClusterCount = 0;
		Box<Integer> tempRef_FinalClusterCount = new Box<>(FinalClusterCount);
		DoSynchronization.PipelineDistributedBroadcast(TemporaryLocalClustersforSynchronization.totalTransportedY_t_i, StorageforTransportedClusters.TotalTransportedY_t_i, TemporaryLocalClustersforSynchronization.totalTransportedStatus, StorageforTransportedClusters.TotalTransportedStatus, TemporaryLocalClustersforSynchronization.totalTransportedCreatedIndex, StorageforTransportedClusters.TotalTransportedCreatedIndex, TemporaryLocalClustersforSynchronization.totalTransportedOriginalHost, StorageforTransportedClusters.TotalTransportedOriginalHost, tempRef_FinalClusterCount);
		FinalClusterCount = tempRef_FinalClusterCount.content;
		StorageforTransportedClusters.SizeOfTransportedArray = FinalClusterCount;
		Program.ActualMaxTransportedClusterStorage = Math.max(FinalClusterCount + 1, Program.ActualMaxTransportedClusterStorage);

		// Set P_t and sigmas which were not transported
		// Add Created Index if new
		for (int Clustersfromafar = 0; Clustersfromafar < StorageforTransportedClusters.SizeOfTransportedArray; Clustersfromafar++)
		{
			int CreatedIndex = StorageforTransportedClusters.TotalTransportedCreatedIndex[Clustersfromafar];
			if (ClusteringSolution.UniversalMapping[CreatedIndex] == null)
			{
				ClusteringSolution.UniversalMapping[CreatedIndex] = new ClusterIndirection(ClusteringSolution.CurrentIteration, -1 - Clustersfromafar);
			}
			ClusteringSolution.UniversalMapping[CreatedIndex].PackedHost = StorageforTransportedClusters.TotalTransportedOriginalHost[Clustersfromafar];
			StorageforTransportedClusters.TotalTransportedP_t[Clustersfromafar] = StorageforTransportedClusters.TotalTransportedY_t_i[Clustersfromafar][Program.ParameterVectorDimension];
			if (StorageforTransportedClusters.TotalTransportedStatus[Clustersfromafar][0] == 3)
			{
				Box<double[]> tempRef_Object = new Box<>(StorageforTransportedClusters.TotalTransportedSigma_t_i[Clustersfromafar]);
				Program.CalculateSigma(StorageforTransportedClusters.TotalTransportedY_t_i[Clustersfromafar],
                        tempRef_Object);
				StorageforTransportedClusters.TotalTransportedSigma_t_i[Clustersfromafar] = tempRef_Object.content;
			}
		}

	} // End MajorSynchronizationTransportDistributedClusterCenters()

	// public static int PrintOuts = 0;
	public static void SetClustersforaPoint(ClusteringSolution Solution) throws MPIException {
		DAVectorUtility.StartSubTimer(13);
		GlobalReductions.FindVectorDoubleSum FindDiagnosticSums_Points = new GlobalReductions.FindVectorDoubleSum(DAVectorUtility.ThreadCount, 8);

        // Note - parallel for
        launchHabaneroApp(() -> {
            forallChunked(0, DAVectorUtility.ThreadCount - 1,
                    (threadIndex) -> // End loop over Point dependent quantities
                    {
                        int[] ProposedClusters = new int[Program.maxNcentperPoint];
                        double[] ProposedDistances = new double[Program.maxNcentperPoint];
                        int indexlen = DAVectorUtility.PointsperThread[threadIndex];
                        int beginpoint = DAVectorUtility.StartPointperThread[threadIndex] - DAVectorUtility.PointStart_Process;
                        for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++) {
                            int NumberofClustersTobeUsed = 0;
                            // double SmallDistce = -1.0;
                            FindDiagnosticSums_Points.addapoint(threadIndex, 1.0, 0);
                            // int AssignmentStatus = DATriangleInequality.SetAssociatedCenters(alpha, out NumberofClustersTobeUsed, ProposedClusters, out SmallDistce, ProposedDistances );
                            Box<Integer> tempRef_NumberofClustersTobeUsed = new Box<>(NumberofClustersTobeUsed);
                            int AssignmentStatus = DATriangleInequality.SetAssociatedCenters(alpha,
                                    tempRef_NumberofClustersTobeUsed, ProposedClusters);
                            NumberofClustersTobeUsed = tempRef_NumberofClustersTobeUsed.content;
                            FindDiagnosticSums_Points.addapoint(threadIndex, 1.0, 2 + AssignmentStatus);
                            if (AssignmentStatus < 0) {
                                DAVectorUtility.printAndThrowRuntimeException(
                                        "Error in Triangle Inequality at Point " + (alpha + DAVectorUtility.PointStart_Process) + " Old Number " + Solution.NumClusters_alpha_[alpha]);

                            }
                            if (NumberofClustersTobeUsed <= 0) {
                                DAVectorUtility.printAndThrowRuntimeException(
                                        "Error due to zero New Number of Clusters Point " + (alpha + DAVectorUtility.PointStart_Process) + " Old Number " + Solution.NumClusters_alpha_[alpha]);

                            }
                            //  Now make list for storing back
                            double[] Mvalues = new double[NumberofClustersTobeUsed];
                            double Msum = 0.0;
                            for (int IndirectClusterIndex = 0; IndirectClusterIndex < NumberofClustersTobeUsed; IndirectClusterIndex++) {
                                Mvalues[IndirectClusterIndex] = 0.0;
                                if (Solution.NumClusters_alpha_[alpha] == 0) {
                                    continue;
                                }
                                int ActiveClusterIndex = ProposedClusters[IndirectClusterIndex];
                                int OldIndirectClusterIndex = Solution.MapClusterToIndirect(alpha, ActiveClusterIndex);
                                if (OldIndirectClusterIndex < 0) {
                                    continue;
                                }
                                Mvalues[IndirectClusterIndex] = Solution.M_alpha_kpointer_[alpha][OldIndirectClusterIndex];
                                Msum += Mvalues[IndirectClusterIndex];
                            }
                            if (Msum < 0.2) {
                                for (int IndirectClusterIndex = 0; IndirectClusterIndex < NumberofClustersTobeUsed; IndirectClusterIndex++) {
                                    Mvalues[IndirectClusterIndex] = 1.0 / (double) NumberofClustersTobeUsed;
                                }
                                Msum = 1.0;
                                FindDiagnosticSums_Points.addapoint(threadIndex, 1.0, 6);
                            }

                            //  Reset List
                            int NumberChange = Solution.NumClusters_alpha_[alpha] - NumberofClustersTobeUsed;
                    /*
                    if ( (Solution.Ncent_Global == 140) && (Solution.Temperature < 0.1) && (PrintOuts < 100) && ( NumberChange != 0) && (DAVectorUtility.MPI_Rank == 0) )
                    {
                        int bestindex = ClusteringSolution.RealClusterIndices[ProposedClusters[0]];
                        double deltaY = DAVectorParallelism.getNOTSquaredScaledDistancebetweenVectors(DATriangleInequality.CenterY_k_i_Current[bestindex],
                                DATriangleInequality.PointPosition[alpha], DATriangleInequality.CenterSigma_k_i_Current[bestindex]);

                        string message = " alpha " + alpha.ToString() + " size " + NumberofClustersTobeUsed.ToString() + " # Chg " + NumberChange.ToString()
                            + " Small " + SmallDistce.ToString("E4") + " Best " +  bestindex.ToString() + " " + deltaY.ToString("E4") + " new ";
                        for (int IndirectClusterIndex = 0; IndirectClusterIndex < NumberofClustersTobeUsed; IndirectClusterIndex++)
                        {
                            int ActiveClusterIndex = ProposedClusters[IndirectClusterIndex];
                            int OldIndirectClusterIndex = Solution.MapClusterToIndirect(alpha, ActiveClusterIndex);
                            if (OldIndirectClusterIndex >= 0)
                                continue;
                            int RealClusterIndex = ClusteringSolution.RealClusterIndices[ActiveClusterIndex];
                            deltaY = DAVectorParallelism.getNOTSquaredScaledDistancebetweenVectors(DATriangleInequality.CenterY_k_i_Current[RealClusterIndex],
                                DATriangleInequality.PointPosition[alpha], DATriangleInequality.CenterSigma_k_i_Current[RealClusterIndex]);
                            message += RealClusterIndex.ToString() + " " + deltaY.ToString("E4") + " " + ProposedDistances[IndirectClusterIndex].ToString("E4") + " ";
                        }
                        message += " old ";

                        for (int IndirectClusterIndex = 0; IndirectClusterIndex < Solution.NumClusters_alpha_[alpha]; IndirectClusterIndex++)
                        {
                            int RealClusterIndex1 = -1;
                            int ActiveClusterIndex1 = -1;
                            int RemoteIndex = 0;
                            VectorAnnealIterate.ClusterPointersforaPoint(alpha, IndirectClusterIndex, ref RealClusterIndex1, ref ActiveClusterIndex1, ref RemoteIndex);
                            for (int IndirectClusterIndex1 = 0; IndirectClusterIndex1 < NumberofClustersTobeUsed; IndirectClusterIndex1++)
                            {
                                if (ActiveClusterIndex1 == ProposedClusters[IndirectClusterIndex1])
                                {
                                    ActiveClusterIndex1 = -1;
                                    break;
                                }
                            }
                            if (ActiveClusterIndex1 < 0)
                                continue;
                            deltaY = DAVectorParallelism.getNOTSquaredScaledDistancebetweenVectors(DATriangleInequality.CenterY_k_i_Current[RealClusterIndex1],
                                DATriangleInequality.PointPosition[alpha], DATriangleInequality.CenterSigma_k_i_Current[RealClusterIndex1]);
                            message += RealClusterIndex1.ToString() + " " + deltaY.ToString("E4") + " " ;
                        }
                        DAVectorUtility.SALSAPrint(0, message);
                        ++PrintOuts;
                    }
                    */
                            FindDiagnosticSums_Points.addapoint(threadIndex, (double) NumberChange, 7);
                            Solution.NumClusters_alpha_[alpha] = NumberofClustersTobeUsed;
                            for (int IndirectClusterIndex = 0; IndirectClusterIndex < NumberofClustersTobeUsed; IndirectClusterIndex++) {
                                Solution.M_alpha_kpointer_[alpha][IndirectClusterIndex] = Mvalues[IndirectClusterIndex] / Msum;
                                Solution.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex] = ClusteringSolution.MapActivetoCreatedIndex(
                                        ProposedClusters[IndirectClusterIndex], Solution);
                            }
                            for (int IndirectClusterIndex = 0; IndirectClusterIndex < NumberofClustersTobeUsed; IndirectClusterIndex++) {
                                ParallelClustering.runningSolution.LegalCluster(alpha, IndirectClusterIndex);
                            }
                        }

                    });
        });

        Solution.DiffMalpha_k_Set = -1;

		FindDiagnosticSums_Points.sumoverthreadsandmpi();
		for (int DiagnosticLoop = 0; DiagnosticLoop < 8; DiagnosticLoop++)
		{
			Program.ClustersperCenterDiagnostics[DiagnosticLoop] += FindDiagnosticSums_Points.TotalVectorSum[DiagnosticLoop];
		}
		VectorAnnealIterate.NumberCountChanges = Program.ClustersperCenterDiagnostics[7] / Program.ClustersperCenterDiagnostics[0];
		DAVectorUtility.StopSubTimer(13);

    } // End SetClustersforaPoint(ClusteringSolution Solution)

	// Control Major Synchronization for distributed and global case
	//  ongoing = true is normal case and works whether distributed or global mode
	//  ongoing = false switches INTO distributed mode
	public static void ManageMajorSynchronization(boolean ongoing) throws MPIException {
		//  Set up Active Cluster Lists
		DAVectorUtility.StartSubTimer(6);
		ParallelClustering.runningSolution.SetActiveClusters();
		VectorAnnealIterate.NumberCountChanges = 0.0;

		//  Case of non distributed operation
		if (ongoing && (!ParallelClustering.runningSolution.DistributedExecutionMode))
		{
			++ClusteringSolution.CurrentIteration; // Increment Iteration Number
			DAVectorUtility.SALSAPrint(0, "Debug:ClusteringSolution.NumberLocalActiveClusters  " + ClusteringSolution.NumberLocalActiveClusters);
			DAVectorUtility.SALSAPrint(0, "Debug:ClusteringSolution.UniversalMapping Length " + ClusteringSolution.UniversalMapping.length);
			for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++)
			{
				int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
				int CreatedIndex = ParallelClustering.runningSolution.LocalCreatedIndex[RealClusterIndex];
				DAVectorUtility.SALSAPrint(0, "Debug:CreatedIndex " + CreatedIndex);
				DAVectorUtility.SALSAPrint(0, "Debug:ClusteringSolution.UniversalMapping " + ClusteringSolution.UniversalMapping[CreatedIndex]);
				ClusteringSolution.UniversalMapping[CreatedIndex].Availability = 1 + RealClusterIndex;
				ClusteringSolution.UniversalMapping[CreatedIndex].IterationSet = ClusteringSolution.CurrentIteration;
				ClusteringSolution.LocalHost[LocalActiveClusterIndex] = 0;
			}
			if (Program.UseTriangleInequality_DA > 0)
			{
				ParallelClustering.runningSolution.SetClusterWidths();
				DATriangleInequality.NextIteration();
				DistributedClusteringSolution.SetClustersforaPoint(ParallelClustering.runningSolution);
			}
			ParallelClustering.runningSolution.DiffMalpha_k_Set = -1;

			ParallelClustering.runningSolution.SetClusterSizes();
			ParallelClustering.runningSolution.SetClusterWidths();
			++Program.NumberMajorSynchs1;
			DAVectorUtility.StopSubTimer(6);
			return;
		}

		++Program.NumberMajorSynchs2;
		Program.CountClusters2 += ParallelClustering.runningSolution.Ncent_ThisNode;

		// Reset Host Information
		boolean UseCluster1 = true;
		if (ongoing)
		{
			UseCluster1 = false;
		}
		ClusterHostlinkage(UseCluster1);

		// If update, Transport Cluster Information and populate local storage in StorageforTransportedClusters
		if (ongoing)
		{
			MajorSynchronizationTransportDistributedClusterCenters();
		}

		//  Set up Distributed Mode for first time setting StorageforTransportedClusters from local data
		else
		{
			SetupDistributedMode();
		}

		//  Process Clusters moved between nodes
		ProcessChangedClusters(ongoing);
		DAVectorUtility.StopSubTimer(6);

		//  Set up a new iteration
		SetupNewIteration();

		//  Initialize storage pointers
		DAVectorUtility.StartSubTimer(11);
		final int LocalNumberClusters = ClusteringSolution.NumberLocalActiveClusters;
		int TotalNumberClusters = StorageforTransportedClusters.SizeOfTransportedArray + LocalNumberClusters;
		ClusteringSolution.NumberAvailableActiveClusters = TotalNumberClusters;
		if (TotalNumberClusters <= 0)
		{
			DAVectorUtility.printAndThrowRuntimeException("No Clusters in this Process " + DAVectorUtility.MPI_Rank + " Remote " + StorageforTransportedClusters.SizeOfTransportedArray + " Local " + LocalNumberClusters);

		}


		final Range[] TotalClusterRanges = RangePartitioner.Partition(TotalNumberClusters, DAVectorUtility.ThreadCount);

		//  Initialize Mapping of Local and Remote Clusters to Node and Thread Accumulations
		//  Initialize ClusteringSolution.UniversalMapping for this Node at this NEW Iteration number
		//  This assumes that any child cluster is created fully but values of Malpha are not set any where
        // Note - parallel for
        launchHabaneroApp(() -> {
            forallChunked(0, DAVectorUtility.ThreadCount - 1,
                    (threadIndex) -> // End Parallel Section over Total Number of Clusters
                    {
                        //  It may not have any thread indices if points for this cluster are outside this node
                        int BeginClusterIndex = TotalClusterRanges[threadIndex].getStartIndex();
                        int ClusterIndexLength = TotalClusterRanges[threadIndex].getLength();
                        for (int LocalActiveClusterIndex = BeginClusterIndex; LocalActiveClusterIndex < BeginClusterIndex + ClusterIndexLength; LocalActiveClusterIndex++) {
                            int RemoteIndex = -1;
                            int RealClusterIndex = -1;
                            if (LocalActiveClusterIndex >= LocalNumberClusters) {
                                RemoteIndex = LocalActiveClusterIndex - LocalNumberClusters;
                            } else {
                                RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
                            }
                            if (RemoteIndex == -1) {
                                ClusteringSolution.LocalNodeAccPosition[LocalActiveClusterIndex] = 0;
                                int CreatedIndex = ParallelClustering.runningSolution.LocalCreatedIndex[RealClusterIndex];
                                ClusteringSolution.UniversalMapping[CreatedIndex].Availability = 1 + RealClusterIndex;
                                ClusteringSolution.UniversalMapping[CreatedIndex].IterationSet = ClusteringSolution.CurrentIteration;
                            } else {
                                StorageforTransportedClusters.TransportedNodeAccPosition[RemoteIndex] = -1;
                                int CreatedIndex = StorageforTransportedClusters.TotalTransportedCreatedIndex[RemoteIndex];
                                ClusteringSolution.UniversalMapping[CreatedIndex].Availability = -1 - RemoteIndex;
                                ClusteringSolution.UniversalMapping[CreatedIndex].IterationSet = ClusteringSolution.CurrentIteration;
                            }
                            for (int ThreadNo = 0; ThreadNo < DAVectorUtility.ThreadCount; ThreadNo++) {
                                if (RemoteIndex == -1) {
                                    ClusteringSolution.LocalThreadAccPosition[LocalActiveClusterIndex][ThreadNo] = -1;
                                } else {
                                    StorageforTransportedClusters.TransportedThreadAccPosition[RemoteIndex][ThreadNo] = -1;
                                }
                            }
                        }

                    }
            ); // End Parallel Section over Total Number of Clusters
        });

        //  Loop over points in this node
		//  Now find which clusters are actually used and remove clusters that have disappeared
		//  Process Pending Cluster Splits
		//  Set up Accumulation arrays by setting any used to 0 (initialized to -1 above)

        // Note - parallel for
        launchHabaneroApp(() -> {
            forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) ->
            {
                double[] Malpha_perpoint = new double[ClusteringSolution.NumberAvailableActiveClusters];
                int[] CreatedIndex_perpoint = new int[ClusteringSolution.NumberAvailableActiveClusters];
                int indexlen = DAVectorUtility.PointsperThread[threadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[threadIndex] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++) {
                    int IndirectSize = ParallelClustering.runningSolution.NumClusters_alpha_[alpha];
                    int NewIndirectSize = 0;
                    double LostSize = 0.0;
                    boolean changed = false;
                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++) {
                        int CreatedIndex = ParallelClustering.runningSolution.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex];
                        if ((CreatedIndex < 0) || (ClusteringSolution.UniversalMapping[CreatedIndex].IterationSet != ClusteringSolution.CurrentIteration)) {
                            LostSize += ParallelClustering.runningSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex];
                            changed = true;
                            /* if ((alpha + DAVectorUtility.PointStart_Process == 22625) && (ClusteringSolution.CurrentIteration > 6500))
                            {
                                DAVectorUtility.SALSAFullPrint(0, " Lost " + IndirectClusterIndex.ToString() + " " + ParallelClustering.runningSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex].ToString("F4")
                                    + " Mapped " + ClusteringSolution.UniversalMapping[CreatedIndex].Availability);
                            } */
                            continue;
                        }
                        /* if ((alpha + DAVectorUtility.PointStart_Process == 22625) && (ClusteringSolution.CurrentIteration > 6500))
                        {
                            DAVectorUtility.SALSAFullPrint(0, " OK " + IndirectClusterIndex.ToString() + " " + ParallelClustering.runningSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex].ToString("F4")
                                + " Mapped " + ClusteringSolution.UniversalMapping[CreatedIndex].Availability);
                        } */
                        Malpha_perpoint[NewIndirectSize] = ParallelClustering.runningSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex];
                        CreatedIndex_perpoint[NewIndirectSize] = CreatedIndex;
                        int AvailabilityIndex = ClusteringSolution.UniversalMapping[CreatedIndex].Availability;
                        int RemoteIndex = -1;
                        int SplitInfo = 0;
                        if (AvailabilityIndex > 0) {
                            --AvailabilityIndex;
                            if (ParallelClustering.runningSolution.LocalStatus[AvailabilityIndex] < 0) {
                                LostSize += ParallelClustering.runningSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex];
                                changed = true;
                                continue;
                            }
                            int ActiveClusterIndex = ClusteringSolution.ActiveClusterIndices[AvailabilityIndex];
                            SplitInfo = ParallelClustering.runningSolution.LocalSplitCreatedIndex[AvailabilityIndex];
                            ClusteringSolution.LocalNodeAccPosition[ActiveClusterIndex] = 0;
                            ClusteringSolution.LocalThreadAccPosition[ActiveClusterIndex][threadIndex] = 0;
                        } else {
                            RemoteIndex = -AvailabilityIndex - 1;
                            SplitInfo = StorageforTransportedClusters.TotalTransportedStatus[RemoteIndex][1];
                            StorageforTransportedClusters.TransportedNodeAccPosition[RemoteIndex] = 0;
                            StorageforTransportedClusters.TransportedThreadAccPosition[RemoteIndex][threadIndex] = 0;
                        }
                        /* if ((alpha + DAVectorUtility.PointStart_Process == 22625) && (ClusteringSolution.CurrentIteration > 6500))
                        {
                            DAVectorUtility.SALSAFullPrint(0, " OK Split " + IndirectClusterIndex.ToString() + " " + RemoteIndex.ToString() + " " + SplitInfo.ToString());
                        } */
                        if (SplitInfo > 0) {
                            DAVectorUtility.printAndThrowRuntimeException(
                                    " Created Index " + CreatedIndex + " has illegal split flag " + (SplitInfo - 1));

                        }
                        if (SplitInfo < 0) {
                            Malpha_perpoint[NewIndirectSize] *= 0.5;
                            ++NewIndirectSize;
                            Malpha_perpoint[NewIndirectSize] = Malpha_perpoint[NewIndirectSize - 1];
                            int ChildCreatedIndex = -SplitInfo - 1;
                            CreatedIndex_perpoint[NewIndirectSize] = ChildCreatedIndex;
                            changed = true;
                            if (ClusteringSolution.UniversalMapping[ChildCreatedIndex].IterationSet != ClusteringSolution.CurrentIteration) {
                                DAVectorUtility.printAndThrowRuntimeException(
                                        " Child Created Index " + ChildCreatedIndex + " is set up wrong -- parent " + CreatedIndex);

                            }
                            int ChildClusterIndex = ClusteringSolution.UniversalMapping[ChildCreatedIndex].Availability;
                            if (ChildClusterIndex > 0) {
                                int ActiveChildClusterIndex = ClusteringSolution.ActiveClusterIndices[ChildClusterIndex - 1];
                                ClusteringSolution.LocalNodeAccPosition[ActiveChildClusterIndex] = 0;
                                ClusteringSolution.LocalThreadAccPosition[ActiveChildClusterIndex][threadIndex] = 0;
                            } else {
                                StorageforTransportedClusters.TransportedNodeAccPosition[-ChildClusterIndex - 1] = 0;
                                StorageforTransportedClusters.TransportedThreadAccPosition[-ChildClusterIndex - 1][threadIndex] = 0;
                            }
                        }
                        ++NewIndirectSize;
                    }
                    //  Reset Cluster Information for this point
                    if (changed) {
                        if (NewIndirectSize > ClusteringSolution.TargetMinimumClustersperPoint) {
                            ParallelClustering.runningSolution.NumClusters_alpha_[alpha] = NewIndirectSize;
                            boolean multiplication = true;
                            double fudge;
                            if (LostSize >= 0.99) {
                                multiplication = false;
                                fudge = 0.0;
                            } else {
                                fudge = 1.0 / (1.0 - LostSize);
                            }
                            for (int IndirectClusterIndex = 0; IndirectClusterIndex < NewIndirectSize; IndirectClusterIndex++) {
                                ParallelClustering.runningSolution.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex] = CreatedIndex_perpoint[IndirectClusterIndex];
                                if (multiplication) {
                                    ParallelClustering.runningSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex] = Malpha_perpoint[IndirectClusterIndex] * fudge;
                                } else {
                                    ParallelClustering.runningSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex] = 1.0 / NewIndirectSize;
                                }
                            }
                        } else {
                            ParallelClustering.runningSolution.SetClustersforaPoint(alpha);
                            NewIndirectSize = ParallelClustering.runningSolution.NumClusters_alpha_[alpha];
                            for (int IndirectClusterIndex = 0; IndirectClusterIndex < NewIndirectSize; IndirectClusterIndex++) {
                                int CreatedIndex = ParallelClustering.runningSolution.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex];
                                int AvailabilityIndex = ClusteringSolution.UniversalMapping[CreatedIndex].Availability;
                                /* if ((alpha + DAVectorUtility.PointStart_Process == 22625) && (ClusteringSolution.CurrentIteration > 6500))
                                {
                                    DAVectorUtility.SALSAFullPrint(0, " Redone " + IndirectClusterIndex.ToString() + " " + ParallelClustering.runningSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex].ToString("F4")
                                        + " Created Index " + CreatedIndex.ToString() + " Mapped " + AvailabilityIndex.ToString() );
                                } */
                                if (AvailabilityIndex > 0) {
                                    --AvailabilityIndex;
                                    int ActiveClusterIndex = ClusteringSolution.ActiveClusterIndices[AvailabilityIndex];
                                    ClusteringSolution.LocalNodeAccPosition[ActiveClusterIndex] = 0;
                                    ClusteringSolution.LocalThreadAccPosition[ActiveClusterIndex][threadIndex] = 0;
                                } else {
                                    int RemoteIndex = -AvailabilityIndex - 1;
                                    StorageforTransportedClusters.TransportedNodeAccPosition[RemoteIndex] = 0;
                                    StorageforTransportedClusters.TransportedThreadAccPosition[RemoteIndex][threadIndex] = 0;
                                }
                            }
                        }
                    }
                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < NewIndirectSize; IndirectClusterIndex++) {
                        ParallelClustering.runningSolution.LegalCluster(alpha, IndirectClusterIndex);
                    }
                }

            }); // End loop Processing Point dependent quantities
        });

        //  Unset any split information and Set Accumulation Orders and start setting NodeAccMetaData
		//  Done Sequentially in each node
		int NumberNodeClusterAccIndices = 0;
		int[] NumberThreadClusterAccIndices = new int[DAVectorUtility.ThreadCount];
		for (int ThreadNo = 0; ThreadNo < DAVectorUtility.ThreadCount; ThreadNo++)
		{
			NumberThreadClusterAccIndices[ThreadNo] = 0;
		}

		for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++)
		{
			int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
			ParallelClustering.runningSolution.LocalSplitCreatedIndex[RealClusterIndex] = 0;
			if (ClusteringSolution.LocalNodeAccPosition[LocalActiveClusterIndex] < 0)
			{
				continue;
			}
			ClusteringSolution.LocalNodeAccPosition[LocalActiveClusterIndex] = NumberNodeClusterAccIndices;
			NodeAccMetaData.NodeAccumulationCreatedIndices[NumberNodeClusterAccIndices] = ParallelClustering.runningSolution.LocalCreatedIndex[RealClusterIndex];
			++NumberNodeClusterAccIndices;

			for (int ThreadNo = 0; ThreadNo < DAVectorUtility.ThreadCount; ThreadNo++)
			{
				if (ClusteringSolution.LocalThreadAccPosition[LocalActiveClusterIndex][ThreadNo] < 0)
				{
					continue;
				}
				ClusteringSolution.LocalThreadAccPosition[LocalActiveClusterIndex][ThreadNo] = NumberThreadClusterAccIndices[ThreadNo];
				++NumberThreadClusterAccIndices[ThreadNo];
			}
		}

		for (int RemoteIndex = 0; RemoteIndex < StorageforTransportedClusters.SizeOfTransportedArray; RemoteIndex++)
		{
			StorageforTransportedClusters.TotalTransportedStatus[RemoteIndex][1] = 0;
			if (StorageforTransportedClusters.TransportedNodeAccPosition[RemoteIndex] < 0)
			{
				continue;
			}
			StorageforTransportedClusters.TransportedNodeAccPosition[RemoteIndex] = NumberNodeClusterAccIndices;
			NodeAccMetaData.NodeAccumulationCreatedIndices[NumberNodeClusterAccIndices] = StorageforTransportedClusters.TotalTransportedCreatedIndex[RemoteIndex];
			++NumberNodeClusterAccIndices;

			for (int ThreadNo = 0; ThreadNo < DAVectorUtility.ThreadCount; ThreadNo++)
			{
				if (StorageforTransportedClusters.TransportedThreadAccPosition[RemoteIndex][ThreadNo] < 0)
				{
					continue;
				}
				StorageforTransportedClusters.TransportedThreadAccPosition[RemoteIndex][ThreadNo] = NumberThreadClusterAccIndices[ThreadNo];
				++NumberThreadClusterAccIndices[ThreadNo];
			}
		}

		//  Node Accumulation Book Keeping
		NodeAccMetaData.NumberofPointsperNode = NumberNodeClusterAccIndices;
		Program.ActualMaxNumberAccumulationsperNode = Math.max(Program.ActualMaxNumberAccumulationsperNode, NumberNodeClusterAccIndices);
		ParallelNodeAccumulationRanges = RangePartitioner.Partition(NumberNodeClusterAccIndices, DAVectorUtility.ThreadCount);
        System.arraycopy(NumberThreadClusterAccIndices, 0, NodeAccMetaData.NumberofPointsperThread, 0,
                DAVectorUtility.ThreadCount);

		//  Associate Thread accumulation to Node Accumulation arrays
		AssociateThreadtoNodeAccumulation();
		DAVectorUtility.StopSubTimer(11);

		DAVectorUtility.StartSubTimer(12);
		ParallelClustering.runningSolution.SetClusterSizes();
		ParallelClustering.runningSolution.SetClusterWidths();
		DAVectorUtility.StopSubTimer(12);

	} // End ManageMajorSynchronization()

	//  Associate Thread Accumnulation Positions to Node Accumulation Positions
	//  Also set Status and Hosts
	//  Parallel over Node Accumulation Positions
	public static void AssociateThreadtoNodeAccumulation()
	{
        // Note - parallel for
        launchHabaneroApp(() -> {
            forallChunked(0, DAVectorUtility.ThreadCount - 1,
                    (threadIndex) ->  //  End Parallel Section over NodeAccumulationIndex
                    {
                        //  Set role and message structure of Node Accumulation Points
                        int beginindex = ParallelNodeAccumulationRanges[threadIndex].getStartIndex();
                        int indexlength = ParallelNodeAccumulationRanges[threadIndex].getLength();
                        int NodeStoragePosition = -1;
                        int TransportedStoragePosition = -1;
                        int NodeAccumulationPosition = -1;
                        int ThreadAccumulationPosition = -1;
                        for (int NodeAccumulationIndex = beginindex; NodeAccumulationIndex < beginindex + indexlength; NodeAccumulationIndex++) {
                            int createdindex = NodeAccMetaData.NodeAccumulationCreatedIndices[NodeAccumulationIndex];
                            for (int ThreadNo = 0; ThreadNo < DAVectorUtility.ThreadCount; ThreadNo++) {
                                Box<Integer> tempRef_NodeStoragePosition = new Box<>(
                                        NodeStoragePosition);
                                Box<Integer> tempRef_TransportedStoragePosition = new Box<>(
                                        TransportedStoragePosition);
                                Box<Integer> tempRef_NodeAccumulationPosition = new Box<>(
                                        NodeAccumulationPosition);
                                Box<Integer> tempRef_ThreadAccumulationPosition = new Box<>(
                                        ThreadAccumulationPosition);
                                int isitOK = IndicesperCluster(createdindex, ThreadNo, tempRef_NodeStoragePosition,
                                        tempRef_TransportedStoragePosition, tempRef_NodeAccumulationPosition,
                                        tempRef_ThreadAccumulationPosition);
                                NodeStoragePosition = tempRef_NodeStoragePosition.content;
                                TransportedStoragePosition = tempRef_TransportedStoragePosition.content;
                                NodeAccumulationPosition = tempRef_NodeAccumulationPosition.content;
                                ThreadAccumulationPosition = tempRef_ThreadAccumulationPosition.content;
                                if (isitOK < 0) {
                                    DAVectorUtility.printAndThrowRuntimeException(
                                            " Inconsistent Created Index " + createdindex + " in Node  Accumulation " + NodeAccumulationIndex);

                                }
                                if (NodeAccumulationPosition != NodeAccumulationIndex) {
                                    DAVectorUtility.printAndThrowRuntimeException(
                                            " Inconsistent Created Index " + createdindex + " " + DAVectorUtility.MPI_Rank + " " + ClusteringSolution.NumberLocalActiveClusters + " in Node  Accumulation 1 " + NodeAccumulationIndex + " and Node  Accumulation 2 " + NodeAccumulationPosition);

                                }
                                NodeAccMetaData.AccumulationNodetoThreadClusterAssociations[NodeAccumulationIndex][ThreadNo] = ThreadAccumulationPosition;
                            }
                            if (NodeStoragePosition >= 0) {
                                int ActiveClusterIndex = ClusteringSolution.ActiveClusterIndices[NodeStoragePosition];
                                int LocalClusterStatus = ParallelClustering.runningSolution.LocalStatus[NodeStoragePosition];
                                if ((LocalClusterStatus <= -1) || (LocalClusterStatus == 3)) {
                                    DAVectorUtility.printAndThrowRuntimeException(
                                            " Created Index " + createdindex + " in Node  Accumulation " + NodeAccumulationIndex + " and Node  Storage " + NodeStoragePosition + " Deleted");

                                }
                                NodeAccMetaData.NodeAccumulationClusterStatus[NodeAccumulationIndex] = LocalClusterStatus;
                                NodeAccMetaData.NodeAccumulationClusterHosts[NodeAccumulationIndex] = ClusteringSolution.LocalHost[ActiveClusterIndex];
                            } else {
                                if (TransportedStoragePosition < 0) {
                                    DAVectorUtility.printAndThrowRuntimeException(
                                            " Created Index " + createdindex + " in Node Accumulation " + NodeAccumulationIndex + " in neither category");

                                }
                                NodeAccMetaData.NodeAccumulationClusterStatus[NodeAccumulationIndex] = 3;
                                NodeAccMetaData.NodeAccumulationClusterHosts[NodeAccumulationIndex] = StorageforTransportedClusters.TotalTransportedOriginalHost[TransportedStoragePosition];
                            }
                        }

                    }
            );
        });
    } // End AssociateThreadtoNodeAccumulation()

	public static void SetupNewIteration() throws MPIException {
		++ClusteringSolution.CurrentIteration; // Increment Iteration Number
		int Testiteration = ClusteringSolution.CurrentIteration;

        // Note - MPI Call - Broadcast - int
        if (DAVectorUtility.MPI_Size > 1){
//            DAVectorUtility.MPI_communicator.<Integer>Broadcast(tempRef_Testiteration, 0);
            Testiteration = DAVectorUtility.mpiOps.broadcast(Testiteration,0);
        }
		if (Testiteration != ClusteringSolution.CurrentIteration)
		{
			DAVectorUtility.printAndThrowRuntimeException(" Inconsistent Iteration " + Testiteration + " " + ClusteringSolution.CurrentIteration + " with rank " + DAVectorUtility.MPI_Rank);

		}
		int NumberTransported = 0;
		for (int ActiveClusterIndex = 0; ActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; ActiveClusterIndex++)
		{
			int RealClusterIndex = ClusteringSolution.RealClusterIndices[ActiveClusterIndex];
			if (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] < 2)
			{
				continue;
			}
			int PackedHost = ClusteringSolution.LocalHost[ActiveClusterIndex];
//C# TO JAVA CONVERTER WARNING: The right shift operator was not replaced by Java's logical right shift operator since the left operand was not confirmed to be of an unsigned type, but you should review whether the logical right shift operator (>>>) is more appropriate:
			int H1 = (PackedHost >> ClusteringSolution.PACKINGSHIFT) & ClusteringSolution.PACKINGMASK;
//C# TO JAVA CONVERTER WARNING: The right shift operator was not replaced by Java's logical right shift operator since the left operand was not confirmed to be of an unsigned type, but you should review whether the logical right shift operator (>>>) is more appropriate:
			int H2 = PackedHost >> (2 * ClusteringSolution.PACKINGSHIFT);
			if ((H1 == DAVectorUtility.MPI_Rank) && (H2 == DAVectorUtility.MPI_Rank))
			{
				continue;
			}
			++NumberTransported;
		}
		if (VectorAnnealIterate.EMIterationCount % Program.PrintInterval == 0)
		{
			String endinfo = "";
			if (!Program.CalculateIndividualWidths)
			{
				endinfo = " Average Width " + String.format("%1$4.3E", ParallelClustering.runningSolution.TotaloverVectorIndicesAverageWidth);
			}
			else
			{
				endinfo = " Average Widths ";
				for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
				{
					endinfo += String.format("%1$4.3E", ParallelClustering.runningSolution.AverageWidth[VectorIndex]) + " ";
				}
			}
			String LocalMessage = "Clstrs " + ParallelClustering.runningSolution.Ncent_ThisNode + " Dltd " + ClusteringSolution.ClustersDeleted + " Splt " + ClusteringSolution.ClustersSplit + " Mvd " + ClusteringSolution.ClustersMoved + " FromRmt " + StorageforTransportedClusters.SizeOfTransportedArray + " ToRmt " + NumberTransported;
			String OverallMessage = "Major Synch Iteration " + ClusteringSolution.CurrentIteration + " Total Clstrs " + ParallelClustering.runningSolution.Ncent_Global + " Temperature " + String.format("%1$5.4f", ParallelClustering.runningSolution.Temperature) + " " + endinfo;
			if (ParallelClustering.runningSolution.SpongeCluster >= 0)
			{
				OverallMessage = "Sponge Count " + String.format("%1$3.2f", ParallelClustering.runningSolution.C_k_[ParallelClustering.runningSolution.SpongeCluster]) + " Factor " + String.format("%1$3.2f", Program.SpongeFactor) + " " + OverallMessage;
			}
			DAVectorUtility.SALSASyncPrint(1, OverallMessage, LocalMessage);
			DAVectorUtility.SALSASyncPrint(2, "Host", HostLinkageMessage);
		}

		ClusteringSolution.ClustersDeleted = 0;
		ClusteringSolution.ClustersMoved = 0;
		ClusteringSolution.ClustersSplit = 0;

    } // End SetupNewIteration()

	public static void ProcessChangedClusters(boolean ongoing) throws MPIException {
		// Deal with deleted clusters removed as zero size
		// Process Any Clusters Moved from afar to here
		// This involves moving clusters from transported to local storagr and deleting moved clusters from local
		boolean changedclusters = false;
		ClusteringSolution.ClustersDeleted = 0;
		ClusteringSolution.ClustersMoved = 0;

		//  Delete Clusters transported elsewhere or those deleted
		for (int RealClusterIndex = 0; RealClusterIndex < ParallelClustering.runningSolution.Ncent_ThisNode; RealClusterIndex++)
		{
			if (RealClusterIndex == ParallelClustering.runningSolution.SpongeCluster)
			{
				continue;
			}
			int CreatedIndex = ParallelClustering.runningSolution.LocalCreatedIndex[RealClusterIndex];
			if (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] < 0)
			{
				changedclusters = true;
				continue;
			}
			int NewHost = ClusteringSolution.LocalHost[RealClusterIndex] & ClusteringSolution.PACKINGMASK;
			if (NewHost == DAVectorUtility.MPI_Rank)
			{
				ClusteringSolution.UniversalMapping[CreatedIndex].Ymapping = ParallelClustering.runningSolution.Y_k_i_[RealClusterIndex][0];
				continue;
			}

			ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] = -1;
			++ClusteringSolution.ClustersMoved;
			changedclusters = true;
		}
		boolean GlobalChangedClusters = changedclusters;
		DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming5);
        if (DAVectorUtility.MPI_Size > 1){
            // Note - MPI Call - Allreduce - boolean - logical OR
//            GlobalChangedClusters = DAVectorUtility.MPI_communicator.<Boolean>Allreduce(GlobalChangedClusters, Operation<Boolean>.LogicalOr);
            GlobalChangedClusters = DAVectorUtility.mpiOps.allReduce(GlobalChangedClusters, MPI.LOR);
        }
		DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming5);

		if (changedclusters)
		{
			ParallelClustering.runningSolution.CleanupClusters();
		}
		if (GlobalChangedClusters)
		{
			ParallelClustering.runningSolution.SetActiveClusters();
			for (int ActiveClusterIndex = 0; ActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; ActiveClusterIndex++)
			{
				int RealClusterIndex = ClusteringSolution.RealClusterIndices[ActiveClusterIndex];
				ClusteringSolution.LocalHost[ActiveClusterIndex] = ClusteringSolution.UniversalMapping[ParallelClustering.runningSolution.LocalCreatedIndex[RealClusterIndex]].PackedHost;
			}
		}

		//  Move Clusters from afar into local storage with status 2
		if (DAVectorUtility.MPI_Size == 1)
		{
			return;
		}

		int NewTransportedCount = 0;
		for (int RemoteIndex = 0; RemoteIndex < StorageforTransportedClusters.SizeOfTransportedArray; RemoteIndex++)
		{
			int RemotePackedHost = StorageforTransportedClusters.TotalTransportedOriginalHost[RemoteIndex];
			if ((RemotePackedHost & ClusteringSolution.PACKINGMASK) == DAVectorUtility.MPI_Rank)
			{ // Move this cluster to be a local distributed cluster
				//  C_k_ must be set later
				int NewCenterIndex = ParallelClustering.runningSolution.Ncent_ThisNode;
				int CreatedIndex = StorageforTransportedClusters.TotalTransportedCreatedIndex[RemoteIndex];
				ParallelClustering.runningSolution.LocalCreatedIndex[NewCenterIndex] = CreatedIndex;
				ParallelClustering.runningSolution.P_k_[NewCenterIndex] = StorageforTransportedClusters.TotalTransportedP_t[RemoteIndex];
				ParallelClustering.runningSolution.LocalStatus[NewCenterIndex] = 2;
				ParallelClustering.runningSolution.LocalSplitCreatedIndex[NewCenterIndex] = StorageforTransportedClusters.TotalTransportedStatus[RemoteIndex][1];
				ParallelClustering.runningSolution.Splittable_k_[NewCenterIndex] = -1;
				ParallelClustering.runningSolution.SplitPriority_k_[NewCenterIndex] = 2;
				ClusteringSolution.UniversalMapping[CreatedIndex].PackedHost = StorageforTransportedClusters.TotalTransportedOriginalHost[RemoteIndex];

				for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
				{
					ParallelClustering.runningSolution.Y_k_i_[NewCenterIndex][VectorIndex] = StorageforTransportedClusters.TotalTransportedY_t_i[RemoteIndex][VectorIndex];
					ParallelClustering.runningSolution.Sigma_k_i_[NewCenterIndex][VectorIndex] = StorageforTransportedClusters.TotalTransportedSigma_t_i[RemoteIndex][VectorIndex];
				}
				++ParallelClustering.runningSolution.Ncent_ThisNode;
				continue;
			}
			if (NewTransportedCount < RemoteIndex)
			{
				StorageforTransportedClusters.TotalTransportedCreatedIndex[NewTransportedCount] = StorageforTransportedClusters.TotalTransportedCreatedIndex[RemoteIndex];
				StorageforTransportedClusters.TotalTransportedOriginalHost[NewTransportedCount] = StorageforTransportedClusters.TotalTransportedOriginalHost[RemoteIndex];
				StorageforTransportedClusters.TotalTransportedP_t[NewTransportedCount] = StorageforTransportedClusters.TotalTransportedP_t[RemoteIndex];
				for (int VectorIndex = 0; VectorIndex <= Program.ParameterVectorDimension; VectorIndex++)
				{
					StorageforTransportedClusters.TotalTransportedY_t_i[NewTransportedCount][VectorIndex] = StorageforTransportedClusters.TotalTransportedY_t_i[RemoteIndex][VectorIndex];
					if (VectorIndex < Program.ParameterVectorDimension)
					{
						StorageforTransportedClusters.TotalTransportedSigma_t_i[NewTransportedCount][VectorIndex] = StorageforTransportedClusters.TotalTransportedSigma_t_i[RemoteIndex][VectorIndex];
					}
				}
				StorageforTransportedClusters.TotalTransportedY_t_i[NewTransportedCount][Program.ParameterVectorDimension] = StorageforTransportedClusters.TotalTransportedY_t_i[RemoteIndex][Program.ParameterVectorDimension];
                System.arraycopy(StorageforTransportedClusters.TotalTransportedStatus[RemoteIndex], 0,
                        StorageforTransportedClusters.TotalTransportedStatus[NewTransportedCount], 0, 2);
			}
			++NewTransportedCount;
		}

		// Set new counts
		if (NewTransportedCount < StorageforTransportedClusters.SizeOfTransportedArray)
		{
			changedclusters = true;
			StorageforTransportedClusters.SizeOfTransportedArray = NewTransportedCount;
			Program.ActualMaxTransportedClusterStorage = Math.max(NewTransportedCount + 1, Program.ActualMaxTransportedClusterStorage);
		}
		GlobalChangedClusters = changedclusters;
		DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming5);
        if (DAVectorUtility.MPI_Size > 1){
            // Note - MPI Call - Allreduce - boolean - logical OR
//            GlobalChangedClusters = DAVectorUtility.MPI_communicator.<Boolean>Allreduce(GlobalChangedClusters, Operation<Boolean>.LogicalOr);
            GlobalChangedClusters = DAVectorUtility.mpiOps.allReduce(GlobalChangedClusters, MPI.LOR);
        }
		DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming5);

		if (GlobalChangedClusters)
		{
			ParallelClustering.runningSolution.SetActiveClusters();
			for (int ActiveClusterIndex = 0; ActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; ActiveClusterIndex++)
			{
				int RealClusterIndex = ClusteringSolution.RealClusterIndices[ActiveClusterIndex];
				ClusteringSolution.LocalHost[ActiveClusterIndex] = ClusteringSolution.UniversalMapping[ParallelClustering.runningSolution.LocalCreatedIndex[RealClusterIndex]].PackedHost;
			}
		}

	} // ProcessMovedClusters()

	//  Switch from global to distributed mode
	//  Note later ProcessMovedClusters() will delete clusters not hosted on each node
	//  This routine sets in StorageforTransportedClusters array all clusters NOT hosted here but relevant
	public static void SetupDistributedMode() throws MPIException {
		int RemoteIndex = 0;
		int LocalType2 = 0;
		int SendUp = 0;
		int SendDown = 0;
		int Type01 = 0;
		for (int ActiveClusterIndex = 0; ActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; ActiveClusterIndex++)
		{
			int RealClusterIndex = ClusteringSolution.RealClusterIndices[ActiveClusterIndex];
			if (RealClusterIndex == ParallelClustering.runningSolution.SpongeCluster)
			{
				Type01++;
				continue;
			}
			int PackedHost = ClusteringSolution.LocalHost[ActiveClusterIndex];
			int H = PackedHost & ClusteringSolution.PACKINGMASK;
//C# TO JAVA CONVERTER WARNING: The right shift operator was not replaced by Java's logical right shift operator since the left operand was not confirmed to be of an unsigned type, but you should review whether the logical right shift operator (>>>) is more appropriate:
			int H1 = (PackedHost >> ClusteringSolution.PACKINGSHIFT) & ClusteringSolution.PACKINGMASK;
//C# TO JAVA CONVERTER WARNING: The right shift operator was not replaced by Java's logical right shift operator since the left operand was not confirmed to be of an unsigned type, but you should review whether the logical right shift operator (>>>) is more appropriate:
			int H2 = PackedHost >> (2 * ClusteringSolution.PACKINGSHIFT);
			if (H == DAVectorUtility.MPI_Rank)
			{
				ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] = 2;
				++LocalType2;
				if (H2 > H)
				{
					++SendUp;
				}
				if (H1 < H)
				{
					++SendDown;
				}
				continue;
			}

			if (DAVectorUtility.MPI_Rank < H1)
			{
				continue;
			}
			if (DAVectorUtility.MPI_Rank > H2)
			{
				continue;
			}

			// A Cluster that affects this node but not hosted here
			if (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] < 0)
			{
				continue;
			}
			if (RemoteIndex >= Program.MaxMPITransportBuffer)
			{
				DAVectorUtility.printAndThrowRuntimeException("MPI Buffer over limit " + RemoteIndex + " Limit " + Program.MaxMPITransportBuffer + " Clusters " + ParallelClustering.runningSolution.Ncent_Global + " Local " + ParallelClustering.runningSolution.Ncent_ThisNode);

			}
			StorageforTransportedClusters.TotalTransportedCreatedIndex[RemoteIndex] = ParallelClustering.runningSolution.LocalCreatedIndex[RealClusterIndex];
			StorageforTransportedClusters.TotalTransportedP_t[RemoteIndex] = ParallelClustering.runningSolution.P_k_[RealClusterIndex];
			StorageforTransportedClusters.TotalTransportedOriginalHost[RemoteIndex] = PackedHost;
			for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
			{
				StorageforTransportedClusters.TotalTransportedY_t_i[RemoteIndex][VectorIndex] = ParallelClustering.runningSolution.Y_k_i_[RealClusterIndex][VectorIndex];
				StorageforTransportedClusters.TotalTransportedSigma_t_i[RemoteIndex][VectorIndex] = ParallelClustering.runningSolution.Sigma_k_i_[RealClusterIndex][VectorIndex];
			}
			StorageforTransportedClusters.TotalTransportedY_t_i[RemoteIndex][Program.ParameterVectorDimension] = ParallelClustering.runningSolution.P_k_[RealClusterIndex];

			StorageforTransportedClusters.TotalTransportedStatus[RemoteIndex][0] = 3;
			StorageforTransportedClusters.TotalTransportedStatus[RemoteIndex][1] = ParallelClustering.runningSolution.LocalSplitCreatedIndex[RealClusterIndex];
			++RemoteIndex;
		} // End Loop over Clusters being dispersed

		StorageforTransportedClusters.SizeOfTransportedArray = RemoteIndex;
		Program.ActualMaxTransportedClusterStorage = Math.max(RemoteIndex + 1, Program.ActualMaxTransportedClusterStorage);

		int[] ClusterCountsperNode = new int[DAVectorUtility.MPI_Size];
		int[] UpClusterCountsperNode = new int[DAVectorUtility.MPI_Size];
		int[] DownClusterCountsperNode = new int[DAVectorUtility.MPI_Size];
		DAVectorUtility.StartSubTimer(DAVectorUtility.MPIGATHERTiming);

        if (DAVectorUtility.MPI_Size > 1){
            // Note - MPI Call - Allgather - int
//            ClusterCountsperNode = DAVectorUtility.MPI_communicator.<Integer>Allgather(LocalType2);
            DAVectorUtility.mpiOps.allGather(LocalType2,ClusterCountsperNode);
            // Note - MPI Call - Allgather - int
//            UpClusterCountsperNode = DAVectorUtility.MPI_communicator.<Integer>Allgather(SendUp);
            DAVectorUtility.mpiOps.allGather(SendUp,UpClusterCountsperNode);
            // Note - MPI Call - Allgather - int
//            DownClusterCountsperNode = DAVectorUtility.MPI_communicator.<Integer>Allgather(SendDown);
            DAVectorUtility.mpiOps.allGather(SendDown,DownClusterCountsperNode);
        } else {
            ClusterCountsperNode[0] = LocalType2;
            UpClusterCountsperNode[0] = SendUp;
            DownClusterCountsperNode[0] = SendDown;
        }

		DAVectorUtility.StopSubTimer(DAVectorUtility.MPIGATHERTiming);
		String message = "Global " + Type01;
		int TotalClusterCount = Type01;
		for (int nodeindex = 0; nodeindex < DAVectorUtility.MPI_Size; nodeindex++)
		{
			TotalClusterCount += ClusterCountsperNode[nodeindex];
			message += " Node " + nodeindex + ":" + ClusterCountsperNode[nodeindex] + "(Up:" + UpClusterCountsperNode[nodeindex] + " Down:" + DownClusterCountsperNode[nodeindex] + ")";
		}
		DAVectorUtility.SALSAPrint(1, "Clusters Distributed " + message);
		if (TotalClusterCount != ParallelClustering.runningSolution.Ncent_Global)
		{
			DAVectorUtility.printAndThrowRuntimeException("Inconsistent Cluster Counts " + TotalClusterCount + " " + ParallelClustering.runningSolution.Ncent_Global);

		}

	} //  End SetupDistributedMode()

} // End DistributedClusteringSolution
 // End Namespace edu.indiana.soic.spidal.davs