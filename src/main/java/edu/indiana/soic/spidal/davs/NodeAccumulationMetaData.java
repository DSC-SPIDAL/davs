package edu.indiana.soic.spidal.davs;

// Hybrid Point and Cluster Parallelism
// We exploit a natural one dimensional decomposition in target problem but most ideas in method are general
//  One dimensional decomposition simplifies MPI communication and determination of clusters of interest
//
//  On input assume data (in Y0) is sorted with Processor 0 holding lowest values in sorted decomposed quantity
//  Generate cuts in Y0 based on equal number of points per hread
//  Clusters are assigned based on these point generated cuts. One processor makes these decisions and so ambiguity at cuts is not important
//  A processor can control a cluster stored just outside it quite correctly
//
//  Clusters are labelled by Created Index which is 1 + CreatingProcessor + Creation Number in Processor * MULTIPLIER
//  There is a look up table in every node that maps CreatedIndex to 
//  Global Cluster Numbers
//  Clusters controlled by THIS node and 
//  Clusters controlled by other nodes but transported to this node
//  Nothing for rest of clusters
//  Created Indices are Unique and information set labelled by iteration count to identify stale data
//
//  Set up all key initial parameters and arrays
//  Start in Global mode
//  When clusters reach some limit, switch to distributed mode
//  Do Major Synchronization involving
//      Transmission of Cluster values and Changes
//      Iteration of Point--Cluster assignments
//      Initialize distributed reductions
//
//  Iterate Little loops updating Y T using distributed reductions with "minor synchronization"
//
//  Examine quality of clusters and splitting (involves distributed computation of correlations)
//  Do Major Synchronization if clusters deleted or split or just have changed a lot

//  Data structures used to store ALL (local and remote)Clusters in Node in Reductions
public class NodeAccumulationMetaData
{
	public int NumberofPointsperNode = 0; // Count of Positions in Node Accumulation
	public int[] NumberofPointsperThread; // Count of Positions in Thread Accumulations
	public int[] NodeAccumulationCreatedIndices; // Current list of CreatedIndices for Node Accumulations
	public int[] NodeAccumulationClusterStatus; //  -2 -1 Deleted 0 Sponge 1 Global 2 Local Distribututed 3 From Another Node
	public int[] NodeAccumulationClusterHosts; // Original Hosts in form H + PACKINGMULTIPLIER (H1 + PACKINGMULTIPLIER * H2 )
	public int[][] AccumulationNodetoThreadClusterAssociations; // Associated Thread Locations for Node Accumulations

	public NodeAccumulationMetaData(int MaxAccumulations)
	{
		NumberofPointsperThread = new int[DAVectorUtility.ThreadCount];
		NodeAccumulationCreatedIndices = new int[MaxAccumulations];
		NodeAccumulationClusterStatus = new int[MaxAccumulations];
		NodeAccumulationClusterHosts = new int[MaxAccumulations];
		AccumulationNodetoThreadClusterAssociations = new int[MaxAccumulations][];
		for (int AccumulationIndex = 0; AccumulationIndex < MaxAccumulations; AccumulationIndex++)
		{
			AccumulationNodetoThreadClusterAssociations[AccumulationIndex] = new int[DAVectorUtility.ThreadCount];
		}
	}

} // End NodeAccumulationMetaData