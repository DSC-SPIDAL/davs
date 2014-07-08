package edu.indiana.soic.spidal.davs;

//  Used to store data on Clusters transported from other nodes
public class TransportedClusterStoredInformation
{
	public int SizeOfTransportedArray; // Current Size of Transported Array
	public double[][] TotalTransportedY_t_i; // Center Y_t_i of Transported Cluster (with one extra index to load P_k)
	public double[] TotalTransportedP_t; // Multiplier P_t of Transported Cluster
	public double[][] TotalTransportedSigma_t_i; // Sigma of Transported Cluster (calculated NOT transported)
	public int[] TotalTransportedCreatedIndex; // Created Index of Transported Cluster
	public int[] TotalTransportedOriginalHost; // Host (i.e. node or MPI process) where cluster home in form H + MULTIPLIER (H1 + MULTIPLIER * H2 )
	public int[][] TotalTransportedStatus; // Status of Transported Cluster [0] -2 Moved -1 deleted 0 Sponge 1 Global 2 Distributed Stored here (Not possible) 3 Distributed Stored Elsewhere (default)
											 //  Status[1] 1 + CreatedIndex of Parent or -1 -CreatedIndex of child if Split UNTIL SPLIT PROPAGATED
	public int[] TransportedNodeAccPosition; // Position in Node Accumulation of this remote cluster; -1 NOT used
	public int[][] TransportedThreadAccPosition; // Position in Thread Accumulation of this remote cluster; -1 NOT used

	public TransportedClusterStoredInformation(int MaxStorage)
	{
		SizeOfTransportedArray = 0;
		TotalTransportedOriginalHost = new int[MaxStorage];
		TotalTransportedCreatedIndex = new int[MaxStorage];
		TotalTransportedStatus = new int[MaxStorage][];
		TotalTransportedY_t_i = new double[MaxStorage][];
		TotalTransportedP_t = new double[MaxStorage];
		TotalTransportedSigma_t_i = new double[MaxStorage][];
		TransportedNodeAccPosition = new int[MaxStorage];
		TransportedThreadAccPosition = new int[MaxStorage][];

		for (int TransportIndex = 0; TransportIndex < MaxStorage; TransportIndex++)
		{
			TotalTransportedY_t_i[TransportIndex] = new double[Program.ParameterVectorDimension + 1];
			TotalTransportedSigma_t_i[TransportIndex] = new double[Program.ParameterVectorDimension];
			TotalTransportedStatus[TransportIndex] = new int[2];
			TransportedThreadAccPosition[TransportIndex] = new int[DAVectorUtility.ThreadCount];
		}

	}
} // End TransportedClusterStoredInformation