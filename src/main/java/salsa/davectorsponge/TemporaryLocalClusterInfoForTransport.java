package salsa.davectorsponge;

/**
 * Used to temporarily store data on Distributed Clusters hosted on this node for transportation to other nodes
 */
public class TemporaryLocalClusterInfoForTransport {
    public int sizeOfTransportedArray; // Current Size of Transported Array
    public double[][] totalTransportedY_t_i; // Center Y_t_i of Transported Cluster with last index P[k]so dimension
    public int[] totalTransportedCreatedIndex; // Created Index of Transported Cluster
    // Host (i.e. node or MPI process) where cluster home in form H + MULTIPLIER (H1 + MULTIPLIER * H2 )
    public int[] totalTransportedOriginalHost;
    // Status of Transported Cluster [0] -2 Moved -1 deleted 0 Sponge 1 Global 2 Distributed Stored here (Not
    // possible) 3 Distributed Stored Elsewhere (default)
    // Status[1] 1 + CreatedIndex of Parent or -1 -CreatedIndex of child if Split UNTIL SPLIT PROPAGATED
    public int[][] totalTransportedStatus;

    public TemporaryLocalClusterInfoForTransport(int maxStorage) {
        sizeOfTransportedArray = 0;
        totalTransportedOriginalHost = new int[maxStorage];
        totalTransportedCreatedIndex = new int[maxStorage];
        totalTransportedStatus = new int[maxStorage][];
        totalTransportedY_t_i = new double[maxStorage][];
        for (int TransportIndex = 0; TransportIndex < maxStorage; TransportIndex++) {
            totalTransportedY_t_i[TransportIndex] = new double[DAVectorSponge.ParameterVectorDimension + 1];
            totalTransportedStatus[TransportIndex] = new int[2];
        }
    }

}