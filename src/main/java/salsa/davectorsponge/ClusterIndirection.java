package salsa.davectorsponge;

// For fast look up of cluster locations
//  IterationSet allows one to designate which information is really up to date
public class ClusterIndirection
{
	public int IterationSet; // Iteration on which this cluster information set
	public int Availability; // >0 1 + Local Cluster Storage location (Real Cluster Not Active Cluster)  =0 NOT available; < 0 -1 - TransportedPosition
	public int GlobalClusterNumber; // Occasionally set to Global Cluster Number
	public int PackedHost; // =0 unless Distributed when Packed Host index in form H + MULTIPLIER (H1 + MULTIPLIER * H2 ) where H must be current host labelled by ActiveCluster Index
	public double Ymapping; // Value of mapped Y last major iteration

	public ClusterIndirection(int IterationSetINPUT, int AvailabilityINPUT)
	{
		IterationSet = IterationSetINPUT;
		Availability = AvailabilityINPUT;
		GlobalClusterNumber = -1;
		PackedHost = 0;
		Ymapping = 0.0;
	}

} // End ClusterIndirection