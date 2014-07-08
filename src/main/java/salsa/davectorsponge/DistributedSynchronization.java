package salsa.davectorsponge;

import mpi.MPI;
import mpi.MPIException;
import salsa.general.Box;

public class DistributedSynchronization
{ // Do the pipeline distributed broadcast with several different options

	public static MPITransportComponentPacket TransportComponent; // Place to receive transported components

	public static class TransportviaPipeline
	{
        // Note - Logging
        /*static Logger log = LogManager.getLogger(TransportviaPipeline.class.getName());*/

		private int InitialArraySize;
		private int NumberofDoubleComponents;
		private int NumberofIntegerComponents;
		private int UpdateMode; // = 0 Replace in place; used to send round new values of cluster parameters;
		// UpdateMode = 1 Increment in place; used in distributed reduction;
		// UpdateMode = 2 Gather targeted clusters when setting up remote clusters and propagating deleted, moved and split(new) clusters
		private int StorageMode; // = 0 use with UpdateMode=2;
		// StorageMode = 1 Local Cluster Storage (NodeStoragePosition) NOT USED
		// StorageMode = 2 Remote Clusters Stored Locally (TransportedStoragePosition)
		// StorageMode = 3 Node Accumulation Array (NodeAccumulationPosition)
		private boolean HostRangeProcessing; // True use Host Range for destination; False Send only to host

		public TransportviaPipeline(int _UpdateMode, boolean _HostRangeProcessing, int _NumberofDoubleComponents,
                                    int _NumberofIntegerComponents, int _InitialArraySize, int _StorageMode)
		{
			NumberofDoubleComponents = _NumberofDoubleComponents;
			NumberofIntegerComponents = _NumberofIntegerComponents;
			UpdateMode = _UpdateMode;
			InitialArraySize = _InitialArraySize;
			HostRangeProcessing = _HostRangeProcessing;
			StorageMode = _StorageMode;
		}

		//  Note FinalClusterCount, FinalCreatedIndex and FinalHostSpecification ONLY used in UpdateMode 2 and only set in this case
        public final void PipelineDistributedBroadcast(double[][] InitialDoubleComponents, double[][] FinalDoubleComponents,
                                                       int[][] InitialIntegerComponents, int[][] FinalIntegerComponents,
                                                       int[] InitialCreatedIndex, int[] FinalCreatedIndex,
                                                       int[] InitialHostSpecification, int[] FinalHostSpecification,
                                                       Box<Integer> FinalClusterCount) throws MPIException {
            FinalClusterCount.content = 0;
            if (DAVectorUtility.MPI_Size <= 1) {
                return;
            }

            // Now process distributed clusters
            // Variables for processing createdindex
            int NodeStoragePosition = -1;
            int TransportedStoragePosition = -1;
            int NodeAccumulationPosition = -1;
            int ThreadAccumulationPosition = -1;

            //  Place where received data stored
            int FinalDataLocationIndex = -1;

            ++DAVectorSponge.NumberPipelineGroups; // Increment calls of this routine

            int[] DownbySteps = new int[DAVectorUtility.MPI_Size];
            int[] UpbySteps = new int[DAVectorUtility.MPI_Size];
            int[] DownbyStepsTotal = new int[DAVectorUtility.MPI_Size];
            int[] UpbyStepsTotal = new int[2 * DAVectorUtility.MPI_Size];
            for (int PipelineSteps = 0; PipelineSteps < DAVectorUtility.MPI_Size; PipelineSteps++) {
                DownbySteps[PipelineSteps] = 0;
                UpbySteps[PipelineSteps] = 0;
            }

            //  Set NumberUp and NumberDown
            for (int ClusterIndirectIndex = 0; ClusterIndirectIndex < InitialArraySize; ClusterIndirectIndex++) {
                int PackedHost = InitialHostSpecification[ClusterIndirectIndex];
                for (int PipelineSteps = 1; PipelineSteps < DAVectorUtility.MPI_Size; PipelineSteps++) {
                    if (HostRangeProcessing) {
//C# TO JAVA CONVERTER WARNING: The right shift operator was not replaced by Java's logical right shift operator since the left operand was not confirmed to be of an unsigned type, but you should review whether the logical right shift operator (>>>) is more appropriate:
                        int H1 = PackedHost >>> ClusteringSolution.PACKINGSHIFT;
//C# TO JAVA CONVERTER WARNING: The right shift operator was not replaced by Java's logical right shift operator since the left operand was not confirmed to be of an unsigned type, but you should review whether the logical right shift operator (>>>) is more appropriate:
                        int H2 = H1 >>> ClusteringSolution.PACKINGSHIFT;
                        H1 = H1 & ClusteringSolution.PACKINGMASK;
                        if (H2 > (DAVectorUtility.MPI_Rank + PipelineSteps - 1)) {
                            ++UpbySteps[PipelineSteps];
                        }
                        if (H1 < (DAVectorUtility.MPI_Rank - PipelineSteps + 1)) {
                            ++DownbySteps[PipelineSteps];
                        }
                    } else {
                        int H = PackedHost & ClusteringSolution.PACKINGMASK;
                        if (H > (DAVectorUtility.MPI_Rank + PipelineSteps - 1)) {
                            ++UpbySteps[PipelineSteps];
                        }
                        if (H < (DAVectorUtility.MPI_Rank - PipelineSteps + 1)) {
                            ++DownbySteps[PipelineSteps];
                        }
                    }
                }
            }
            for (int PipelineSteps = 0; PipelineSteps < DAVectorUtility.MPI_Size; PipelineSteps++) {
                UpbyStepsTotal[PipelineSteps] = UpbySteps[PipelineSteps];
                UpbyStepsTotal[PipelineSteps + DAVectorUtility.MPI_Size] = DownbySteps[PipelineSteps];
            }
            DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming4);

            if (DAVectorUtility.MPI_Size > 1) {
                // Note - MPI Call - Allreduce - int[] - sum
//                UpbyStepsTotal = DAVectorUtility.MPI_communicator.<Integer>Allreduce(UpbyStepsTotal, Operation<Integer>.Add);
                DAVectorUtility.mpiOps.allReduce(UpbyStepsTotal, MPI.SUM);
            }

            DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming4);
            System.arraycopy(UpbyStepsTotal, DAVectorUtility.MPI_Size, DownbyStepsTotal, 0, DAVectorUtility.MPI_Size);


            // Variables Used for Up and Down Sections
            boolean Initialstep;
            int CurrentNode = DAVectorUtility.MPI_Rank;
            int ReceivedTotal = 0;
            int NumberClustertoSend = 0;
            int NumberDoubletoSend = 0;
            int NumberIntegertoSend = 0;
            int IsItOK;

            // Process Clusters going Up the Chain
            Initialstep = true;
            int LocalTotal = UpbySteps[1];
            int StepsUp = 0;

            while (true) {
                // Decide if ANY node needs to communicate Up
                ++StepsUp;
                if (StepsUp >= DAVectorUtility.MPI_Size){
                    break;
                }

                int JobTotal = UpbyStepsTotal[StepsUp];
                if (JobTotal == 0) {
                    break;
                }


                // Some Nodes want to go up the line
                int SourceProc;
                int DestProc;
                int SourceTag = 0; // Random Number
                int DestTag = 0;
                SourceProc = DAVectorUtility.MPI_Size - 1;
                DestProc = 0;
                if (CurrentNode != 0) {
                    SourceProc = CurrentNode - 1;
                }
                if (CurrentNode != (DAVectorUtility.MPI_Size - 1)) {
                    DestProc = CurrentNode + 1;
                } else {
                    LocalTotal = 0;
                }
                MPITransportComponentPacket SendBuffer = new MPITransportComponentPacket(LocalTotal, NumberofDoubleComponents, NumberofIntegerComponents); // Sent Buffer is EXACT size
                NumberClustertoSend = 0;
                NumberDoubletoSend = 0;
                NumberIntegertoSend = 0;

                if (LocalTotal > 0) { // If no data here, just send dummy packet
                    if (Initialstep) { // Construct message to send from Initial Arrays
                        for (int ClusterSendPointer = 0; ClusterSendPointer < InitialArraySize; ClusterSendPointer++) {
                            int PackedHost = InitialHostSpecification[ClusterSendPointer];
                            if (HostRangeProcessing) {
//C# TO JAVA CONVERTER WARNING: The right shift operator was not replaced by Java's logical right shift operator since the left operand was not confirmed to be of an unsigned type, but you should review whether the logical right shift operator (>>>) is more appropriate:
                                int H2 = PackedHost >>> (2 * ClusteringSolution.PACKINGSHIFT);
                                if (H2 <= DAVectorUtility.MPI_Rank) {
                                    continue;
                                }
                            } else {
                                int H = PackedHost & ClusteringSolution.PACKINGMASK;
                                if (H <= DAVectorUtility.MPI_Rank) {
                                    continue;
                                }
                            }
                            SendBuffer.setAssociatedCreatedIndexAt(NumberClustertoSend, InitialCreatedIndex[ClusterSendPointer]);
                            SendBuffer.setClusterHostRangeAt(NumberClustertoSend, InitialHostSpecification[ClusterSendPointer]);
                            if (NumberofDoubleComponents > 0) {
                                for (int ComponentIndex = 0; ComponentIndex < NumberofDoubleComponents; ComponentIndex++) {
                                    SendBuffer.setClusterDoubleComponentAt(NumberDoubletoSend, InitialDoubleComponents[ClusterSendPointer][ComponentIndex]);
                                    ++NumberDoubletoSend;
                                }
                            }
                            if (NumberofIntegerComponents > 0) {
                                for (int ComponentIndex = 0; ComponentIndex < NumberofIntegerComponents; ComponentIndex++) {
                                    SendBuffer.setClusterIntegerComponentAt(NumberIntegertoSend, InitialIntegerComponents[ClusterSendPointer][ComponentIndex]);
                                    ++NumberIntegertoSend;
                                }
                            }
                            ++NumberClustertoSend;
                            if (NumberClustertoSend >= LocalTotal) {
                                break;
                            }
                        }
                    } else { // Construct message to send from en passant data
                        for (int ReceivedClusterIndex = 0; ReceivedClusterIndex < ReceivedTotal; ReceivedClusterIndex++) {
                            int PackedHost = TransportComponent.getClusterHostRangeAt(ReceivedClusterIndex);
                            if (HostRangeProcessing) {
//C# TO JAVA CONVERTER WARNING: The right shift operator was not replaced by Java's logical right shift operator since the left operand was not confirmed to be of an unsigned type, but you should review whether the logical right shift operator (>>>) is more appropriate:
                                int H2 = PackedHost >>> (2 * ClusteringSolution.PACKINGSHIFT);
                                if (H2 <= DAVectorUtility.MPI_Rank) {
                                    continue;
                                }
                            } else {
                                int H = PackedHost & ClusteringSolution.PACKINGMASK;
                                if (H <= DAVectorUtility.MPI_Rank) {
                                    continue;
                                }
                            }
                            SendBuffer.setAssociatedCreatedIndexAt(NumberClustertoSend, TransportComponent.getAssociatedCreatedIndexAt(ReceivedClusterIndex));
                            SendBuffer.setClusterHostRangeAt(NumberClustertoSend, TransportComponent.getClusterHostRangeAt(ReceivedClusterIndex));
                            if (NumberofDoubleComponents > 0) {
                                int OverallIndex = NumberofDoubleComponents * ReceivedClusterIndex;
                                for (int ComponentIndex = 0; ComponentIndex < NumberofDoubleComponents; ComponentIndex++) {
                                    SendBuffer.setClusterDoubleComponentAt(NumberDoubletoSend, TransportComponent.getClusterDoubleComponentAt(OverallIndex));
                                    ++NumberDoubletoSend;
                                    ++OverallIndex;
                                }
                            }
                            if (NumberofIntegerComponents > 0) {
                                int OverallIndex = NumberofIntegerComponents * ReceivedClusterIndex;
                                for (int ComponentIndex = 0; ComponentIndex < NumberofIntegerComponents; ComponentIndex++) {
                                    SendBuffer.setClusterIntegerComponentAt(NumberIntegertoSend, TransportComponent.getClusterIntegerComponentAt(OverallIndex));
                                    ++NumberIntegertoSend;
                                    ++OverallIndex;
                                }
                            }
                            ++NumberClustertoSend;
                            if (NumberClustertoSend >= LocalTotal) {
                                break;
                            }
                        }

                    }
                } // End Case where there is Local Data to Send

                // Send data in a pipeline forward
                DAVectorUtility.StartSubTimer(DAVectorUtility.MPISENDRECEIVETiming);

                if (DAVectorUtility.MPI_Size > 1) {
                    // Note - MPI Call - Sendreceive - MPITransportComponentPacket
                    TransportComponent = DAVectorUtility.mpiOps.sendReceive(SendBuffer, DestProc, DestTag, SourceProc, SourceTag);
                } else {
                    TransportComponent = SendBuffer;
                }

                DAVectorUtility.StopSubTimer(DAVectorUtility.MPISENDRECEIVETiming);

                ++DAVectorSponge.NumberPipelineSteps;
                DAVectorSponge.NumberofPipelineClusters += SendBuffer.getNumberOfClusters();

                //  Examine Data passed from lower ranked processor
                //  Set new LocalTotal and Store Data
                ReceivedTotal = TransportComponent.getNumberOfClusters();
                DAVectorSponge.ActualMaxMPITransportBuffer = Math.max(DAVectorSponge.ActualMaxMPITransportBuffer, ReceivedTotal);
                LocalTotal = 0; // Count Number of Clusters on next step

                if (NumberofDoubleComponents != TransportComponent.getNumberOfDoubleComponents()) {
                    DAVectorUtility.printAndThrowRuntimeException(" Double Components Inconsistent " + NumberofDoubleComponents + " " + TransportComponent.getNumberOfDoubleComponents() + " in Rank " + DAVectorUtility.MPI_Rank + " Bad");

                }
                if (NumberofIntegerComponents != TransportComponent.getNumberOfIntegerComponents()) {
                    DAVectorUtility.printAndThrowRuntimeException(" Integer Components Inconsistent " + NumberofIntegerComponents + " " + TransportComponent.getNumberOfIntegerComponents() + " in Rank " + DAVectorUtility.MPI_Rank + " Bad");

                }

                if (ReceivedTotal > 0) {
                    for (int ReceivedClusterIndex = 0; ReceivedClusterIndex < ReceivedTotal; ReceivedClusterIndex++) {
                        int PackedHost = TransportComponent.getClusterHostRangeAt(ReceivedClusterIndex);
                        if (HostRangeProcessing) {
//C# TO JAVA CONVERTER WARNING: The right shift operator was not replaced by Java's logical right shift operator since the left operand was not confirmed to be of an unsigned type, but you should review whether the logical right shift operator (>>>) is more appropriate:
                            int H2 = PackedHost >>> (2 * ClusteringSolution.PACKINGSHIFT);
                            if (H2 < DAVectorUtility.MPI_Rank) {
                                DAVectorUtility.printAndThrowRuntimeException(" Transported host " + PackedHost + " in Rank " + DAVectorUtility.MPI_Rank + " Bad Up Range");

                            }
                            if (H2 > DAVectorUtility.MPI_Rank) {
                                ++LocalTotal;
                            }
                        } else {
                            int H = PackedHost & ClusteringSolution.PACKINGMASK;
                            if (H < DAVectorUtility.MPI_Rank) {
                                DAVectorUtility.printAndThrowRuntimeException(" Transported host " + PackedHost + " in Rank " + DAVectorUtility.MPI_Rank + " Bad Up Host");

                            }
                            if (H > DAVectorUtility.MPI_Rank) {
                                ++LocalTotal;
                                continue;
                            }
                        }
                        int host = PackedHost & ClusteringSolution.PACKINGMASK;
                        int CreatedIndex = TransportComponent.getAssociatedCreatedIndexAt(ReceivedClusterIndex);
                        if (UpdateMode < 2) {
                            FinalDataLocationIndex = -1;
                            Box<Integer> tempRef_NodeStoragePosition = new Box<>(NodeStoragePosition);
                            Box<Integer> tempRef_TransportedStoragePosition = new Box<>(TransportedStoragePosition);
                            Box<Integer> tempRef_NodeAccumulationPosition = new Box<>(NodeAccumulationPosition);
                            Box<Integer> tempRef_ThreadAccumulationPosition = new Box<>(ThreadAccumulationPosition);
                            IsItOK = DistributedClusteringSolution.IndicesperCluster(CreatedIndex, -1, tempRef_NodeStoragePosition, tempRef_TransportedStoragePosition, tempRef_NodeAccumulationPosition, tempRef_ThreadAccumulationPosition);
                            NodeStoragePosition = tempRef_NodeStoragePosition.content;
                            TransportedStoragePosition = tempRef_TransportedStoragePosition.content;
                            NodeAccumulationPosition = tempRef_NodeAccumulationPosition.content;
                            ThreadAccumulationPosition = tempRef_ThreadAccumulationPosition.content;
                            if (StorageMode == 1) {
                                FinalDataLocationIndex = NodeStoragePosition;
                            }
                            if (StorageMode == 2) {
                                FinalDataLocationIndex = TransportedStoragePosition;
                            }
                            if (StorageMode == 3) {
                                FinalDataLocationIndex = NodeAccumulationPosition;
                            }
                            if ((host == DAVectorUtility.MPI_Rank) && (IsItOK != 0)) {
                                DAVectorUtility.printAndThrowRuntimeException(" Transported Created Index " + CreatedIndex + " in Rank " + DAVectorUtility.MPI_Rank + " Bad with code " + IsItOK + " host " + host + " Update mode " + UpdateMode);

                            }
                        } else { // UpdateMode 2
                            IsItOK = 0;
                            FinalDataLocationIndex = FinalClusterCount.content;
                            ++FinalClusterCount.content;
                            FinalCreatedIndex[FinalDataLocationIndex] = CreatedIndex;
                            FinalHostSpecification[FinalDataLocationIndex] = PackedHost;
                        }
                        if (IsItOK >= 0) {
                            if (FinalDataLocationIndex == -1) {
                                DAVectorUtility.printAndThrowRuntimeException(" Transported Created Index " + CreatedIndex + " in Rank " + DAVectorUtility.MPI_Rank + " Bad with Storage Mode " + StorageMode + " host " + host + " Update mode " + UpdateMode);

                            }

                            if (NumberofDoubleComponents > 0) {
                                String message = "";
                                int OverallDoubleIndex = NumberofDoubleComponents * ReceivedClusterIndex;
                                for (int ComponentIndex = 0; ComponentIndex < NumberofDoubleComponents; ComponentIndex++) {
                                    if ((UpdateMode == 0) || (UpdateMode == 2)) {
                                        FinalDoubleComponents[FinalDataLocationIndex][ComponentIndex] = TransportComponent.getClusterDoubleComponentAt(OverallDoubleIndex);
                                    }
                                    if (UpdateMode == 1) {
                                        FinalDoubleComponents[FinalDataLocationIndex][ComponentIndex] += TransportComponent.getClusterDoubleComponentAt(OverallDoubleIndex);
                                    }
                                    message += " * " + String.format("%1$4.3E", TransportComponent.getClusterDoubleComponentAt(OverallDoubleIndex)) + " " + String.format("%1$4.3f", FinalDoubleComponents[FinalDataLocationIndex][ComponentIndex]);
                                    ++OverallDoubleIndex;
                                }
                                /* if (CreatedIndex == 901)
								    DAVectorUtility.SALSAFullPrint(1, "Up901 Transport " + UpdateMode.ToString() + " " + FinalDataLocationIndex.ToString()  + message); */
                            }

                            if (NumberofIntegerComponents > 0) {
                                int OverallIntegerIndex = NumberofIntegerComponents * ReceivedClusterIndex;
                                for (int ComponentIndex = 0; ComponentIndex < NumberofIntegerComponents; ComponentIndex++) {
                                    if ((UpdateMode == 0) || (UpdateMode == 2)) {
                                        FinalIntegerComponents[FinalDataLocationIndex][ComponentIndex] = TransportComponent.getClusterIntegerComponentAt(OverallIntegerIndex);
                                    }
                                    if (UpdateMode == 1) {
                                        FinalIntegerComponents[FinalDataLocationIndex][ComponentIndex] += TransportComponent.getClusterIntegerComponentAt(OverallIntegerIndex);
                                    }
                                    ++OverallIntegerIndex;
                                }
                            }
                        } // End case where location found IsItOK >= 0
                    } // end Loop over ReceivedClusterIndex
                } // End case where ReceivedTotal > 0
                Initialstep = false;

            } // End While over MPI pipeline steps for pipeline going UP the chain

            // Process Clusters going Down the Chain
            Initialstep = true;
            int StepsDown = 0;
            LocalTotal = DownbySteps[1];

            while (true) {
                StepsDown++;
                if (StepsDown >= DAVectorUtility.MPI_Size){
                    break;
                }

                int JobTotal = DownbyStepsTotal[StepsDown];
                if (JobTotal == 0) {
                    break;
                }


                // Some Nodes want to go up the line
                int SourceProc;
                int DestProc;
                DestProc = DAVectorUtility.MPI_Size - 1;
                SourceProc = 0;
                int SourceTag = 22; // Random Number
                int DestTag = 22;
                if (CurrentNode != 0) {
                    DestProc = CurrentNode - 1;
                } else {
                    LocalTotal = 0;
                }
                if (CurrentNode != (DAVectorUtility.MPI_Size - 1)) {
                    SourceProc = CurrentNode + 1;
                }
                MPITransportComponentPacket SendBuffer = new MPITransportComponentPacket(LocalTotal, NumberofDoubleComponents, NumberofIntegerComponents); // Sent Buffer is EXACT size
                NumberClustertoSend = 0;
                NumberDoubletoSend = 0;
                NumberIntegertoSend = 0;

                if (LocalTotal > 0) { // If no data here, just send dummy packet

                    if (Initialstep) { // Construct message to send from local accumulation arrays
                        for (int ClusterIndirectIndex = 0; ClusterIndirectIndex < InitialArraySize; ClusterIndirectIndex++) {
                            int PackedHost = InitialHostSpecification[ClusterIndirectIndex];
                            if (HostRangeProcessing) {
//C# TO JAVA CONVERTER WARNING: The right shift operator was not replaced by Java's logical right shift operator since the left operand was not confirmed to be of an unsigned type, but you should review whether the logical right shift operator (>>>) is more appropriate:
                                int H1 = (PackedHost >>> ClusteringSolution.PACKINGSHIFT) & ClusteringSolution.PACKINGMASK;
                                if (H1 >= DAVectorUtility.MPI_Rank) {
                                    continue;
                                }
                            } else {
                                int H = PackedHost & ClusteringSolution.PACKINGMASK;
                                if (H >= DAVectorUtility.MPI_Rank) {
                                    continue;
                                }
                            }
                            SendBuffer.setAssociatedCreatedIndexAt(NumberClustertoSend, InitialCreatedIndex[ClusterIndirectIndex]);
                            SendBuffer.setClusterHostRangeAt(NumberClustertoSend, InitialHostSpecification[ClusterIndirectIndex]);
                            if (NumberofDoubleComponents > 0) {
                                for (int ComponentIndex = 0; ComponentIndex < NumberofDoubleComponents; ComponentIndex++) {
                                    SendBuffer.setClusterDoubleComponentAt(NumberDoubletoSend, InitialDoubleComponents[ClusterIndirectIndex][ComponentIndex]);
                                    ++NumberDoubletoSend;
                                }
                            }
                            if (NumberofIntegerComponents > 0) {
                                for (int ComponentIndex = 0; ComponentIndex < NumberofIntegerComponents; ComponentIndex++) {
                                    SendBuffer.setClusterIntegerComponentAt(NumberIntegertoSend, InitialIntegerComponents[ClusterIndirectIndex][ComponentIndex]);
                                    ++NumberIntegertoSend;
                                }
                            }
                            ++NumberClustertoSend;
                            if (NumberClustertoSend >= LocalTotal) {
                                break;
                            }
                        }
                    } else { // Construct message to send from en passant data
                        for (int ReceivedClusterIndex = 0; ReceivedClusterIndex < ReceivedTotal; ReceivedClusterIndex++) {
                            int PackedHost = TransportComponent.getClusterHostRangeAt(ReceivedClusterIndex);
                            if (HostRangeProcessing) {
//C# TO JAVA CONVERTER WARNING: The right shift operator was not replaced by Java's logical right shift operator since the left operand was not confirmed to be of an unsigned type, but you should review whether the logical right shift operator (>>>) is more appropriate:
                                int H1 = (PackedHost >>> ClusteringSolution.PACKINGSHIFT) & ClusteringSolution.PACKINGMASK;
                                if (H1 >= DAVectorUtility.MPI_Rank) {
                                    continue;
                                }
                            } else {
                                int H = PackedHost & ClusteringSolution.PACKINGMASK;
                                if (H >= DAVectorUtility.MPI_Rank) {
                                    continue;
                                }
                            }
                            SendBuffer.setAssociatedCreatedIndexAt(NumberClustertoSend, TransportComponent.getAssociatedCreatedIndexAt(ReceivedClusterIndex));
                            SendBuffer.setClusterHostRangeAt(NumberClustertoSend, TransportComponent.getClusterHostRangeAt(ReceivedClusterIndex));
                            if (NumberofDoubleComponents > 0) {
                                int OverallIndex = NumberofDoubleComponents * ReceivedClusterIndex;
                                for (int ComponentIndex = 0; ComponentIndex < NumberofDoubleComponents; ComponentIndex++) {
                                    SendBuffer.setClusterDoubleComponentAt(NumberDoubletoSend, TransportComponent.getClusterDoubleComponentAt(OverallIndex));
                                    ++NumberDoubletoSend;
                                    ++OverallIndex;
                                }
                            }
                            if (NumberofIntegerComponents > 0) {
                                int OverallIndex = NumberofIntegerComponents * ReceivedClusterIndex;
                                for (int ComponentIndex = 0; ComponentIndex < NumberofIntegerComponents; ComponentIndex++) {
                                    SendBuffer.setClusterIntegerComponentAt(NumberIntegertoSend, TransportComponent.getClusterIntegerComponentAt(OverallIndex));
                                    ++NumberIntegertoSend;
                                    ++OverallIndex;
                                }
                            }
                            ++NumberClustertoSend;
                            if (NumberClustertoSend >= LocalTotal) {
                                break;
                            }
                        }
                    } // end en passant data is source of information

                } // End Case where there is Local Data to Send

                // Send data in a pipeline backwards
                DAVectorUtility.StartSubTimer(DAVectorUtility.MPISENDRECEIVETiming);

                if (DAVectorUtility.MPI_Size > 1) {
                    // Note - MPI Call - Sendreceive - MPITransportComponentPacket
                    TransportComponent = DAVectorUtility.mpiOps.sendReceive(SendBuffer, DestProc, DestTag, SourceProc, SourceTag);
                } else {
                    TransportComponent = SendBuffer;
                }
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPISENDRECEIVETiming);
                ++DAVectorSponge.NumberPipelineSteps;
                DAVectorSponge.NumberofPipelineClusters += SendBuffer.getNumberOfClusters();

                //  Examine Data passed from higher ranked processor
                ReceivedTotal = TransportComponent.getNumberOfClusters();
                DAVectorSponge.ActualMaxMPITransportBuffer = Math.max(DAVectorSponge.ActualMaxMPITransportBuffer, ReceivedTotal);
                LocalTotal = 0;

                if (NumberofDoubleComponents != TransportComponent.getNumberOfDoubleComponents()) {
                    DAVectorUtility.printAndThrowRuntimeException(" Double Components Inconsistent " + NumberofDoubleComponents + " " + TransportComponent.getNumberOfDoubleComponents() + " in Rank " + DAVectorUtility.MPI_Rank + " Bad");

                }
                if (NumberofIntegerComponents != TransportComponent.getNumberOfIntegerComponents()) {
                    DAVectorUtility.printAndThrowRuntimeException(" Integer Components Inconsistent " + NumberofIntegerComponents + " " + TransportComponent.getNumberOfIntegerComponents() + " in Rank " + DAVectorUtility.MPI_Rank + " Bad");

                }

                if (ReceivedTotal > 0) {
                    for (int ReceivedClusterIndex = 0; ReceivedClusterIndex < ReceivedTotal; ReceivedClusterIndex++) {
                        int PackedHost = TransportComponent.getClusterHostRangeAt(ReceivedClusterIndex);
                        if (HostRangeProcessing) {
//C# TO JAVA CONVERTER WARNING: The right shift operator was not replaced by Java's logical right shift operator since the left operand was not confirmed to be of an unsigned type, but you should review whether the logical right shift operator (>>>) is more appropriate:
                            int H1 = (PackedHost >>> ClusteringSolution.PACKINGSHIFT) & ClusteringSolution.PACKINGMASK;
                            if (H1 > DAVectorUtility.MPI_Rank) {
                                DAVectorUtility.printAndThrowRuntimeException(" Transported host " + PackedHost + " in Rank " + DAVectorUtility.MPI_Rank + " Bad Down Range");

                            }
                            if (H1 < DAVectorUtility.MPI_Rank) {
                                ++LocalTotal;
                            }
                        } else {
                            int H = PackedHost & ClusteringSolution.PACKINGMASK;
                            if (H > DAVectorUtility.MPI_Rank) {
                                DAVectorUtility.printAndThrowRuntimeException(" Transported host " + PackedHost + " in Rank " + DAVectorUtility.MPI_Rank + " Bad Down Not Range");

                            }
                            if (H < DAVectorUtility.MPI_Rank) {
                                ++LocalTotal;
                                continue;
                            }
                        }
                        int host = PackedHost & ClusteringSolution.PACKINGMASK;
                        int CreatedIndex = TransportComponent.getAssociatedCreatedIndexAt(ReceivedClusterIndex);
                        if (UpdateMode < 2) {
                            FinalDataLocationIndex = -1;
                            Box<Integer> tempRef_NodeStoragePosition2 = new Box<>(NodeStoragePosition);
                            Box<Integer> tempRef_TransportedStoragePosition2 = new Box<>(TransportedStoragePosition);
                            Box<Integer> tempRef_NodeAccumulationPosition2 = new Box<>(NodeAccumulationPosition);
                            Box<Integer> tempRef_ThreadAccumulationPosition2 = new Box<>(ThreadAccumulationPosition);
                            IsItOK = DistributedClusteringSolution.IndicesperCluster(CreatedIndex, -1, tempRef_NodeStoragePosition2, tempRef_TransportedStoragePosition2, tempRef_NodeAccumulationPosition2, tempRef_ThreadAccumulationPosition2);
                            NodeStoragePosition = tempRef_NodeStoragePosition2.content;
                            TransportedStoragePosition = tempRef_TransportedStoragePosition2.content;
                            NodeAccumulationPosition = tempRef_NodeAccumulationPosition2.content;
                            ThreadAccumulationPosition = tempRef_ThreadAccumulationPosition2.content;
                            if (StorageMode == 1) {
                                FinalDataLocationIndex = NodeStoragePosition;
                            }
                            if (StorageMode == 2) {
                                FinalDataLocationIndex = TransportedStoragePosition;
                            }
                            if (StorageMode == 3) {
                                FinalDataLocationIndex = NodeAccumulationPosition;
                            }
                            if ((host == DAVectorUtility.MPI_Rank) && (IsItOK != 0)) {
                                DAVectorUtility.printAndThrowRuntimeException(" Transported Created Index " + CreatedIndex + " in Rank " + DAVectorUtility.MPI_Rank + " Bad with code " + IsItOK + " host " + host + " Update mode " + UpdateMode);

                            }
                        } else {
                            IsItOK = 0;
                            FinalDataLocationIndex = FinalClusterCount.content;
                            ++FinalClusterCount.content;
                            FinalCreatedIndex[FinalDataLocationIndex] = CreatedIndex;
                            FinalHostSpecification[FinalDataLocationIndex] = PackedHost;
                        }
                        if (IsItOK >= 0) {
                            if (FinalDataLocationIndex == -1) {
                                DAVectorUtility.printAndThrowRuntimeException(" Transported Created Index " + CreatedIndex + " in Rank " + DAVectorUtility.MPI_Rank + " Bad with Storage Mode " + StorageMode + " host " + host + " Update mode " + UpdateMode);

                            }

                            if (NumberofDoubleComponents > 0) {
                                String message = "";
                                int OverallDoubleIndex = NumberofDoubleComponents * ReceivedClusterIndex;
                                for (int ComponentIndex = 0; ComponentIndex < NumberofDoubleComponents; ComponentIndex++) {
                                    if ((UpdateMode == 0) || (UpdateMode == 2)) {
                                        FinalDoubleComponents[FinalDataLocationIndex][ComponentIndex] = TransportComponent.getClusterDoubleComponentAt(OverallDoubleIndex);
                                    }
                                    if (UpdateMode == 1) {
                                        FinalDoubleComponents[FinalDataLocationIndex][ComponentIndex] += TransportComponent.getClusterDoubleComponentAt(OverallDoubleIndex);
                                    }
                                    message += " * " + String.format("%1$4.3E", TransportComponent.getClusterDoubleComponentAt(OverallDoubleIndex)) + " " + String.format("%1$4.3f", FinalDoubleComponents[FinalDataLocationIndex][ComponentIndex]);
                                    ++OverallDoubleIndex;
                                }
								/* if (CreatedIndex == 901)
								    DAVectorUtility.SALSAFullPrint(1, "Dn901 Transport " + UpdateMode.ToString() + " " + FinalDataLocationIndex.ToString() + message); */
                            }

                            if (NumberofIntegerComponents > 0) {
                                int OverallIntegerIndex = NumberofIntegerComponents * ReceivedClusterIndex;
                                for (int ComponentIndex = 0; ComponentIndex < NumberofIntegerComponents; ComponentIndex++) {
                                    if ((UpdateMode == 0) || (UpdateMode == 2)) {
                                        FinalIntegerComponents[FinalDataLocationIndex][ComponentIndex] = TransportComponent.getClusterIntegerComponentAt(OverallIntegerIndex);
                                    }
                                    if (UpdateMode == 1) {
                                        FinalIntegerComponents[FinalDataLocationIndex][ComponentIndex] += TransportComponent.getClusterIntegerComponentAt(OverallIntegerIndex);
                                    }
                                    ++OverallIntegerIndex;
                                }
                            }
                        }
                    } // End Loop over Received Clusters
                } // End case when data received ReceivedTotal > 0

                Initialstep = false;

            } // End While over MPI pipeline steps going DOWN the chain

        } // End PipelineDistributedBroadcast

	} // End TransportviaPipeline

} // End class DistributedSynchronization
 // End Namespace SALSALIbrary
