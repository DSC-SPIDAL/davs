package edu.indiana.soic.spidal.davs;

//	Class VectorClass **************************************************

import edu.rice.hj.api.SuspendableException;
import mpi.MPIException;
import edu.indiana.soic.spidal.general.Box;

import java.util.Random;

import static edu.rice.hj.Module1.forallChunked;

public class VectorClass
{

	public double Eigenvalue; // Current eigenvalue
	public int EigenStatus; // Indicator of status of eigenvalue
	public double[] Eigenvector; // Eigenvector

	public double[][] CenterEigenvector;
	public double[] CenterEigenvalue;
	public double[] InitVector;
	public double[] FirstTerm;
	public int[] CenterEigenstatus;
	public int[] CenterEigenconvergence;

	public ClusteringSolution CurrentSolution;

	// Find the minimum eigenvalue of second derivative matrix -- called from shouldSplit
	public final void getEigenvaluefromMatrix(double[][] SecondDerivMatrix)
	{

		Eigenvector = new double[Program.ParameterVectorDimension];
		if (!Program.CalculateEigenvaluesfromMatrix)
		{ // Re-use Earlier Calculation of results for all clusters
		}

		//  Calculate Eigenvalues from Matrix
		if (Program.ParameterVectorDimension != 2)
		{
			DAVectorUtility.printAndThrowRuntimeException(" Illegal Vector Dimension " + Program.ParameterVectorDimension);

		}

		//  Case of Two Dimensions
		EigenStatus = 1;
		double tmp = SecondDerivMatrix[0][0] - SecondDerivMatrix[1][1];
		tmp = tmp * tmp + 4.0 * SecondDerivMatrix[0][1] * SecondDerivMatrix[0][1];
		Eigenvalue = 0.5 * (SecondDerivMatrix[0][0] + SecondDerivMatrix[1][1] - Math.sqrt(tmp));
		Eigenvector[0] = -SecondDerivMatrix[0][1];
		Eigenvector[1] = SecondDerivMatrix[1][1] - Eigenvalue;

		// Normalize
		tmp = 0.0;
		for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
		{
			tmp += Eigenvector[VectorIndex] * Eigenvector[VectorIndex];
		}
		tmp = 1.0 / Math.sqrt(tmp);
		for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
		{
			Eigenvector[VectorIndex] *= tmp;
		}

	} // End getEigenvalue(double[,] SecondDerivMatrix)

	public final void SetAllEigenvaluesIteratively(ClusteringSolution Solution) throws MPIException {
		if (Solution.DistributedExecutionMode)
		{
			DAVectorUtility.printAndThrowRuntimeException(" Illegal Eigenvalue and Parallelization Combination ");

		}
		if (Program.SigmaMethod > 0)
		{
			DAVectorUtility.printAndThrowRuntimeException(" Illegal Eigenvalue and Sigma Method Combination " + Program.SigmaMethod);

		}
		this.CurrentSolution = Solution;
		this.CenterEigenvector = this.CurrentSolution.Eigenvector_k_i;
		this.CenterEigenvalue = this.CurrentSolution.Eigenvalue_k;
		this.InitVector = new double[Program.ParameterVectorDimension];
		this.FirstTerm = new double[this.CurrentSolution.Ncent_Global];
		this.CenterEigenstatus = new int[this.CurrentSolution.Ncent_Global];
		this.CenterEigenconvergence = new int[this.CurrentSolution.Ncent_Global];

		Random random = new Random();
		double InitNorm = 0.0;
		for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
		{
			InitVector[VectorIndex] = -0.5 + random.nextDouble();
			InitNorm += InitVector[VectorIndex] * InitVector[VectorIndex];
		}
		InitNorm = 1.0 / Math.sqrt(InitNorm);
		for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
		{
			InitVector[VectorIndex] *= InitNorm;
		}

		//  Initialization Loop over Clusters
		int somethingtodo = 0;
		for (int ClusterIndex = 0; ClusterIndex < this.CurrentSolution.Ncent_Global; ClusterIndex++)
		{
			this.CenterEigenconvergence[ClusterIndex] = 0;
			this.CenterEigenstatus[ClusterIndex] = 0;
			this.FirstTerm[ClusterIndex] = 0;
			if (this.CurrentSolution.Splittable_k_[ClusterIndex] != 1)
			{
				continue;
			}
			++somethingtodo;
            System.arraycopy(InitVector, 0, this.CenterEigenvector[ClusterIndex], 0,
                    Program.ParameterVectorDimension);
		} // End Loop over Clusters
		if (somethingtodo == 0)
		{
			return;
		}

		final GlobalReductions.FindVectorDoubleSum FindClusterFirstTerm = new GlobalReductions.FindVectorDoubleSum(DAVectorUtility.ThreadCount, this.CurrentSolution.Ncent_Global);
		final GlobalReductions.FindDoubleSum FindNumberScalarProducts = new GlobalReductions.FindDoubleSum(DAVectorUtility.ThreadCount);

		for (int NumPowerIterations = 0; NumPowerIterations < Program.PowerIterationLimit; NumPowerIterations++)
		{
			somethingtodo = 0;
			for (int ClusterIndex = 0; ClusterIndex < this.CurrentSolution.Ncent_Global; ClusterIndex++)
			{
				if (this.CurrentSolution.LocalStatus[ClusterIndex] != 1)
				{
					continue;
				}
				if (this.CurrentSolution.Splittable_k_[ClusterIndex] != 1)
				{
					continue;
				}
				if (this.CenterEigenconvergence[ClusterIndex] == 0)
				{
					++somethingtodo;
				}
			}
			if (somethingtodo == 0)
			{
				break;
			}

			final GlobalReductions.FindVectorDoubleSum3 FindNewPowerVectors = new GlobalReductions.FindVectorDoubleSum3(DAVectorUtility.ThreadCount, Program.ParameterVectorDimension, this.CurrentSolution.Ncent_Global);

            final int NumPowerIterationsLoopVar = NumPowerIterations;
            // Note - parallel for
            try {
                forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) -> {
                    FindNewPowerVectors.startthread(threadIndex);
                    double[] PartVector = new double[Program.ParameterVectorDimension];
                    int indexlen = DAVectorUtility.PointsperThread[threadIndex];
                    int beginpoint = DAVectorUtility.StartPointperThread[threadIndex] - DAVectorUtility.PointStart_Process;
                    for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++) {
                        int IndirectSize = this.CurrentSolution.NumClusters_alpha_[alpha];
                        for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++) {
                            int RealClusterIndex = -1;
                            int RemoteIndex = -1;
                            int ActiveClusterIndex = -1;
                            Box<Integer> tempRef_RealClusterIndex = new Box<>(RealClusterIndex);
                            Box<Integer> tempRef_ActiveClusterIndex = new Box<>(ActiveClusterIndex);
                            Box<Integer> tempRef_RemoteIndex = new Box<>(RemoteIndex);
                            VectorAnnealIterate
                                    .ClusterPointersforaPoint(alpha, IndirectClusterIndex, tempRef_RealClusterIndex,
                                            tempRef_ActiveClusterIndex, tempRef_RemoteIndex);
                            RealClusterIndex = tempRef_RealClusterIndex.content;
                            ActiveClusterIndex = tempRef_ActiveClusterIndex.content;
                            RemoteIndex = tempRef_RemoteIndex.content;
                            if (this.CurrentSolution.Splittable_k_[RealClusterIndex] != 1) {
                                continue;
                            }
                            double Mvalue = this.CurrentSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex];
                            if (NumPowerIterationsLoopVar == 0) {
                                FindClusterFirstTerm.addapoint(threadIndex, Mvalue, RealClusterIndex);
                            }
                            double multiplier = 0.0;
                            for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension;
                                 VectorIndex++) {
                                PartVector[VectorIndex] = this.CurrentSolution.Y_k_i_[RealClusterIndex][VectorIndex] -
                                        Program.PointPosition[alpha][VectorIndex];
                                multiplier += PartVector[VectorIndex] * CenterEigenvector[RealClusterIndex][VectorIndex];
                            }
                            FindNumberScalarProducts.addapoint(threadIndex, 1.0);
                            double wgt = Mvalue * multiplier;
                            for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension;
                                 VectorIndex++) {
                                PartVector[VectorIndex] *= wgt;
                            }
                            FindNewPowerVectors.addapoint(threadIndex, PartVector, RealClusterIndex);
                        }
                    } // End Loop over points

                }); // End loop initialing Point dependent quantities
            } catch (SuspendableException e) {
                DAVectorUtility.printAndThrowRuntimeException(e.getMessage());
            }

            FindNewPowerVectors.sumoverthreadsandmpi();
			for (int ClusterIndex = 0; ClusterIndex < this.CurrentSolution.Ncent_Global; ClusterIndex++)
			{
				if (this.CurrentSolution.LocalStatus[ClusterIndex] != 1)
				{
					continue;
				}
				if ((this.CurrentSolution.Splittable_k_[ClusterIndex] != 1) || (this.CenterEigenconvergence[ClusterIndex] != 0))
				{
					continue;
				}
				double[] sums = new double[3]; // Old.New Old.Old New.New
				for (int loop = 0; loop < 3; loop++)
				{
					sums[loop] = 0.0;
				}
				for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
				{
					int TotalIndex = VectorIndex + ClusterIndex * Program.ParameterVectorDimension;
					double newvalue = FindNewPowerVectors.TotalVectorSum[TotalIndex];
					double oldvalue = CenterEigenvector[ClusterIndex][VectorIndex];
					sums[0] += oldvalue * newvalue;
					sums[1] += oldvalue * oldvalue;
					sums[2] += newvalue * newvalue;
					CenterEigenvector[ClusterIndex][VectorIndex] = newvalue;
				}

				//  Decide if finished and set eigenvalue
				double CandidateEigenvalue = sums[0] / sums[1];
				boolean LegalEigenvalue = (CandidateEigenvalue > 0.0);
				LegalEigenvalue = DAVectorUtility.SynchronizeMPIBoolean(LegalEigenvalue);

				//	Check if converged
				//	Do this in one process ONLY 
				if ((NumPowerIterations > 5) && LegalEigenvalue)
				{ // Arbitrary choice for Number of Power Iterations Cut

					int EigenvalueDone = 0;
					if (DAVectorUtility.MPI_Rank == 0)
					{ // Decisions can only be made in one process
						if (Math.abs(CandidateEigenvalue - this.CenterEigenvalue[ClusterIndex]) > CandidateEigenvalue * Program.eigenvaluechange)
						{
							++EigenvalueDone;
						}
						double delta = sums[2] - 2.0 * sums[0] * CandidateEigenvalue + sums[1] * CandidateEigenvalue * CandidateEigenvalue; // (Ax- Eigenvalue*Axold)**2
						if (Math.abs(delta) > CandidateEigenvalue * CandidateEigenvalue * Program.eigenvectorchange)
						{
							++EigenvalueDone;
						}
					} // End Test on Convergence
					EigenvalueDone = DAVectorUtility.SynchronizeMPIInteger(EigenvalueDone);

					if (EigenvalueDone == 0)
					{
						this.CenterEigenconvergence[ClusterIndex] = 1 + NumPowerIterations;
					}
				}
				this.CenterEigenvalue[ClusterIndex] = CandidateEigenvalue;

				//  Normalize current Power Vector to 1
				double wgt = 1.0 / Math.sqrt(sums[2]);
				for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
				{
					CenterEigenvector[ClusterIndex][VectorIndex] *= wgt;
				}
			} // End Loop over Clusters

		} //  End Loop over NumPowerIterations

		FindClusterFirstTerm.sumoverthreadsandmpi();
		FindNumberScalarProducts.sumoverthreadsandmpi();
		Program.SumEigenSPCalcs += FindNumberScalarProducts.Total;

		for (int ClusterIndex = 0; ClusterIndex < this.CurrentSolution.Ncent_Global; ClusterIndex++)
		{
			this.CenterEigenstatus[ClusterIndex] = 0;
			if (this.CurrentSolution.LocalStatus[ClusterIndex] != 1)
			{
				continue;
			}
			if ((this.CurrentSolution.Splittable_k_[ClusterIndex] != 1) || (this.CenterEigenconvergence[ClusterIndex] <= 0))
			{
				continue;
			}
			this.CenterEigenstatus[ClusterIndex] = 1;
			this.FirstTerm[ClusterIndex] = FindClusterFirstTerm.TotalVectorSum[ClusterIndex];
			double tmp = this.CenterEigenvalue[ClusterIndex] / this.CurrentSolution.Temperature;
			this.CenterEigenvalue[ClusterIndex] = this.FirstTerm[ClusterIndex] - tmp;
		}


	} // End SetEigenvaluesIteratively(ClusteringSolution Solution)

	public final void getEigenvaluefromIteration(int ClusterIndex)
	{ // Eigenvector stored in Solution

		if (this.CenterEigenconvergence[ClusterIndex] <= 0)
		{
			this.CenterEigenstatus[ClusterIndex] = 0;
		}
		this.EigenStatus = this.CenterEigenstatus[ClusterIndex];
		if (this.EigenStatus <= 0)
		{
			this.Eigenvalue = 1.0;
		}
		else
		{
			this.Eigenvalue = this.CenterEigenvalue[ClusterIndex];
		}

	} // End getEigenvaluefromIteration(int ClusterIndex)

} // End VectorClass
 // End namespace edu.indiana.soic.spidal.davs
