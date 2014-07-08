package edu.indiana.soic.spidal.configuration;

import edu.indiana.soic.spidal.configuration.sections.DAVectorSpongeSection;

public class ConfigurationMgr {
    private String configurationFilePath;
    public DAVectorSpongeSection daVectorSpongeSection;

    public ConfigurationMgr(String configurationFilePath) {
        this.configurationFilePath = configurationFilePath;
        daVectorSpongeSection = new DAVectorSpongeSection(configurationFilePath);


    }

    public static ConfigurationMgr LoadConfiguration(String configurationFilePath){
        // TODO - Fix configuration management
        return new ConfigurationMgr(configurationFilePath);
    }

    // TODO - Tester method - remove when releasing
    /*public static void main(String[] args) {
        ConfigurationMgr mgr = ConfigurationMgr.LoadConfiguration("E:\\Sali\\git\\iu\\Main\\Java\\SalsaTPL_r4085_co7.17.13_DAVectorSponge_Serial\\Samples\\DAVectorSpong_Sample.properties");
        DAVectorSpongeSection da = mgr.daVectorSpongeSection;

        System.out.println(da.ClusterFile);
        System.out.println(da.DistanceMatrixFile);
        System.out.println(da.LabelFile);
        System.out.println(da.TimingFile);
        System.out.println(da.SummaryFile);

        System.out.println(da.UseSponge);
        System.out.println(da.SpongeFactor1);
        System.out.println(da.SpongeFactor2);
        System.out.println(da.SpongePOption);
        System.out.println(da.SpongePWeight);
        System.out.println(da.CreateSpongeScaledSquaredWidth);
        System.out.println(da.ContinuousClustering);
        System.out.println(da.ParameterVectorDimension);

        System.out.println(da.SpongeTemperature1);
        System.out.println(da.SpongeTemperature2);
        System.out.println(da.RestartTemperature);

        System.out.println(da.NumberDataPoints);
        System.out.println(da.SelectedInputLabel);
        System.out.println(da.OutputFileType);
        System.out.println(da.Replicate);

        System.out.println(da.SigmaMethod);
        System.out.println(da.FinalTargetTemperature);
        System.out.println(da.FinalTargetSigma0);
        System.out.println(da.InitialSigma0);

        System.out.println(da.ClusterCountOutput);
        System.out.println(da.NumberNearbyClusters);
        System.out.println(da.NearbySpongePointLimit);

        System.out.println(da.ProcessingOption);

        System.out.println(da.CacheLineSize);
        System.out.println(da.ClusterPrintNumber);
        System.out.println(da.PrintInterval);
        System.out.println(da.RemovalDiagnosticPrint);

        for(double t : da.MagicTemperatures){
            System.out.print(t + ",");
        }

        System.out.println(da.MagicIndex);

        System.out.println(da.MaxNcentPerNode);
        System.out.println(da.MaxNcentTotal);
        System.out.println(da.TargetNcentPerPoint);
        System.out.println(da.TargetMinimumNcentPerPoint);
        System.out.println(da.MaxNcentPerPoint);

        System.out.println(da.MaxIntegerComponents);
        System.out.println(da.MaxDoubleComponents);
        System.out.println(da.MaxMPITransportBuffer);
        System.out.println(da.MaxNumberAccumulationsPerNode);
        System.out.println(da.MaxTransportedClusterStorage);

        System.out.println(da.ExpArgumentCut1);
        System.out.println(da.ExpArgumentCut2);
        System.out.println(da.ExpArgumentCut3);
        System.out.println(da.Tminimum);

        System.out.println(da.InitialNcent);
        System.out.println(da.MinimumCountForClusterCk);
        System.out.println(da.MinimumCountForClusterCkWithSponge);
        System.out.println(da.MinimuCountForClusterPoints);
        System.out.println(da.CountForClusterCkToBeZero);
        System.out.println(da.AddSpongeScaledWidthSquared);

        System.out.println(da.InitialCoolingFactor);
        System.out.println(da.FineCoolingFactor);
        System.out.println(da.WaitIterations);

        System.out.println(da.IterationAtEnd);
        System.out.println(da.ConvergenceLoopLimit);

        System.out.println(da.FreezingLimit);
        System.out.println(da.MalphaMaxChange);
        System.out.println(da.MaxNumberSplitClusters);
        System.out.println(da.ConvergeIntermediateClusters);
        System.out.println(da.TooSmallToSplit);
        System.out.println(da.ScaledWidthSquaredToSplit);

        System.out.println(da.ClusterLimitForDistribution);
        System.out.println(da.TemperatureLimitForDistribution);

        System.out.println(da.DebugPrintOption);
        System.out.println(da.ConsoleDebugOutput);
    }*/

}
