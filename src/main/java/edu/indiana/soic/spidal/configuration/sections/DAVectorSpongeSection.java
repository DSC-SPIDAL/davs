package edu.indiana.soic.spidal.configuration.sections;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.Properties;

public class DAVectorSpongeSection {

    public DAVectorSpongeSection(String configurationFilePath)
    {
        Properties p = new Properties();
        try {
            p.load(new FileInputStream(configurationFilePath));

            ClusterFile = getProperty(p,"ClusterFile","cluster.txt");
            DistanceMatrixFile = getProperty(p,"DistanceMatrixFile", "distance.bin");
            LabelFile = getProperty(p,"LabelFile", "labels.txt");
            ComparisonClusterFile = getProperty(p,"ComparisonClusterFile", "ComparisonClusterFile.txt");
            RestartClusterFile = getProperty(p,"RestartClusterFile", "RestartClusterFile.txt");
            TimingFile = getProperty(p,"TimingFile", "timings.txt");
            SummaryFile = getProperty(p,"SummaryFile","summary.txt");


            UseSponge = Boolean.parseBoolean(getProperty(p,"UseSponge", "false"));
            SpongeFactor1 = Double.parseDouble(getProperty(p,"SpongeFactor1", "3.0"));
            SpongeFactor2 = Double.parseDouble(getProperty(p,"SpongeFactor2", "3.0"));
            SpongePOption = Integer.parseInt(getProperty(p,"SpongePOption", "1"));
            SpongePWeight = Double.parseDouble(getProperty(p,"SpongePWeight", "0.1"));
            CreateSpongeScaledSquaredWidth = Double.parseDouble(getProperty(p,"CreateSpongeScaledSquaredWidth","-1.0" ));
            ContinuousClustering = Boolean.parseBoolean(getProperty(p,"ContinuousClustering","true"));
            ParameterVectorDimension = Integer.parseInt(getProperty(p,"ParameterVectorDimension","2"));

            SpongeTemperature1 = Double.parseDouble(getProperty(p,"SpongeTemperature1","-1.0" ));
            SpongeTemperature2 = Double.parseDouble(getProperty(p,"SpongeTemperature2","-1.0" ));
            RestartTemperature = Double.parseDouble(getProperty(p,"RestartTemperature","-1.0" ));

            NumberDataPoints = Integer.parseInt(getProperty(p,"NumberDataPoints","-1"));
            SelectedInputLabel = Integer.parseInt(getProperty(p,"SelectedInputLabel","6"));
            InputFileType = Integer.parseInt(getProperty(p,"InputFileType","0"));
            Replicate = Integer.parseInt(getProperty(p,"Replicate","1"));


            CompareSolution = Integer.parseInt(getProperty(p,"CompareSolution","-1"));
            ComparisonInputFileType = Integer.parseInt(getProperty(p,"ComparisonInputFileType","0"));
            ComparisonSelectedInputLabel = Integer.parseInt(getProperty(p,"ComparisonSelectedInputLabel","2"));
            RestartSelectedInputLabel = Integer.parseInt(getProperty(p,"RestartSelectedInputLabel","-100000000"));
            RestartInputFileType = Integer.parseInt(getProperty(p,"RestartInputFileType","1"));

            SigmaMethod = Integer.parseInt(getProperty(p,"SigmaMethod","0"));

            String SigmaVectorParametersIString = getProperty(p,"SigmaVectorParameters_i", "0.000598,2.35");
            String [] splits = SigmaVectorParametersIString.split(",");
            SigmaVectorParameters_i = new double[splits.length];
            for (int i = 0; i < splits.length; ++i){
                SigmaVectorParameters_i[i] = Double.parseDouble(splits[i]);
            }

            FinalTargetTemperature = Double.parseDouble(getProperty(p,"FinalTargetTemperature","3.0" ));
            FinalTargetSigma0 = Double.parseDouble(getProperty(p,"FinalTargetSigma0","0.0" ));
            InitialSigma0 = Double.parseDouble(getProperty(p,"InitialSigma0","0.0" ));

            ClusterCountOutput = Integer.parseInt(getProperty(p,"ClusterCountOutput", "0"));
            NumberNearbyClusters = Integer.parseInt(getProperty(p,"NumberNearbyClusters", "5"));
            NearbySpongePointLimit = Double.parseDouble(getProperty(p,"NearbySpongePointLimit", "-1.0"));

            ProcessingOption = Integer.parseInt(getProperty(p,"ProcessingOption", "0"));

            CacheLineSize = Integer.parseInt(getProperty(p,"CacheLineSize", "0"));
            ClusterPrintNumber = Integer.parseInt(getProperty(p,"ClusterPrintNumber", "5"));
            PrintInterval = Integer.parseInt(getProperty(p,"PrintInterval", "3"));
            RemovalDiagnosticPrint = Boolean.parseBoolean(getProperty(p,"RemovalDiagnosticPrint", "false"));

            String magicTemperaturesString = getProperty(p,"MagicTemperatures", "4.0,3.0,2.0,1.0,0.5");
            splits = magicTemperaturesString.split(",");
            MagicTemperatures = new double[splits.length];
            for (int i = 0; i < splits.length; ++i){
                MagicTemperatures[i] = Double.parseDouble(splits[i]);
            }

            MagicIndex = Integer.parseInt(getProperty(p,"MagicIndex", "0"));

            MaxNcentPerNode = Integer.parseInt(getProperty(p,"MaxNcentPerNode", "0"));
            MaxNcentTotal = Integer.parseInt(getProperty(p,"MaxNcentTotal", "0"));
            maxNcentCreated = Integer.parseInt(getProperty(p,"maxNcentCreated", "0"));
            TargetNcentPerPoint = Integer.parseInt(getProperty(p,"TargetNcentPerPoint", "20"));
            TargetMinimumNcentPerPoint = Integer.parseInt(getProperty(p,"TargetMinimumNcentPerPoint", "1"));
            MaxNcentPerPoint = Integer.parseInt(getProperty(p,"MaxNcentPerPoint", "25"));

            MaxIntegerComponents = Integer.parseInt(getProperty(p,"MaxIntegerComponents", "2"));
            MaxDoubleComponents = Integer.parseInt(getProperty(p,"MaxDoubleComponents", "3"));
            MaxMPITransportBuffer = Integer.parseInt(getProperty(p,"MaxMPITransportBuffer", "500"));
            MaxNumberAccumulationsPerNode = Integer.parseInt(getProperty(p,"MaxNumberAccumulationsPerNode", "30000"));
            MaxTransportedClusterStorage = Integer.parseInt(getProperty(p,"MaxTransportedClusterStorage", "500"));

            ExpArgumentCut1 = Double.parseDouble(getProperty(p,"ExpArgumentCut1", "20.0"));
            ExpArgumentCut2 = Double.parseDouble(getProperty(p,"ExpArgumentCut2", "40.0"));
            ExpArgumentCut3 = Double.parseDouble(getProperty(p,"ExpArgumentCut3", "50.0"));
            Tminimum = Double.parseDouble(getProperty(p,"Tminimum", "-1000.0"));

            InitialNcent = Integer.parseInt(getProperty(p,"InitialNcent", "1"));
            MinimumCountForClusterCk = Double.parseDouble(getProperty(p,"MinimumCountForClusterCk", "1.0"));
            MinimumCountForClusterCkWithSponge = Double.parseDouble(getProperty(p,"MinimumCountForClusterCkWithSponge", "1.5"));
            MinimuCountForClusterPoints = Integer.parseInt(getProperty(p,"MinimuCountForClusterPoints", "2"));
            CountForClusterCkToBeZero = Double.parseDouble(getProperty(p,"CountForClusterCkToBeZero", "0.001"));
            AddSpongeScaledWidthSquared = Double.parseDouble(getProperty(p,"AddSpongeScaledWidthSquared", "-1.0"));

            InitialCoolingFactor = Double.parseDouble(getProperty(p,"InitialCoolingFactor", "0.9875"));
            InitialCoolingFactor1 = Double.parseDouble(getProperty(p,"InitialCoolingFactor1", "0.9875"));
            InitialCoolingFactor2 = Double.parseDouble(getProperty(p,"InitialCoolingFactor2", "0.99375"));
            FineCoolingFactor = Double.parseDouble(getProperty(p,"FineCoolingFactor", "0.9975"));
            FineCoolingFactor1 = Double.parseDouble(getProperty(p,"FineCoolingFactor1", "0.9975"));
            FineCoolingFactor2 = Double.parseDouble(getProperty(p,"FineCoolingFactor2", "0.999375"));
            CoolingTemperatureSwitch = Double.parseDouble(getProperty(p,"CoolingTemperatureSwitch", "12.0"));
            WaitIterations = Integer.parseInt(getProperty(p,"WaitIterations", "10"));
            WaitIterationsConverge = Integer.parseInt(getProperty(p,"WaitIterationsConverge", "4"));

            IterationAtEnd = Integer.parseInt(getProperty(p,"IterationAtEnd", "2000"));
            ConvergenceLoopLimit = Integer.parseInt(getProperty(p,"ConvergenceLoopLimit", "20"));

            FreezingLimit = Double.parseDouble(getProperty(p,"FreezingLimit", "0.002"));
            MalphaMaxChange = Double.parseDouble(getProperty(p,"MalphaMaxChange", "0.005"));
            MalphaMaxChange1 = Double.parseDouble(getProperty(p,"MalphaMaxChange1", "0.005"));
            MaxNumberSplitClusters = Integer.parseInt(getProperty(p,"MaxNumberSplitClusters", "3"));
            ConvergeIntermediateClusters = Boolean.parseBoolean(getProperty(p,"ConvergeIntermediateClusters", "false"));
            TooSmallToSplit = Double.parseDouble(getProperty(p,"TooSmallToSplit", "4.0"));
            MinimumScaledWidthSquaredToSplit = Double.parseDouble(getProperty(p,"MinimumScaledWidthSquaredToSplit", "1.5"));
            ScaledSquaredDistanceAtClosenessTest = Double.parseDouble(getProperty(p,"ScaledSquaredDistanceAtClosenessTest", "0.5"));

            ClusterLimitForDistribution = Integer.parseInt(getProperty(p,"ClusterLimitForDistribution", "-1"));
            TemperatureLimitForDistribution = Double.parseDouble(getProperty(p,"TemperatureLimitForDistribution", "-1.0"));
            TemperatureForClosenessTest = Double.parseDouble(getProperty(p,"TemperatureForClosenessTest", "4.0"));

            DebugPrintOption = Integer.parseInt(getProperty(p,"DebugPrintOption", "1"));
            ConsoleDebugOutput = Boolean.parseBoolean(getProperty(p,"ConsoleDebugOutput", "true"));

            LCMSmode = Integer.parseInt(getProperty(p,"LCMSmode","1"));

            //Cluster Quality Parameters

            HistogramOccupationMax = Integer.parseInt(getProperty(p,"HistogramOccupationMax","202"));
            HistogramDistancesMax = Integer.parseInt(getProperty(p,"HistogramDistancesMax","200"));

        } catch (IOException e) {
            e.printStackTrace();
        }


    }

    private static String getProperty(Properties p, String name, String def) {
        String val = System.getProperty(name);
        if (val == null) {
            if (def != null) {
                val = p.getProperty(name, def);
            } else {
                val = p.getProperty(name);
            }
        }
        return val;
    }

    public boolean UseSponge;
    public double SpongeFactor1;
    public double SpongeFactor2;
    public int SpongePOption;
    public double SpongePWeight;
    public double CreateSpongeScaledSquaredWidth;
    public boolean ContinuousClustering;
    public int ParameterVectorDimension;

    public double SpongeTemperature1;
    public double SpongeTemperature2;
    public double RestartTemperature;

    public int NumberDataPoints;
    public int SelectedInputLabel;
    public int InputFileType;
    public int Replicate;

    public int CompareSolution;
    public int ComparisonInputFileType;
    public int ComparisonSelectedInputLabel;
    public int RestartSelectedInputLabel;
    public int RestartInputFileType;

    public int SigmaMethod;
    public double[] SigmaVectorParameters_i;
    public double FinalTargetTemperature;
    public double FinalTargetSigma0;
    public double InitialSigma0;

    public int ClusterCountOutput;
    public int NumberNearbyClusters;
    public double NearbySpongePointLimit;

    public int ProcessingOption;

    public int CacheLineSize;
    public int ClusterPrintNumber;
    public int PrintInterval;
    public boolean RemovalDiagnosticPrint;
    public double[] MagicTemperatures;
    public int MagicIndex;

    public int MaxNcentPerNode;
    public int MaxNcentTotal;
    public int maxNcentCreated;
    public int TargetNcentPerPoint;
    public int TargetMinimumNcentPerPoint;
    public int MaxNcentPerPoint;

    public int MaxIntegerComponents;
    public int MaxDoubleComponents;
    public int MaxMPITransportBuffer;
    public int MaxNumberAccumulationsPerNode;
    public int MaxTransportedClusterStorage;

    public double ExpArgumentCut1;
    public double ExpArgumentCut2;
    public double ExpArgumentCut3;
    public double Tminimum;

    public int InitialNcent;
    public double MinimumCountForClusterCk;
    public double MinimumCountForClusterCkWithSponge;
    public int MinimuCountForClusterPoints;
    public double CountForClusterCkToBeZero;
    public double AddSpongeScaledWidthSquared;

    public double InitialCoolingFactor;
    public double InitialCoolingFactor1;
    public double InitialCoolingFactor2;
    public double FineCoolingFactor;
    public double FineCoolingFactor1;
    public double FineCoolingFactor2;
    public double CoolingTemperatureSwitch;
    public int WaitIterations;
    public int WaitIterationsConverge;

    public int IterationAtEnd;
    public int ConvergenceLoopLimit;

    public double FreezingLimit;
    public double MalphaMaxChange;
    public double MalphaMaxChange1;
    public int MaxNumberSplitClusters;
    public boolean ConvergeIntermediateClusters;
    public double TooSmallToSplit;
    public double MinimumScaledWidthSquaredToSplit;
    public double ScaledSquaredDistanceAtClosenessTest;

    public int ClusterLimitForDistribution;
    public double TemperatureLimitForDistribution;
    public double TemperatureForClosenessTest;

    public int NodeCount;
    public int ThreadCount;

    public int DebugPrintOption;
    public boolean ConsoleDebugOutput;

    public String ClusterFile;
    public String DistanceMatrixFile;
    public String LabelFile;
    public String ComparisonClusterFile;
    public String RestartClusterFile;
    public String TimingFile;
    public String SummaryFile;

    public int LCMSmode;

    //Cluster Quality Parameters
    public int HistogramOccupationMax;
    public int HistogramDistancesMax;

    public boolean isUseSponge() {
        return UseSponge;
    }

    public void setUseSponge(boolean useSponge) {
        UseSponge = useSponge;
    }

    public double getSpongeFactor1() {
        return SpongeFactor1;
    }

    public void setSpongeFactor1(double spongeFactor1) {
        SpongeFactor1 = spongeFactor1;
    }

    public double getSpongeFactor2() {
        return SpongeFactor2;
    }

    public void setSpongeFactor2(double spongeFactor2) {
        SpongeFactor2 = spongeFactor2;
    }

    public int getSpongePOption() {
        return SpongePOption;
    }

    public void setSpongePOption(int spongePOption) {
        SpongePOption = spongePOption;
    }

    public double getSpongePWeight() {
        return SpongePWeight;
    }

    public void setSpongePWeight(double spongePWeight) {
        SpongePWeight = spongePWeight;
    }

    public double getCreateSpongeScaledSquaredWidth() {
        return CreateSpongeScaledSquaredWidth;
    }

    public void setCreateSpongeScaledSquaredWidth(double createSpongeScaledSquaredWidth) {
        CreateSpongeScaledSquaredWidth = createSpongeScaledSquaredWidth;
    }

    public boolean isContinuousClustering() {
        return ContinuousClustering;
    }

    public void setContinuousClustering(boolean continuousClustering) {
        ContinuousClustering = continuousClustering;
    }

    public int getParameterVectorDimension() {
        return ParameterVectorDimension;
    }

    public void setParameterVectorDimension(int parameterVectorDimension) {
        ParameterVectorDimension = parameterVectorDimension;
    }

    public double getSpongeTemperature1() {
        return SpongeTemperature1;
    }

    public void setSpongeTemperature1(double spongeTemperature1) {
        SpongeTemperature1 = spongeTemperature1;
    }

    public double getSpongeTemperature2() {
        return SpongeTemperature2;
    }

    public void setSpongeTemperature2(double spongeTemperature2) {
        SpongeTemperature2 = spongeTemperature2;
    }

    public double getRestartTemperature() {
        return RestartTemperature;
    }

    public void setRestartTemperature(double restartTemperature) {
        RestartTemperature = restartTemperature;
    }

    public int getNumberDataPoints() {
        return NumberDataPoints;
    }

    public void setNumberDataPoints(int numberDataPoints) {
        NumberDataPoints = numberDataPoints;
    }

    public int getSelectedInputLabel() {
        return SelectedInputLabel;
    }

    public void setSelectedInputLabel(int selectedInputLabel) {
        SelectedInputLabel = selectedInputLabel;
    }

    public int getInputFileType() {
        return InputFileType;
    }

    public void setInputFileType(int inputFileType) {
        InputFileType = inputFileType;
    }

    public int getReplicate() {
        return Replicate;
    }

    public void setReplicate(int replicate) {
        Replicate = replicate;
    }

    public int getCompareSolution() {
        return CompareSolution;
    }

    public void setCompareSolution(int compareSolution) {
        CompareSolution = compareSolution;
    }

    public int getComparisonInputFileType() {
        return ComparisonInputFileType;
    }

    public void setComparisonInputFileType(int comparisonInputFileType) {
        ComparisonInputFileType = comparisonInputFileType;
    }

    public int getComparisonSelectedInputLabel() {
        return ComparisonSelectedInputLabel;
    }

    public void setComparisonSelectedInputLabel(int comparisonSelectedInputLabel) {
        ComparisonSelectedInputLabel = comparisonSelectedInputLabel;
    }

    public int getRestartSelectedInputLabel() {
        return RestartSelectedInputLabel;
    }

    public void setRestartSelectedInputLabel(int restartSelectedInputLabel) {
        RestartSelectedInputLabel = restartSelectedInputLabel;
    }

    public int getRestartInputFileType() {
        return RestartInputFileType;
    }

    public void setRestartInputFileType(int restartInputFileType) {
        RestartInputFileType = restartInputFileType;
    }

    public int getSigmaMethod() {
        return SigmaMethod;
    }

    public void setSigmaMethod(int sigmaMethod) {
        SigmaMethod = sigmaMethod;
    }

    public double[] getSigmaVectorParameters_i() {
        return SigmaVectorParameters_i;
    }

    public void setSigmaVectorParameters_i(double[] sigmaVectorParameters_i) {
        SigmaVectorParameters_i = sigmaVectorParameters_i;
    }

    public double getFinalTargetTemperature() {
        return FinalTargetTemperature;
    }

    public void setFinalTargetTemperature(double finalTargetTemperature) {
        FinalTargetTemperature = finalTargetTemperature;
    }

    public double getFinalTargetSigma0() {
        return FinalTargetSigma0;
    }

    public void setFinalTargetSigma0(double finalTargetSigma0) {
        FinalTargetSigma0 = finalTargetSigma0;
    }

    public double getInitialSigma0() {
        return InitialSigma0;
    }

    public void setInitialSigma0(double initialSigma0) {
        InitialSigma0 = initialSigma0;
    }

    public int getClusterCountOutput() {
        return ClusterCountOutput;
    }

    public void setClusterCountOutput(int clusterCountOutput) {
        ClusterCountOutput = clusterCountOutput;
    }

    public int getNumberNearbyClusters() {
        return NumberNearbyClusters;
    }

    public void setNumberNearbyClusters(int numberNearbyClusters) {
        NumberNearbyClusters = numberNearbyClusters;
    }

    public double getNearbySpongePointLimit() {
        return NearbySpongePointLimit;
    }

    public void setNearbySpongePointLimit(double nearbySpongePointLimit) {
        NearbySpongePointLimit = nearbySpongePointLimit;
    }

    public int getProcessingOption() {
        return ProcessingOption;
    }

    public void setProcessingOption(int processingOption) {
        ProcessingOption = processingOption;
    }

    public int getCacheLineSize() {
        return CacheLineSize;
    }

    public void setCacheLineSize(int cacheLineSize) {
        CacheLineSize = cacheLineSize;
    }

    public int getClusterPrintNumber() {
        return ClusterPrintNumber;
    }

    public void setClusterPrintNumber(int clusterPrintNumber) {
        ClusterPrintNumber = clusterPrintNumber;
    }

    public int getPrintInterval() {
        return PrintInterval;
    }

    public void setPrintInterval(int printInterval) {
        PrintInterval = printInterval;
    }

    public boolean isRemovalDiagnosticPrint() {
        return RemovalDiagnosticPrint;
    }

    public void setRemovalDiagnosticPrint(boolean removalDiagnosticPrint) {
        RemovalDiagnosticPrint = removalDiagnosticPrint;
    }

    public double[] getMagicTemperatures() {
        return MagicTemperatures;
    }

    public void setMagicTemperatures(double[] magicTemperatures) {
        MagicTemperatures = magicTemperatures;
    }

    public int getMagicIndex() {
        return MagicIndex;
    }

    public void setMagicIndex(int magicIndex) {
        MagicIndex = magicIndex;
    }

    public int getMaxNcentPerNode() {
        return MaxNcentPerNode;
    }

    public void setMaxNcentPerNode(int maxNcentPerNode) {
        MaxNcentPerNode = maxNcentPerNode;
    }

    public int getMaxNcentTotal() {
        return MaxNcentTotal;
    }

    public void setMaxNcentTotal(int maxNcentTotal) {
        MaxNcentTotal = maxNcentTotal;
    }

    public int getMaxNcentCreated() {
        return maxNcentCreated;
    }

    public void setMaxNcentCreated(int maxNcentCreated) {
        this.maxNcentCreated = maxNcentCreated;
    }

    public int getTargetNcentPerPoint() {
        return TargetNcentPerPoint;
    }

    public void setTargetNcentPerPoint(int targetNcentPerPoint) {
        TargetNcentPerPoint = targetNcentPerPoint;
    }

    public int getTargetMinimumNcentPerPoint() {
        return TargetMinimumNcentPerPoint;
    }

    public void setTargetMinimumNcentPerPoint(int targetMinimumNcentPerPoint) {
        TargetMinimumNcentPerPoint = targetMinimumNcentPerPoint;
    }

    public int getMaxNcentPerPoint() {
        return MaxNcentPerPoint;
    }

    public void setMaxNcentPerPoint(int maxNcentPerPoint) {
        MaxNcentPerPoint = maxNcentPerPoint;
    }

    public int getMaxIntegerComponents() {
        return MaxIntegerComponents;
    }

    public void setMaxIntegerComponents(int maxIntegerComponents) {
        MaxIntegerComponents = maxIntegerComponents;
    }

    public int getMaxDoubleComponents() {
        return MaxDoubleComponents;
    }

    public void setMaxDoubleComponents(int maxDoubleComponents) {
        MaxDoubleComponents = maxDoubleComponents;
    }

    public int getMaxMPITransportBuffer() {
        return MaxMPITransportBuffer;
    }

    public void setMaxMPITransportBuffer(int maxMPITransportBuffer) {
        MaxMPITransportBuffer = maxMPITransportBuffer;
    }

    public int getMaxNumberAccumulationsPerNode() {
        return MaxNumberAccumulationsPerNode;
    }

    public void setMaxNumberAccumulationsPerNode(int maxNumberAccumulationsPerNode) {
        MaxNumberAccumulationsPerNode = maxNumberAccumulationsPerNode;
    }

    public int getMaxTransportedClusterStorage() {
        return MaxTransportedClusterStorage;
    }

    public void setMaxTransportedClusterStorage(int maxTransportedClusterStorage) {
        MaxTransportedClusterStorage = maxTransportedClusterStorage;
    }

    public double getExpArgumentCut1() {
        return ExpArgumentCut1;
    }

    public void setExpArgumentCut1(double expArgumentCut1) {
        ExpArgumentCut1 = expArgumentCut1;
    }

    public double getExpArgumentCut2() {
        return ExpArgumentCut2;
    }

    public void setExpArgumentCut2(double expArgumentCut2) {
        ExpArgumentCut2 = expArgumentCut2;
    }

    public double getExpArgumentCut3() {
        return ExpArgumentCut3;
    }

    public void setExpArgumentCut3(double expArgumentCut3) {
        ExpArgumentCut3 = expArgumentCut3;
    }

    public double getTminimum() {
        return Tminimum;
    }

    public void setTminimum(double tminimum) {
        Tminimum = tminimum;
    }

    public int getInitialNcent() {
        return InitialNcent;
    }

    public void setInitialNcent(int initialNcent) {
        InitialNcent = initialNcent;
    }

    public double getMinimumCountForClusterCk() {
        return MinimumCountForClusterCk;
    }

    public void setMinimumCountForClusterCk(double minimumCountForClusterCk) {
        MinimumCountForClusterCk = minimumCountForClusterCk;
    }

    public double getMinimumCountForClusterCkWithSponge() {
        return MinimumCountForClusterCkWithSponge;
    }

    public void setMinimumCountForClusterCkWithSponge(double minimumCountForClusterCkWithSponge) {
        MinimumCountForClusterCkWithSponge = minimumCountForClusterCkWithSponge;
    }

    public int getMinimuCountForClusterPoints() {
        return MinimuCountForClusterPoints;
    }

    public void setMinimuCountForClusterPoints(int minimuCountForClusterPoints) {
        MinimuCountForClusterPoints = minimuCountForClusterPoints;
    }

    public double getCountForClusterCkToBeZero() {
        return CountForClusterCkToBeZero;
    }

    public void setCountForClusterCkToBeZero(double countForClusterCkToBeZero) {
        CountForClusterCkToBeZero = countForClusterCkToBeZero;
    }

    public double getAddSpongeScaledWidthSquared() {
        return AddSpongeScaledWidthSquared;
    }

    public void setAddSpongeScaledWidthSquared(double addSpongeScaledWidthSquared) {
        AddSpongeScaledWidthSquared = addSpongeScaledWidthSquared;
    }

    public double getInitialCoolingFactor() {
        return InitialCoolingFactor;
    }

    public void setInitialCoolingFactor(double initialCoolingFactor) {
        InitialCoolingFactor = initialCoolingFactor;
    }

    public double getInitialCoolingFactor1() {
        return InitialCoolingFactor1;
    }

    public void setInitialCoolingFactor1(double initialCoolingFactor1) {
        InitialCoolingFactor1 = initialCoolingFactor1;
    }

    public double getInitialCoolingFactor2() {
        return InitialCoolingFactor2;
    }

    public void setInitialCoolingFactor2(double initialCoolingFactor2) {
        InitialCoolingFactor2 = initialCoolingFactor2;
    }

    public double getFineCoolingFactor() {
        return FineCoolingFactor;
    }

    public void setFineCoolingFactor(double fineCoolingFactor) {
        FineCoolingFactor = fineCoolingFactor;
    }

    public double getFineCoolingFactor1() {
        return FineCoolingFactor1;
    }

    public void setFineCoolingFactor1(double fineCoolingFactor1) {
        FineCoolingFactor1 = fineCoolingFactor1;
    }

    public double getFineCoolingFactor2() {
        return FineCoolingFactor2;
    }

    public void setFineCoolingFactor2(double fineCoolingFactor2) {
        FineCoolingFactor2 = fineCoolingFactor2;
    }

    public double getCoolingTemperatureSwitch() {
        return CoolingTemperatureSwitch;
    }

    public void setCoolingTemperatureSwitch(double coolingTemperatureSwitch) {
        CoolingTemperatureSwitch = coolingTemperatureSwitch;
    }

    public int getWaitIterations() {
        return WaitIterations;
    }

    public void setWaitIterations(int waitIterations) {
        WaitIterations = waitIterations;
    }

    public int getWaitIterationsConverge() {
        return WaitIterationsConverge;
    }

    public void setWaitIterationsConverge(int waitIterationsConverge) {
        WaitIterationsConverge = waitIterationsConverge;
    }

    public int getIterationAtEnd() {
        return IterationAtEnd;
    }

    public void setIterationAtEnd(int iterationAtEnd) {
        IterationAtEnd = iterationAtEnd;
    }

    public int getConvergenceLoopLimit() {
        return ConvergenceLoopLimit;
    }

    public void setConvergenceLoopLimit(int convergenceLoopLimit) {
        ConvergenceLoopLimit = convergenceLoopLimit;
    }

    public double getFreezingLimit() {
        return FreezingLimit;
    }

    public void setFreezingLimit(double freezingLimit) {
        FreezingLimit = freezingLimit;
    }

    public double getMalphaMaxChange() {
        return MalphaMaxChange;
    }

    public void setMalphaMaxChange(double malphaMaxChange) {
        MalphaMaxChange = malphaMaxChange;
    }

    public double getMalphaMaxChange1() {
        return MalphaMaxChange1;
    }

    public void setMalphaMaxChange1(double malphaMaxChange1) {
        MalphaMaxChange1 = malphaMaxChange1;
    }

    public int getMaxNumberSplitClusters() {
        return MaxNumberSplitClusters;
    }

    public void setMaxNumberSplitClusters(int maxNumberSplitClusters) {
        MaxNumberSplitClusters = maxNumberSplitClusters;
    }

    public boolean isConvergeIntermediateClusters() {
        return ConvergeIntermediateClusters;
    }

    public void setConvergeIntermediateClusters(boolean convergeIntermediateClusters) {
        ConvergeIntermediateClusters = convergeIntermediateClusters;
    }

    public double getTooSmallToSplit() {
        return TooSmallToSplit;
    }

    public void setTooSmallToSplit(double tooSmallToSplit) {
        TooSmallToSplit = tooSmallToSplit;
    }

    public double getMinimumScaledWidthSquaredToSplit() {
        return MinimumScaledWidthSquaredToSplit;
    }

    public void setMinimumScaledWidthSquaredToSplit(double minimumScaledWidthSquaredToSplit) {
        MinimumScaledWidthSquaredToSplit = minimumScaledWidthSquaredToSplit;
    }

    public double getScaledSquaredDistanceAtClosenessTest() {
        return ScaledSquaredDistanceAtClosenessTest;
    }

    public void setScaledSquaredDistanceAtClosenessTest(double scaledSquaredDistanceAtClosenessTest) {
        ScaledSquaredDistanceAtClosenessTest = scaledSquaredDistanceAtClosenessTest;
    }

    public int getClusterLimitForDistribution() {
        return ClusterLimitForDistribution;
    }

    public void setClusterLimitForDistribution(int clusterLimitForDistribution) {
        ClusterLimitForDistribution = clusterLimitForDistribution;
    }

    public double getTemperatureLimitForDistribution() {
        return TemperatureLimitForDistribution;
    }

    public void setTemperatureLimitForDistribution(double temperatureLimitForDistribution) {
        TemperatureLimitForDistribution = temperatureLimitForDistribution;
    }

    public double getTemperatureForClosenessTest() {
        return TemperatureForClosenessTest;
    }

    public void setTemperatureForClosenessTest(double temperatureForClosenessTest) {
        TemperatureForClosenessTest = temperatureForClosenessTest;
    }

    public int getNodeCount() {
        return NodeCount;
    }

    public void setNodeCount(int nodeCount) {
        NodeCount = nodeCount;
    }

    public int getThreadCount() {
        return ThreadCount;
    }

    public void setThreadCount(int threadCount) {
        ThreadCount = threadCount;
    }

    public int getDebugPrintOption() {
        return DebugPrintOption;
    }

    public void setDebugPrintOption(int debugPrintOption) {
        DebugPrintOption = debugPrintOption;
    }

    public boolean isConsoleDebugOutput() {
        return ConsoleDebugOutput;
    }

    public void setConsoleDebugOutput(boolean consoleDebugOutput) {
        ConsoleDebugOutput = consoleDebugOutput;
    }

    public String getClusterFile() {
        return ClusterFile;
    }

    public void setClusterFile(String clusterFile) {
        ClusterFile = clusterFile;
    }

    public String getDistanceMatrixFile() {
        return DistanceMatrixFile;
    }

    public void setDistanceMatrixFile(String distanceMatrixFile) {
        DistanceMatrixFile = distanceMatrixFile;
    }

    public String getLabelFile() {
        return LabelFile;
    }

    public void setLabelFile(String labelFile) {
        LabelFile = labelFile;
    }

    public String getComparisonClusterFile() {
        return ComparisonClusterFile;
    }

    public void setComparisonClusterFile(String comparisonClusterFile) {
        ComparisonClusterFile = comparisonClusterFile;
    }

    public String getRestartClusterFile() {
        return RestartClusterFile;
    }

    public void setRestartClusterFile(String restartClusterFile) {
        RestartClusterFile = restartClusterFile;
    }

    public String getTimingFile() {
        return TimingFile;
    }

    public void setTimingFile(String timingFile) {
        TimingFile = timingFile;
    }

    public String getSummaryFile() {
        return SummaryFile;
    }

    public void setSummaryFile(String summaryFile) {
        SummaryFile = summaryFile;
    }
}




